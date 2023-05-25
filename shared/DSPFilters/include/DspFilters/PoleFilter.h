/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_POLEFILTER_H
#define DSPFILTERS_POLEFILTER_H

#include "DspFilters/Cascade.h"
#include "DspFilters/Common.h"
#include "DspFilters/MathSupplement.h"

namespace Dsp {

/*
 * Base for filters designed via algorithmic placement of poles and zeros.
 *
 * Typically, the filter is first designed as a half-band low pass or
 * low shelf analog filter (s-plane). Then, using a transformation such
 * as the ones from Constantinides, the poles and zeros of the analog filter
 * are calculated in the z-plane.
 *
 */

// Factored implementations to reduce template instantiations

// Not sure we need this anymore.
class DigitalPoleFilterBase : public Cascade {
 public:
  // This gets the poles/zeros directly from the digital
  // prototype. It is used to double check the correctness
  // of the recovery of pole/zeros from biquad coefficients.
  //
  // It can also be used to accelerate the interpolation
  // of pole/zeros for parameter modulation, since a pole
  // filter already has them calculated

#if 1
  // Commenting this out will pass the call to the Cascade,
  // which tries to compute the poles and zeros from the biquad
  // coefficients.
  const PoleZero<double>& getPoleZeros() const {
    return static_cast<PoleZero<double>>(m_digitalProto);
  }
#endif

 protected:
  Layout m_digitalProto;
};

// Serves a container to hold the analog prototype
// and the digital pole/zero layout.
template <class AnalogPrototype>
class AnalogPoleFilterBase : public DigitalPoleFilterBase {
 protected:
  void setPrototypeStorage(const LayoutBase<double>& analogStorage,
                           const LayoutBase<double>& digitalStorage) {
    m_analogProto.reserve(analogStorage);
    this->m_digitalProto = digitalStorage;
  }

 protected:
  AnalogPrototype m_analogProto;
};

//------------------------------------------------------------------------------

// Storage for pole filters
template <class BaseClass, int MaxAnalogPoles,
          int MaxDigitalPoles = MaxAnalogPoles>
struct PoleFilter : BaseClass, CascadeStages<(MaxDigitalPoles + 1) / 2> {
  PoleFilter() {
    // This glues together the factored base classes
    // with the templatized storage classes.
    BaseClass::setCascadeStorage(this->getCascadeStorage());
    BaseClass::setPrototypeStorage(m_analogStorage, m_digitalStorage);
  }

 private:
  Layout m_analogStorage;
  Layout m_digitalStorage;
};

//------------------------------------------------------------------------------

/*
 * s-plane to z-plane transforms
 *
 * For pole filters, an analog prototype is created via placement of
 * poles and zeros in the s-plane. The analog prototype is either
 * a halfband low pass or a halfband low shelf. The poles, zeros,
 * and normalization parameters are transformed into the z-plane
 * using variants of the bilinear transformation.
 *
 */

class C2DTransformBase {
 public:
  // transforming a layout
  void c2d(const Layout& analog, Layout& digital) {
    for (const auto& pole : analog.poles()) {
      auto transformedPole = transform(pole.first);
      if (pole.order() == 1)
        digital.append(transformedPole, false);
      else if (pole.order() == 2){
        if (pole.isReal()){
          digital.append(transformedPole, false);
          transformedPole = transform(pole.second);
          digital.append(transformedPole, false);
        } else {
          digital.appendConjugatePairs(transformedPole, false);
        }
      }
    }
    for (const auto& zero : analog.zeros()) {
      auto transformedZero = transform(zero.first);
      if (zero.order() == 1)
        digital.append(transformedZero, true);
      else if (zero.order() == 2){
        if (zero.isReal()){
          digital.append(transformedZero, true);
          transformedZero = transform(zero.second);
          digital.append(transformedZero, true);
        } else {
          digital.appendConjugatePairs(transformedZero, true);
        }
      }
    }
  }

 protected:
  // transforming a single pole or zero
  virtual std::vector<std::complex<double>> transform(const std::complex<double>& x) = 0;
};

// low pass to low pass
class LowPassTransform : public C2DTransformBase {
 public:
  LowPassTransform(double fc, Layout& digital, const Layout& analog) {
    digital.clear();
    // prewarp
    f = std::tan(doublePi * fc);
    c2d(analog, digital);
    digital.setNormal(analog.getNormalW(), analog.getNormalGain());
  }

 protected:
  std::vector<std::complex<double>> transform(const std::complex<double>& c) override {
    std::vector<std::complex<double>> ret{};
    if (c == infinity<double>()){
      ret.emplace_back(-1., 0.);
      return ret;
    }
    // frequency transform
    auto c1 = f * c;
    // bilinear low pass transform
    ret.emplace_back((1. + c1) / (1. - c1));
    return ret;
  }

 private:
  double f;
};

//------------------------------------------------------------------------------

// low pass to high pass
class HighPassTransform : public C2DTransformBase {
 public:
  HighPassTransform(double fc, Layout& digital, const Layout& analog) {
    digital.clear();
    // prewarp
    f = 1. / std::tan(doublePi * fc);
    c2d(analog, digital);
    digital.setNormal(doublePi - analog.getNormalW(), analog.getNormalGain());
  }

 protected:
  std::vector<std::complex<double>> transform(const std::complex<double>& c) override {
    std::vector<std::complex<double>> ret{};
    if (c == infinity<double>()){
      ret.emplace_back(1., 0.);
      return ret;
    }
    // frequency transform
    auto c1 = f * c;
    // bilinear high pass transform
    ret.emplace_back(-(1. + c1) / (1. - c1));
    return ret;
  }

 private:
  double f;
};

//------------------------------------------------------------------------------

// low pass to band pass transform
class BandPassTransform : DenormalPrevention<double>, public C2DTransformBase {
 public:
  BandPassTransform(double fc, double fw, Layout& digital,
                    const Layout& analog) {
    // handle degenerate cases efficiently
    // THIS DOESNT WORK because the cascade states won't match
#if 0
    const double fw_2 = fw / 2;
    if (fc - fw_2 < 0) {
      LowPassTransform::LowPassTransform (fc + fw_2, digital, analog);
    } else if (fc + fw_2 >= 0.5) {
      HighPassTransform::HighPassTransform (fc - fw_2, digital, analog);
    } else
#endif
    digital.clear();
    const double ww = 2 * doublePi * fw;
    // pre-calcs
    wc2 = 2 * doublePi * fc - (ww / 2);
    wc = wc2 + ww;

    // what is this crap?
    wc2 = std::clamp(wc2, dc(), doublePi - dc());

    a = std::cos((wc + wc2) * 0.5) / std::cos((wc - wc2) * 0.5);
    b = 1. / std::tan((wc - wc2) * 0.5);
    a2 = a * a;
    b2 = b * b;
    ab = a * b;
    ab_2 = 2 * ab;

    c2d(analog, digital);

    double wn = analog.getNormalW();
    digital.setNormal(
        2 * std::atan(std::sqrt(std::tan((wc + wn) * 0.5) * std::tan((wc2 + wn) * 0.5))),
        analog.getNormalGain());
  }

 protected:
  std::vector<std::complex<double>> transform(const std::complex<double>& c) override {
    if (c == infinity<double>()) return std::vector<std::complex<double>>{-1., 1.};

    auto c1 = (1. + c) / (1. - c);  // bilinear

    std::complex<double> v = 0;
    v = addmul(v, 4 * (b2 * (a2 - 1) + 1), c1);
    v += 8 * (b2 * (a2 - 1) - 1);
    v *= c1;
    v += 4 * (b2 * (a2 - 1) + 1);
    v = std::sqrt(v);

    std::complex<double> u = -v;
    u = addmul(u, ab_2, c1);
    u += ab_2;

    v = addmul(v, ab_2, c1);
    v += ab_2;

    std::complex<double> d = 0;
    d = addmul(d, 2 * (b - 1), c1) + 2 * (1 + b);

    return std::vector<std::complex<double>>{u / d, v / d};
  }

 private:
  double wc, wc2, a, b, a2, b2, ab, ab_2;
};

//------------------------------------------------------------------------------

// low pass to band stop transform
class BandStopTransform : DenormalPrevention<double>, public C2DTransformBase {
 public:
  BandStopTransform(double fc, double fw, Layout& digital,
                    const Layout& analog) {
    digital.clear();

    const double ww = 2 * doublePi * fw;

    wc2 = 2 * doublePi * fc - (ww / 2);
    wc = wc2 + ww;

    // this is crap
    if (wc2 < dc()) wc2 = dc();
    if (wc > doublePi - dc()) wc = doublePi - dc();

    a = std::cos((wc + wc2) * .5) / std::cos((wc - wc2) * .5);
    b = std::tan((wc - wc2) * .5);
    a2 = a * a;
    b2 = b * b;

    c2d(analog, digital);

    if (fc < 0.25)
      digital.setNormal(doublePi, analog.getNormalGain());
    else
      digital.setNormal(0, analog.getNormalGain());
  }

 private:
  std::vector<std::complex<double>> transform(const std::complex<double>& c) override {
    std::complex<double> c1 = -1.;
    if (c != infinity<double>())
      c1 = (1. + c) / (1. - c);  // bilinear

    std::complex<double> u(0);
    u = addmul(u, 4 * (b2 + a2 - 1), c);
    u += 8 * (b2 - a2 + 1);
    u *= c;
    u += 4 * (a2 + b2 - 1);
    u = std::sqrt(u);

    std::complex<double> v = u * -.5;
    v += a;
    v = addmul(v, -a, c);

    u *= .5;
    u += a;
    u = addmul(u, -a, c);

    std::complex<double> d(b + 1);
    d = addmul(d, b - 1, c);

    return std::vector<std::complex<double>>{u / d, v / d};
  }

  double wc, wc2, a, b, a2, b2;
};

}  // namespace Dsp

#endif
