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

#ifndef DSPFILTERS_BIQUAD_H
#define DSPFILTERS_BIQUAD_H

#include "DspFilters/Common.h"
#include "DspFilters/MathSupplement.h"
#include "DspFilters/Types.h"

namespace Dsp {

/*
 * Holds coefficients for a second order Infinite Impulse Response
 * digital filter. This is the building block for all IIR filters.
 *
 */

// Factored interface to prevent outsiders from fiddling
class BiquadBase : protected Layout {
 public:
  // Calculate filter response at the given normalized frequency.
  complex_t<double> response(const double& normalizedFrequency) const {
    const double a0 = getA0();
    const double a1 = getA1();
    const double a2 = getA2();
    const double b0 = getB0();
    const double b1 = getB1();
    const double b2 = getB2();

    const double w = 2 * doublePi * normalizedFrequency;
    const complex_t<double> czn1 = std::polar(1., -w);
    const complex_t<double> czn2 = std::polar(1., -2 * w);
    complex_t<double> ch(1);
    complex_t<double> cbot(1);

    complex_t<double> ct(b0 / a0);
    complex_t<double> cb(1);
    ct = addmul(ct, b1 / a0, czn1);
    ct = addmul(ct, b2 / a0, czn2);
    cb = addmul(cb, a1 / a0, czn1);
    cb = addmul(cb, a2 / a0, czn2);
    ch *= ct;
    cbot *= cb;

    return ch / cbot;
  }

  std::vector<PoleZero<double>> getPoleZeros() const {
    std::vector<PoleZero<double>> vpz;
    vpz.emplace_back(*this);
    return vpz;
  }

  double getA0() const { return m_a0; }
  double getA1() const { return m_a[1] * m_a0; }
  double getA2() const { return m_a[2] * m_a0; }
  double getB0() const { return m_b[0] * m_a0; }
  double getB1() const { return m_b[1] * m_a0; }
  double getB2() const { return m_b[2] * m_a0; }

  // Process a block of samples in the given form
  template <class StateType, typename Sample>
  void process(int numSamples, Sample* dest, StateType& state) const {
    while (--numSamples >= 0) {
      *dest = state.process(*dest, *this);
      dest++;
    }
  }

 protected:
  //
  // These are protected so you can't mess with RBJ biquads
  //

  void setCoefficients(const double& a0, const double& a1, const double& a2,
                       const double& b0, const double& b1, const double& b2) {
    // #ifndef NDEBUG
    //     assert (!Dsp::is_nan(a0));
    //     assert (!Dsp::is_nan(b0));
    //     assert (!Dsp::is_nan(a1));
    //     assert (!Dsp::is_nan(a2));
    //     assert (!Dsp::is_nan(b1));
    //     assert (!Dsp::is_nan(b2));
    // #endif
    m_a.resize(3, 1.);
    m_b.resize(3, 1.);

    if (!(Dsp::isnan(a0) || Dsp::isnan(a1) || Dsp::isnan(a2) ||
          Dsp::isnan(b0) || Dsp::isnan(b1) || Dsp::isnan(b2))) {
      m_a0 = a0;
      m_a[1] = a1 / a0;
      m_a[2] = a2 / a0;
      m_b[0] = b0 / a0;
      m_b[1] = b1 / a0;
      m_b[2] = b2 / a0;
    }
  }

  void setOnePole(const complex_t<double>& pole,
                  const complex_t<double>& zero) {
#if 0
  pole = adjust_imag (pole);
  zero = adjust_imag (zero);
#else
    assert(pole.imag() == 0);
    assert(zero.imag() == 0);
#endif

    const double a0 = 1;
    const double a1 = -pole.real();
    const double a2 = 0;
    const double b0 = 1;
    const double b1 = -zero.real();
    const double b2 = 0;

    setCoefficients(a0, a1, a2, b0, b1, b2);
  }

  void setTwoPole(const complex_t<double>& pole1,
                  const complex_t<double>& zero1,
                  const complex_t<double>& pole2,
                  const complex_t<double>& zero2) {
#if 0
  pole1 = adjust_imag (pole1);
  pole2 = adjust_imag (pole2);
  zero1 = adjust_imag (zero1);
  zero2 = adjust_imag (zero2);
#endif

    const double a0 = 1.;
    double a1;
    double a2;

    if (pole1.imag() != 0) {
      // assert (pole2 == std::conj (pole1));

      a1 = -2 * pole1.real();
      a2 = std::norm(pole1);
    } else {
      assert(pole2.imag() == 0);

      a1 = -(pole1.real() + pole2.real());
      a2 = pole1.real() * pole2.real();
    }

    const double b0 = 1;
    double b1;
    double b2;

    if (zero1.imag() != 0) {
      // assert (zero2 == std::conj (zero1));

      b1 = -2 * zero1.real();
      b2 = std::norm(zero1);
    } else {
      assert(zero2.imag() == 0);

      b1 = -(zero1.real() + zero2.real());
      b2 = zero1.real() * zero2.real();
    }

    setCoefficients(a0, a1, a2, b0, b1, b2);
  }

  void setPoleZero(const PoleZero<double>& pair) {
    if (pair.isSinglePole())
      setOnePole(pair.poles.first, pair.zeros.first);
    else
      setTwoPole(pair.poles.first, pair.zeros.first, pair.poles.second,
                 pair.zeros.second);
  }

  void setPoleZeroForm(const BiquadPoleState& bps) {
    setPoleZero(bps);
    applyScale(bps.gain);
  }

  void setIdentity() { setCoefficients(1., 0., 0., 1., 0., 0.); };

  void applyScale(const double& scale) {
    setCoefficients(m_a0, m_a[1], m_a[2], m_b[0] * scale, m_b[1] * scale,
                    m_b[2] * scale);
  };

 public:
  double m_a0;
  std::vector<double> m_a, m_b;
};

//------------------------------------------------------------------------------

// Expresses a biquad as a pair of pole/zeros, with gain
// values so that the coefficients can be reconstructed precisely.
struct BiquadPoleState : PoleZero<double> {
  BiquadPoleState() {}

  explicit BiquadPoleState(const BiquadBase& s) {
    const double a0 = s.getA0();
    const double a1 = s.getA1();
    const double a2 = s.getA2();
    const double b0 = s.getB0();
    const double b1 = s.getB1();
    const double b2 = s.getB2();

    if (a2 == 0 && b2 == 0) {
      // single pole
      poles.first = -a1;
      zeros.first = -b0 / b1;
      poles.second = 0;
      zeros.second = 0;
    } else {
      {
        const complex_t<double> c =
            std::sqrt(complex_t<double>(a1 * a1 - 4 * a0 * a2, 0));
        double d = 1. / (2. * a0);
        poles.first = -(a1 + c) * d;
        poles.second = (c - a1) * d;
        assert(!poles.isnan());
      }

      {
        const complex_t<double> c =
            sqrt(complex_t<double>(b1 * b1 - 4 * b0 * b2, 0));
        double d = 1. / (2. * b0);
        zeros.first = -(b1 + c) * d;
        zeros.second = (c - b1) * d;
        assert(!zeros.isnan());
      }
    }

    gain = b0 / a0;
  }

  double gain;
};

// More permissive interface for fooling around
class Biquad : public BiquadBase {
 public:
  Biquad() { setIdentity(); };

  explicit Biquad(const BiquadPoleState& bps) { setPoleZeroForm(bps); };

 public:
  // Process a block of samples, interpolating from the old section's
  // coefficients to this section's coefficients, over numSamples. This
  // implements smooth parameter changes.

  template <class StateType, typename Sample>
  void smoothProcess1(int numSamples, Sample* dest, StateType& state,
                      Biquad sectionPrev) const {
    double t = 1. / numSamples;
    double da1 = (m_a[1] - sectionPrev.m_a[1]) * t;
    double da2 = (m_a[2] - sectionPrev.m_a[2]) * t;
    double db0 = (m_b[0] - sectionPrev.m_b[0]) * t;
    double db1 = (m_b[1] - sectionPrev.m_b[1]) * t;
    double db2 = (m_b[2] - sectionPrev.m_b[2]) * t;

    while (--numSamples >= 0) {
      sectionPrev.m_a[1] += da1;
      sectionPrev.m_a[2] += da2;
      sectionPrev.m_b[0] += db0;
      sectionPrev.m_b[1] += db1;
      sectionPrev.m_b[2] += db2;

      *dest = state.process(*dest, sectionPrev);
      dest++;
    }
  }

  // Process a block of samples, interpolating from the old section's pole/zeros
  // to this section's pole/zeros, over numSamples. The interpolation is done
  // in the z-plane using polar coordinates.
  template <class StateType, typename Sample>
  void smoothProcess2(int numSamples, Sample* dest, StateType& state,
                      BiquadPoleState zPrev) const {
    BiquadPoleState z(*this);
    double t = 1. / numSamples;
    complex_t<double> dp0 = (z.poles.first - zPrev.poles.first) * t;
    complex_t<double> dp1 = (z.poles.second - zPrev.poles.second) * t;
    complex_t<double> dz0 = (z.zeros.first - zPrev.zeros.first) * t;
    complex_t<double> dz1 = (z.zeros.second - zPrev.zeros.second) * t;
    double dg = (z.gain - zPrev.gain) * t;

    while (--numSamples >= 0) {
      zPrev.poles.first += dp0;
      zPrev.poles.second += dp1;
      zPrev.zeros.first += dz0;
      zPrev.zeros.second += dz1;
      zPrev.gain += dg;

      *dest = state.process(*dest, Biquad(zPrev));
      dest++;
    }
  }

 public:
  // Export these as public

  void setOnePole(const complex_t<double>& pole,
                  const complex_t<double>& zero) {
    BiquadBase::setOnePole(pole, zero);
  }

  void setTwoPole(const complex_t<double>& pole1,
                  const complex_t<double>& zero1,
                  const complex_t<double>& pole2,
                  const complex_t<double>& zero2) {
    BiquadBase::setTwoPole(pole1, zero1, pole2, zero2);
  }

  void setPoleZero(const PoleZero<double>& pair) {
    BiquadBase::setPoleZero(pair);
  }

  void applyScale(const double& scale) { BiquadBase::applyScale(scale); }
};

}  // namespace Dsp

#endif
