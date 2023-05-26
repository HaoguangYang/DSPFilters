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

#ifndef DSPFILTERS_BUTTERWORTH_H
#define DSPFILTERS_BUTTERWORTH_H

#include "DspFilters/Cascade.h"
#include "DspFilters/Common.h"
#include "DspFilters/Design.h"
#include "DspFilters/Filter.h"
#include "DspFilters/PoleFilter.h"

namespace Dsp {

/*
 * Filters with Butterworth response characteristics
 *
 */

namespace Butterworth {

// Half-band analog prototypes (s-plane)
class AnalogLowPass : public Layout {
 public:
  AnalogLowPass() {
    setNormal(0, 1);
  }

  void design(const int numPoles) {
    if (getNumPoles() != numPoles) {
      clear();

      const double n2 = 2 * numPoles;
      const int pairs = numPoles / 2;
      for (int i = 0; i < pairs; ++i) {
        std::complex<double> c =
            std::polar(1., doublePi_2 + (2 * i + 1) * doublePi / n2);
        assert(!std::isinf(c.real()));
        assert(!std::isinf(c.imag()));
        appendConjugatePairs(c, false);
      }

      if (numPoles % 2) append(-1, false);
    }
  }
};

//------------------------------------------------------------------------------
class AnalogLowShelf : public Layout {
 public:
  AnalogLowShelf() {
    setNormal(doublePi, 1);
  }

  void design(const int& numPoles, const double& gainDb) {
    if (getNumPoles() != numPoles || m_gainDb != gainDb) {
      m_gainDb = gainDb;

      clear();

      const double n2 = numPoles * 2.;
      const double g = std::pow(std::pow(10., gainDb / 20), 1. / n2);
      const double gp = -1. / g;
      const double gz = -g;

      const int pairs = numPoles / 2;
      for (int i = 1; i <= pairs; ++i) {
        const double theta = doublePi * (0.5 - (2 * i - 1) / n2);
        appendConjugatePairs(std::polar(std::abs(gp), theta),
                                  std::polar(std::abs(gz), theta));
      }

      if (numPoles % 2) append(gp, gz);
    }
  }

 private:
  double m_gainDb;
};

//------------------------------------------------------------------------------

// Factored implementations to reduce template instantiations
struct LowPassBase : AnalogPoleFilterBase<AnalogLowPass> {
  void setup(int order, double sampleRate, double cutoffFrequency) {
    m_analogProto.design(order);
    LowPassTransform(cutoffFrequency / sampleRate, m_digitalProto,
                     m_analogProto);
    Cascade::setLayout(m_digitalProto);
  }
};

struct HighPassBase : AnalogPoleFilterBase<AnalogLowPass> {
  void setup(int order, double sampleRate, double cutoffFrequency) {
    m_analogProto.design(order);
    HighPassTransform(cutoffFrequency / sampleRate, m_digitalProto,
                      m_analogProto);
    Cascade::setLayout(m_digitalProto);
  }
};

struct BandPassBase : AnalogPoleFilterBase<AnalogLowPass> {
  void setup(int order, double sampleRate, double centerFrequency,
             double widthFrequency) {
    m_analogProto.design(order);
    BandPassTransform(centerFrequency / sampleRate, widthFrequency / sampleRate,
                      m_digitalProto, m_analogProto);
    Cascade::setLayout(m_digitalProto);
  }
};

struct BandStopBase : AnalogPoleFilterBase<AnalogLowPass> {
  void setup(int order, double sampleRate, double centerFrequency,
             double widthFrequency) {
    this->m_analogProto.design(order);
    BandStopTransform(centerFrequency / sampleRate, widthFrequency / sampleRate,
                      m_digitalProto, m_analogProto);
    Cascade::setLayout(m_digitalProto);
  }
};

struct LowShelfBase : AnalogPoleFilterBase<AnalogLowShelf> {
  void setup(int order, double sampleRate, double cutoffFrequency,
             double gainDb) {
    this->m_analogProto.design(order, gainDb);
    LowPassTransform(cutoffFrequency / sampleRate, m_digitalProto,
                     m_analogProto);
    Cascade::setLayout(m_digitalProto);
  }
};

struct HighShelfBase : AnalogPoleFilterBase<AnalogLowShelf> {
  void setup(int order, double sampleRate, double cutoffFrequency,
             double gainDb) {
    this->m_analogProto.design(order, gainDb);
    HighPassTransform(cutoffFrequency / sampleRate, m_digitalProto,
                      m_analogProto);
    Cascade::setLayout(m_digitalProto);
  }
};

struct BandShelfBase : AnalogPoleFilterBase<AnalogLowShelf> {
  void setup(int order, double sampleRate, double centerFrequency,
             double widthFrequency, double gainDb) {
    this->m_analogProto.design(order, gainDb);
    BandPassTransform(centerFrequency / sampleRate, widthFrequency / sampleRate,
                      m_digitalProto, m_analogProto);
    // HACK!
    // m_digitalProto.setNormal (((centerFrequency/sampleRate) < 0.25) ?
    // doublePi : 0, 1);
    Cascade::setLayout(m_digitalProto);
  }
};

//------------------------------------------------------------------------------

//
// Raw filters
//

template <int MaxOrder>
struct LowPass : PoleFilter<LowPassBase, MaxOrder> {};

template <int MaxOrder>
struct HighPass : PoleFilter<HighPassBase, MaxOrder> {};

template <int MaxOrder>
struct BandPass : PoleFilter<BandPassBase, MaxOrder, MaxOrder * 2> {};

template <int MaxOrder>
struct BandStop : PoleFilter<BandStopBase, MaxOrder, MaxOrder * 2> {};

template <int MaxOrder>
struct LowShelf : PoleFilter<LowShelfBase, MaxOrder> {};

template <int MaxOrder>
struct HighShelf : PoleFilter<HighShelfBase, MaxOrder> {};

template <int MaxOrder>
struct BandShelf : PoleFilter<BandShelfBase, MaxOrder, MaxOrder * 2> {};

//------------------------------------------------------------------------------

//
// Gui-friendly Design layer
//

namespace Design {

struct TypeIBase : DesignBase {
  enum { NumParams = 3 };

  static int getNumParams() { return 3; }

  static const ParamInfo getParamInfo_2() {
    return ParamInfo::defaultCutoffFrequencyParam();
  }
};

template <class FilterClass>
struct TypeI : TypeIBase, FilterClass {
  void setParams(const Params& params) {
    FilterClass::setup(int(params[1]), params[0], params[2]);
  }
};

struct TypeIIBase : DesignBase {
  enum { NumParams = 4 };

  static int getNumParams() { return 4; }

  static const ParamInfo getParamInfo_2() {
    return ParamInfo::defaultCenterFrequencyParam();
  }

  static const ParamInfo getParamInfo_3() {
    return ParamInfo::defaultBandwidthHzParam();
  }
};

template <class FilterClass>
struct TypeII : TypeIIBase, FilterClass {
  void setParams(const Params& params) {
    FilterClass::setup(int(params[1]), params[0], params[2], params[3]);
  }
};

struct TypeIIIBase : DesignBase {
  enum { NumParams = 4 };

  static int getNumParams() { return 4; }

  static const ParamInfo getParamInfo_2() {
    return ParamInfo::defaultCutoffFrequencyParam();
  }

  static const ParamInfo getParamInfo_3() {
    return ParamInfo::defaultGainParam();
  }
};

template <class FilterClass>
struct TypeIII : TypeIIIBase, FilterClass {
  void setParams(const Params& params) {
    FilterClass::setup(int(params[1]), params[0], params[2], params[3]);
  }
};

struct TypeIVBase : DesignBase {
  enum { NumParams = 5 };

  static int getNumParams() { return 5; }

  static const ParamInfo getParamInfo_2() {
    return ParamInfo::defaultCenterFrequencyParam();
  }

  static const ParamInfo getParamInfo_3() {
    return ParamInfo::defaultBandwidthHzParam();
  }

  static const ParamInfo getParamInfo_4() {
    return ParamInfo::defaultGainParam();
  }
};

template <class FilterClass>
struct TypeIV : TypeIVBase, FilterClass {
  void setParams(const Params& params) {
    FilterClass::setup(int(params[1]), params[0], params[2], params[3],
                       params[4]);
  }
};

// Factored kind and name

struct LowPassDescription {
  static Kind getKind() { return kindLowPass; }
  static const char* getName() { return "Butterworth Low Pass"; }
};

struct HighPassDescription {
  static Kind getKind() { return kindHighPass; }
  static const char* getName() { return "Butterworth High Pass"; }
};

struct BandPassDescription {
  static Kind getKind() { return kindBandPass; }
  static const char* getName() { return "Butterworth Band Pass"; }
};

struct BandStopDescription {
  static Kind getKind() { return kindBandStop; }
  static const char* getName() { return "Butterworth Band Stop"; }
};

struct LowShelfDescription {
  static Kind getKind() { return kindLowShelf; }
  static const char* getName() { return "Butterworth Low Shelf"; }
};

struct HighShelfDescription {
  static Kind getKind() { return kindHighShelf; }
  static const char* getName() { return "Butterworth High Shelf"; }
};

struct BandShelfDescription {
  static Kind getKind() { return kindBandShelf; }
  static const char* getName() { return "Butterworth Band Shelf"; }
};

// This glues on the Order parameter
template <int MaxOrder, template <class> class TypeClass,
          template <int> class FilterClass>
struct OrderBase : TypeClass<FilterClass<MaxOrder>> {
  const ParamInfo getParamInfo_1() const {
    return ParamInfo(idOrder, "Order", "Order", 1, MaxOrder, 2,
                     &ParamInfo::Int_toControlValue,
                     &ParamInfo::Int_toNativeValue, &ParamInfo::Int_toString);
  }
};

//------------------------------------------------------------------------------

//
// Design filters
//

template <int MaxOrder>
struct LowPass : OrderBase<MaxOrder, TypeI, Butterworth::LowPass>,
                 LowPassDescription {};

template <int MaxOrder>
struct HighPass : OrderBase<MaxOrder, TypeI, Butterworth::HighPass>,
                  HighPassDescription {};

template <int MaxOrder>
struct BandPass : OrderBase<MaxOrder, TypeII, Butterworth::BandPass>,
                  BandPassDescription {};

template <int MaxOrder>
struct BandStop : OrderBase<MaxOrder, TypeII, Butterworth::BandStop>,
                  BandStopDescription {};

template <int MaxOrder>
struct LowShelf : OrderBase<MaxOrder, TypeIII, Butterworth::LowShelf>,
                  LowShelfDescription {};

template <int MaxOrder>
struct HighShelf : OrderBase<MaxOrder, TypeIII, Butterworth::HighShelf>,
                   HighShelfDescription {};

template <int MaxOrder>
struct BandShelf : OrderBase<MaxOrder, TypeIV, Butterworth::BandShelf>,
                   BandShelfDescription {};

}  // namespace Design

}  // namespace Butterworth

}  // namespace Dsp

#endif

