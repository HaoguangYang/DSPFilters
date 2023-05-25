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

#ifndef DSPFILTERS_BESSEL_H
#define DSPFILTERS_BESSEL_H

#include "DspFilters/Cascade.h"
#include "DspFilters/Common.h"
#include "DspFilters/Design.h"
#include "DspFilters/Filter.h"
#include "DspFilters/PoleFilter.h"
#include "DspFilters/RootFinder.h"

namespace Dsp {

/*
 * Filters with Bessel response characteristics
 *
 */

namespace Bessel {

// A Workspace is necessary to find roots
struct WorkspaceBase {
  WorkspaceBase(RootFinderBase* rootsBase) : roots(*rootsBase) {}

  RootFinderBase& roots;

 private:
  WorkspaceBase(WorkspaceBase&);
  WorkspaceBase& operator=(WorkspaceBase&) { return *this; }
};

template <int MaxOrder>
struct Workspace : WorkspaceBase {
  Workspace() : WorkspaceBase(&m_roots) {}

 private:
  RootFinder<MaxOrder> m_roots;
};

//------------------------------------------------------------------------------
// returns the k-th zero based coefficient of the reverse bessel polynomial of
// degree n
static double reversebessel(const int& k, const int& n) {
  if (n < 0 || n < k) throw std::runtime_error("Unreasonable Bessel expansion");
  return static_cast<double>(factorial(2 * n - k)) /
         (static_cast<double>(factorial(n - k) * factorial(k)) *
          static_cast<double>(pow(2, n - k)));
}

// Half-band analog prototypes (s-plane)
class AnalogLowPass : public Layout {
 public:
  AnalogLowPass() { setNormal(0, 1); }

  void design(const int numPoles, WorkspaceBase* w) {
    if (getNumPoles() != numPoles) {
      clear();
      RootFinderBase& solver(w->roots);
      for (int i = 0; i < numPoles + 1; ++i)
        solver.coef()[i] = reversebessel(i, numPoles);
      solver.solve(numPoles);
      
      // FIXME: wrong formulation. refer to pole_filter.
      const int pairs = numPoles / 2;
      for (int i = 0; i < pairs; ++i) {
        complex_t<double> c = solver.root()[i];
        appendConjugatePairs(c, false);
      }

      if (numPoles % 2)
        append(solver.root()[pairs].real(), false);
    }
  }
};

//------------------------------------------------------------------------------
class AnalogLowShelf : public Layout {
 public:
  AnalogLowShelf() { setNormal(doublePi, 1); }

  void design(int numPoles, double gainDb, WorkspaceBase* w) {
    if (getNumPoles() != numPoles || m_gainDb != gainDb) {
      m_gainDb = gainDb;
      clear();

      const double G = std::pow(10., gainDb / 20) - 1.;
      RootFinderBase& poles(w->roots);
      for (int i = 0; i < numPoles + 1; ++i)
        poles.coef()[i] = reversebessel(i, numPoles);
      poles.solve(numPoles);

      RootFinder<50> zeros;
      for (int i = 0; i < numPoles + 1; ++i)
        zeros.coef()[i] = reversebessel(i, numPoles);
      double a0 = reversebessel(0, numPoles);
      zeros.coef()[0] += G * a0;
      zeros.solve(numPoles);

      // FIXME: wrong formulation. refer to pole_filter.
      const int pairs = numPoles / 2;
      for (int i = 0; i < pairs; ++i) {
        complex_t<double> p = poles.root()[i];
        complex_t<double> z = zeros.root()[i];
        appendConjugatePairs(p, z);
      }

      if (numPoles % 2)
        append(poles.root()[pairs].real(), zeros.root()[pairs].real());
    }
  }

 private:
  double m_gainDb;
};

//------------------------------------------------------------------------------

// Factored implementations to reduce template instantiations
struct LowPassBase : AnalogPoleFilterBase<AnalogLowPass> {
  void setup(int order, double sampleRate, double cutoffFrequency,
             WorkspaceBase* w) {
    m_analogProto.design(order, w);

    LowPassTransform(cutoffFrequency / sampleRate, this->m_digitalProto,
                         this->m_analogProto);

    Cascade::setLayout(this->m_digitalProto);
  }
};

struct HighPassBase : AnalogPoleFilterBase<AnalogLowPass> {
  void setup(int order, double sampleRate, double cutoffFrequency,
             WorkspaceBase* w) {
    m_analogProto.design(order, w);

    HighPassTransform(cutoffFrequency / sampleRate, this->m_digitalProto,
                          this->m_analogProto);

    Cascade::setLayout(this->m_digitalProto);
  }
};

struct BandPassBase : AnalogPoleFilterBase<AnalogLowPass> {
  void setup(int order, double sampleRate, double centerFrequency,
             double widthFrequency, WorkspaceBase* w) {
    m_analogProto.design(order, w);

    BandPassTransform(centerFrequency / sampleRate,
                          widthFrequency / sampleRate, this->m_digitalProto,
                          this->m_analogProto);

    Cascade::setLayout(this->m_digitalProto);
  }
};

struct BandStopBase : AnalogPoleFilterBase<AnalogLowPass> {
  void setup(int order, double sampleRate, double centerFrequency,
             double widthFrequency, WorkspaceBase* w) {
    m_analogProto.design(order, w);

    BandStopTransform(centerFrequency / sampleRate,
                          widthFrequency / sampleRate, this->m_digitalProto,
                          this->m_analogProto);

    Cascade::setLayout(this->m_digitalProto);
  }
};

struct LowShelfBase : AnalogPoleFilterBase<AnalogLowShelf> {
  void setup(int order, double sampleRate, double cutoffFrequency,
             double gainDb, WorkspaceBase* w) {
    m_analogProto.design(order, gainDb, w);

    LowPassTransform(cutoffFrequency / sampleRate, this->m_digitalProto,
                         this->m_analogProto);

    Cascade::setLayout(this->m_digitalProto);
  }
};

//------------------------------------------------------------------------------

//
// Raw filters
//

template <int MaxOrder>
struct LowPass : PoleFilter<LowPassBase, MaxOrder> {
  void setup(int order, double sampleRate, double cutoffFrequency) {
    Workspace<MaxOrder> w;
    LowPassBase::setup(order, sampleRate, cutoffFrequency, &w);
  }
};

template <int MaxOrder>
struct HighPass : PoleFilter<HighPassBase, MaxOrder> {
  void setup(int order, double sampleRate, double cutoffFrequency) {
    Workspace<MaxOrder> w;
    HighPassBase::setup(order, sampleRate, cutoffFrequency, &w);
  }
};

template <int MaxOrder>
struct BandPass : PoleFilter<BandPassBase, MaxOrder, MaxOrder * 2> {
  void setup(int order, double sampleRate, double centerFrequency,
             double widthFrequency) {
    Workspace<MaxOrder> w;
    BandPassBase::setup(order, sampleRate, centerFrequency, widthFrequency, &w);
  }
};

template <int MaxOrder>
struct BandStop : PoleFilter<BandStopBase, MaxOrder, MaxOrder * 2> {
  void setup(int order, double sampleRate, double centerFrequency,
             double widthFrequency) {
    Workspace<MaxOrder> w;
    BandStopBase::setup(order, sampleRate, centerFrequency, widthFrequency, &w);
  }
};

template <int MaxOrder>
struct LowShelf : PoleFilter<LowShelfBase, MaxOrder, MaxOrder * 2> {
  void setup(int order, double sampleRate, double cutoffFrequency,
             double gainDb) {
    Workspace<MaxOrder> w;
    LowShelfBase::setup(order, sampleRate, cutoffFrequency, gainDb, &w);
  }
};

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
  static const char* getName() { return "Bessel Low Pass"; }
};

struct HighPassDescription {
  static Kind getKind() { return kindHighPass; }
  static const char* getName() { return "Bessel High Pass"; }
};

struct BandPassDescription {
  static Kind getKind() { return kindBandPass; }
  static const char* getName() { return "Bessel Band Pass"; }
};

struct BandStopDescription {
  static Kind getKind() { return kindBandStop; }
  static const char* getName() { return "Bessel Band Stop"; }
};

struct LowShelfDescription {
  static Kind getKind() { return kindLowShelf; }
  static const char* getName() { return "Bessel Low Shelf"; }
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
// Gui-friendly Design layer
//

template <int MaxOrder>
struct LowPass : OrderBase<MaxOrder, TypeI, Bessel::LowPass>,
                 LowPassDescription {};

template <int MaxOrder>
struct HighPass : OrderBase<MaxOrder, TypeI, Bessel::HighPass>,
                  HighPassDescription {};

template <int MaxOrder>
struct BandPass : OrderBase<MaxOrder, TypeII, Bessel::BandPass>,
                  BandPassDescription {};

template <int MaxOrder>
struct BandStop : OrderBase<MaxOrder, TypeII, Bessel::BandStop>,
                  BandStopDescription {};

/*
 * NOT IMPLEMENTED
 *
 */
template <int MaxOrder>
struct LowShelf : OrderBase<MaxOrder, TypeIII, Bessel::LowShelf>,
                  LowShelfDescription {};

}  // namespace Design

}  // namespace Bessel

}  // namespace Dsp

#endif

/* This is a test of svn:external */
