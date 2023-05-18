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

#include "DspFilters/Common.h"
#include "DspFilters/Cascade.h"
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
template <typename FP>
struct WorkspaceBase
{
  WorkspaceBase (RootFinderBase<FP>* rootsBase)
    : roots (*rootsBase)
  {
  }

  RootFinderBase<FP>& roots;

private:
  WorkspaceBase (WorkspaceBase<FP>&);
  WorkspaceBase<FP>& operator= (WorkspaceBase<FP>&);
};

template <int MaxOrder, typename FP>
struct Workspace : WorkspaceBase<FP>
{
  Workspace ()
    : WorkspaceBase (&m_roots)
  {
  }

private:
  RootFinder<MaxOrder, FP> m_roots;
};

//------------------------------------------------------------------------------

// Half-band analog prototypes (s-plane)
template<typename FP>
class AnalogLowPass : public LayoutBase
{
public:
  AnalogLowPass ();

  void design (const int numPoles,
               WorkspaceBase<FP>* w);
};

//------------------------------------------------------------------------------
template <typename FP>
class AnalogLowShelf : public LayoutBase
{
public:
  AnalogLowShelf ();

  void design (int numPoles,
               FP gainDb,
               WorkspaceBase<FP>* w);

private:
  double m_gainDb;
};

//------------------------------------------------------------------------------

// Factored implementations to reduce template instantiations
template<typename FP>
struct LowPassBase : AnalogPoleFilterBase <AnalogLowPass<FP>>
{
  void setup (int order,
              FP sampleRate,
              FP cutoffFrequency,
              WorkspaceBase<FP>* w);
};

template<typename FP>
struct HighPassBase : AnalogPoleFilterBase <AnalogLowPass<FP>>
{
  void setup (int order,
              FP sampleRate,
              FP cutoffFrequency,
              WorkspaceBase<FP>* w);
};

template<typename FP>
struct BandPassBase : AnalogPoleFilterBase <AnalogLowPass<FP>>
{
  void setup (int order,
              FP sampleRate,
              FP centerFrequency,
              FP widthFrequency,
              WorkspaceBase* w);
};

template<typename FP>
struct BandStopBase : AnalogPoleFilterBase <AnalogLowPass<FP>>
{
  void setup (int order,
              FP sampleRate,
              FP centerFrequency,
              FP widthFrequency,
              WorkspaceBase* w);
};

template<typename FP>
struct LowShelfBase : AnalogPoleFilterBase <AnalogLowShelf<FP>>
{
  void setup (int order,
              FP sampleRate,
              FP cutoffFrequency,
              FP gainDb,
              WorkspaceBase* w);
};

//------------------------------------------------------------------------------

//
// Raw filters
//

template <int MaxOrder>
struct LowPass : PoleFilter <LowPassBase, MaxOrder>
{
  void setup (int order,
              double sampleRate,
              double cutoffFrequency)
  {
    Workspace <MaxOrder> w;
    LowPassBase::setup (order,
                        sampleRate,
                        cutoffFrequency,
                        &w);
  }
};

template <int MaxOrder>
struct HighPass : PoleFilter <HighPassBase, MaxOrder>
{
  void setup (int order,
              double sampleRate,
              double cutoffFrequency)
  {
    Workspace <MaxOrder> w;
    HighPassBase::setup (order,
                         sampleRate,
                         cutoffFrequency,
                         &w);
  }
};

template <int MaxOrder>
struct BandPass : PoleFilter <BandPassBase, MaxOrder, MaxOrder*2>
{
  void setup (int order,
              double sampleRate,
              double centerFrequency,
              double widthFrequency)
  {
    Workspace <MaxOrder> w;
    BandPassBase::setup (order,
                         sampleRate,
                         centerFrequency,
                         widthFrequency,
                         &w);
  }
};

template <int MaxOrder>
struct BandStop : PoleFilter <BandStopBase, MaxOrder, MaxOrder*2>
{
  void setup (int order,
              double sampleRate,
              double centerFrequency,
              double widthFrequency)
  {
    Workspace <MaxOrder> w;
    BandStopBase::setup (order,
                         sampleRate,
                         centerFrequency,
                         widthFrequency,
                         &w);
  }
};

template <int MaxOrder>
struct LowShelf : PoleFilter <LowShelfBase, MaxOrder, MaxOrder*2>
{
  void setup (int order,
              double sampleRate,
              double cutoffFrequency,
              double gainDb)
  {
    Workspace <MaxOrder> w;
    LowShelfBase::setup (order,
                         sampleRate,
                         cutoffFrequency,
                         gainDb,
                         &w);
  }
};

//------------------------------------------------------------------------------

//
// Gui-friendly Design layer
//

namespace Design {

struct TypeIBase : DesignBase
{
  enum
  {
    NumParams = 3
  };

  static int getNumParams ()
  {
    return 3;
  }

  static const ParamInfo getParamInfo_2 ()
  {
    return ParamInfo::defaultCutoffFrequencyParam ();
  }
};

template <class FilterClass>
struct TypeI : TypeIBase, FilterClass
{
  void setParams (const Params& params)
  {
    FilterClass::setup (int(params[1]), params[0], params[2]);
  }
};

struct TypeIIBase : DesignBase
{
  enum
  {
    NumParams = 4
  };

  static int getNumParams ()
  {
    return 4;
  }

  static const ParamInfo getParamInfo_2 ()
  {
    return ParamInfo::defaultCenterFrequencyParam ();
  }

  static const ParamInfo getParamInfo_3 ()
  {
    return ParamInfo::defaultBandwidthHzParam ();
  }
};

template <class FilterClass>
struct TypeII : TypeIIBase, FilterClass
{
  void setParams (const Params& params)
  {
    FilterClass::setup (int(params[1]), params[0], params[2], params[3]);
  }
};

struct TypeIIIBase : DesignBase
{
  enum
  {
    NumParams = 4
  };

  static int getNumParams ()
  {
    return 4;
  }

  static const ParamInfo getParamInfo_2 ()
  {
    return ParamInfo::defaultCutoffFrequencyParam ();
  }

  static const ParamInfo getParamInfo_3 ()
  {
    return ParamInfo::defaultGainParam ();
  }
};

template <class FilterClass>
struct TypeIII : TypeIIIBase, FilterClass
{
  void setParams (const Params& params)
  {
    FilterClass::setup (int(params[1]),
                        params[0],
                        params[2],
                        params[3]);
  }
};

struct TypeIVBase : DesignBase
{
  enum
  {
    NumParams = 5
  };

  static int getNumParams ()
  {
    return 5;
  }

  static const ParamInfo getParamInfo_2 ()
  {
    return ParamInfo::defaultCenterFrequencyParam ();
  }

  static const ParamInfo getParamInfo_3 ()
  {
    return ParamInfo::defaultBandwidthHzParam ();
  }

  static const ParamInfo getParamInfo_4 ()
  {
    return ParamInfo::defaultGainParam ();
  }
};

template <class FilterClass>
struct TypeIV : TypeIVBase, FilterClass
{
  void setParams (const Params& params)
  {
    FilterClass::setup (int(params[1]), params[0], params[2], params[3], params[4]);
  }
};

// Factored kind and name

struct LowPassDescription
{
  static Kind getKind () { return kindLowPass; }
  static const char* getName() { return "Bessel Low Pass"; }
};

struct HighPassDescription
{
  static Kind getKind () { return kindHighPass; }
  static const char* getName() { return "Bessel High Pass"; }
};

struct BandPassDescription
{
  static Kind getKind () { return kindBandPass; }
  static const char* getName() { return "Bessel Band Pass"; }
};

struct BandStopDescription
{
  static Kind getKind () { return kindBandStop; }
  static const char* getName() { return "Bessel Band Stop"; }
};

struct LowShelfDescription
{
  static Kind getKind () { return kindLowShelf; }
  static const char* getName() { return "Bessel Low Shelf"; }
};

// This glues on the Order parameter
template <int MaxOrder,
          template <class> class TypeClass,
          template <int> class FilterClass>
struct OrderBase : TypeClass <FilterClass <MaxOrder> >
{
  const ParamInfo getParamInfo_1 () const
  {
    return ParamInfo (idOrder, "Order", "Order",
                       1, MaxOrder, 2,
                       &ParamInfo::Int_toControlValue,
                       &ParamInfo::Int_toNativeValue,
                       &ParamInfo::Int_toString);

  }
};

//------------------------------------------------------------------------------

//
// Gui-friendly Design layer
//

template <int MaxOrder>
struct LowPass : OrderBase <MaxOrder, TypeI, Bessel::LowPass>,
                 LowPassDescription
{
};

template <int MaxOrder>
struct HighPass : OrderBase <MaxOrder, TypeI, Bessel::HighPass>,
                  HighPassDescription
{
};

template <int MaxOrder>
struct BandPass : OrderBase <MaxOrder, TypeII, Bessel::BandPass>,
                  BandPassDescription
{
};

template <int MaxOrder>
struct BandStop : OrderBase <MaxOrder, TypeII, Bessel::BandStop>,
                  BandStopDescription
{
};

/*
 * NOT IMPLEMENTED
 *
 */
template <int MaxOrder>
struct LowShelf : OrderBase <MaxOrder, TypeIII, Bessel::LowShelf>,
                  LowShelfDescription
{
};

}

}

}

#endif

/* This is a test of svn:external */
