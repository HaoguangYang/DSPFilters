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

#ifndef DSPFILTERS_FILTER_H
#define DSPFILTERS_FILTER_H

#include "DspFilters/Common.h"
#include "DspFilters/MathSupplement.h"
#include "DspFilters/Params.h"
#include "DspFilters/State.h"
#include "DspFilters/Types.h"

namespace Dsp {

/*
 * Filter
 *
 * Full abstraction of a digital IIR filter.
 * Supports run-time introspection and modulation of filter
 * parameters.
 *
 */

class Filter {
 public:
  virtual ~Filter(){};

  virtual Kind getKind() const = 0;

  virtual const std::string getName() const = 0;

  virtual int getNumParams() const = 0;

  virtual ParamInfo getParamInfo(int index) const = 0;

  Params getDefaultParams() const {
    Params params;
    params.clear();
    for (int i = 0; i < getNumParams(); ++i)
      params[i] = getParamInfo(i).getDefaultValue();
    return params;
  }

  const Params& getParams() const { return m_params; }

  double getParam(int paramIndex) const {
    assert(paramIndex >= 0 && paramIndex <= getNumParams());
    return m_params[paramIndex];
  }

  void setParam(int paramIndex, double nativeValue) {
    assert(paramIndex >= 0 && paramIndex <= getNumParams());
    m_params[paramIndex] = nativeValue;
    doSetParams(m_params);
  }

  int findParamId(int paramId) {
    int index = -1;
    for (int i = getNumParams(); --i >= 0;) {
      if (getParamInfo(i).getId() == paramId) {
        index = i;
        break;
      }
    }
    return index;
  }

  void setParamById(int paramId, double nativeValue) {
    for (int i = getNumParams(); --i >= 0;) {
      if (getParamInfo(i).getId() == paramId) {
        setParam(i, nativeValue);
        return;
      }
    }
    assert(0);
  }

  void setParams(const Params& parameters) {
    m_params = parameters;
    doSetParams(parameters);
  }

  // This makes a best-effort to pick up the values
  // of matching parameters from another set. It uses
  // the ParamID information to make the match.
  void copyParamsFrom(Dsp::Filter const* other) {
    // first, set reasonable defaults
    m_params = getDefaultParams();
    if (other) {
      // now loop
      for (int i = 0; i < getNumParams(); ++i) {
        const ParamInfo& paramInfo = getParamInfo(i);
        // find a match
        for (int j = 0; j < other->getNumParams(); ++j) {
          const ParamInfo& otherParamInfo = other->getParamInfo(j);
          if (paramInfo.getId() == otherParamInfo.getId()) {
            // match!
            m_params[i] = paramInfo.clamp(other->getParam(j));
            break;
          }
        }
      }
    }

    doSetParams(m_params);
  };

  virtual std::complex<double> response(double normalizedFrequency) const = 0;

  virtual int getNumChannels() = 0;
  virtual void reset() = 0;

  // virtual void process(int numSamples, double* const* arrayOfChannels) = 0;

 protected:
  virtual void doSetParams(const Params& parameters) = 0;

 private:
  Params m_params;
};

//------------------------------------------------------------------------------

/*
 * FilterDesign
 *
 * This container holds a filter Design (Gui-friendly layer) and
 * optionally combines it with the necessary state information to
 * process channel data.
 *
 */

template <typename FP>
struct DiscreteTransferFcn {
  std::vector<FP> numerator;
  std::vector<FP> denominator;
};

// Factored to reduce template instantiations
// DesignClass must inherit from Layout.
template <class DesignClass>
class FilterDesignBase : public Filter {
 public:
  Kind getKind() const { return m_design.getKind(); }

  const std::string getName() const { return m_design.getName(); }

  int getNumParams() const { return DesignClass::NumParams; }

  Params getDefaultParams() const { return m_design.getDefaultParams(); }

  ParamInfo getParamInfo(int index) const {
    switch (index) {
      case 0:
        return m_design.getParamInfo_0();
      case 1:
        return m_design.getParamInfo_1();
      case 2:
        return m_design.getParamInfo_2();
      case 3:
        return m_design.getParamInfo_3();
      case 4:
        return m_design.getParamInfo_4();
      case 5:
        return m_design.getParamInfo_5();
      case 6:
        return m_design.getParamInfo_6();
      case 7:
        return m_design.getParamInfo_7();
    };

    return ParamInfo();
  }

  // This method gets the designed filter's digital poles and zeros (in Z
  // plane), and cast its type from double to demanded type.
  const PoleZero<double>& getPoleZeros(const bool& sPlane = false) const {
    return m_design.getPoleZeros(sPlane);
  }

  template <typename FP>
  DiscreteTransferFcn<FP> getCoeffs() const {
    DiscreteTransferFcn<FP> ret;
    ret.numerator = m_design.getNumerator<FP>();
    ret.denominator = m_design.getDenominator<FP>();
  }

  // FIXME: not implemented yet.
  template <typename FP>
  std::pair<DiscreteTransferFcn<FP>, bool> smoothGetCoeffs() {
    DiscreteTransferFcn<FP> ret;
    ret.numerator = m_design.getNumerator<FP>();
    ret.denominator = m_design.getDenominator<FP>();
    return std::make_pair(ret, false);
  }

  complex_t response(double normalizedFrequency) const {
    return m_design.response(normalizedFrequency);
  }

 protected:
  void doSetParams(const Params& parameters) { m_design.setParams(parameters); }

 protected:
  DesignClass m_design;
};

// This class does barely anything now
template <class DesignClass, int Channels = 0,
          class StateType = DirectFormII<float>>
class FilterDesign : public FilterDesignBase<DesignClass> {
 public:
  FilterDesign() {}

  int getNumChannels() { return Channels; }

  void reset() { m_state.reset(); }

  void process(int numSamples, StateType::SignalType* const* arrayOfChannels) {
    m_state.process(numSamples, arrayOfChannels,
                    FilterDesignBase<DesignClass>::m_design);
  }

 protected:
  ChannelsState<Channels, typename DesignClass::template State<StateType>>
      m_state;
};

//------------------------------------------------------------------------------

/*
 * This container combines a raw filter with state information
 * so it can process channels. In order to set up the filter you
 * must call a setup function directly. Smooth changes are
 * not supported, but this class has a smaller footprint.
 *
 */

template <class FilterDesignClass, class StateType = DirectFormI<float>,
          size_t Channels = 0>
class FilterImpl {
 public:
  typedef StateType::SignalType SignalType;

  FilterImpl() {}

  template <typename... ArgsT>
  FilterImpl(ArgsT... params) {
    rawFilterInitImpl_({params...});
  }

  void setParam(const size_t& channelIndex, const size_t& paramIndex,
                const double& nativeValue) {
    if (paramIndex >= FilterDesignClass::getNumParams())
      throw std::runtime_error(
          "parameter index for raw filter exceeds RawFilter's parameters");
    rawFilters_[channelIndex].setParam(paramIndex, nativeValue);
  }

  void setParam(const size_t& channelIndex,
                const std::initializer_list<double>& params) {
    rawFilters_[channelIndex].setParam(params);
  }

  size_t getNumChannels() { return Channels; }

  void reset() {
    for (auto& i : state_) i.reset();
  }

  void resetChannel(size_t index) { state_[index].reset(); }

  void process(int numSamples, SignalType* const* arrayOfChannels) {
    // in-place processing of the input signal
    for (size_t n = 0; n < Channels; n++)
      state_[n].process(numSamples, arrayOfChannels[n], dtfs_[n]);
  }

  void smoothProcess(int numSamples, SignalType* const* arrayOfChannels) {
    // in-place processing of the input signal
    for (size_t n = 0; n < Channels; n++) {
      if (needsSmoothing_[n]) {
        for (int j = 0; j < numSamples; j++) {
          // Break up the input array to individual samples and process.
          // FIXME: implement this API in FilterDesign class
          auto res = rawFilters_[n].smoothGetCoeffs<SignalType>();
          dtfs_[n] = res.first;
          needsSmoothing_[n] = res.second;
          // TODO: perform a dtfs pre-update and post-update process to obtain
          // the smoothing delta, and add exponential decay to the smoothing
          // delta.
          state_[n].process(1, arrayOfChannels[n][j], dtfs_[n]);
        }
      } else {
        state_[n].process(numSamples, arrayOfChannels[n], dtfs_[n]);
        // TODO: add exponential decay to smoothing delta
      }
    }
  }

 protected:
  void rawFilterInitImpl_(const std::initializer_list<double>& params) {
    if (params.size() != Channels * FilterDesignClass::getNumParams()) {
      if (param.size() != FilterDesignClass::getNumParams())
        throw std::runtime_error(
            "Filter initialization parameters are not aligned! Expected %d or "
            "%d params, got %d.",
            FilterDesignClass::getNumParams(),
            Channels * FilterDesignClass::getNumParams(), params.size());
      else {
        // Initialize all filters to have the same parameters.
        for (size_t n = 0; n < Channels; n++)
          rawFilters_[n] = FilterDesignClass(params);
      }
    } else {
      for (size_t n = 0; n < Channels; n++) {
        rawFilters_[n] = FilterDesignClass(
            {params.begin() + n * FilterDesignClass::getNumParams(),
             (params.begin() + 1) + n * FilterDesignClass::getNumParams()});
      }
    }

    for (size_t n = 0; n < Channels; n++) {
      dtfs_[n] = rawFilters_[n].getCoeffs<SignalType>();
    }
  }

  std::array<StateType, Channels> state_;
  std::array<DiscreteTransferFcn<SignalType>, Channels> dtfs_;
  std::array<bool, Channels> needsSmoothing_{false};
  std::array<SignalType, Chanels> smoothDelta_{SignalType(0.)};
  std::array<FilterDesignClass, Channels> rawFilters_;
};

// Filter class for analysis only
template <class FilterDesignClass, class StateType>
class FilterImpl<FilterDesignClass, StateType, 0> {
 public:
  typedef StateType::SignalType SignalType;

  FilterImpl() {}

  template <typename... ArgsT>
  FilterImpl(ArgsT... params) {
    rawFilter_ = FilterDesignClass({params...});
    dtf_ = rawFilter_.getCoeffs<SignalType>()
  }

  size_t getNumChannels() { return 0; }

  void reset() {
    throw std::logic_error("attempt to reset empty ChannelState");
  }

  void resetChannel(size_t index) {
    (void)index;
    throw std::logic_error("attempt to reset empty ChannelState");
  }

  void process(int numSamples, SignalType* const* arrayOfChannels) {
    (void)numSamples;
    (void)arrayOfChannels;
    throw std::logic_error("attempt to process empty ChannelState");
  }

  void smoothProcess(int numSamples, SignalType* const* arrayOfChannels) {
    (void)numSamples;
    (void)arrayOfChannels;
    throw std::logic_error("attempt to process empty ChannelState");
  }

 protected:
  DiscreteTransferFcn<SignalType> dtf_;
  FilterDesignClass rawFilter_;
};

}  // namespace Dsp

#endif
