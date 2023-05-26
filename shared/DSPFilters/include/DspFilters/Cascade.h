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

#ifndef DSPFILTERS_CASCADE_H
#define DSPFILTERS_CASCADE_H

#include <type_traits>

#include "DspFilters/Biquad.h"
#include "DspFilters/Common.h"
#include "DspFilters/Filter.h"
#include "DspFilters/Layout.h"
#include "DspFilters/MathSupplement.h"

namespace Dsp {

/*
 * Holds coefficients for a cascade of second order sections.
 *
 */

// Factored implementation to reduce template instantiations
class Cascade {
 public:
  template <class StateType>
  class StateBase : private DenormalPrevention<StateType::SignalType> {
   public:
    template <typename Sample>
    inline Sample process(
        const Sample& in, const Cascade& c) {
      StateType* state = m_stateArray;
      const Biquad* stage = c.m_stageArray;
      const typename StateType::SignalType vsa = ac();
      typename StateType::SignalType out = in;
      out = (state++)->process1(out, *stage++, vsa);
      for (int i = 0; i < c.m_numStages - 1; i++) {
        out = (state++)->process1(out, *stage++, typename StateType::SignalType());
      }
      return static_cast<Sample>(out);
    }

   protected:
    StateBase(StateType* stateArray) : m_stateArray(stateArray) {}

   protected:
    StateType* m_stateArray;
  };

  struct Stage : Biquad {};

  struct Storage {
    Storage(int maxStages_, Stage* stageArray_)
        : maxStages(maxStages_), stageArray(stageArray_) {}

    int maxStages;
    Stage* stageArray;
  };

  int getNumStages() const { return m_numStages; }

  const Stage& operator[](int index) {
    assert(index >= 0 && index <= m_numStages);
    return m_stageArray[index];
  }

 public:
  virtual ~Cascade() {}

  // Calculate filter response at the given normalized frequency.
  complex_t<double> response(double normalizedFrequency) const {
    double w = 2 * doublePi * normalizedFrequency;
    const complex_t<double> czn1 = std::polar(1., -w);
    const complex_t<double> czn2 = std::polar(1., -2 * w);
    complex_t<double> ch(1);
    complex_t<double> cbot(1);

    const Biquad* stage = m_stageArray;
    for (int i = m_numStages; --i >= 0; ++stage) {
      complex_t<double> cb(1);
      complex_t<double> ct(stage->getB0() / stage->getA0());
      ct = addmul(ct, stage->getB1() / stage->getA0(), czn1);
      ct = addmul(ct, stage->getB2() / stage->getA0(), czn2);
      cb = addmul(cb, stage->getA1() / stage->getA0(), czn1);
      cb = addmul(cb, stage->getA2() / stage->getA0(), czn2);
      ch *= ct;
      cbot *= cb;
    }

    return ch / cbot;
  }

  std::vector<PoleZero<double>> getPoleZeros() const {
    std::vector<PoleZero<double>> vpz;
    vpz.reserve(m_numStages);

    const Stage* stage = m_stageArray;
    for (int i = m_numStages; --i >= 0;) {
      BiquadPoleState<double> bps(*stage++);
      assert(!bps.isSinglePole() || i == 0);
      vpz.push_back(bps);
    }

    return vpz;
  }

  // Process a block of samples in the given form
  template <class StateType, typename Sample>
  void process(int numSamples, Sample* dest, StateType& state) const {
    while (--numSamples >= 0) {
      *dest = state.process(*dest, *this);
      dest++;
    }
  }

  template <class StateType, typename It>
  void process(It first, It last, StateType& state) const {
    for (; first != last; ++first) {
      *first = state.process(*first, *this);
    }
  }

 protected:
  Cascade() : m_numStages(0), m_maxStages(0), m_stageArray(0){};

  void setCascadeStorage(const Storage& storage) {
    m_numStages = 0;
    m_maxStages = storage.maxStages;
    m_stageArray = storage.stageArray;
  }

  virtual void applyScale(double scale) {
    // For higher order filters it might be helpful
    // to spread this factor between all the stages.
    assert(m_numStages > 0);
    m_stageArray->applyScale(scale);
  }

  void setLayout(const Layout& proto) {
    const int numPoles = proto.getNumPoles();

    Biquad* stage = m_stageArray;
    for (int i = 0; i < m_numStages; ++i, ++stage)
      stage->setPoleZero(proto[i]);

    applyScale(proto.getNormalGain() /
               std::abs(response(proto.getNormalW() / (2 * doublePi))));
  }

 private:
  int m_numStages;
  int m_maxStages;
  Stage* m_stageArray;
};

//------------------------------------------------------------------------------

// Storage for Cascade
template <int MaxStages>
class CascadeStages {
 public:
  template <class StateType>
  class State : public Cascade::StateBase<StateType> {
   public:
    State() : Cascade::StateBase<StateType>(m_states) {
      Cascade::StateBase<StateType>::m_stateArray = m_states;
      reset();
    }

    void reset() {
      StateType* state = m_states;
      for (int i = MaxStages; --i >= 0; ++state) state->reset();
    }

   private:
    StateType m_states[MaxStages];
  };

  /*@Internal*/
  Cascade::Storage getCascadeStorage() {
    return Cascade<FP>::Storage(MaxStages, m_stages);
  }

 private:
  Cascade::Stage m_stages[MaxStages];
};

}  // namespace Dsp

#endif
