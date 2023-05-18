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
template <typename FP>
class BiquadBase
{
public:
  template <class StateType>
  struct State : StateType, private DenormalPrevention
  {
    template <typename Sample>
    inline Sample process (const Sample& in, const BiquadBase<FP>& b)
    {
      return static_cast<Sample> (StateType::process1 (in, b, ac()));
    }
  };

public:
  // Calculate filter response at the given normalized frequency.
  std::complex<FP> response (const FP& normalizedFrequency) const;

  std::vector<PoleZeroPair<FP>> getPoleZeros () const;

  double getA0 () const { return m_a0; }
  double getA1 () const { return m_a[1]*m_a0; }
  double getA2 () const { return m_a[2]*m_a0; }
  double getB0 () const { return m_b[0]*m_a0; }
  double getB1 () const { return m_b[1]*m_a0; }
  double getB2 () const { return m_b[2]*m_a0; }

  // Process a block of samples in the given form
  template <class StateType, typename Sample>
  void process (int numSamples, Sample* dest, StateType& state) const
  {
    while (--numSamples >= 0) {
      *dest = state.process (*dest, *this);
      dest++;
    }
  }

protected:
  //
  // These are protected so you can't mess with RBJ biquads
  //

  void setCoefficients (const FP& a0, const FP& a1, const FP& a2,
                        const FP& b0, const FP& b1, const FP& b2);

  void setOnePole (const std::complex<FP>& pole, const std::complex<FP>& zero);

  void setTwoPole (const std::complex<FP>& pole1, const std::complex<FP>& zero1,
                   const std::complex<FP>& pole2, const std::complex<FP>& zero2);

  void setPoleZeroPair (const PoleZeroPair<FP>& pair)
  {
    if (pair.isSinglePole ())
      setOnePole (pair.poles.first, pair.zeros.first);
    else
      setTwoPole (pair.poles.first, pair.zeros.first,
                  pair.poles.second, pair.zeros.second);
  }

  void setPoleZeroForm (const BiquadPoleState<FP>& bps)
  {
    setPoleZeroPair (bps);
    applyScale (bps.gain);
  }

  void setIdentity () { setCoefficients (1., 0., 0., 1., 0., 0.); };

  void applyScale (const FP& scale) { setCoefficients(m_a0, m_a[1], m_a[2], m_b[0] * scale, m_b[1] * scale, m_b[2] * scale); };

public:
  FP m_a0;
  std::vector<FP> m_a, m_b;
};

//------------------------------------------------------------------------------

// Expresses a biquad as a pair of pole/zeros, with gain
// values so that the coefficients can be reconstructed precisely.
template <typename FP>
struct BiquadPoleState : PoleZeroPair<FP>
{
  BiquadPoleState () { }

  explicit BiquadPoleState (const BiquadBase<FP>& s);

  FP gain;
};

// More permissive interface for fooling around
template <typename FP>
class Biquad : public BiquadBase<FP>
{
public:
  Biquad () { setIdentity (); };

  explicit Biquad (const BiquadPoleState<FP>& bps) { setPoleZeroForm (bps); };

public:
  // Process a block of samples, interpolating from the old section's coefficients
  // to this section's coefficients, over numSamples. This implements smooth
  // parameter changes.

  template <class StateType, typename Sample>
  void smoothProcess1 (int numSamples,
                       Sample* dest,
                       StateType& state,
                       Biquad<FP> sectionPrev) const 
  {
    FP t = 1. / numSamples;
    FP da1 = (m_a[1] - sectionPrev.m_a[1]) * t;
    FP da2 = (m_a[2] - sectionPrev.m_a[2]) * t;
    FP db0 = (m_b[0] - sectionPrev.m_b[0]) * t;
    FP db1 = (m_b[1] - sectionPrev.m_b[1]) * t;
    FP db2 = (m_b[2] - sectionPrev.m_b[2]) * t;

    while (--numSamples >= 0)
    {
      sectionPrev.m_a[1] += da1;
      sectionPrev.m_a[2] += da2;
      sectionPrev.m_b[0] += db0;
      sectionPrev.m_b[1] += db1;
      sectionPrev.m_b[2] += db2;

      *dest = state.process (*dest, sectionPrev);
      dest++;
    }
  }

  // Process a block of samples, interpolating from the old section's pole/zeros
  // to this section's pole/zeros, over numSamples. The interpolation is done
  // in the z-plane using polar coordinates.
  template <class StateType, typename Sample>
  void smoothProcess2 (int numSamples,
                       Sample* dest,
                       StateType& state,
                       BiquadPoleState<FP> zPrev) const 
  {
    BiquadPoleState<FP> z (*this);
    FP t = 1. / numSamples;
    std::complex<FP> dp0 = (z.poles.first  - zPrev.poles.first) * t;
    std::complex<FP> dp1 = (z.poles.second - zPrev.poles.second) * t;
    std::complex<FP> dz0 = (z.zeros.first  - zPrev.zeros.first) * t;
    std::complex<FP> dz1 = (z.zeros.second - zPrev.zeros.second) * t;
    FP dg = (z.gain - zPrev.gain) * t;

    while (--numSamples >= 0)
    {
      zPrev.poles.first += dp0;
      zPrev.poles.second += dp1;
      zPrev.zeros.first += dz0;
      zPrev.zeros.second += dz1;
      zPrev.gain += dg;

      *dest = state.process (*dest, Biquad<FP> (zPrev));
      dest++;
    }
  }

public:
  // Export these as public

  void setOnePole (const std::complex<FP>& pole, const std::complex<FP>& zero)
  {
    BiquadBase<FP>::setOnePole (pole, zero);
  }

  void setTwoPole (const std::complex<FP>& pole1, const std::complex<FP>& zero1,
                   const std::complex<FP>& pole2, const std::complex<FP>& zero2)
  {
    BiquadBase<FP>::setTwoPole (pole1, zero1, pole2, zero2);
  }

  void setPoleZeroPair (const PoleZeroPair<FP>& pair)
  {
    BiquadBase<FP>::setPoleZeroPair (pair);
  }

  void applyScale (const FP& scale)
  {
    BiquadBase<FP>::applyScale (scale);
  }
};

}

#endif
