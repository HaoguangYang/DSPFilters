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

#ifndef DSPFILTERS_TYPES_H
#define DSPFILTERS_TYPES_H

#include "DspFilters/Common.h"
#include "DspFilters/MathSupplement.h"

namespace Dsp {

// A conjugate or real pair
template<typename FP>
struct ComplexPair : complex_pair_t<FP>
{
  ComplexPair ()
  {
  }

  explicit ComplexPair (const std::complex<FP>& c1)
    : complex_pair_t<FP> (c1, 0.)
  {
    assert (isReal());
  }

  ComplexPair (const std::complex<FP>& c1,
               const std::complex<FP>& c2)
    : complex_pair_t<FP> (c1, c2) {}

  bool isConjugate () const
  {
    return second == std::conj (first);
  }

  bool isReal () const
  {
    return first.imag() == 0 && second.imag() == 0;
  }

  // Returns true if this is either a conjugate pair,
  // or a pair of reals where neither is zero.
  bool isMatchedPair () const
  {
    if (first.imag() != 0)
      return second == std::conj (first);
    else
      return second.imag () == 0 &&
             second.real () != 0 &&
             first.real () != 0;
  }

  bool isnan () const
  {
    return Dsp::isnan (first) || Dsp::isnan (second);
  }
};

// A pair of pole/zeros. This fits in a biquad (but is missing the gain)
template<typename FP>
struct PoleZeroPair
{
  ComplexPair<FP> poles;
  ComplexPair<FP> zeros;

  PoleZeroPair () { }

  // single pole/zero
  PoleZeroPair (const std::complex<FP>& p, const std::complex<FP>& z)
    : poles (p), zeros (z)
  {
  }

  // pole/zero pair
  PoleZeroPair (const std::complex<FP>& p1, const std::complex<FP>& z1,
                const std::complex<FP>& p2, const std::complex<FP>& z2)
    : poles (p1, p2)
    , zeros (z1, z2)
  {
  }

  bool isSinglePole () const
  {
    return poles.second == 0. && zeros.second == 0.;
  }

  bool isnan () const
  {
    return poles.isnan() || zeros.isnan();
  }
};

// Identifies the general class of filter
enum Kind
{
  kindLowPass,
  kindHighPass,
  kindBandPass,
  kindBandStop,
  kindLowShelf,
  kindHighShelf,
  kindBandShelf,
  kindOther
};

}

#endif
