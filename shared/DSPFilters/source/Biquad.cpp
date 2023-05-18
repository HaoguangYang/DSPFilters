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

#include "DspFilters/Common.h"
#include "DspFilters/MathSupplement.h"
#include "DspFilters/Biquad.h"

namespace Dsp {

template <typename FP>
BiquadPoleState<FP>::BiquadPoleState (const BiquadBase<FP>& s)
{
  const FP a0 = s.getA0 ();
  const FP a1 = s.getA1 ();
  const FP a2 = s.getA2 ();
  const FP b0 = s.getB0 ();
  const FP b1 = s.getB1 ();
  const FP b2 = s.getB2 ();

  if (a2 == 0 && b2 == 0)
  {
    // single pole
    poles.first = -a1;
    zeros.first = -b0 / b1;
    poles.second = 0;
    zeros.second = 0;
  }
  else
  {
    {
      const std::complex<FP> c = std::sqrt (std::complex<FP> (a1 * a1 - 4 * a0 * a2, 0));
      FP d = 1. / (2. * a0);
      poles.first = -(a1 + c) * d;
      poles.second =  (c - a1) * d;
      assert (!poles.isnan());
    }

    {
      const std::complex<FP> c = sqrt (std::complex<FP> (
        b1 * b1 - 4 * b0 * b2, 0));
      FP d = 1. / (2. * b0);
      zeros.first = -(b1 + c) * d;
      zeros.second =  (c - b1) * d;
      assert (!zeros.isnan());
    }
  }

  gain = b0 / a0;
}

//------------------------------------------------------------------------------
template <typename FP>
std::complex<FP> BiquadBase<FP>::response (const FP& normalizedFrequency) const
{
  const FP a0 = getA0 ();
  const FP a1 = getA1 ();
  const FP a2 = getA2 ();
  const FP b0 = getB0 ();
  const FP b1 = getB1 ();
  const FP b2 = getB2 ();

  const FP w = 2 * doublePi * normalizedFrequency;
  const std::complex<FP> czn1 = std::polar (1., -w);
  const std::complex<FP> czn2 = std::polar (1., -2 * w);
  std::complex<FP> ch (1);
  std::complex<FP> cbot (1);

  std::complex<FP> ct (b0/a0);
  std::complex<FP> cb (1);
  ct = addmul (ct, b1/a0, czn1);
  ct = addmul (ct, b2/a0, czn2);
  cb = addmul (cb, a1/a0, czn1);
  cb = addmul (cb, a2/a0, czn2);
  ch   *= ct;
  cbot *= cb;

  return ch / cbot;
}

template <typename FP>
std::vector<PoleZeroPair<FP>> BiquadBase<FP>::getPoleZeros () const
{
  std::vector<PoleZeroPair<FP>> vpz;
  BiquadPoleState<FP> bps (*this);
  vpz.emplace_back (bps);
  return vpz;
}

template <typename FP>
void BiquadBase<FP>::setCoefficients (const FP& a0, const FP& a1, const FP& a2,
                                  const FP& b0, const FP& b1, const FP& b2)
{
//#ifndef NDEBUG
//    assert (!Dsp::is_nan(a0));
//    assert (!Dsp::is_nan(b0));
//    assert (!Dsp::is_nan(a1));
//    assert (!Dsp::is_nan(a2));
//    assert (!Dsp::is_nan(b1));
//    assert (!Dsp::is_nan(b2));
//#endif
    m_a.resize(3, 1.);
    m_b.resize(3, 1.);
    
    if (!(Dsp::isnan (a0) || Dsp::isnan (a1) || Dsp::isnan (a2) ||
        Dsp::isnan (b0) || Dsp::isnan (b1) || Dsp::isnan (b2)))
    {
        m_a0 = a0;
        m_a[1] = a1/a0;
        m_a[2] = a2/a0;
        m_b[0] = b0/a0;
        m_b[1] = b1/a0;
        m_b[2] = b2/a0;
    }
}

template <typename FP>
void BiquadBase<FP>::setOnePole (const std::complex<FP>& pole, const std::complex<FP>& zero)
{
#if 0
  pole = adjust_imag (pole);
  zero = adjust_imag (zero);
#else
  assert (pole.imag() == 0); 
  assert (zero.imag() == 0);
#endif
  
  const double a0 = 1;
  const double a1 = -pole.real();
  const double a2 = 0;
  const double b0 = 1;
  const double b1 = -zero.real();
  const double b2 = 0;

  setCoefficients (a0, a1, a2, b0, b1, b2);
}

template <typename FP>
void BiquadBase<FP>::setTwoPole (const std::complex<FP>& pole1, const std::complex<FP>& zero1,
                             const std::complex<FP>& pole2, const std::complex<FP>& zero2)
{
#if 0
  pole1 = adjust_imag (pole1);
  pole2 = adjust_imag (pole2);
  zero1 = adjust_imag (zero1);
  zero2 = adjust_imag (zero2);
#endif

  const FP a0 = 1.;
  FP a1;
  FP a2;

  if (pole1.imag() != 0)
  {
    //assert (pole2 == std::conj (pole1));

    a1 = -2 * pole1.real();
    a2 = std::norm (pole1);
  }
  else
  {
    assert (pole2.imag() == 0);

    a1 = -(pole1.real() + pole2.real());
    a2 =   pole1.real() * pole2.real();
  }

  const double b0 = 1;
  double b1;
  double b2;

  if (zero1.imag() != 0)
  {
    //assert (zero2 == std::conj (zero1));

    b1 = -2 * zero1.real();
    b2 = std::norm (zero1);
  }
  else
  {
    assert (zero2.imag() == 0);

    b1 = -(zero1.real() + zero2.real());
    b2 =   zero1.real() * zero2.real();
  }

  setCoefficients (a0, a1, a2, b0, b1, b2);
}

}
