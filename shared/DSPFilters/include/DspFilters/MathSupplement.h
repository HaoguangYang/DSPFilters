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

#ifndef DSPFILTERS_MATHSUPPLEMENT_H
#define DSPFILTERS_MATHSUPPLEMENT_H

#include "DspFilters/Common.h"

namespace Dsp {

const double doublePi   = M_PIf64;
const double doublePi_2	= M_PI_2f64;
const double doubleLn2  = std::log(2.);
const double doubleLn10	= std::log(10.);

template <typename Real>
using complex_pair_t = std::pair<std::complex<Real>, std::complex<Real>>;

template<typename Real>
inline std::complex<Real> solve_quadratic_1 (const Real& a, const Real& b, const Real& c)
{
  return (-b + std::sqrt (std::complex<Real> (b*b - 4*a*c, 0))) / (2. * a);
}

template<typename Real>
inline std::complex<Real> solve_quadratic_2 (const Real& a, const Real& b, const Real& c)
{
  return (-b - std::sqrt (std::complex<Real> (b*b - 4*a*c, 0))) / (2. * a);
}

template<typename Real>
inline const std::complex<Real> infinity()
{
  return std::complex<Real>(std::numeric_limits<Real>::infinity());
}

template<typename Real>
inline const std::complex<Real> adjust_imag (const std::complex<Real>& c)
{
  if (std::fabs (c.imag()) < 1e-30)
    return std::complex<Real>(c.real(), 0.);
  else
    return c;
}

template <typename Ty, typename To>
inline std::complex<Ty> addmul (const std::complex<Ty>& c,
                                const Ty& v,
                                const std::complex<To>& c1)
{
  return std::complex <Ty> (
    c.real() + v * c1.real(), c.imag() + v * c1.imag());
}

template <typename Ty>
inline std::complex<Ty> recip (const std::complex<Ty>& c)
{
  Ty n = 1.0 / std::norm (c);
  
  return std::complex<Ty> (n * c.real(), n * c.imag());
}

template <typename Ty>
inline Ty asinh (const Ty& x)
{
  return std::asinh (x);
}

template <typename Ty>
inline Ty acosh (const Ty& x)
{
  return std::acosh (x);
}

template <typename Ty>
inline bool isnan (const Ty& v)
{
  return !(v == v);
}

template <>
inline bool isnan (const float& v) { return std::isnan(v); }

template <>
inline bool isnan (const double& v) { return std::isnan(v); }

template <>
inline bool isnan (const long double& v) { return std::isnan(v); }

template <typename FP>
inline bool isnan (const std::complex<FP>& v)
{
  return Dsp::isnan (v.real()) || Dsp::isnan (v.imag());
}

//------------------------------------------------------------------------------

/*
 * Hack to prevent denormals
 *
 */

//const double anti_denormal_vsa = 1e-16; // doesn't prevent denormals
//const double anti_denormal_vsa = 0;
const double anti_denormal_vsa = 1e-8;

class DenormalPrevention
{
public:
  DenormalPrevention ()
    : m_v (anti_denormal_vsa)
  {
  }

  // small alternating current
  inline double ac ()
  {
    return m_v = -m_v;
  }

  // small direct current
  static inline double dc ()
  {
    return anti_denormal_vsa;
  }

private:
  double m_v;
};

}

#endif
