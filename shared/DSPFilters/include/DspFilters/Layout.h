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

#ifndef DSPFILTERS_LAYOUT_H
#define DSPFILTERS_LAYOUT_H

#include <type_traits>

#include "DspFilters/Common.h"
#include "DspFilters/MathSupplement.h"

namespace Dsp {

//
// Describes a filter as a collection of poles and zeros along with
// normalization information to achieve a specified gain at a specified
// frequency. The poles and zeros may lie either in the s or the z plane.
//

// Base uses pointers to reduce template instantiations
class Layout : public PoleZero<double>{
 public:
  // setter and getter for normalization gain and frequency
  void setNormal(const double& w, const double& g) {
    normalW_ = w;
    normalGain_ = g;
  }

  template <typename FP = double>
  FP getNormalW() const {
    return FP(normalW_);
  }

  template <typename FP = double>
  FP getNormalGain() const {
    return FP(normalGain_);
  }

  template <typename FP>
  std::vector<FP> getDenominator() const {
    return Layout::rootPairsToCoeffs<FP>(poles());
  }

  template <typename FP>
  std::vector<FP> getNumerator() const {
    return Layout::rootPairsToCoeffs<FP>(zeros());
  }

 protected:
  template <typename FP>
  static std::vector<FP> rootPairsToCoeffs(
      const std::vector<RootPair<double>>& rootPairs, const size_t& order = 0) {
    if (!rootPairs.size()) return std::vector<FP>{FP(1)};
    size_t coeffSize = order;
    // recalculate size to reserve for coefficients
    if (order < rootPairs.size()) {
      coeffSize = std::reduce(rootPairs.cbegin(), rootPairs.cend(), 0,
                              [](const size_t& x, const RootPair<double>& y) {
                                return x + y.order();
                              });
    }
    std::vector<double> res;
    res.reserve(coeffSize + 1);
    for (const auto& rp : rootPairs) {
      // convolute over each pairs of roots to get resultant polynomial
      // coefficients
      switch (rp.order()) {
        case 0:
          throw std::runtime_error("Unitialized root pair not allowed!");
        case 1:
          // roots contained in rp is of first order. Convolute with [1, -r]
          res = conv(res, {1., -rp.first.real()});
          continue;
        case 2:
          // roots contained in rp is of second order and is a conjugate pair.
          // Convolute with [1,-2*r.real(),r.real()*r.real()+r.imag()*r.imag()]
          res = conv(res, {1., -2. * rp.first.real(),
                           rp.first.real() * rp.first.real() +
                               rp.first.imag() * rp.first.imag()});
          continue;
        default:
          continue;
      }
    }
    // coefficients need to be casted to desired type.
    if (std::is_same<T, double>::value) return res;
    std::vector<FP> coeffs;
    coeffs.reserve(res.size());
    std::transform(res.cbegin(), res.cend(), coeffs.begin(),
                   [](const double& x) { return FP(x); });
    return coeffs;
  }

 private:
  double normalW_;
  double normalGain_;
};

}  // namespace Dsp

#endif
