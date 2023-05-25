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
template <typename FP>
struct RootPair : complex_pair_t<FP> {
  RootPair() {}

  explicit RootPair(const std::complex<FP>& c1)
      : complex_pair_t<FP>(c1, FP(0)) {
    if (!isReal())
      throw std::runtime_error(
          "RootPair: unpaired constructor from std::complex -- imaginary part "
          "must be zero.");
    nRoots_ = 1;
  }

  RootPair(const std::complex<FP>& c1, const std::complex<FP>& c2)
      : complex_pair_t<FP>(c1, c2) {
    nRoots_ = 2;
    if (!isMatchedPair())
      throw std::runtime_error(
          "RootPair: paired constructor from std::complex -- roots are not "
          "paired.");
  }

  bool isConjugate() const { return second == std::conj(first); }

  bool isReal() const { return first.imag() == 0 && second.imag() == 0; }

  // Returns true if this is either a conjugate pair,
  // or a pair of reals where neither is zero.
  bool isMatchedPair() const {
    if (nRoots_ < 2)
      return false;
    else if (first.imag() != 0)
      return second == std::conj(first);
    else
      return second.imag() == 0 && second.real() == first.real();
  }

  bool isnan() const { return Dsp::isnan(first) || Dsp::isnan(second); }

  const uint8_t& order() const { return nRoots_; }

 private:
  uint8_t nRoots_ = 0;
};

// A pair of pole/zeros. This fits in a biquad (but is missing the gain)
template <typename FP>
class PoleZero {
 public:
  PoleZero() : poles_{}, zeros_{}, nPoles_(0), nZeros_(0) {}

  PoleZero(const size_t& np, const size_t& nz)
      : poles_{}, zeros_{}, nPoles_(0), nZeros_(0) {
    reserve(np, nz);
  }

  PoleZero(const std::vector<RootPair<FP>>& p,
           const std::vector<RootPair<FP>>& z)
      : poles_(p), zeros_(z) {
    auto npz = PoleZero::sumOrders(p, z);
    nPoles_ = npz.first;
    nZeros_ = npz.second;
  }

  PoleZero(const std::complex<FP>& p, const std::complex<FP>& z)
      : poles_{RootPair<FP>(p)},
        zeros_{RootPair<FP>(z)},
        nPoles_(1),
        nZeros_(1) {}

  PoleZero(const FP& p, const FP& z)
      : poles_{RootPair<FP>(std::complex<FP>(p, FP(0)))},
        zeros_{RootPair<FP>(std::complex<FP>(z, FP(0)))},
        nPoles_(1),
        nZeros_(1) {}

  PoleZero(const std::complex<FP>& p1, const std::complex<FP>& z1,
           const std::complex<FP>& p2, const std::complex<FP>& z2)
      : poles_{RootPair<FP>(p1, p2)},
        zeros_{RootPair<FP>(z1, z2)},
        nPoles_(2),
        nZeros_(2) {}

  PoleZero(const RootPair<FP>& p, const RootPair<FP>& z)
      : poles_{p}, zeros_{z}, nPoles_(p.order()), nZeros_(z.order()) {}

  template <typename FP2>
  void reserve(const PoleZero<FP2>& other) {
    poles_.reserve(other.poles().size());
    zeros_.reserve(other.zeros().size());
  }

  void reserve(const size_t& np, const size_t& nz) {
    poles_.reserve(np);
    zeros_.reserve(nz);
  }

  // modifiers
  // paired append
  template <typename T>
  void append(const T& p, const T& z) {
    append(p, false);
    append(z, true);
  }

  template <typename T>
  void appendConjugatePairs(const T& p, const T& z) {
    appendConjugatePairs(p, false);
    appendConjugatePairs(z, true);
  }

  // pole/zero pair
  void append(const std::complex<FP>& p1, const std::complex<FP>& z1,
              const std::complex<FP>& p2, const std::complex<FP>& z2) {
    poles_.emplace_back(p1, p2);
    zeros_.emplace_back(z1, z2);
    nPoles_ += 2;
    nZeros_ += 2;
  }

  // single (unpaired) append
  void append(const FP& x, constexpr bool isZero = false) {
    if (isZero) {
      zeros_.emplace_back(std::complex<FP>(x, FP(0)));
      nZeros_++;
    } else {
      poles_.emplace_back(std::complex<FP>(x, FP(0)));
      nPoles_++;
    }
  }

  void append(const std::complex<FP>& x, constexpr bool isZero = false) {
    if (isZero) {
      zeros_.emplace_back(x);
      nZeros_++;
    } else {
      poles_.emplace_back(x);
      nPoles_++;
    }
  }

  void append(const RootPair<FP>& x, constexpr bool isZero = false) {
    if (isZero) {
      zeros_.emplace_back(x);
      nZeros_++;
    } else {
      poles_.emplace_back(x);
      nPoles_++;
    }
  }

  void append(const std::vector<std::complex<FP>>& x,
              constexpr bool isZero = false) {
    auto xpairs = PoleZero::makeRootPair(x);
    if (isZero) {
      zeros_.insert(zeros_.end(), xpairs.begin(), xpairs.end());
      nZeros_ += x.size();
    } else {
      poles_.insert(poles_.end(), xpairs.begin(), xpairs.end());
      nPoles_ += x.size();
    }
  }

  void append(const std::vector<RootPair<FP>>& x,
              constexpr bool isZero = false) {
    size_t ordersToAdd = PoleZero::sumOrders(x);
    if (isZero) {
      zeros_.insert(zeros.end(), x.begin(), x.end());
      nZeros_ += ordersToAdd;
    } else {
      poles_.insert(poles.end(), x.begin(), x.end());
      nPoles_ += ordersToAdd;
    }
  }

  void appendConjugatePairs(const std::complex<FP>& x,
                            constexpr bool isZero = false) {
    if (isZero) {
      zeros_.emplace_back(x, std::conj(x));
      nZeros_ += 2;
    } else {
      poles_.emplace_back(x, std::conj(x));
      nPoles_ += 2;
    }
  }

  void appendConjugatePairs(const std::vector<std::complex<FP>>& x,
                            constexpr bool isZero = false) {
    if (isZero) {
      for (const auto& i : x){
        zeros_.emplace_back(i, std::conj(i));
        nZeros_ += 2;
      }
    } else {
      for (const auto& i : x){
        poles_.emplace_back(i, std::conj(i));
        nPoles_ += 2;
      }
    }
  }

  void sort(constexpr bool sPlane = true, constexpr bool descending = false) {
    PoleZero::sort(poles_, sPlane, descending);
    PoleZero::sort(zeros_, sPlane, descending);
  }

  void inverse() {
    std::swap(poles_, zeros_);
    std::swap(nPoles_, nZeros_);
  }

  void clear() {
    poles_.clear();
    zeros_.clear();
    nPoles_ = 0;
    nZeros_ = 0;
  }

  // properties
  bool isSinglePole(const size_t& index = 0) const {
    return poles_[index].order() == 1;
  }

  bool isSingleZero(const size_t& index = 0) const {
    return zeros_[index].order() == 1;
  }

  bool isnan() const {
    for (const auto& i : poles_) {
      if (i.isnan()) return true;
    }
    for (const auto& i : zeros_) {
      if (i.isnan()) return true;
    }
    return false;
  }

  // getters
  const size_t& getNumPoles() const { return nPoles_; }

  const size_t& getNumZeros() const { return nZeros_; }

  const std::vector<RootPair<FP>>& poles() const { return poles_; }

  const std::vector<RootPair<FP>>& zeros() const { return zeros_; }

  // allow indexing poles and zeros together
  std::pair<RootPair<FP>, RootPair<FP>> getPair(const size_t& pairIndex) const {
    if (pairIndex >= this->poles_.size() || pairIndex >= this->zeros_.size())
      throw std::runtime_error("pairIndex out of range");
    return std::make_pair(this->poles_[pairIndex], this->zeros_[pairIndex]);
  }

  std::pair<RootPair<FP>, RootPair<FP>> operator[](
      const size_t& pairIndex) const {
    return getPair(pairIndex);
  }

 protected:
  static std::pair<size_t, size_t> sumOrders(
      const std::vector<RootPair<FP>>& p, const std::vector<RootPair<FP>>& z) {
    return std::make_pair(PoleZero::sumOrders(p), PoleZero::sumOrders(z));
  }

  static size_t sumOrders(const std::vector<RootPair<FP>>& x) {
    return std::reduce(p.cbegin(), p.cend(), size_t(0),
                       [](const size_t& prev, const RootPair<FP>& p) -> size_t {
                         return prev + p.order();
                       });
  }

  static std::vector<RootPair<FP>> makeRootPair(
      const std::vector<std::complex<FP>>& v) {
    std::vector<RootPair<FP>> ret{};
    if (v.empty()) return ret;
    if (v.size() < 2) {
      ret.emplace_back(v.front());
      return ret;
    }
    std::unordered_map<std::pair<FP, FP>, unsigned int> hashMap;
    for (const auto& v1 : v) {
      auto mapItem = std::make_pair(std::real(v1), std::abs(std::imag(v1)));
      hashMap[mapItem] += 1;
      if (!hashMap[mapItem] % 2) ret.emplace_back(v1, std::conj(v1));
    }
    for (const auto& kvpair : hashMap) {
      // complex root with imaginary parts must come in pairs
      if (!(kvpair.second % 2)) continue;
      if (kvpair.first.second != 0)
        throw std::runtime_error("Complex roots must come in pairs.");
      else
        ret.emplace_back(std::complex<FP>(kvpair.first.first, FP(0)));
    }
  }

  static void sort(std::vector<RootPair<FP>>& roots,
                   constexpr bool sPlane = true,
                   constexpr bool descending = false) {
    auto comp = [](const RootPair<FP>& p1, const RootPair<FP>& p2) -> bool {
      return (p1.first.real() < p2.first.real());
    };
    if (descending && sPlane) {
      comp = [](const RootPair<FP>& p1, const RootPair<FP>& p2) -> bool {
        return (p1.first.real() > p2.first.real());
      };
    } else if (!sPlane) {
      if (!descending) {
        comp = [](const RootPair<FP>& p1, const RootPair<FP>& p2) -> bool {
          return (p1.first.norm() < p2.first.norm());
        };
      } else {
        comp = [](const RootPair<FP>& p1, const RootPair<FP>& p2) -> bool {
          return (p1.first.norm() > p2.first.norm());
        };
      }
    }
    std::sort(roots.begin(), roots.end(), comp);
  }

 private:
  std::vector<RootPair<FP>> poles_, zeros_;
  size_t nPoles_, nZeros_;
};

// Identifies the general class of filter
enum Kind {
  kindLowPass,
  kindHighPass,
  kindBandPass,
  kindBandStop,
  kindLowShelf,
  kindHighShelf,
  kindBandShelf,
  kindOther
};

}  // namespace Dsp

#endif
