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

#ifndef DSPFILTERS_STATE_H
#define DSPFILTERS_STATE_H

#ifdef __SSE3__
#include <pmmintrin.h>
#if (__GNUC__ == 4 && __GNUC_MINOR__ < 8)
#ifndef __m128_buggy_gxx_up_to_4_7_type
#define __m128_buggy_gxx_up_to_4_7_type
union __m128_buggy_gxx_up_to_4_7 {
    __m128 v;
    float e[4];
};
#endif
#endif
#elif defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

#include "DspFilters/Common.h"
#include "DspFilters/Biquad.h"

#include <stdexcept>
#include <iterator>

#include <deque>
#include <initializer_list>
#include <vector>
#include <numeric>
#include <memory>

namespace Dsp {

/*
 * Various forms of state information required to
 * process channels of actual sample data.
 *
 */

/**
 * @brief A ring buffer data structure to store state data of each channel
 * 
 * @tparam T data type
 */
template <typename T>
class RingBuffer {
 public:
  RingBuffer(const size_t& capacity) : buffer_(capacity), capacity_(std::min(capacity, buffer_.max_size())) {
    buffer_.clear();
  }

  void enqueue(const T& item) {
    if (isFull()) {
      buffer_.pop_back(); // Remove the tail item to make space
    }
    buffer_.push_front(item);
  }

  void enqueue(const std::initializer_list<T>& items) {
    for (const auto& item : items) {
      enqueue(item);
    }
  }

  T dequeue() {
    if (isEmpty()) {
      return T();
    }

    T item = buffer_.back();
    buffer_.pop_back();
    return item;
  }

  std::vector<T> dequeue(const size_t& N) {
    std::vector<T> lastElements;
    for (size_t n = 0; n < N && !isEmpty(); n++) {
      lastElements.push_back(buffer_.back());
      buffer_.pop_back();
    }
    return lastElements;
  }

  void clear() { buffer_.clear(); }

  void resize(const size_t& newCapacity) {
    capacity_ = newCapacity;
    if (buffer_.size() > newCapacity)
      buffer_.resize(newCapacity);
  }

  bool isEmpty() const {
    return buffer_.empty();
  }

  bool isFull() const {
    return buffer_.size() == capacity_;
  }

  size_t getSize() const {
    return buffer_.size();
  }

  std::vector<T> getAllElements() const {
    return std::vector<T>(buffer_.begin(), buffer_.end());
  }

  std::deque<T> getRawBuffer() const { return buffer_; }

  T& operator[](size_t index) {
    return buffer_[index];
  }

  const T& operator[](size_t index) const {
    return buffer_[index];
  }

  T* begin() {
    return buffer_.begin();
  }

  T* end() {
    return buffer_.end();
  }

 private:
  std::deque<T> buffer_;
  size_t capacity_;
};

//------------------------------------------------------------------------------

/*
 * State for applying a second order section to a sample using Direct Form I
 *
 * Difference equation:
 *
 *  y[n] = (b0/a0)*x[n] + (b1/a0)*x[n-1] + (b2/a0)*x[n-2]
 *                      - (a1/a0)*y[n-1] - (a2/a0)*y[n-2]  
 */
template <class FP>
class DirectFormI
{
public:
  typedef FP FPType;

  DirectFormI(const BiquadBase<FP>& s) : m_coeffs_(s) {
    size_t N1 = s.m_b.size();
    size_t M1 = s.m_a.size();
    assert(M1 > 0 && s.m_a[0] == 1. && "ERROR: First term of the denominator should be normalized to 1.");
    m_x_ = std::make_unique<RingBuffer<FP>>(N1);
    m_y_ = std::make_unique<RingBuffer<FP>>(M1-1);
    reset();
  }

  DirectFormI ()
  {
    m_x_ = std::make_unique<RingBuffer<FP>>(3);
    m_y_ = std::make_unique<RingBuffer<FP>>(2);
    reset();
  }

  void reset ()
  {
    m_x_->clear();
    m_y_->clear();
  }

  inline FP process1 (const FP& in,
                      const BiquadBase<FP>& s,
                      const FP& vsa) // very small amount
  {
    m_x_->enqueue(in);
    FP out = s.m_b[0]*m_x_[0] + s.m_b[1]*m_x_[1] + s.m_b[2]*m_x_[2]
                       - s.m_a[1]*m_y_[0] - s.m_a[2]*m_y_[1];
             + vsa;
    m_y_->enqueue(out);
    return out;
  }

  inline FP process1 (const FP& in, const FP& vsa) // very small amount
  {
    m_x_->enqueue(in);
    FP out = std::transform_reduce(m_coeffs_.m_b.begin(), m_coeffs_.m_b.end(), m_x_->begin(), FP(0),
                                 std::plus<FP>(), std::multiplies<FP>());
    out -= std::transform_reduce(m_coeffs_.m_a.begin()+1, m_coeffs_.m_a.end(), m_y_->begin(), -vsa,
                                 std::plus<FP>(), std::multiplies<FP>());
    m_y_->enqueue(out);
    return out;
  }

protected:
  std::unique_ptr<RingBuffer<FP>> m_x_; // Ring buffer for x[n]...x[n-N]
  std::unique_ptr<RingBuffer<FP>> m_y_; // Ring buffer for y[n-1]...y[n-M]

  BiquadBase<FP> m_coeffs_;
};

//------------------------------------------------------------------------------

/*
 * State for applying a second order section to a sample using Direct Form II
 *
 * Difference equation:
 *
 *  v[n] =         x[n] - (a1/a0)*v[n-1] - (a2/a0)*v[n-2]
 *  y(n) = (b0/a0)*v[n] + (b1/a0)*v[n-1] + (b2/a0)*v[n-2]
 *
 */
template <class FP>
class DirectFormII
{
public:
  typedef FP FPType;
  static constexpr bool HasSimd = false;

  DirectFormII ()
  {
    reset ();
  }

  void reset ()
  {
    m_v1 = FP();
    m_v2 = FP();
  }

  FP process1 (const FP& in,
               const BiquadBase<FP>& s,
               const FP& vsa)
  {
    FP w   = in - s.m_a[1]*m_v1 - s.m_a[2]*m_v2 + vsa;
    FP out =      s.m_b[0]*w    + s.m_b[1]*m_v1 + s.m_b[2]*m_v2;

    m_v2 = m_v1;
    m_v1 = w;

    return out;
  }

private:
  FP m_v1; // v[-1]
  FP m_v2; // v[-2]
};

//------------------------------------------------------------------------------


/*
 * Transposed Direct Form I and II
 * by lubomir i. ivanov (neolit123 [at] gmail)
 *
 * Reference:
 * http://www.kvraudio.com/forum/viewtopic.php?p=4430351
 *
 */

// I think this one is broken
template <class FP>
class TransposedDirectFormI
{
public:
  typedef FP FPType;
  static constexpr bool HasSimd = false;

  TransposedDirectFormI ()
  {
    reset ();
  }

  void reset ()
  {
    m_v = FP();
    m_s1 = FP();
    m_s1_1 = FP();
    m_s2 = FP();
    m_s2_1 = FP();
    m_s3 = FP();
    m_s3_1 = FP();
    m_s4 = FP();
    m_s4_1 = FP();
  }

  template <typename Sample, typename FP>
  inline Sample process1 (const Sample in,
                          const BiquadBase<FP>& s,
                          const double vsa)
  {
    FP out;

    // can be: in += m_s1_1;
    m_v = in + m_s1_1;
    out = s.m_b[0]*m_v + m_s3_1 + vsa;
    m_s1 = m_s2_1 - s.m_a[1]*m_v;
    m_s2 = -s.m_a[2]*m_v;
    m_s3 = s.m_b[1]*m_v + m_s4_1;
    m_s4 = s.m_b[2]*m_v; 

    m_s4_1 = m_s4;
    m_s3_1 = m_s3;
    m_s2_1 = m_s2;
    m_s1_1 = m_s1;

    return out;
  }

private:
  FP m_v;
  FP m_s1;
  FP m_s1_1;
  FP m_s2;
  FP m_s2_1;
  FP m_s3;
  FP m_s3_1;
  FP m_s4;
  FP m_s4_1;
};

//------------------------------------------------------------------------------
template <class FP>
class TransposedDirectFormII
{
public:
  typedef FP FPType;
  static constexpr bool HasSimd = false;

  TransposedDirectFormII ()
  {
    reset ();
  }

  void reset ()
  {
    m_s1 = FP();
    m_s1_1 = FP();
    m_s2 = FP();
    m_s2_1 = FP();
  }

  template <typename Sample, typename FP>
  inline Sample process1 (const Sample in,
                          const BiquadBase<FP>& s,
                          const FP vsa)
  {
    FP out;

    out = m_s1_1 + s.m_b[0]*in + vsa;
    m_s1 = m_s2_1 + s.m_b[1]*in - s.m_a[1]*out;
    m_s2 = s.m_b[2]*in - s.m_a[2]*out;
    m_s1_1 = m_s1;
    m_s2_1 = m_s2;

    return static_cast<Sample> (out);
  }

private:
  FP m_s1;
  FP m_s1_1;
  FP m_s2;
  FP m_s2_1;
};

//------------------------------------------------------------------------------

// Holds an array of states suitable for multi-channel processing
template <int Channels, class StateType>
class ChannelsState
{
public:
  ChannelsState ()
  {
  }

  int getNumChannels() const
  {
    return Channels;
  }

  void reset ()
  {
    for (int i = 0; i < Channels; ++i)
      m_state[i].reset();
  }

  StateType& operator[] (int index)
  {
    assert (index >= 0 && index < Channels);
    return m_state[index];
  }

  template <class Filter, typename Sample>
  void process (int numSamples,
                Sample* const* arrayOfChannels,
                Filter& filter)
  {
    for (int i = 0; i < Channels; ++i)
      filter.process (numSamples, arrayOfChannels[i], m_state[i]);
  }

  template <class Filter, typename... Iterators>
  void process (Filter& filter,
                Iterators... its
                )
  {
    static_assert(sizeof...(its) == 2 * Channels, "The number of iterator pairs must match the number of channels.");
    process_impl(filter, its...);
  }

private:
  template <class Filter, typename Iterator, typename... Iterators>
  void process_impl (Filter& filter,
                Iterator firstA, Iterator lastA,
                Iterators... its
                )
  {
    static_assert(sizeof...(its) % 2 == 0, "Must pass iterators in pairs");
    filter.process (firstA, lastA, m_state[Channels - sizeof...(its) / 2 - 1]);
    process_impl(filter, its...);
  }

  template <class Filter>
  void process_impl (Filter& filter) {}

  StateType m_state[Channels];
};

// Empty state, can't process anything
template <class StateType>
class ChannelsState <0, StateType>
{
public:
  int getNumChannels() const
  {
    return 0;
  }

  void reset ()
  {
    throw std::logic_error ("attempt to reset empty ChannelState");
  }

  template <class FilterDesign, typename Sample>
  void process (int /*numSamples*/,
                Sample* const* /*arrayOfChannels*/,
                FilterDesign& /*filter*/)
  {
    throw std::logic_error ("attempt to process empty ChannelState");
  }

  template <class Filter, typename It>
  void process (It first,
                Filter& filter)
  {
    throw std::logic_error ("attempt to process empty ChannelState");
  }
};

//------------------------------------------------------------------------------

}

#endif
