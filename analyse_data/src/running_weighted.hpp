#ifndef STATISTICS_RUNNING_WEIGHTED_HPP
#define STATISTICS_RUNNING_WEIGHTED_HPP

#include <cmath>

namespace statistics {

template <typename T = double>
class running_weighted
/// @brief calculates statistics on the fly
///
/// see also:
/// [1] Pebay et al. Comput Stat 31 (2016) 1305.
/// DOI: 10.1007/s00180-015-0637-z
{
  private:
  double n;
  T X, M, W, W2;  // helper variables
  // X <-> weighted mean
  // M <-> sum over squared deviation to the mean
  // W <-> sum over weights
  // W2<-> sum over squared weights
	constexpr T squared(T x) { return x*x; }

  public:
  /// default constructor
  running_weighted() : n(0.), X(0.), M(0.), W(0.), W2(0.) {}
  running_weighted(const running_weighted<T>& other)
      : n(other.n), X(other.X), M(other.M), W(other.W), W2(other.W2){};

  /// add a new data point via operator()
  inline void operator()(const T& x, const T w = 1)
  {
    // running average;
    n += 1;
    const T dX = x - X;
    const T wA = W;
    W  += w;  // new total weight
    X  += w / W * dX;
    M  += wA * squared(w / W * dX) + w * squared(wA / W * dX);
    W2 += squared(w);
  }

  inline void operator()(const running_weighted<T>& other)
  {
    if (W > 0.0) {
      const T XA = X;
      const T XB = other.X;
      const T MA = M;
      const T MB = other.M;
      const T wA = W;
      const T wB = other.W;
      const T dX = XB - XA;
      W  = wA + wB;  // new total weight
      X  = XA + wB / W * dX;
      M  = MA + MB + wA * squared(wB / W * dX) + wB * squared(wA / W * dX);
      n  += other.n;
      W2 += other.W2;
    }
    else {
      X  = other.X;
      M  = other.M;
      W  = other.W;
      n  = other.n;
      W2 = other.W2;
    }
  }

  /// clears all information ( if you want to reuse the same instance
  /// of this class for another dataset )
  void clear()
  {
    n = 0.;
    X = M = W = W2 = 0.;
  }

  /// current average
  inline T mean() const { return X; }
  inline T weight() const { return W; }
  inline int nr_samples() const { return static_cast<int>(n + 0.5); }
  /// unbiased weighted variance;
  /// with
  /// \f$ W := \sum_{i=0}^n w_i \f$ and \f$ W_2 := \sum_{i=0}^n (w_i)^2 \f$
  /// it is defined as
  /// \f$ \sigma^2 =
  ///     \frac{1}{W-W_2/W}\sum_{i=0}^n w_i (x_i - \left<x\right>)^2 \f$
  inline T variance() const
  {
    if (W != 0.) {
      return M * 1.0 / (W - W2/W);
    }
    else {
      return 0.;
    }
  }

  /// current variance
  /// \f$ \sigma^2_{n} =
  ///     \frac{1}{n}\sum_{i=0}^n(x_i - \left<x\right>)^2 \f$
  inline T variance_n() const
  {
    if (W == 0) {
      return 0.;
    }
    else {
      return M / W;
    }
  }

  /// current standard deviation
  /// \f$ \epsilon = \sqrt{\sigma^2_{n-1}} =
  ///    \sqrt{\frac{1}{n-1}\sum_{i=0}^n(x_i - \left<x\right>)^2} \f$
  inline T standard_deviation() const { return sqrt(variance()); }
  inline T standard_deviation_of_mean() const { return sqrt(variance() / W); }
};
}  // namespace statistics

#endif  // STATISTICS_RUNNING_WEIGHTED_HPP
