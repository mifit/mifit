#ifndef mi_math_Random_h
#define mi_math_Random_h

#include <cmath>

namespace mi {
namespace math {

class Random {
  double nextNextGaussian;
  bool haveNextNextGaussian;

protected:


public:

  Random()
    : haveNextNextGaussian(false) {
  }

  virtual ~Random() {
  }

  void resetGaussian() {
    haveNextNextGaussian = false;
  }

  /**
   * Sets the seed of this random number generator.
   *
   * @param seed the seed.
   */
  virtual void setSeed(unsigned int seed) = 0;

  /**
   * Returns the next random, uniformly
   * distributed <code>int</code> value.
   */
  virtual int nextInt() = 0;

  /**
   * Returns the maximum value of a random int.
   */
  virtual int getIntMaximum() = 0;

  /**
   * Returns the minimum value of a random int.
   */
  virtual int getIntMinimum() = 0;

  /**
   * Returns the next random, uniformly
   * distributed <code>boolean</code> value.
   */
  virtual bool nextBoolean() = 0;

  /**
   * Returns the next random, uniformly distributed
   * <code>float</code> value between <code>0.0</code>
   * and <code>1.0</code>. Whether the range is
   * inclusive or exclusive will depend upon the
   * implementation.
   */
  virtual float nextFloat() = 0;

  /**
   * Returns the next random, uniformly
   * distributed <code>float</code> value between
   * <code>0.0</code> and <code>1.0</code> exclusive.
   */
  virtual float nextFloatExclusive() = 0;

  /**
   * Returns the next random, uniformly
   * distributed <code>float</code> value between
   * <code>0.0</code> and <code>1.0</code> inclusive.
   */
  virtual float nextFloatInclusive() = 0;

  /**
   * Returns the next random, uniformly distributed
   * <code>double</code> value between <code>0.0</code>
   * and <code>1.0</code>. Whether the range is
   * inclusive or exclusive will depend upon the
   * implementation.
   */
  virtual double nextDouble() = 0;

  /**
   * Returns the next random, uniformly distributed
   * <code>double</code> value between <code>0.0</code>
   * and <code>1.0</code> exclusive.
   */
  virtual double nextDoubleExclusive() = 0;

  /**
   * Returns the next random, uniformly distributed
   * <code>double</code> value between <code>0.0</code>
   * and <code>1.0</code> inclusive.
   */
  virtual double nextDoubleInclusive() = 0;

  /**
   * Returns the next random, Gaussian
   * ("normally") distributed <code>double</code> value
   * with mean <code>0.0</code> and standard deviation
   * <code>1.0</code>.
   */
  virtual double nextGaussian() {
    if (haveNextNextGaussian) {
      haveNextNextGaussian = false;
      return nextNextGaussian;
    } else {
      double v1;
      double v2;
      double s;
      do {
        v1 = 2.0 * nextDouble() - 1.0;
        v2 = 2.0 * nextDouble() - 1.0;
        s = v1*v1 + v2*v2;
      } while (s >= 1.0 || s == 0.0);
      double multiplier = std::sqrt(-2.0 * std::log(s)/s);
      nextNextGaussian = v1 * multiplier;
      haveNextNextGaussian = true;
      return v2 * multiplier;
    }
  }

  /**
   * Returns the next random, Cauchy distributed
   * <code>double</code> value between <code>0.0</code>
   * <code>1.0</code>.
   */
  virtual double nextCauchy() {
    double v1;
    double v2;
    do {
      v1 = 2.0 * nextDouble() - 1.0;
      v2 = 2.0 * nextDouble() - 1.0;
    } while (v1*v1 + v2*v2 > 1.0 || v1 == 0.0);
    return v2/v1;
  }

};
}
}


#endif
