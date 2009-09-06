#include "fft.h"

int MIMapFactor(int ntest, int prime, int even, int inc) {
  int n;
  int j;
  int times = 0;
  /*
   *  find out if number is factorizable by numbers <= prime.
   *  return next larger number that is.
   *  even is 1 if numbers can be odd or 2 if must be even;
   */
  if (ntest == 0) {
    return (even);
  }
  if (even == EVEN && inc >= 0 && inc < even) {
    inc = even;
  }
  if (even == EVEN && inc < 0 && inc > even) {
    inc = -even;
  }
  if (even == EVEN && inc > 0) {
    ntest = ((ntest+1)/even)*even;
  }
  if (even == EVEN && inc < 0) {
    ntest = ((ntest-1)/even)*even;
  }
  n = ntest;
  while (times < 100) {
    for (j = 2; j <= prime; j++) {
      while (n == (n/j)*j) {
        n = n/j;
      }
    }
    if (n == 1) {
      break;
    }
    ntest += inc;
    n = ntest;
    times++;
  }
  return ntest;
}

