/*
 *   This is an implementation of RGBM algorithm for Gene Regulatory Network
 *   inference from any type of expression data built on top of ennet package, in form of an R package.
 *   Copyright (C) 2016  Raghvendra Mall
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program, see LICENSE.
 */

#include <iostream>
#include <cstdlib>
#include "regression_stump.h"
#include "pcg_basic.h"

#define N 100
#define P 5
#define lf 1

using namespace std;

int main() {

  // In this version of the code, we'll use the global rng, rather than a
  // local one.
  
  // You should *always* seed the RNG.  The usual time to do it is the
  // point in time when you create RNG (typically at the beginning of the
  // program).
  //
  // pcg32_srandom_r takes two 64-bit constants (the initial state, and the
  // rng sequence selector; rngs with different sequence selectors will
  // *never* have random sequences that coincide, at all) - the code below
  // shows three possible ways to do so.
  bool seed_given = false;
  uint64_t seed = 123;
  if (seed_given) {
    pcg32_srandom(42u, seed);
  } else {
    // Seed with a fixed constant
    pcg32_srandom(42u, 54u);
  }
  
	// dummy test
	double x[N * P];
	for (int i = 0; i < N * P; i++) {
		//x[i] = rand() % 1000;
		x[i] = (int)pcg32_boundedrand(1000 - .00000000001);
	}
	double y[N];
	for (int i = 0; i < N; i++) {
		//y[i] = rand() % 1000;
	  y[i] = (int)pcg32_boundedrand(1000 - .00000000001);
	}

	train_regression_stump(N, P, x, y, 0.5, 0.5, lf, 1, 0.01);
//	cout << "done" << endl;
}
