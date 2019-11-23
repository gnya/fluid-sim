/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 23
 */

#ifndef FLUID_MATH_H
#define FLUID_MATH_H

#include "util.h"
#include "advect.h"
#include "jacobi.h"
#include "divergence.h"
#include "gradient.h"
#include "noise.h"

namespace fluid::math {
  void brush2f(float *x, float p_idx[2], float v[2], float r, float dt, int m, int n) {
    using namespace util;

    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        auto n2 = (float) (pow(i - p_idx[0], 2) + pow(j - p_idx[1], 2));
        auto s = exp(-n2 / r);

        x[at2_x(m, i, j)] += v[0] * dt * s;
        x[at2_y(m, i, j)] += v[1] * dt * s;
      }
    }
  }

  void left_boundary2f(float *p, float *u, int m, int n) {
    using namespace util;

    for (int j = 1; j < n - 1; j++) {
      p[at(m, 0, j)] = p[at(m, 0 + 1, j)];
      u[at2_x(m, 0, j)] = -u[at2_x(m, 0 + 1, j)];
      u[at2_y(m, 0, j)] = -u[at2_y(m, 0 + 1, j)];
    }
  }

  void right_boundary2f(float *p, float *u, int m, int n) {
    using namespace util;

    for (int j = 1; j < n - 1; j++) {
      p[at(m, m - 1, j)] = p[at(m, m - 1 - 1, j)];
      u[at2_x(m, m - 1, j)] = -u[at2_x(m, m - 1 - 1, j)];
      u[at2_y(m, m - 1, j)] = -u[at2_y(m, m - 1 - 1, j)];
    }
  }

  void bottom_boundary2f(float *p, float *u, int m, int n) {
    using namespace util;

    for (int i = 1; i < m - 1; i++) {
      p[at(m, i, 0)] = p[at(m, i, 0 + 1)];
      u[at2_x(m, i, 0)] = -u[at2_x(m, i, 0 + 1)];
      u[at2_y(m, i, 0)] = -u[at2_y(m, i, 0 + 1)];
    }
  }

  void top_boundary2f(float *p, float *u, int m, int n) {
    using namespace util;

    for (int i = 1; i < m - 1; i++) {
      p[at(m, i, n - 1)] = p[at(m, i, n - 1 - 1)];
      u[at2_x(m, i, n - 1)] = -u[at2_x(m, i, n - 1 - 1)];
      u[at2_y(m, i, n - 1)] = -u[at2_y(m, i, n - 1 - 1)];
    }
  }
}

#endif //FLUID_MATH_H
