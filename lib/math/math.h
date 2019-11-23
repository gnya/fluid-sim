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
  inline void brush_step2f(float *x, __m256 i, __m256 j, __m256 p_idx_x, __m256 p_idx_y,
                           __m256 r_inv, __m256 vdt) {
    using namespace util::avx;

    __m256 d_i = _mm256_sub_ps(i, p_idx_x);
    __m256 d_j = _mm256_sub_ps(j, p_idx_y);

    __m256 _n2 = _mm256_add_ps(_mm256_mul_ps(d_i, d_i), _mm256_mul_ps(d_j, d_j));
    __m256 _s = exp_m256(_mm256_mul_ps(_n2, r_inv));

    __m256 tmp = _mm256_loadu_ps(x);
    tmp = _mm256_add_ps(tmp, _mm256_mul_ps(vdt, _s));
    _mm256_storeu_ps(x, tmp);
  }

  void brush2f(float *x, float p_idx[2], float v[2], float r, float dt, int m, int n) {
    using namespace util;

    const __m256 _v = _mm256_set_ps(v[1], v[0], v[1], v[0], v[1], v[0], v[1], v[0]);
    const __m256 _p_idx_x = _mm256_set1_ps(p_idx[0]);
    const __m256 _p_idx_y = _mm256_set1_ps(p_idx[1]);
    const __m256 _r_inv = _mm256_set1_ps(-1 / r);
    const __m256 _dt = _mm256_set1_ps(dt);
    const __m256 _vdt = _mm256_mul_ps(_v, _dt);
    const __m256 _idx = _mm256_set_ps(3, 3, 2, 2, 1, 1, 0, 0);

    auto i_parallel_end = m - m % 4;

    for (int j = 0; j < n; ++j) {
      __m256 _j = _mm256_set1_ps((float) j);

      for (int i = 0; i < i_parallel_end; i += 4) {
        __m256 _i = _mm256_add_ps(_mm256_set1_ps((float) i), _idx);

        brush_step2f(&x[at2_x(m, i, j)], _i, _j, _p_idx_x, _p_idx_y, _r_inv, _vdt);
      }

      for (int i = i_parallel_end; i < m; ++i) {
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
