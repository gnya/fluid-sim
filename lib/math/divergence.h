/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 23
 */

#ifndef FLUID_MATH_DIVERGENCE_H
#define FLUID_MATH_DIVERGENCE_H

namespace fluid::math {
  inline void divergence_avx_step2f(float *x, float *x_div, __m256 dx_half, int m) {
    using namespace util::avx;

    __m256 x0 = load_m256_2f(x + 1 * 2 + 0); // x[i + 1, j].x
    __m256 x1 = load_m256_2f(x - 1 * 2 + 0); // x[i - 1, j].x
    __m256 y0 = load_m256_2f(x + m * 2 + 1); // x[i, j + 1].y
    __m256 y1 = load_m256_2f(x - m * 2 + 1); // x[i, j - 1].y
    __m256 tmp = _mm256_add_ps(_mm256_sub_ps(x0, x1), _mm256_sub_ps(y0, y1));

    _mm256_storeu_ps(x_div, _mm256_mul_ps(tmp, dx_half));
  }

  inline void divergence_step2f(float *x, float *x_div, float dx_half, int m,
                                int c_l = 1, int c_r = 1, int c_b = 1, int c_t = 1) {
    float x0 = *(x + 1 * 2 * c_r + 0); // x[i + 1, j].x
    float x1 = *(x - 1 * 2 * c_l + 0); // x[i - 1, j].x
    float y0 = *(x + m * 2 * c_t + 1); // x[i, j + 1].y
    float y1 = *(x - m * 2 * c_b + 1); // x[i, j - 1].y

    // x_div[i, j] = ...
    *x_div = ((x0 - x1) + (y0 - y1)) * dx_half;
  }

  void divergence2f(float *x, float *x_div, float dx, int m, int n) {
    const float dx_half = dx / 2;
    const __m256 _dx_half = _mm256_set1_ps(dx_half);
    auto i_parallel_end = (m - 1) - (m - 2) % 8;

    // inner of boundary
    for (int j = 1; j < n - 1; ++j) {
      for (int i = 1; i < i_parallel_end; i += 8) {
        int idx = util::at(m, i, j);
        divergence_avx_step2f(&x[idx * 2], &x_div[idx], _dx_half, m);
      }

      for (int i = i_parallel_end; i < m - 1; ++i) {
        int idx = util::at(m, i, j);
        divergence_step2f(&x[idx * 2], &x_div[idx], dx_half, m);
      }
    }

    for (int j = 1; j < n - 1; ++j) {
      int idx_l = util::at(m, 0, j); // left boundary
      int idx_r = util::at(m, m - 1, j); // right boundary

      divergence_step2f(&x[idx_l * 2], &x_div[idx_l], dx_half, m, 0, 1, 1, 1);
      divergence_step2f(&x[idx_r * 2], &x_div[idx_r], dx_half, m, 1, 0, 1, 1);
    }

    for (int i = 1; i < m - 1; ++i) {
      int idx_b = util::at(m, i, 0); // bottom boundary
      int idx_t = util::at(m, i, n - 1); // top boundary

      divergence_step2f(&x[idx_b * 2], &x_div[idx_b], dx_half, m, 1, 1, 0, 1);
      divergence_step2f(&x[idx_t * 2], &x_div[idx_t], dx_half, m, 1, 1, 1, 0);
    }
  }
}

#endif //FLUID_MATH_DIVERGENCE_H
