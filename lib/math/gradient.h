/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 23
 */

#ifndef FLUID_MATH_GRADIENT_H
#define FLUID_MATH_GRADIENT_H

namespace fluid::math {
  inline void gradient_avx_step2f(const float *p, const float *w, float *u, __m128 dx_half, int m) {
    using namespace util::avx;

    __m256 w_term = _mm256_loadu_ps(w); // w[i, j]

    __m128 p_x0 = _mm_loadu_ps(p + 1); // p[i + 1, j]
    __m128 p_x1 = _mm_loadu_ps(p - 1); // p[i - 1, j]
    __m128 p_y0 = _mm_loadu_ps(p + m); // p[i, j + 1]
    __m128 p_y1 = _mm_loadu_ps(p - m); // p[i, j - 1]
    p_x0 = _mm_mul_ps(_mm_sub_ps(p_x0, p_x1), dx_half);
    p_y0 = _mm_mul_ps(_mm_sub_ps(p_y0, p_y1), dx_half);
    __m256 p_term = marge_m256_2f(p_x0, p_y0);

    _mm256_storeu_ps(u, _mm256_sub_ps(w_term, p_term));
  }

  inline void gradient_step2f(const float *p, const float *w, float *u, float dx_half, int m,
                              int c_l = 1, int c_r = 1, int c_b = 1, int c_t = 1) {
    float w_x = *(w + 0); // w[i, j].x
    float w_y = *(w + 1); // w[i, j].y

    float p_x0 = *(p + 1 * c_r); // p[i + 1, j]
    float p_x1 = *(p - 1 * c_l); // p[i - 1, j]
    float p_y0 = *(p + m * c_t); // p[i, j + 1]
    float p_y1 = *(p - m * c_b); // p[i, j - 1]

    // u[i, j].x = ...
    // u[i, j].y = ...
    *(u + 0) = w_x - (p_x0 - p_x1) * dx_half;
    *(u + 1) = w_y - (p_y0 - p_y1) * dx_half;
  }

  void gradient2f(const float *p, const float *w, float *u, float dx, int m, int n) {
    const float dx_half = dx / 2;
    const __m128 _dx_half = _mm_set1_ps(dx_half);
    auto i_parallel_end = (m - 1) - (m - 2) % 4;

    // inner of boundary
    for (int j = 1; j < n - 1; ++j) {
      for (int i = 1; i < i_parallel_end; i += 4) {
        int idx = util::at(m, i, j);
        gradient_avx_step2f(&p[idx], &w[idx * 2], &u[idx * 2], _dx_half, m);
      }

      for (int i = i_parallel_end; i < m - 1; ++i) {
        int idx = util::at(m, i, j);
        gradient_step2f(&p[idx], &w[idx * 2], &u[idx * 2], dx_half, m);
      }
    }

    for (int j = 1; j < n - 1; ++j) {
      int idx_l = util::at(m, 0, j); // left boundary
      int idx_r = util::at(m, m - 1, j); // right boundary

      gradient_step2f(&p[idx_l], &w[idx_l * 2], &u[idx_l * 2], dx_half, m, 0, 1, 1, 1);
      gradient_step2f(&p[idx_r], &w[idx_r * 2], &u[idx_r * 2], dx_half, m, 1, 0, 1, 1);
    }

    for (int i = 1; i < m - 1; ++i) {
      int idx_b = util::at(m, i, 0); // bottom boundary
      int idx_t = util::at(m, i, n - 1); // top boundary

      gradient_step2f(&p[idx_b], &w[idx_b * 2], &u[idx_b * 2], dx_half, m, 1, 1, 0, 1);
      gradient_step2f(&p[idx_t], &w[idx_t * 2], &u[idx_t * 2], dx_half, m, 1, 1, 1, 0);
    }
  }
}

#endif //FLUID_MATH_GRADIENT_H
