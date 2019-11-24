/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 23
 */

#ifndef FLUID_MATH_JACOBI_H
#define FLUID_MATH_JACOBI_H

namespace fluid::math {
  inline void jacobi_avx_step(const float *x, const float *b, float *x_new,
                              __m256 a, __m256 b_inv, int m, int dim) {
    __m256 tmp;

    tmp = _mm256_mul_ps(a  , _mm256_loadu_ps(b)); // b[i, j]
    tmp = _mm256_add_ps(tmp, _mm256_loadu_ps(x - dim * m)); // x[i, j - 1]
    tmp = _mm256_add_ps(tmp, _mm256_loadu_ps(x + dim * m)); // x[i, j + 1]
    tmp = _mm256_add_ps(tmp, _mm256_loadu_ps(x - dim * 1)); // x[i - 1, j]
    tmp = _mm256_add_ps(tmp, _mm256_loadu_ps(x + dim * 1)); // x[i + 1, j]
    tmp = _mm256_mul_ps(tmp, b_inv);

    _mm256_storeu_ps(x_new, tmp);
  }

  inline void jacobi_step1f(const float *x, const float *b, float *x_new,
                            float a, float b_inv, int m,
                            int c_l = 1, int c_r = 1, int c_b = 1, int c_t = 1) {
    float x00 = *(x - m * c_b); // x[i, j - 1]
    float x01 = *(x + m * c_t); // x[i, j + 1]
    float x10 = *(x - 1 * c_l); // x[i - 1, j]
    float x11 = *(x + 1 * c_r); // x[i + 1, j]
    float b_a = *b * a; // b[i, j] * a

    // x_new[i, j] = ...
    *x_new = (x00 + x01 + x10 + x11 + b_a) * b_inv;
  }

  void jacobi1f(const float *x, const float *b, float *x_new,
                float a, float b_inv, int m, int n) {
    const __m256 _a = _mm256_set1_ps(a);
    const __m256 _b_inv = _mm256_set1_ps(b_inv);
    auto i_parallel_end = (m - 1) - (m - 2) % 8;

    // inner of boundary
    for (int j = 1; j < n - 1; ++j) {
      for (int i = 1; i < i_parallel_end; i += 8) {
        int idx = util::at(m, i, j);
        jacobi_avx_step(&x[idx], &b[idx], &x_new[idx], _a, _b_inv, m, 1);
      }

      for (int i = i_parallel_end; i < m - 1; ++i) {
        int idx = util::at(m, i, j);
        jacobi_step1f(&x[idx], &b[idx], &x_new[idx], a, b_inv, m);
      }
    }

    for (int j = 1; j < n - 1; ++j) {
      int idx_l = util::at(m, 0, j); // left boundary
      int idx_r = util::at(m, m - 1, j); // right boundary

      jacobi_step1f(&x[idx_l], &b[idx_l], &x_new[idx_l], a, b_inv, m, 0, 1, 1, 1);
      jacobi_step1f(&x[idx_r], &b[idx_r], &x_new[idx_r], a, b_inv, m, 1, 0, 1, 1);
    }

    for (int i = 1; i < m - 1; ++i) {
      int idx_b = util::at(m, i, 0); // bottom boundary
      int idx_t = util::at(m, i, n - 1); // top boundary

      jacobi_step1f(&x[idx_b], &b[idx_b], &x_new[idx_b], a, b_inv, m, 1, 1, 0, 1);
      jacobi_step1f(&x[idx_t], &b[idx_t], &x_new[idx_t], a, b_inv, m, 1, 1, 1, 0);
    }
  }

  inline void jacobi_step2f(const float *x, const float *b, float *x_new,
                            float a, float b_inv, int m,
                            int c_l = 1, int c_r = 1, int c_b = 1, int c_t = 1) {
    float x00_x = *(x - 2 * m * c_b + 0); // x[i, j - 1].x
    float x00_y = *(x - 2 * m * c_b + 1); // x[i, j - 1].y
    float x01_x = *(x + 2 * m * c_t + 0); // x[i, j + 1].x
    float x01_y = *(x + 2 * m * c_t + 1); // x[i, j + 1].y
    float x10_x = *(x - 2 * 1 * c_l + 0); // x[i - 1, j].x
    float x10_y = *(x - 2 * 1 * c_l + 1); // x[i - 1, j].y
    float x11_x = *(x + 2 * 1 * c_r + 0); // x[i + 1, j].x
    float x11_y = *(x + 2 * 1 * c_r + 1); // x[i + 1, j].y
    float b_x_a = *(b + 0) * a; // b[i, j].x * a
    float b_y_a = *(b + 1) * a; // b[i, j].y * a

    // x_new[i, j].x = ...
    // x_new[i, j].y = ...
    *(x_new + 0) = (x00_x + x01_x + x10_x + x11_x + b_x_a) * b_inv;
    *(x_new + 1) = (x00_y + x01_y + x10_y + x11_y + b_y_a) * b_inv;
  }

  void jacobi2f(const float *x, const float *b, float *x_new,
                float a, float b_inv, int m, int n) {
    const __m256 _a = _mm256_set1_ps(a);
    const __m256 _b_inv = _mm256_set1_ps(b_inv);
    auto i_parallel_end = (m - 1) - (m - 2) % 4;

    // inner of boundary
    for (int j = 1; j < n - 1; ++j) {
      for (int i = 1; i < i_parallel_end; i += 4) {
        int idx = util::at2_x(m, i, j);
        jacobi_avx_step(&x[idx], &b[idx], &x_new[idx], _a, _b_inv, m, 2);
      }

      for (int i = i_parallel_end; i < m - 1; ++i) {
        int idx = util::at2_x(m, i, j);
        jacobi_step2f(&x[idx], &b[idx], &x_new[idx], a, b_inv, m);
      }
    }

    for (int j = 1; j < n - 1; ++j) {
      int idx_l = util::at2_x(m, 0, j); // left boundary
      int idx_r = util::at2_x(m, m - 1, j); // right boundary

      jacobi_step2f(&x[idx_l], &b[idx_l], &x_new[idx_l], a, b_inv, m, 0, 1, 1, 1);
      jacobi_step2f(&x[idx_r], &b[idx_r], &x_new[idx_r], a, b_inv, m, 1, 0, 1, 1);
    }

    for (int i = 1; i < m - 1; ++i) {
      int idx_b = util::at2_x(m, i, 0); // bottom boundary
      int idx_t = util::at2_x(m, i, n - 1); // top boundary

      jacobi_step2f(&x[idx_b], &b[idx_b], &x_new[idx_b], a, b_inv, m, 1, 1, 0, 1);
      jacobi_step2f(&x[idx_t], &b[idx_t], &x_new[idx_t], a, b_inv, m, 1, 1, 1, 0);
    }
  }
}

#endif //FLUID_MATH_JACOBI_H
