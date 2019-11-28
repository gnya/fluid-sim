/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 23
 */

#ifndef FLUID_MATH_ADVECT_H
#define FLUID_MATH_ADVECT_H

#include <ctgmath>

namespace fluid::math {
  template <typename Func>
  void advect_scaled(const float *q, const float *u, float *q_new, float dt, float dx,
                     int q_m, int q_n, int u_m, int u_n, int dim, Func lerp_func) {
    using namespace util;

    float scale_m = (float) q_m / (float) u_m;
    float scale_n = (float) q_n / (float) u_n;

    for (int u_j = 0; u_j < u_n; u_j++) {
      for (int u_i = 0; u_i < u_m; u_i++) {
        int i_min = (int) ((float) (u_i + 0) * scale_m);
        int i_max = (int) ((float) (u_i + 1) * scale_m);
        int j_min = (int) ((float) (u_j + 0) * scale_n);
        int j_max = (int) ((float) (u_j + 1) * scale_n);

        // dq_x = -u[u_i, u_j].x * dt * dx * scale_m
        // dq_y = -u[u_i, u_j].y * dt * dx * scale_n
        float dq_x = -(*u++) * dt * dx * scale_m;
        float dq_y = -(*u++) * dt * dx * scale_n;

        int dq_i = (int) std::floor(dq_x);
        int dq_j = (int) std::floor(dq_y);

        float dec_x = dq_x - (float) dq_i;
        float dec_y = dq_y - (float) dq_j;

        bool c_i = i_min + dq_i >= 0 && i_max + dq_i < q_m - 1;
        bool c_j = j_min + dq_j >= 0 && j_max + dq_j < q_n - 1;

        if (c_i && c_j) {
          const float *q_ptr = &q[at(dim, q_m, 0, i_min + dq_i, j_min + dq_j)];
          float *q_new_ptr = &q_new[at(dim, q_m, 0, i_min, j_min)];

          for (int j = 0; j < j_max - j_min; j++) {
            for (int i = 0; i < i_max - i_min; i++) {
              lerp_func(q_ptr + i * dim, q_new_ptr + i * dim, dec_x, dec_y, q_m);
            }

            q_ptr += q_m * dim;
            q_new_ptr += q_m * dim;
          }
        } else {
          for (int j = j_min; j < j_max; j++) {
            int q_j = std::max(0, std::min(j + dq_j, q_n - 2));

            for (int i = i_min; i < i_max; i++) {
              int q_i = std::max(0, std::min(i + dq_i, q_m - 2));

              const float *q_ptr = &q[at(dim, q_m, 0, q_i, q_j)];
              float *q_new_ptr = &q_new[at(dim, q_m, 0, i, j)];

              lerp_func(q_ptr, q_new_ptr, dec_x, dec_y, q_m);
            }
          }
        }
      }
    }
  }

  void advect_scaled2f(const float *q, const float *u, float *q_new,
                          float dt, float dx, int q_m, int q_n, int u_m, int u_n) {
    return advect_scaled(q, u, q_new, dt, dx, q_m, q_n, u_m, u_n, 2, util::lerp2d2f);
  }

  void advect_scaled4f(const float *q, const float *u, float *q_new,
                       float dt, float dx, int q_m, int q_n, int u_m, int u_n) {
    return advect_scaled(q, u, q_new, dt, dx, q_m, q_n, u_m, u_n, 4, util::lerp2d4f);
  }

  void advect2f(const float *q, const float *u, float *q_new,
                    float dt, float dx, int m, int n) {
    using namespace util;

    const __m256i _idx_i = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    const __m256  _idx_f = _mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0);
    const __m256i _m     = _mm256_set1_epi32(m);
    const __m256i _min_i = _mm256_set1_epi32(0);
    const __m256i _min_j = _mm256_set1_epi32(0);
    const __m256i _max_i = _mm256_set1_epi32(m - 2);
    const __m256i _max_j = _mm256_set1_epi32(n - 2);
    const __m256  _dtdx  = _mm256_set1_ps(dt * dx);
    auto i_parallel_end = m - m % 8;

    for (int j = 0; j < n; ++j) {
      __m256 _j = _mm256_set1_ps((float) j);

      for (int i = 0; i < i_parallel_end; i += 8) {
        __m256 _i = _mm256_add_ps(_mm256_set1_ps((float) i), _idx_f);

        // calc dec_x and dec_y
        __m256 x = _mm256_i32gather_ps(u + 0, _idx_i, 2 * sizeof(float));
        __m256 y = _mm256_i32gather_ps(u + 1, _idx_i, 2 * sizeof(float));
        x = _mm256_sub_ps(_i, _mm256_mul_ps(x, _dtdx));
        y = _mm256_sub_ps(_j, _mm256_mul_ps(y, _dtdx));
        __m256 x_floor = _mm256_floor_ps(x);
        __m256 y_floor = _mm256_floor_ps(y);
        __m256 dec_x = _mm256_sub_ps(x, x_floor);
        __m256 dec_y = _mm256_sub_ps(y, y_floor);

        // calc q_i and q_j
        __m256i q_i = _mm256_cvtps_epi32(x_floor);
        __m256i q_j = _mm256_cvtps_epi32(y_floor);
        q_i = avx::clip_m256i(q_i, _min_i, _max_i);
        q_j = avx::clip_m256i(q_j, _min_j, _max_j);
        __m256i idx = _mm256_add_epi32(q_i, _mm256_mullo_epi32(q_j, _m));
        idx = _mm256_slli_epi32(idx, 1);

        avx::lerp2d2f(q, q_new, idx, dec_x, dec_y, m);
        q_new += 16; u += 16;
      }

      for (int i = i_parallel_end; i < m; ++i) {
        // calc dec_x and dec_y
        float x = (float) i - *(u + 0) * dt * dx;
        float y = (float) j - *(u + 1) * dt * dx;
        float x_floor = std::floor(x);
        float y_floor = std::floor(y);
        float dec_x = x - x_floor;
        float dec_y = y - y_floor;

        // calc q_i and q_j
        int q_i = clip((int) x_floor, 0, m - 2);
        int q_j = clip((int) y_floor, 0, n - 2);
        int idx = at2_x(m, q_i, q_j);

        lerp2d2f(q + idx, q_new, dec_x, dec_y, m);
        q_new += 2; u += 2;
      }
    }
  }
}

#endif //FLUID_MATH_ADVECT_H
