/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 23
 */

#ifndef FLUID_MATH_ADVECT_H
#define FLUID_MATH_ADVECT_H

#include <ctgmath>

namespace fluid::math {
  void advect_scaled2f(const float *q, const float *u, float *q_new,
                       float dt, float dx, int q_m, int q_n, int u_m, int u_n) {
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
          const float *q_ptr = &q[at2_x(q_m, i_min + dq_i, j_min + dq_j)];
          float *q_new_ptr = &q_new[at2_x(q_m, i_min, j_min)];

          for (int j = 0; j < j_max - j_min; j++) {
            for (int i = 0; i < i_max - i_min; i++) {
              lerp2d2f(q_ptr + i * 2, q_new_ptr + i * 2, dec_x, dec_y, q_m);
            }

            q_ptr += q_m * 2;
            q_new_ptr += q_m * 2;
          }
        } else {
          for (int j = j_min; j < j_max; j++) {
            int q_j = std::max(0, std::min(j + dq_j, q_n - 2));

            for (int i = i_min; i < i_max; i++) {
              int q_i = std::max(0, std::min(i + dq_i, q_m - 2));

              const float *q_ptr = &q[at2_x(q_m, q_i, q_j)];
              float *q_new_ptr = &q_new[at2_x(q_m, i, j)];

              lerp2d2f(q_ptr, q_new_ptr, dec_x, dec_y, q_m);
            }
          }
        }
      }
    }
  }

  void advect2f(const float *q, const float *u, float *q_new,
                float dt, float dx, int m, int n) {
    using namespace util;

    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        // x = i - u[i, j].x * dt * dx;
        // y = j - u[i, j].y * dt * dx;
        float x = (float) i - *u++ * dt * dx;
        float y = (float) j - *u++ * dt * dx;

        int q_i = (int) std::floor(x);
        int q_j = (int) std::floor(y);

        float dec_x = x - (float) q_i;
        float dec_y = y - (float) q_j;

        q_i = std::max(0, std::min(q_i, m - 2));
        q_j = std::max(0, std::min(q_j, n - 2));

        lerp2d2f(&q[at2_x(m, q_i, q_j)], q_new, dec_x, dec_y, m);

        q_new += 2;
      }
    }
  }
}

#endif //FLUID_MATH_ADVECT_H
