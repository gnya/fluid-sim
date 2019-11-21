/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 09
 */

#ifndef FLUID_SIM_FLUID_H
#define FLUID_SIM_FLUID_H

#include <ctgmath>
#include <immintrin.h>

#include "noise.h"

namespace fluid {
  using namespace std;

  namespace util {
    inline int at(int m, int i, int j) {
      return i + j * m;
    }

    inline int at2_x(int m, int i, int j) {
      return 0 + i * 2 + j * 2 * m;
    }

    inline int at2_y(int m, int i, int j) {
      return 1 + i * 2 + j * 2 * m;
    }

    inline float lerp(float a, float b, float t) {
      return a + (b - a) * t;
    }

    inline void lerp2d2f(const float *q, float *q_new, float t_x, float t_y, int m) {
      // q0_x = lerp(q[q_i + 0, q_j + 0].x, q[q_i + 1, q_j + 0].x, t_x)
      float q0_x = lerp(*(q + 0), *(q + 2), t_x);
      // q0_y = lerp(q[q_i + 0, q_j + 0].y, q[q_i + 1, q_j + 0].y, t_x)
      float q0_y = lerp(*(q + 1), *(q + 3), t_x);

      q += m * 2;

      // q1_x = lerp(q[q_i + 0, q_j + 1].x, q[q_i + 1, q_j + 1].x, t_x)
      float q1_x = lerp(*(q + 0), *(q + 2), t_x);
      // q1_y = lerp(q[q_i + 0, q_j + 1].y, q[q_i + 1, q_j + 1].y, t_x)
      float q1_y = lerp(*(q + 1), *(q + 3), t_x);

      // q_new[i, j].x = lerp(q0_x, q1_x, t_y);
      *(q_new + 0) = lerp(q0_x, q1_x, t_y);
      // q_new[i, j].y = lerp(q0_y, q1_y, t_y);
      *(q_new + 1) = lerp(q0_y, q1_y, t_y);
    }

    void advect_scaled2f(float *q, float *q_new, float *u,
                         float dt, float dx,
                         int q_m, int q_n,
                         int u_m, int u_n) {
      float scale_m = (float) q_m / (float) u_m;
      float scale_n = (float) q_n / (float) u_n;

      for (int u_j = 0; u_j < u_n; u_j++) {
        for (int u_i = 0; u_i < u_m; u_i++) {
          // du_x = -u[u_i, u_j].x * dt * dx
          float du_x = -(*u++) * dt * dx;
          // du_y = -u[u_i, u_j].y * dt * dx
          float du_y = -(*u++) * dt * dx;

          int du_i = (int) floor(du_x);
          int du_j = (int) floor(du_y);

          float dec_x = du_x - du_i;
          float dec_y = du_y - du_j;

          int i_bgn = u_i * scale_m, i_end = (u_i + 1) * scale_m;
          int j_bgn = u_j * scale_n, j_end = (u_j + 1) * scale_n;

          bool c_i = i_bgn + du_i >= 0 && i_end + du_i < q_m - 1;
          bool c_j = j_bgn + du_j >= 0 && j_end + du_j < q_n - 1;

          if (c_i && c_j) {
            float *q_ptr = &q[at2_x(q_m, i_bgn + du_i, j_bgn + du_j)];
            float *q_new_ptr = &q_new[at2_x(q_m, i_bgn, j_bgn)];

            for (int j = 0; j < j_end - j_bgn; j++) {
              for (int i = 0; i < i_end - i_bgn; i++) {
                lerp2d2f(q_ptr + i * 2, q_new_ptr + i * 2, dec_x, dec_y, q_m);
              }

              q_ptr += q_m * 2;
              q_new_ptr += q_m * 2;
            }
          } else {
            for (int j = j_bgn; j < j_end; j++) {
              int q_j = std::max(0, std::min(j + du_j, q_n - 2));

              for (int i = i_bgn; i < i_end; i++) {
                int q_i = std::max(0, std::min(i + du_i, q_m - 2));

                float *q_ptr = &q[at2_x(q_m, q_i, q_j)];
                float *q_new_ptr = &q_new[at2_x(q_m, i, j)];

                lerp2d2f(q_ptr, q_new_ptr, dec_x, dec_y, q_m);
              }
            }
          }
        }
      }
    }

    void advect2f(float *q, float *q_new, float *u,
                  float dt, float dx, int m, int n) {
      for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
          // x = i - u[i, j].x * dt * dx;
          float x = i - *u++ * dt * dx;
          // y = j - u[i, j].y * dt * dx;
          float y = j - *u++ * dt * dx;

          int q_i = (int) floor(x);
          int q_j = (int) floor(y);

          float dec_x = x - q_i;
          float dec_y = y - q_j;

          q_i = std::max(0, std::min(q_i, m - 2));
          q_j = std::max(0, std::min(q_j, n - 2));

          lerp2d2f(&q[at2_x(m, q_i, q_j)], q_new, dec_x, dec_y, m);

          q_new += 2;
        }
      }
    }

    inline void jacobi_step2f(float *x, float *x_new, float *b,
                              float a, float r_b, int m,
                              int c_l = 1, int c_r = 1,
                              int c_b = 1, int c_t = 1) {
      float x00_x = *(x - 2 * m * c_b + 0); // x[i, j - 1].x
      float x00_y = *(x - 2 * m * c_b + 1); // x[i, j - 1].y
      float x01_x = *(x + 2 * m * c_t + 0); // x[i, j + 1].x
      float x01_y = *(x + 2 * m * c_t + 1); // x[i, j + 1].y
      float x10_x = *(x - 2 * 1 * c_l + 0); // x[i - 1, j].x
      float x10_y = *(x - 2 * 1 * c_l + 1); // x[i - 1, j].y
      float x11_x = *(x + 2 * 1 * c_r + 0); // x[i + 1, j].x
      float x11_y = *(x + 2 * 1 * c_r + 1); // x[i + 1, j].y
      float b_x_a = *(b + 0) * a; // b[i, j].x * a;
      float b_y_a = *(b + 1) * a; // b[i, j].y * a;

      // x_new[i, j].x
      *(x_new + 0) = (x00_x + x01_x + x10_x + x11_x + b_x_a) * r_b;
      // x_new[i, j].y
      *(x_new + 1) = (x00_y + x01_y + x10_y + x11_y + b_y_a) * r_b;
    }

    void jacobi2f(float *x, float *x_new, float *b,
                  float a, float r_b, int m, int n) {
      __m256 _a = _mm256_set1_ps(a);
      __m256 _r_b = _mm256_set1_ps(r_b);
      auto i_parallel_end = (m - 1) - (m - 2) % 4;

      // inner of boundary
      for (int j = 1; j < n - 1; ++j) {
        for (int i = 1; i < i_parallel_end; i += 4) {
          __m256 x00 = _mm256_load_ps(&x[at2_x(m, i, j - 1)]);
          __m256 x01 = _mm256_load_ps(&x[at2_x(m, i, j + 1)]);
          __m256 x10 = _mm256_load_ps(&x[at2_x(m, i - 1, j)]);
          __m256 x11 = _mm256_load_ps(&x[at2_x(m, i + 1, j)]);
          __m256 _b  = _mm256_load_ps(&b[at2_x(m, i, j)]);

          __m256 tmp;

          tmp = _mm256_mul_ps(_a, _b);
          tmp = _mm256_add_ps(tmp, x00);
          tmp = _mm256_add_ps(tmp, x01);
          tmp = _mm256_add_ps(tmp, x10);
          tmp = _mm256_add_ps(tmp, x11);
          tmp = _mm256_mul_ps(tmp, _r_b);

          _mm256_storeu_ps(&x_new[at2_x(m, i, j)], tmp);
        }

        for (int i = i_parallel_end; i < m - 1; ++i) {
          int idx = at2_x(m, i, j);

          jacobi_step2f(&x[idx], &x_new[idx], &b[idx], a, r_b, m);
        }
      }

      for (int j = 1; j < n - 1; ++j) {
        // left boundary
        {
          int idx = at2_x(m, 0, j);

          jacobi_step2f(&x[idx], &x_new[idx], &b[idx], a, r_b, m, 0, 1, 1, 1);
        }

        // right boundary
        {
          int idx = at2_x(m, m - 1, j);

          jacobi_step2f(&x[idx], &x_new[idx], &b[idx], a, r_b, m, 1, 0, 1, 1);
        }
      }

      for (int i = 1; i < m - 1; ++i) {
        // bottom boundary
        {
          int idx = at2_x(m, i, 0);

          jacobi_step2f(&x[idx], &x_new[idx], &b[idx], a, r_b, m, 1, 1, 0, 1);
        }

        // top boundary
        {
          int idx = at2_x(m, i, n - 1);

          jacobi_step2f(&x[idx], &x_new[idx], &b[idx], a, r_b, m, 1, 1, 1, 0);
        }
      }
    }

    inline void jacobi_step1f(float *x, float *x_new, float *b,
                              float a, float r_b, int m,
                              int c_l = 1, int c_r = 1,
                              int c_b = 1, int c_t = 1) {
      float x00 = *(x - m * c_b); // x[i, j - 1]
      float x01 = *(x + m * c_t); // x[i, j + 1]
      float x10 = *(x - 1 * c_l); // x[i - 1, j]
      float x11 = *(x + 1 * c_r); // x[i + 1, j]
      float b_a = *b * a; // b[i, j] * a;

      // x_new[i, j]
      *x_new = (x00 + x01 + x10 + x11 + b_a) * r_b;
    }

    void jacobi1f(float *x, float *x_new, float *b,
                  float a, float r_b, int m, int n) {
      __m256 _a = _mm256_set1_ps(a);
      __m256 _r_b = _mm256_set1_ps(r_b);
      auto i_parallel_end = (m - 1) - (m - 2) % 8;

      // inner of boundary
      for (int j = 1; j < n - 1; ++j) {
        for (int i = 1; i < i_parallel_end; i += 8) {
          __m256 x00 = _mm256_load_ps(&x[at(m, i, j - 1)]);
          __m256 x01 = _mm256_load_ps(&x[at(m, i, j + 1)]);
          __m256 x10 = _mm256_load_ps(&x[at(m, i - 1, j)]);
          __m256 x11 = _mm256_load_ps(&x[at(m, i + 1, j)]);
          __m256 _b  = _mm256_load_ps(&b[at(m, i, j)]);

          __m256 tmp;

          tmp = _mm256_mul_ps(_a, _b);
          tmp = _mm256_add_ps(tmp, x00);
          tmp = _mm256_add_ps(tmp, x01);
          tmp = _mm256_add_ps(tmp, x10);
          tmp = _mm256_add_ps(tmp, x11);
          tmp = _mm256_mul_ps(tmp, _r_b);

          _mm256_storeu_ps(&x_new[at(m, i, j)], tmp);
        }

        for (int i = i_parallel_end; i < m - 1; ++i) {
          int idx = at(m, i, j);

          jacobi_step1f(&x[idx], &x_new[idx], &b[idx], a, r_b, m);
        }
      }

      for (int j = 1; j < n - 1; ++j) {
        // left boundary
        {
          int idx = at(m, 0, j);

          jacobi_step1f(&x[idx], &x_new[idx], &b[idx], a, r_b, m, 0, 1, 1, 1);
        }

        // right boundary
        {
          int idx = at(m, m - 1, j);

          jacobi_step1f(&x[idx], &x_new[idx], &b[idx], a, r_b, m, 1, 0, 1, 1);
        }
      }

      for (int i = 1; i < m - 1; ++i) {
        // bottom boundary
        {
          int idx = at(m, i, 0);

          jacobi_step1f(&x[idx], &x_new[idx], &b[idx], a, r_b, m, 1, 1, 0, 1);
        }

        // top boundary
        {
          int idx = at(m, i, n - 1);

          jacobi_step1f(&x[idx], &x_new[idx], &b[idx], a, r_b, m, 1, 1, 1, 0);
        }
      }
    }

    void brush2f(float *x, float p_idx[2], float v[2],
                 float r, float dt, int m, int n) {
      for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
          auto n2 = (float) (pow(i - p_idx[0], 2) + pow(j - p_idx[1], 2));
          auto s = exp(-n2 / r);

          x[at2_x(m, i, j)] += v[0] * dt * s;
          x[at2_y(m, i, j)] += v[1] * dt * s;
        }
      }
    }

    void divergence2f(float *x, float *x_div,
                      float dx, int m, int n) {
      // inner of boundary
      for (int j = 1; j < n - 1; j++) {
        for (int i = 1; i < m - 1; i++) {
          x_div[at(m, i, j)] = dx * ((x[at2_x(m, i + 1, j)] - x[at2_x(m, i - 1, j)]) +
                                     (x[at2_y(m, i, j + 1)] - x[at2_y(m, i, j - 1)])) / 2;
        }
      }

      for (int j = 1; j < n - 1; j++) {
        // left boundary
        {
          x_div[at(m, 0, j)] = dx * ((x[at2_x(m, 0 + 1, j)] - x[at2_x(m, 0 - 0, j)]) +
                                     (x[at2_y(m, 0, j + 1)] - x[at2_y(m, 0, j - 1)])) / 2;
        }

        // right boundary
        {
          x_div[at(m, m - 1, j)] = dx * ((x[at2_x(m, m - 1 + 0, j)] - x[at2_x(m, m - 1 - 1, j)]) +
                                         (x[at2_y(m, m - 1, j + 1)] - x[at2_y(m, m - 1, j - 1)])) / 2;
        }
      }

      for (int i = 1; i < m - 1; i++) {
        // bottom boundary
        {
          x_div[at(m, i, 0)] = dx * ((x[at2_x(m, i + 1, 0)] - x[at2_x(m, i - 1, 0)]) +
                                     (x[at2_y(m, i, 0 + 1)] - x[at2_y(m, i, 0 - 0)])) / 2;
        }

        // top boundary
        {
          x_div[at(m, i, n - 1)] = dx * ((x[at2_x(m, i + 1, n - 1)] - x[at2_x(m, i - 1, n - 1)]) +
                                         (x[at2_y(m, i, n - 1 + 0)] - x[at2_y(m, i, n - 1 - 1)])) / 2;
        }
      }
    }

    void gradient2f(float *p, float *w, float *u,
                    float dx, int m, int n) {
      // inner of boundary
      for (int j = 1; j < n - 1; j++) {
        for (int i = 1; i < m - 1; i++) {
          u[at2_x(m, i, j)] = w[at2_x(m, i, j)];
          u[at2_y(m, i, j)] = w[at2_y(m, i, j)];

          u[at2_x(m, i, j)] -= dx * (p[at(m, i + 1, j)] - p[at(m, i - 1, j)]) / 2;
          u[at2_y(m, i, j)] -= dx * (p[at(m, i, j + 1)] - p[at(m, i, j - 1)]) / 2;
        }
      }

      for (int j = 1; j < n - 1; j++) {
        // left boundary
        {
          u[at2_x(m, 0, j)] = w[at2_x(m, 0, j)];
          u[at2_y(m, 0, j)] = w[at2_y(m, 0, j)];

          u[at2_x(m, 0, j)] -= dx * (p[at(m, 0 + 1, j)] - p[at(m, 0 - 0, j)]) / 2;
          u[at2_y(m, 0, j)] -= dx * (p[at(m, 0, j + 1)] - p[at(m, 0, j - 1)]) / 2;
        }

        // right boundary
        {
          u[at2_x(m, m - 1, j)] = w[at2_x(m, m - 1, j)];
          u[at2_y(m, m - 1, j)] = w[at2_y(m, m - 1, j)];

          u[at2_x(m, m - 1, j)] -= dx * (p[at(m, m - 1 + 0, j)] - p[at(m, m - 1 - 1, j)]) / 2;
          u[at2_y(m, m - 1, j)] -= dx * (p[at(m, m - 1, j + 1)] - p[at(m, m - 1, j - 1)]) / 2;
        }
      }

      for (int i = 1; i < m - 1; i++) {
        // bottom boundary
        {
          u[at2_x(m, i, 0)] = w[at2_x(m, i, 0)];
          u[at2_y(m, i, 0)] = w[at2_y(m, i, 0)];

          u[at2_x(m, i, 0)] -= dx * (p[at(m, i + 1, 0)] - p[at(m, i - 1, 0)]) / 2;
          u[at2_y(m, i, 0)] -= dx * (p[at(m, i, 0 + 1)] - p[at(m, i, 0 - 0)]) / 2;
        }

        // top boundary
        {
          u[at2_x(m, i, n - 1)] = w[at2_x(m, i, n - 1)];
          u[at2_y(m, i, n - 1)] = w[at2_y(m, i, n - 1)];

          u[at2_x(m, i, n - 1)] -= dx * (p[at(m, i + 1, n - 1)] - p[at(m, i - 1, n - 1)]) / 2;
          u[at2_y(m, i, n - 1)] -= dx * (p[at(m, i, n - 1 + 0)] - p[at(m, i, n - 1 - 1)]) / 2;
        }
      }
    }

    void left_boundary2f(float *p, float *u, int m, int n) {
      for (int j = 1; j < n - 1; j++) {
        p[at(m, 0, j)] = p[at(m, 0 + 1, j)];
        u[at2_x(m, 0, j)] = -u[at2_x(m, 0 + 1, j)];
        u[at2_y(m, 0, j)] = -u[at2_y(m, 0 + 1, j)];
      }
    }

    void right_boundary2f(float *p, float *u, int m, int n) {
      for (int j = 1; j < n - 1; j++) {
        p[at(m, m - 1, j)] = p[at(m, m - 1 - 1, j)];
        u[at2_x(m, m - 1, j)] = -u[at2_x(m, m - 1 - 1, j)];
        u[at2_y(m, m - 1, j)] = -u[at2_y(m, m - 1 - 1, j)];
      }
    }

    void bottom_boundary2f(float *p, float *u, int m, int n) {
      for (int i = 1; i < m - 1; i++) {
        p[at(m, i, 0)] = p[at(m, i, 0 + 1)];
        u[at2_x(m, i, 0)] = -u[at2_x(m, i, 0 + 1)];
        u[at2_y(m, i, 0)] = -u[at2_y(m, i, 0 + 1)];
      }
    }

    void top_boundary2f(float *p, float *u, int m, int n) {
      for (int i = 1; i < m - 1; i++) {
        p[at(m, i, n - 1)] = p[at(m, i, n - 1 - 1)];
        u[at2_x(m, i, n - 1)] = -u[at2_x(m, i, n - 1 - 1)];
        u[at2_y(m, i, n - 1)] = -u[at2_y(m, i, n - 1 - 1)];
      }
    }
  }

  class Fluid {
  protected:
    int _m{}, _n{};

    float _dt{}, _dx{}, _v{};

    float *_u{};
    float *_w{};
    float *_w_tmp{};
    float *_w_div{};
    float *_p{};
    float *_p_tmp{};
  public:
    Fluid() = default;

    Fluid(int m, int n, float dt, float dx, float v) {
      _m = m;
      _n = n;

      _dt = dt;
      _dx = dx;
      _v  = v;

      _u     = new (std::align_val_t{32}) float[m * n * 2];
      _w     = new (std::align_val_t{32}) float[m * n * 2];
      _w_tmp = new (std::align_val_t{32}) float[m * n * 2];
      _w_div = new (std::align_val_t{32}) float[m * n * 1];
      _p     = new (std::align_val_t{32}) float[m * n * 1];
      _p_tmp = new (std::align_val_t{32}) float[m * n * 1];
    }

    virtual ~Fluid() {
      std::cout << "[MEM RELEASE] fluid" << std::endl;
      delete[] _u;
      delete[] _w;
      delete[] _w_tmp;
      delete[] _w_div;
      delete[] _p;
      delete[] _p_tmp;
    }

    void advect_velocity() {
      util::advect2f(_u, _w, _u, _dt, _dx, _m, _n);
    }

    void diffuse(int n_jacob) {
      float a = (float) pow(_dx, 2) / (_v * _dt);
      float r_b = 1.0f / (4 + a);

      for (int i = 0; i < n_jacob / 2; i++) {
        util::jacobi2f(_w, _w_tmp, _w, a, r_b, _m, _n);
        util::jacobi2f(_w_tmp, _w, _w_tmp, a, r_b, _m, _n);
      }
    }

    void add_force(float pos[2], float f[2], float r) {
      util::brush2f(_w, pos, f, r, _dt, _m, _n);
    }

    void projection(int n_jacob) {
      float a = - (float) pow(_dx, 2);
      float r_b = 1.0f / 4;

      util::divergence2f(_w, _w_div, _dx, _m, _n);

      for (int i = 0; i < n_jacob / 2; i++) {
        util::jacobi1f(_p, _p_tmp, _w_div, a, r_b, _m, _n);
        util::jacobi1f(_p_tmp, _p, _w_div, a, r_b, _m, _n);
      }

      util::gradient2f(_p, _w, _u, _dx, _m, _n);
    }

    void left_boundary() {
      util::left_boundary2f(_p, _u, _m, _n);
    }

    void right_boundary() {
      util::right_boundary2f(_p, _u, _m, _n);
    }

    void bottom_boundary() {
      util::bottom_boundary2f(_p, _u, _m, _n);
    }

    void top_boundary() {
      util::top_boundary2f(_p, _u, _m, _n);
    }

    void boundary() {
      left_boundary(); right_boundary();
      bottom_boundary(); top_boundary();
    }

    static void accelerate_by_single_vector(Fluid &f, float v[2]);

    static void accelerate_by_perlin_noise(Fluid &f, int x_seed = 0, int y_seed = 1,
                                           float amplitude = 1, float scale = 0.01f);
  };

  void Fluid::accelerate_by_single_vector(Fluid &f, float v[2]) {
    for (int j = 1; j < f._n; j++) {
      for (int i = 1; i < f._m; i++) {
        f._u[util::at2_x(f._m, i, j)] += v[0];
        f._u[util::at2_y(f._m, i, j)] += v[1];
      }
    }
  }

  void Fluid::accelerate_by_perlin_noise(Fluid &f, int x_seed, int y_seed,
                                         float amplitude, float scale) {
    noise::OctavePerlinNoise x_noise(x_seed);
    noise::OctavePerlinNoise y_noise(y_seed);

    for (int j = 1; j < f._n; j++) {
      for (int i = 1; i < f._m; i++) {
        float x_n = x_noise(i * scale, j * scale);
        float y_n = y_noise(i * scale, j * scale);

        f._u[util::at2_x(f._m, i, j)] += (x_n * 2 - 1) * amplitude;
        f._u[util::at2_y(f._m, i, j)] += (y_n * 2 - 1) * amplitude;
      }
    }
  }
}

#endif //FLUID_SIM_FLUID_H
