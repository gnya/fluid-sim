/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 09
 */

#ifndef FLUID_SIM_FLUID_H
#define FLUID_SIM_FLUID_H

#include <ctgmath>

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

    void advect_scaled2f(float *q, float *q_new, float *u,
                         float dt, float dx,
                         int q_m, int q_n,
                         int u_m, int u_n) {
      float scale_m = (float) q_m / (float) u_m;
      float scale_n = (float) q_n / (float) u_n;

      for (int u_j = 0; u_j < u_n; u_j++) {
        for (int u_i = 0; u_i < u_m; u_i++) {
          float du_x = -u[at2_x(u_m, u_i, u_j)] * dt * dx;
          float du_y = -u[at2_y(u_m, u_i, u_j)] * dt * dx;

          int du_i = (int) floor(du_x);
          int du_j = (int) floor(du_y);

          float dec_x = du_x - du_i;
          float dec_y = du_y - du_j;

          for (int j = u_j * scale_n; j < (u_j + 1) * scale_n; j++) {
            int q_j = std::max(0, std::min(j + du_j, q_n - 2));

            for (int i = u_i * scale_m; i < (u_i + 1) * scale_m; i++) {
              int q_i = std::max(0, std::min(i + du_i, q_m - 2));

              float q00_x = q[at2_x(q_m, q_i + 0, q_j + 0)];
              float q00_y = q[at2_y(q_m, q_i + 0, q_j + 0)];

              float q01_x = q[at2_x(q_m, q_i + 0, q_j + 1)];
              float q01_y = q[at2_y(q_m, q_i + 0, q_j + 1)];

              float q10_x = q[at2_x(q_m, q_i + 1, q_j + 0)];
              float q10_y = q[at2_y(q_m, q_i + 1, q_j + 0)];

              float q11_x = q[at2_x(q_m, q_i + 1, q_j + 1)];
              float q11_y = q[at2_y(q_m, q_i + 1, q_j + 1)];

              float q0_x = lerp(q00_x, q10_x, dec_x);
              float q0_y = lerp(q00_y, q10_y, dec_x);

              float q1_x = lerp(q01_x, q11_x, dec_x);
              float q1_y = lerp(q01_y, q11_y, dec_x);

              q_new[at2_x(q_m, i, j)] = lerp(q0_x, q1_x, dec_y);
              q_new[at2_y(q_m, i, j)] = lerp(q0_y, q1_y, dec_y);
            }
          }
        }
      }
    }

    void advect2f(float *q, float *q_new, float *u,
                  float dt, float dx, int m, int n) {
      for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
          float x = i - u[at2_x(m, i, j)] * dt * dx;
          float y = j - u[at2_y(m, i, j)] * dt * dx;

          int q_i = (int) floor(x);
          int q_j = (int) floor(y);

          float dec_x = x - q_i;
          float dec_y = y - q_j;

          q_i = std::max(0, std::min(q_i, m - 2));
          q_j = std::max(0, std::min(q_j, n - 2));

          float q00_x = q[at2_x(m, q_i + 0, q_j + 0)];
          float q00_y = q[at2_y(m, q_i + 0, q_j + 0)];

          float q01_x = q[at2_x(m, q_i + 0, q_j + 1)];
          float q01_y = q[at2_y(m, q_i + 0, q_j + 1)];

          float q10_x = q[at2_x(m, q_i + 1, q_j + 0)];
          float q10_y = q[at2_y(m, q_i + 1, q_j + 0)];

          float q11_x = q[at2_x(m, q_i + 1, q_j + 1)];
          float q11_y = q[at2_y(m, q_i + 1, q_j + 1)];

          float q0_x = lerp(q00_x, q10_x, dec_x);
          float q0_y = lerp(q00_y, q10_y, dec_x);

          float q1_x = lerp(q01_x, q11_x, dec_x);
          float q1_y = lerp(q01_y, q11_y, dec_x);

          q_new[at2_x(m, i, j)] = lerp(q0_x, q1_x, dec_y);
          q_new[at2_y(m, i, j)] = lerp(q0_y, q1_y, dec_y);
        }
      }
    }

    void jacobi2f(float *x, float *x_new, float *b,
                  float a, float r_b, int m, int n) {
      // inner of boundary
      for (int j = 1; j < n - 1; j++) {
        for (int i = 1; i < m - 1; i++) {
          float x00_x = x[at2_x(m, i, j - 1)];
          float x00_y = x[at2_y(m, i, j - 1)];
          float x01_x = x[at2_x(m, i, j + 1)];
          float x01_y = x[at2_y(m, i, j + 1)];
          float x10_x = x[at2_x(m, i - 1, j)];
          float x10_y = x[at2_y(m, i - 1, j)];
          float x11_x = x[at2_x(m, i + 1, j)];
          float x11_y = x[at2_y(m, i + 1, j)];
          float b_x_a = b[at2_x(m, i, j)] * a;
          float b_y_a = b[at2_y(m, i, j)] * a;

          x_new[at2_x(m, i, j)] = (x00_x + x01_x + x10_x + x11_x + b_x_a) * r_b;
          x_new[at2_y(m, i, j)] = (x00_y + x01_y + x10_y + x11_y + b_y_a) * r_b;
        }
      }

      for (int j = 1; j < n - 1; j++) {
        // left boundary
        {
          float x00_x = x[at2_x(m, 0, j - 1)];
          float x00_y = x[at2_y(m, 0, j - 1)];
          float x01_x = x[at2_x(m, 0, j + 1)];
          float x01_y = x[at2_y(m, 0, j + 1)];
          float x10_x = x[at2_x(m, 0 - 0, j)];
          float x10_y = x[at2_y(m, 0 - 0, j)];
          float x11_x = x[at2_x(m, 0 + 1, j)];
          float x11_y = x[at2_y(m, 0 + 1, j)];
          float b_x_a = b[at2_x(m, 0, j)] * a;
          float b_y_a = b[at2_y(m, 0, j)] * a;

          x_new[at2_x(m, 0, j)] = (x00_x + x01_x + x10_x + x11_x + b_x_a) * r_b;
          x_new[at2_y(m, 0, j)] = (x00_y + x01_y + x10_y + x11_y + b_y_a) * r_b;
        }

        // right boundary
        {
          float x00_x = x[at2_x(m, m - 1, j - 1)];
          float x00_y = x[at2_y(m, m - 1, j - 1)];
          float x01_x = x[at2_x(m, m - 1, j + 1)];
          float x01_y = x[at2_y(m, m - 1, j + 1)];
          float x10_x = x[at2_x(m, m - 1 - 1, j)];
          float x10_y = x[at2_y(m, m - 1 - 1, j)];
          float x11_x = x[at2_x(m, m - 1 + 0, j)];
          float x11_y = x[at2_y(m, m - 1 + 0, j)];
          float b_x_a = b[at2_x(m, m - 1, j)] * a;
          float b_y_a = b[at2_y(m, m - 1, j)] * a;

          x_new[at2_x(m, m - 1, j)] = (x00_x + x01_x + x10_x + x11_x + b_x_a) * r_b;
          x_new[at2_y(m, m - 1, j)] = (x00_y + x01_y + x10_y + x11_y + b_y_a) * r_b;
        }
      }

      for (int i = 1; i < m - 1; i++) {
        // bottom boundary
        {
          float x00_x = x[at2_x(m, i, 0 - 0)];
          float x00_y = x[at2_y(m, i, 0 - 0)];
          float x01_x = x[at2_x(m, i, 0 + 1)];
          float x01_y = x[at2_y(m, i, 0 + 1)];
          float x10_x = x[at2_x(m, i - 1, 0)];
          float x10_y = x[at2_y(m, i - 1, 0)];
          float x11_x = x[at2_x(m, i + 1, 0)];
          float x11_y = x[at2_y(m, i + 1, 0)];
          float b_x_a = b[at2_x(m, i, 0)] * a;
          float b_y_a = b[at2_y(m, i, 0)] * a;

          x_new[at2_x(m, i, 0)] = (x00_x + x01_x + x10_x + x11_x + b_x_a) * r_b;
          x_new[at2_y(m, i, 0)] = (x00_y + x01_y + x10_y + x11_y + b_y_a) * r_b;
        }

        // top boundary
        {
          float x00_x = x[at2_x(m, i, n - 1 - 1)];
          float x00_y = x[at2_y(m, i, n - 1 - 1)];
          float x01_x = x[at2_x(m, i, n - 1 + 0)];
          float x01_y = x[at2_y(m, i, n - 1 + 0)];
          float x10_x = x[at2_x(m, i - 1, n - 1)];
          float x10_y = x[at2_y(m, i - 1, n - 1)];
          float x11_x = x[at2_x(m, i + 1, n - 1)];
          float x11_y = x[at2_y(m, i + 1, n - 1)];
          float b_x_a = b[at2_x(m, i, n - 1)] * a;
          float b_y_a = b[at2_y(m, i, n - 1)] * a;

          x_new[at2_x(m, i, n - 1)] = (x00_x + x01_x + x10_x + x11_x + b_x_a) * r_b;
          x_new[at2_y(m, i, n - 1)] = (x00_y + x01_y + x10_y + x11_y + b_y_a) * r_b;
        }
      }
    }

    void jacobi1f(float *x, float *x_new, float *b,
                  float a, float r_b, int m, int n) {
      // inner of boundary
      for (int j = 1; j < n - 1; j++) {
        for (int i = 1; i < m - 1; i++) {
          float x00 = x[at(m, i, j - 1)];
          float x01 = x[at(m, i, j + 1)];
          float x10 = x[at(m, i - 1, j)];
          float x11 = x[at(m, i + 1, j)];
          float b_a = b[at(m, i, j)] * a;

          x_new[at(m, i, j)] = (x00 + x01 + x10 + x11 + b_a) * r_b;
        }
      }

      for (int j = 1; j < n - 1; j++) {
        // left boundary
        {
          float x00 = x[at(m, 0, j - 1)];
          float x01 = x[at(m, 0, j + 1)];
          float x10 = x[at(m, 0 - 0, j)];
          float x11 = x[at(m, 0 + 1, j)];
          float b_a = b[at(m, 0, j)] * a;

          x_new[at(m, 0, j)] = (x00 + x01 + x10 + x11 + b_a) * r_b;
        }

        // right boundary
        {
          float x00 = x[at(m, m - 1, j - 1)];
          float x01 = x[at(m, m - 1, j + 1)];
          float x10 = x[at(m, m - 1 - 1, j)];
          float x11 = x[at(m, m - 1 + 0, j)];
          float b_a = b[at(m, m - 1, j)] * a;

          x_new[at(m, m - 1, j)] = (x00 + x01 + x10 + x11 + b_a) * r_b;
        }
      }

      for (int i = 1; i < m - 1; i++) {
        // bottom boundary
        {
          float x00 = x[at(m, i, 0 - 0)];
          float x01 = x[at(m, i, 0 + 1)];
          float x10 = x[at(m, i - 1, 0)];
          float x11 = x[at(m, i + 1, 0)];
          float b_a = b[at(m, i, 0)] * a;

          x_new[at(m, i, 0)] = (x00 + x01 + x10 + x11 + b_a) * r_b;
        }

        // top boundary
        {
          float x00 = x[at(m, i - 1, n - 1)];
          float x01 = x[at(m, i, n - 1 - 1)];
          float x10 = x[at(m, i, n - 1 + 0)];
          float x11 = x[at(m, i + 1, n - 1)];
          float b_a = b[at(m, i, n - 1)] * a;

          x_new[at(m, i, n - 1)] = (x00 + x01 + x10 + x11 + b_a) * r_b;
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

      _u     = new float[m * n * 2];
      _w     = new float[m * n * 2];
      _w_tmp = new float[m * n * 2];
      _w_div = new float[m * n * 1];
      _p     = new float[m * n * 1];
      _p_tmp = new float[m * n * 1];
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
