/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 09
 */

#ifndef FLUID_SIM_FLUID_H
#define FLUID_SIM_FLUID_H

#include <cstring>
#include <ctgmath>

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

      for (int j = 0; j < q_n; j++) {
        int u_j = (int) ((float) j / scale_n);

        for (int i = 0; i < q_m; i++) {
          int u_i = (int) ((float) i / scale_m);

          float x = i - dt * dx * u[at2_x(u_m, u_i, u_j)];
          float y = j - dt * dx * u[at2_y(u_m, u_i, u_j)];

          // clip
          if (x < 0) {
            x = 0.5f;
          } else if (x >= q_m - 1) {
            x = q_m - 1 - 0.5f;
          }

          if (y < 0) {
            y = 0.5f;
          } else if (y >= q_n - 1) {
            y = q_n - 1 - 0.5f;
          }

          auto q_i = (int) floor(x);
          auto q_j = (int) floor(y);

          float dec_x = x - floor(x);
          float dec_y = y - floor(y);

          float q0_x = lerp(q[at2_x(q_m, q_i + 0, q_j + 0)], q[at2_x(q_m, q_i + 1, q_j + 0)], dec_x);
          float q0_y = lerp(q[at2_y(q_m, q_i + 0, q_j + 0)], q[at2_y(q_m, q_i + 1, q_j + 0)], dec_x);

          float q1_x = lerp(q[at2_x(q_m, q_i + 0, q_j + 1)], q[at2_x(q_m, q_i + 1, q_j + 1)], dec_x);
          float q1_y = lerp(q[at2_y(q_m, q_i + 0, q_j + 1)], q[at2_y(q_m, q_i + 1, q_j + 1)], dec_x);

          q_new[at2_x(q_m, i, j)] = lerp(q0_x, q1_x, dec_y);
          q_new[at2_y(q_m, i, j)] = lerp(q0_y, q1_y, dec_y);
        }
      }
    }

    void advect2f(float *q, float *q_new, float *u,
                  float dt, float dx, int m, int n) {
      for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
          float x = i - dt * dx * u[at2_x(m, i, j)];
          float y = j - dt * dx * u[at2_y(m, i, j)];

          // clip
          if (x < 0) {
            x = 0.5f;
          } else if (x >= m - 1) {
            x = m - 1 - 0.5f;
          }

          if (y < 0) {
            y = 0.5f;
          } else if (y >= n - 1) {
            y = n - 1 - 0.5f;
          }

          // rectangle index
          auto q_i = (int) floor(x);
          auto q_j = (int) floor(y);

          float dec_x = x - floor(x);
          float dec_y = y - floor(y);

          float q0_x = lerp(q[at2_x(m, q_i + 0, q_j + 0)], q[at2_x(m, q_i + 1, q_j + 0)], dec_x);
          float q0_y = lerp(q[at2_y(m, q_i + 0, q_j + 0)], q[at2_y(m, q_i + 1, q_j + 0)], dec_x);

          float q1_x = lerp(q[at2_x(m, q_i + 0, q_j + 1)], q[at2_x(m, q_i + 1, q_j + 1)], dec_x);
          float q1_y = lerp(q[at2_y(m, q_i + 0, q_j + 1)], q[at2_y(m, q_i + 1, q_j + 1)], dec_x);

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
          x_new[at2_x(m, i, j)] = (x[at2_x(m, i - 1, j)] + x[at2_x(m, i + 1, j)] +
                                   x[at2_x(m, i, j - 1)] + x[at2_x(m, i, j + 1)] +
                                   b[at2_x(m, i, j)] * a) * r_b;
          x_new[at2_y(m, i, j)] = (x[at2_y(m, i - 1, j)] + x[at2_y(m, i + 1, j)] +
                                   x[at2_y(m, i, j - 1)] + x[at2_y(m, i, j + 1)] +
                                   b[at2_y(m, i, j)] * a) * r_b;
        }
      }

      for (int j = 1; j < n - 1; j++) {
        // left boundary
        {
          x_new[at2_x(m, 0, j)] = (x[at2_x(m, 0 - 0, j)] + x[at2_x(m, 0 + 1, j)] +
                                   x[at2_x(m, 0, j - 1)] + x[at2_x(m, 0, j + 1)] +
                                   b[at2_x(m, 0, j)] * a) * r_b;
          x_new[at2_y(m, 0, j)] = (x[at2_y(m, 0 - 0, j)] + x[at2_y(m, 0 + 1, j)] +
                                   x[at2_y(m, 0, j - 1)] + x[at2_y(m, 0, j + 1)] +
                                   b[at2_y(m, 0, j)] * a) * r_b;
        }

        // right boundary
        {
          x_new[at2_x(m, m - 1, j)] = (x[at2_x(m, m - 1 - 1, j)] + x[at2_x(m, m - 1 + 0, j)] +
                                       x[at2_x(m, m - 1, j - 1)] + x[at2_x(m, m - 1, j + 1)] +
                                       b[at2_x(m, m - 1, j)] * a) * r_b;
          x_new[at2_y(m, m - 1, j)] = (x[at2_y(m, m - 1 - 1, j)] + x[at2_y(m, m - 1 + 0, j)] +
                                       x[at2_y(m, m - 1, j - 1)] + x[at2_y(m, m - 1, j + 1)] +
                                       b[at2_y(m, m - 1, j)] * a) * r_b;
        }
      }

      for (int i = 1; i < m - 1; i++) {
        // bottom boundary
        {
          x_new[at2_x(m, i, 0)] = (x[at2_x(m, i - 1, 0)] + x[at2_x(m, i + 1, 0)] +
                                   x[at2_x(m, i, 0 - 0)] + x[at2_x(m, i, 0 + 1)] +
                                   b[at2_x(m, i, 0)] * a) * r_b;
          x_new[at2_y(m, i, 0)] = (x[at2_y(m, i - 1, 0)] + x[at2_y(m, i + 1, 0)] +
                                   x[at2_y(m, i, 0 - 0)] + x[at2_y(m, i, 0 + 1)] +
                                   b[at2_y(m, i, 0)] * a) * r_b;
        }

        // top boundary
        {
          x_new[at2_x(m, i, n - 1)] = (x[at2_x(m, i - 1, n - 1)] + x[at2_x(m, i + 1, n - 1)] +
                                       x[at2_x(m, i, n - 1 - 1)] + x[at2_x(m, i, n - 1 + 0)] +
                                       b[at2_x(m, i, n - 1)] * a) * r_b;
          x_new[at2_y(m, i, n - 1)] = (x[at2_y(m, i - 1, n - 1)] + x[at2_y(m, i + 1, n - 1)] +
                                       x[at2_y(m, i, n - 1 - 1)] + x[at2_y(m, i, n - 1 + 0)] +
                                       b[at2_y(m, i, n - 1)] * a) * r_b;
        }
      }
    }

    void jacobi1f(float *x, float *x_new, float *b,
                  float a, float r_b, int m, int n) {
      // inner of boundary
      for (int j = 1; j < n - 1; j++) {
        for (int i = 1; i < m - 1; i++) {
          x_new[at(m, i, j)] = (x[at(m, i - 1, j)] + x[at(m, i + 1, j)] +
                                x[at(m, i, j - 1)] + x[at(m, i, j + 1)] +
                                b[at(m, i, j)] * a) * r_b;
        }
      }

      for (int j = 1; j < n - 1; j++) {
        // left boundary
        {
          x_new[at(m, 0, j)] = (x[at(m, 0 - 0, j)] + x[at(m, 0 + 1, j)] +
                                x[at(m, 0, j - 1)] + x[at(m, 0, j + 1)] +
                                b[at(m, 0, j)] * a) * r_b;
        }

        // right boundary
        {
          x_new[at(m, m - 1, j)] = (x[at(m, m - 1 - 1, j)] + x[at(m, m - 1 + 0, j)] +
                                    x[at(m, m - 1, j - 1)] + x[at(m, m - 1, j + 1)] +
                                    b[at(m, m - 1, j)] * a) * r_b;
        }
      }

      for (int i = 1; i < m - 1; i++) {
        // bottom boundary
        {
          x_new[at(m, i, 0)] = (x[at(m, i - 1, 0)] + x[at(m, i + 1, 0)] +
                                x[at(m, i, 0 - 0)] + x[at(m, i, 0 + 1)] +
                                b[at(m, i, 0)] * a) * r_b;
        }

        // top boundary
        {
          x_new[at(m, i, n - 1)] = (x[at(m, i - 1, n - 1)] + x[at(m, i + 1, n - 1)] +
                                    x[at(m, i, n - 1 - 1)] + x[at(m, i, n - 1 + 0)] +
                                    b[at(m, i, n - 1)] * a) * r_b;
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

    void boundary2f(float *p, float *u, int m, int n) {
      for (int j = 1; j < n - 1; j++) {
        // left boundary
        {
          p[at(m, 0, j)] = p[at(m, 0 + 1, j)];
          u[at2_x(m, 0, j)] = -u[at2_x(m, 0 + 1, j)];
          u[at2_y(m, 0, j)] = -u[at2_y(m, 0 + 1, j)];
        }

        // right boundary
        {
          p[at(m, m - 1, j)] = p[at(m, m - 1 - 1, j)];
          u[at2_x(m, m - 1, j)] = -u[at2_x(m, m - 1 - 1, j)];
          u[at2_y(m, m - 1, j)] = -u[at2_y(m, m - 1 - 1, j)];
        }
      }

      for (int i = 1; i < m - 1; i++) {
        // bottom boundary
        {
          p[at(m, i, 0)] = p[at(m, i, 0 + 1)];
          u[at2_x(m, i, 0)] = -u[at2_x(m, i, 0 + 1)];
          u[at2_y(m, i, 0)] = -u[at2_y(m, i, 0 + 1)];
        }

        // top boundary
        {
          p[at(m, i, n - 1)] = p[at(m, i, n - 1 - 1)];
          u[at2_x(m, i, n - 1)] = -u[at2_x(m, i, n - 1 - 1)];
          u[at2_y(m, i, n - 1)] = -u[at2_y(m, i, n - 1 - 1)];
        }
      }
    }
  }

  class Fluid {
  protected:
    int _u_m{}, _u_n{};

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
      _u_m = m;
      _u_n = n;

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
      util::advect2f(_u, _w, _u, _dt, _dx, _u_m, _u_n);
    }

    void diffuse(int n_jacob) {
      float a = (float) pow(_dx, 2) / (_v * _dt);
      float r_b = 1.0f / (4 + a);

      for (int i = 0; i < n_jacob; i++) {
        util::jacobi2f(_w, _w_tmp, _w, a, r_b, _u_m, _u_n);
        memcpy(_w, _w_tmp, _u_m * _u_n * 2 * sizeof(float));
      }
    }

    void add_force(float pos[2], float f[2], float r) {
      util::brush2f(_w, pos, f, r, _dt, _u_m, _u_n);
    }

    void projection(int n_jacob) {
      float a = - (float) pow(_dx, 2);
      float r_b = 1.0f / 4;

      util::divergence2f(_w, _w_div, _dx, _u_m, _u_n);

      for (int i = 0; i < n_jacob; i++) {
        util::jacobi1f(_p, _p_tmp, _w_div, a, r_b, _u_m, _u_n);
        memcpy(_p, _p_tmp, _u_m * _u_n * 1 * sizeof(float));
      }

      util::gradient2f(_p, _w, _u, _dx, _u_m, _u_n);
    }

    void boundary() {
      util::boundary2f(_p, _u, _u_m, _u_n);
    }
  };
}

#endif //FLUID_SIM_FLUID_H
