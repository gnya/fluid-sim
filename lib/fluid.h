/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 09
 */

#ifndef FLUID_FLUID_H
#define FLUID_FLUID_H

#include "math/math.h"

namespace fluid {
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

      _u     = new (std::align_val_t{32}) float[m * n * 2]();
      _w     = new (std::align_val_t{32}) float[m * n * 2]();
      _w_tmp = new (std::align_val_t{32}) float[m * n * 2]();
      _w_div = new (std::align_val_t{32}) float[m * n * 1]();
      _p     = new (std::align_val_t{32}) float[m * n * 1]();
      _p_tmp = new (std::align_val_t{32}) float[m * n * 1]();
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
      math::advect2f(_u, _u, _w, _dt, _dx, _m, _n);
    }

    void diffuse(int n_jacob) {
      const float a = (float) pow(_dx, 2) / (_v * _dt);
      const float b_inv = 1.0f / (4 + a);

      for (int i = 0; i < n_jacob / 2; i++) {
        math::jacobi2f(_w, _w, _w_tmp, a, b_inv, _m, _n);
        math::jacobi2f(_w_tmp, _w_tmp, _w, a, b_inv, _m, _n);
      }
    }

    void add_force(float pos[2], float f[2], float r) {
      math::brush2f(_w, pos, f, r, _dt, _m, _n);
    }

    void projection(int n_jacob) {
      const float a = - (float) pow(_dx, 2);
      const float b_inv = 1.0f / 4;

      math::divergence2f(_w, _w_div, _dx, _m, _n);

      for (int i = 0; i < n_jacob / 2; i++) {
        math::jacobi1f(_p, _w_div, _p_tmp, a, b_inv, _m, _n);
        math::jacobi1f(_p_tmp, _w_div, _p, a, b_inv, _m, _n);
      }

      math::gradient2f(_p, _w, _u, _dx, _m, _n);
    }

    void left_boundary() {
      math::left_boundary2f(_p, _u, _m, _n);
    }

    void right_boundary() {
      math::right_boundary2f(_p, _u, _m, _n);
    }

    void bottom_boundary() {
      math::bottom_boundary2f(_p, _u, _m, _n);
    }

    void top_boundary() {
      math::top_boundary2f(_p, _u, _m, _n);
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
    for (int j = 0; j < f._n; j++) {
      for (int i = 0; i < f._m; i++) {
        f._u[util::at2_x(f._m, i, j)] += v[0];
        f._u[util::at2_y(f._m, i, j)] += v[1];
      }
    }
  }

  void Fluid::accelerate_by_perlin_noise(Fluid &f, int seed_x, int seed_y, float amp, float scale) {
    math::noise::add_noise2f(f._u, f._m, f._n, seed_x, seed_y, amp, scale);
  }
}

#endif //FLUID_FLUID_H
