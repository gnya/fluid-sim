/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 15
 */


#ifndef FLUID_INK_FLUID_H
#define FLUID_INK_FLUID_H

#include "fluid.h"

#include <png++/png.hpp>

namespace fluid {
  class InkFluid : public Fluid {
  protected:
    int _width{}, _height{};

    float *_ink{};
    float *_ink_tmp{};
    float *_ink_src{};
  public:
    InkFluid() = default;

    InkFluid(int width, int height, float scale_x, float scale_y, float dt, float dx, float v)
      : Fluid((int) ((float) width * scale_x), (int) ((float) height * scale_y), dt, dx, v) {
      _width = width;
      _height = height;

      _ink     = new (std::align_val_t{32}) float[width * height * 4]();
      _ink_tmp = new (std::align_val_t{32}) float[width * height * 4]();
      _ink_src = new (std::align_val_t{32}) float[width * height * 4]();
    }

    ~InkFluid() override {
      std::cout << "[MEM RELEASE] ink fluid" << std::endl;
      delete[] _ink;
      delete[] _ink_tmp;
      delete[] _ink_src;
    }

    void advect_ink() {
      math::advect_scaled4f(_ink, _u, _ink_tmp, _dt, _dx, _width, _height, _m, _n);
      math::add4f(_ink_tmp, _ink_src, _ink, _width, _height, 0, 255);
    }

    void add_ink(float pos[2], float c[4], float r) {
      math::brush4f(_ink, pos, c, r, _dt, _width, _height);
    }

    void add_ink_src(float pos[2], float c[4], float r) {
      math::brush4f(_ink_src, pos, c, r, _dt, _width, _height);
    }

    template<typename pixel>
    void set_ink(const png::pixel_buffer<pixel> &src);

    template<typename pixel>
    void get_ink(png::solid_pixel_buffer<pixel> &dst);
  };

  template<>
  void InkFluid::set_ink(const png::pixel_buffer<png::rgb_pixel> &src) {
    for (int j = 0; j < _height; j++) {
      for (int i = 0; i < _width; i++) {
        _ink[util::at(4, _width, 0, i, j)] = src[j][i].red;
        _ink[util::at(4, _width, 1, i, j)] = src[j][i].green;
        _ink[util::at(4, _width, 2, i, j)] = src[j][i].blue;
        _ink[util::at(4, _width, 3, i, j)] = 255;
      }
    }
  }

  template<>
  void InkFluid::set_ink(const png::pixel_buffer<png::rgba_pixel> &src) {
    for (int j = 0; j < _height; j++) {
      for (int i = 0; i < _width; i++) {
        _ink[util::at(4, _width, 0, i, j)] = src[j][i].red;
        _ink[util::at(4, _width, 1, i, j)] = src[j][i].green;
        _ink[util::at(4, _width, 2, i, j)] = src[j][i].blue;
        _ink[util::at(4, _width, 3, i, j)] = src[j][i].alpha;
      }
    }
  }

  template <>
  void InkFluid::get_ink(png::solid_pixel_buffer<png::rgb_pixel> &dst) {
    for (int j = 0; j < _height; j++) {
      for (int i = 0; i < _width; i++) {
        int x = (int) _ink[util::at(4, _width, 0, i, j)];
        int y = (int) _ink[util::at(4, _width, 1, i, j)];
        int z = (int) _ink[util::at(4, _width, 2, i, j)];

        dst[j][i] = png::rgb_pixel(x, y, z);
      }
    }
  }

  template <>
  void InkFluid::get_ink(png::solid_pixel_buffer<png::rgba_pixel> &dst) {
    for (int j = 0; j < _height; j++) {
      for (int i = 0; i < _width; i++) {
        int x = (int) _ink[util::at(4, _width, 0, i, j)];
        int y = (int) _ink[util::at(4, _width, 1, i, j)];
        int z = (int) _ink[util::at(4, _width, 2, i, j)];
        int w = (int) _ink[util::at(4, _width, 3, i, j)];

        dst[j][i] = png::rgba_pixel(x, y, z, w);
      }
    }
  }
}

#endif //FLUID_INK_FLUID_H
