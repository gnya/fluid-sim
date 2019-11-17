/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 09
 */

#ifndef FLUID_SIM_MAP_FLUID_H
#define FLUID_SIM_MAP_FLUID_H

#include "fluid.h"

#include <png++/png.hpp>

namespace fluid {
  using namespace png;

  class MapFluid : public Fluid {
  protected:
    int _width{}, _height{};

    float *_map{};
    float *_map_tmp{};
  public:
    MapFluid() = default;

    MapFluid(int width, int height, float scale_x, float scale_y, float dt, float dx, float v)
      : Fluid((int) ((float) width * scale_x), (int) ((float) height * scale_y), dt, dx, v) {
      _width = width;
      _height = height;

      _map     = new float[width * height * 2];
      _map_tmp = new float[width * height * 2];

      // reset map
      for (int j = 0; j < _height; j++) {
        for (int i = 0; i < _width; i++) {
          _map[util::at2_x(_width, i, j)] = (float) i;
          _map[util::at2_y(_width, i, j)] = (float) j;
        }
      }
    }

    ~MapFluid() override {
      std::cout << "[MEM RELEASE] map fluid" << std::endl;
      delete[] _map;
      delete[] _map_tmp;
    }

    void advect_map() {
      util::advect_scaled2f(_map, _map_tmp, _u, _dt, _dx, _width, _height, _m, _n);
      std::swap(_map, _map_tmp);
    }

    template<typename pixel>
    void mapping(const pixel_buffer<pixel> &src, solid_pixel_buffer<pixel> &dst) {
      // mapping
      for (int j = 0; j < _height; j++) {
        for (int i = 0; i < _width; i++) {
          int x = (int) _map[util::at2_x(_width, i, j)];
          int y = (int) _map[util::at2_y(_width, i, j)];

          dst[j][i] = src[_height - y - 1][x];
        }
      }
    }
  };
}

#endif //FLUID_SIM_MAP_FLUID_H
