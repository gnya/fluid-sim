/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 15
 */

#ifndef FLUID_MATH_NOISE_H
#define FLUID_MATH_NOISE_H

#include <random>
#include <algorithm>

namespace fluid::math::noise {
  class PerlinNoise {
  private:
    // hash lookup table size
    const int _perlin_table_size = 256;

    // constant value for convert hash to radian
    const float _hash2rad = 2.0f * M_PI / 255.0f;

    std::vector<unsigned char> _table;

    float (*_func)(float);
  public:
    explicit PerlinNoise(int seed = 0, float (*func)(float) = util::lerp_func::quintic) {
      std::mt19937 engine(seed);

      _table.resize(_perlin_table_size);
      std::iota(_table.begin(), _table.end(), 0);
      std::shuffle(_table.begin(), _table.end(), engine);

      _func = func;
    }

    unsigned char table(int x, int y) {
      return _table[(_table[x % _perlin_table_size] + y) % _perlin_table_size];
    }

    float grad(int x_hash, int y_hash, float x, float y) {
      float theta = (float) table(x_hash, y_hash) * _hash2rad;

      return x * cos(theta) + y * sin(theta);
    }

    float operator ()(float x, float y) {
      // calculate unsigned x, y
      unsigned char x_hash = ((int) floor(x)) & 0xff;
      unsigned char y_hash = ((int) floor(y)) & 0xff;

      float xf = x - floor(x);
      float yf = y - floor(y);

      float g00 = grad(x_hash + 0, y_hash + 0, xf - 0, yf - 0);
      float g01 = grad(x_hash + 1, y_hash + 0, xf - 1, yf - 0);
      float g10 = grad(x_hash + 0, y_hash + 1, xf - 0, yf - 1);
      float g11 = grad(x_hash + 1, y_hash + 1, xf - 1, yf - 1);

      float gx0 = util::lerp(g00, g01, xf, _func);
      float gx1 = util::lerp(g10, g11, xf, _func);

      return (util::lerp(gx0, gx1, yf, _func) + 1) / 2; // [-1, 1] -> [0, 1]
    }
  };

  class OctavePerlinNoise : PerlinNoise {
  public:
    explicit OctavePerlinNoise(int seed = 0) : PerlinNoise(seed) {}

    float operator ()(float x, float y, int octave = 8, float persistance = 0.5f) {
      float n = 0, n_max = 0, f = 1, a = 1;

      for (int i = 0; i < octave; i++) {
        n += PerlinNoise::operator()(x * f, y * f) * a;
        n_max += a;

        f *= 2;
        a *= persistance;
      }

      return n / n_max; // normalize to [0, 1]
    }
  };
}

#endif //FLUID_MATH_NOISE_H
