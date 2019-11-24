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

    int table(int i, int j) const {
      return _table[(_table[i % _perlin_table_size] + j) % _perlin_table_size];
    }

    float grad(int x_i, int y_i, float x, float y) const {
      const float sqrt2_inv = 1 / M_SQRT2;
      unsigned char hash = table(x_i, y_i);

      switch (hash & 0b111) {
        case 0b000:
          return +x;
        case 0b001:
          return +y;
        case 0b010:
          return -y;
        case 0b011:
          return -x;
        case 0b100:
          return sqrt2_inv * x + sqrt2_inv * y;
        case 0b101:
          return -sqrt2_inv * x + sqrt2_inv * y;
        case 0b110:
          return sqrt2_inv * x - sqrt2_inv * y;
        case 0b111:
          return -sqrt2_inv * x - sqrt2_inv * y;
        default:
          return 0;
      }
    }

    float operator()(float x, float y) const {
      using namespace util;

      // calculate unsigned x, y
      unsigned char x_i = ((int) std::floor(x)) & 0xff;
      unsigned char y_i = ((int) std::floor(y)) & 0xff;

      float x_f = x - std::floor(x);
      float y_f = y - std::floor(y);

      float g00 = grad(x_i + 0, y_i + 0, x_f - 0, y_f - 0);
      float g01 = grad(x_i + 1, y_i + 0, x_f - 1, y_f - 0);
      float g10 = grad(x_i + 0, y_i + 1, x_f - 0, y_f - 1);
      float g11 = grad(x_i + 1, y_i + 1, x_f - 1, y_f - 1);

      float gx0 = lerp(g00, g01, x_f, _func);
      float gx1 = lerp(g10, g11, x_f, _func);

      return (lerp(gx0, gx1, y_f, _func) + 1) / 2; // [-1, 1] -> [0, 1]
    }
  };

  class OctavePerlinNoise : PerlinNoise {
  public:
    explicit OctavePerlinNoise(int seed = 0) : PerlinNoise(seed) {}

    float operator ()(float x, float y, int octave = 8, float persistance = 0.5f) const {
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

namespace fluid::math::noise::avx {
  class PerlinNoise {
  private:
    // hash lookup table size
    const int _perlin_table_size = 256;

    std::vector<int> _table;

    __m128 (*_func)(__m128);
  public:
    explicit PerlinNoise(int seed = 0, __m128 (*func)(__m128) = util::avx::lerp_func::quintic) {
      std::mt19937 engine(seed);

      _table.resize(_perlin_table_size);
      std::iota(_table.begin(), _table.end(), 0);
      std::shuffle(_table.begin(), _table.end(), engine);

      _func = func;
    }

    int table(int i, int j) const {
      return _table[(_table[i % _perlin_table_size] + j) % _perlin_table_size];
    }

    __m128i table(__m128i x_i, __m128i y_i) const {
      __attribute__((aligned(32))) int i[4] = {0};
      __attribute__((aligned(32))) int j[4] = {0};
      _mm_store_si128((__m128i*) i, x_i);
      _mm_store_si128((__m128i*) j, y_i);

      return _mm_set_epi32(
          table(i[3], j[3]), table(i[2], j[2]),
          table(i[1], j[1]), table(i[0], j[0])
      );
    }

    __m128 grad(__m128i x_i, __m128i y_i, __m128 x, __m128 y) const {
      const __m128 _one           = _mm_set1_ps(1);
      const __m128 _two           = _mm_set1_ps(2);
      const __m128 _sqrt2_inv     = _mm_set1_ps(1 / M_SQRT2 - 1);
      const __m128i _hash_mask    = _mm_set1_epi32(0b111);
      const __m128i _sign_x_mask  = _mm_set1_epi32(0b001);
      const __m128i _sign_y_mask  = _mm_set1_epi32(0b010);
      const __m128i _group_mask   = _mm_set1_epi32(0b100);

      __m128i hash = _mm_and_si128(table(x_i, y_i), _hash_mask);

      // 0bit 0/1 : x +/-, 1bit 0/1 : y +/-
      __m128i sign_x_i = _mm_and_si128(hash, _sign_x_mask);
      __m128i sign_y_i = _mm_srli_epi32(_mm_and_si128(hash, _sign_y_mask), 1);
      __m128 sign_x_f = _mm_cvtepi32_ps(sign_x_i);
      __m128 sign_y_f = _mm_cvtepi32_ps(sign_y_i);

      // group0 : 0, group1 : 1 / sqrt2
      __m128i group_i = _mm_srli_epi32(_mm_and_si128(hash, _group_mask), 2);
      __m128 group_f = _mm_cvtepi32_ps(group_i);
      __m128 c_xy = _mm_add_ps(group_f, _mm_mul_ps(group_f, _sqrt2_inv));

      // group0 : 0b000, 0b011 -> x, 0b001, 0b010 -> y
      __m128i y_flag = _mm_xor_si128(sign_x_i, sign_y_i);
      __m128 group_f_inv = _mm_sub_ps(_one, group_f);
      __m128 c_x = _mm_sub_ps(_one, _mm_cvtepi32_ps(y_flag));
      __m128 c_y = _mm_sub_ps(_one, c_x);
      c_x = _mm_add_ps(_mm_mul_ps(c_x, group_f_inv), c_xy);
      c_y = _mm_add_ps(_mm_mul_ps(c_y, group_f_inv), c_xy);

      c_x = _mm_mul_ps(c_x, _mm_sub_ps(_one, _mm_mul_ps(sign_x_f, _two)));
      c_y = _mm_mul_ps(c_y, _mm_sub_ps(_one, _mm_mul_ps(sign_y_f, _two)));

      return _mm_add_ps(_mm_mul_ps(x, c_x), _mm_mul_ps(y, c_y));
    }

    __m128 operator ()(__m128 x, __m128 y) const {
      using namespace util::avx;

      const __m128i _one_i = _mm_set1_epi32(1);
      const __m128 _one_f = _mm_set1_ps(1);
      const __m128 _half_f = _mm_set1_ps(0.5f);
      const __m128i _0xff = _mm_set1_epi32(0xff);

      // calculate unsigned x, y
      __m128 x_floor = _mm_floor_ps(x);
      __m128 y_floor = _mm_floor_ps(y);
      __m128i x_i0 = _mm_and_si128(_mm_cvtps_epi32(x_floor), _0xff), x_i1;
      __m128i y_i0 = _mm_and_si128(_mm_cvtps_epi32(y_floor), _0xff), y_i1;
      x_i1 = _mm_add_epi32(x_i0, _one_i);
      y_i1 = _mm_add_epi32(y_i0, _one_i);

      __m128 x_f0 = _mm_sub_ps(x ,x_floor), x_f1;
      __m128 y_f0 = _mm_sub_ps(y, y_floor), y_f1;
      x_f1 = _mm_sub_ps(x_f0, _one_f);
      y_f1 = _mm_sub_ps(y_f0, _one_f);

      __m128 g00 = grad(x_i0, y_i0, x_f0, y_f0);
      __m128 g01 = grad(x_i1, y_i0, x_f1, y_f0);
      __m128 g10 = grad(x_i0, y_i1, x_f0, y_f1);
      __m128 g11 = grad(x_i1, y_i1, x_f1, y_f1);

      __m128 gx0 = lerp(g00, g01, x_f0, _func);
      __m128 gx1 = lerp(g10, g11, x_f0, _func);
      __m128 g   = lerp(gx0, gx1, y_f0, _func);

      return _mm_mul_ps(_mm_add_ps(g, _one_f), _half_f); // [-1, 1] -> [0, 1]
    }
  };

  class OctavePerlinNoise : PerlinNoise {
  public:
    explicit OctavePerlinNoise(int seed = 0) : PerlinNoise(seed) {}

    __m128 operator ()(__m128 x, __m128 y, int octave = 8, float persistance = 0.5f) const {
      const __m128 _two = _mm_set1_ps(2);
      const __m128 _persistance = _mm_set1_ps(persistance);

      __m128 n     = _mm_set1_ps(0);
      __m128 n_max = _mm_set1_ps(0);
      __m128 f     = _mm_set1_ps(1);
      __m128 a     = _mm_set1_ps(1);

      for (int i = 0; i < octave; i++) {
        __m128 p = PerlinNoise::operator()(_mm_mul_ps(x, f), _mm_mul_ps(y, f));
        n = _mm_add_ps(n, _mm_mul_ps(p, a));
        n_max = _mm_add_ps(n_max, a);

        f = _mm_mul_ps(f, _two);
        a = _mm_mul_ps(a, _persistance);
      }

      return _mm_div_ps(n, n_max); // normalize to [0, 1]
    }
  };
}

namespace fluid::math::noise {
  void add_noise2f(float *u, int m, int n, int seed_x, int seed_y, float amp, float scale) {
    using namespace util;
    using namespace math;

    const __m128 _scale = _mm_set1_ps(scale);
    const __m128 _idx = _mm_set_ps(3, 2, 1, 0);
    const __m256 _amp = _mm256_set1_ps(amp);
    const __m256 _one = _mm256_set1_ps(1);
    const __m256 _two = _mm256_set1_ps(2);
    auto i_parallel_end = m - m % 4;

    noise::avx::OctavePerlinNoise x_noise(seed_x);
    noise::avx::OctavePerlinNoise y_noise(seed_y);

    for (int j = 0; j < n; j++) {
      __m128 _j = _mm_set1_ps((float) j);

      for (int i = 0; i < i_parallel_end; i += 4) {
        __m128 _i = _mm_set1_ps((float) i);

        __m128 x = _mm_mul_ps(_mm_add_ps(_i, _idx), _scale);
        __m128 y = _mm_mul_ps(_j, _scale);

        __m256 noise = util::avx::marge_m256_2f(x_noise(x, y), y_noise(x, y));
        __m256 a = _mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(noise, _two), _one), _amp);

        __m256 tmp = _mm256_loadu_ps(&u[at2_x(m, i, j)]);
        _mm256_storeu_ps(&u[at2_x(m, i, j)], _mm256_add_ps(tmp, a));
      }

      if (i_parallel_end < m) {
        __m128 _i = _mm_set1_ps((float) i_parallel_end);

        __m128 x = _mm_mul_ps(_mm_add_ps(_i, _idx), _scale);
        __m128 y = _mm_mul_ps(_j, _scale);

        __m256 noise = util::avx::marge_m256_2f(x_noise(x, y), y_noise(x, y));
        __m256 a = _mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(noise, _two), _one), _amp);

        __attribute__((aligned(32))) float _a[8] = {0};
        _mm256_store_ps(_a, a);

        for (int i = i_parallel_end; i < m; ++i) {
          u[at2_x(m, i, j)] += _a[(i - i_parallel_end) * 2 + 0];
          u[at2_y(m, i, j)] += _a[(i - i_parallel_end) * 2 + 1];
        }
      }
    }
  }
}

#endif //FLUID_MATH_NOISE_H
