/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 23
 */

#ifndef FLUID_MATH_UTIL_H
#define FLUID_MATH_UTIL_H

#include <chrono>

namespace fluid::util {
  template<class Process>
  void timer(Process p, const char* prefix, int loop = 1, bool enable = true) {
    using namespace std::chrono;

    auto bgn = high_resolution_clock::now();
    for (int i = 0; i < loop; i++) p();
    auto end = high_resolution_clock::now();
    auto t = duration_cast<nanoseconds>(end - bgn).count() / 1000.0f;

    if (enable) {
      std::cout << "[debug] (" << prefix << ") time: ";
      if (t < 1000.0f) {
        std::cout << t << "ns" << std::endl;
      } else if (t / 1000.0f < 1000.0f) {
        std::cout << t / 1000.0f << "ms" << std::endl;
      } else {
        std::cout << t / 1000.0f / 1000.0f << "s" << std::endl;
      }
    }
  }

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

  namespace lerp_func {
    // quintic equation
    float quintic(float t) {
      return 6 * pow(t, 5.0f) - 15 * pow(t, 4.0f) + 10 * pow(t, 3.0f);
    }
  }

  inline float lerp(float a, float b, float t, float (*f)(float)) {
    return a + (b - a) * f(t);
  }

  inline void lerp2d2f(const float *q, float *q_new, float t_x, float t_y, int m) {
    // q0_x = lerp(q[q_i + 0, q_j + 0].x, q[q_i + 1, q_j + 0].x, t_x)
    // q0_y = lerp(q[q_i + 0, q_j + 0].y, q[q_i + 1, q_j + 0].y, t_x)
    float q0_x = lerp(*(q + 0), *(q + 2), t_x);
    float q0_y = lerp(*(q + 1), *(q + 3), t_x);

    q += m * 2;

    // q1_x = lerp(q[q_i + 0, q_j + 1].x, q[q_i + 1, q_j + 1].x, t_x)
    // q1_y = lerp(q[q_i + 0, q_j + 1].y, q[q_i + 1, q_j + 1].y, t_x)
    float q1_x = lerp(*(q + 0), *(q + 2), t_x);
    float q1_y = lerp(*(q + 1), *(q + 3), t_x);

    // q_new[i, j].x = lerp(q0_x, q1_x, t_y);
    // q_new[i, j].y = lerp(q0_y, q1_y, t_y);
    *(q_new + 0) = lerp(q0_x, q1_x, t_y);
    *(q_new + 1) = lerp(q0_y, q1_y, t_y);
  }

  namespace avx {
    void show_m128i(__m128i v, const char* prefix) {
      __attribute__((aligned(32))) int t[4] = {0};
      _mm_store_si128((__m128i*) t, v);

      std::cout << "[debug] (" << prefix << ")";
      std::cout << " 0:" << t[0] << " 1:" << t[1];
      std::cout << " 2:" << t[2] << " 3:" << t[3];
      std::cout << std::endl;
    }

    void show_m128(__m128 v, const char* prefix) {
      __attribute__((aligned(32))) float t[4] = {0};
      _mm_store_ps(t, v);

      std::cout << "[debug] (" << prefix << ")";
      std::cout << " 0:" << t[0] << " 1:" << t[1];
      std::cout << " 2:" << t[2] << " 3:" << t[3];
      std::cout << std::endl;
    }

    void show_m256(__m256 v, const char* prefix) {
      __attribute__((aligned(32))) float t[8] = {0};
      _mm256_store_ps(t, v);

      std::cout << "[debug] (" << prefix << ")";
      std::cout << " 0:" << t[0] << " 1:" << t[1];
      std::cout << " 2:" << t[2] << " 3:" << t[3];
      std::cout << " 4:" << t[4] << " 5:" << t[5];
      std::cout << " 6:" << t[6] << " 7:" << t[7];
      std::cout << std::endl;
    }

    inline __m256 set_m256_2f(float v[2]) {
      return _mm256_set_ps(v[1], v[0], v[1], v[0], v[1], v[0], v[1], v[0]);
    }

    inline __m256 clip_m256(__m256 x, __m256 min, __m256 max) {
      return _mm256_min_ps(_mm256_max_ps(x, min), max);
    }

    namespace lerp_func {
      inline __m128 quintic(__m128 x) {
        // y = 6 * x ^ 5 - 15 * t ^ 4 + 10 * x ^ 3;
        __m128 x3 = _mm_mul_ps(_mm_mul_ps(x, x), x);
        __m128 x4 = _mm_mul_ps(x3, x);
        __m128 x5 = _mm_mul_ps(x4, x);
        x3 = _mm_mul_ps(x3, _mm_set1_ps(10));
        x4 = _mm_mul_ps(x4, _mm_set1_ps(15));
        x5 = _mm_mul_ps(x5, _mm_set1_ps(6));

        return _mm_add_ps(x3, _mm_sub_ps(x5, x4));
      }
    }

    inline __m128 lerp(__m128 a, __m128 b, __m128 t, __m128 (*f)(__m128)) {
      // y = a + (b - a) * f(t);
      return _mm_add_ps(a, _mm_mul_ps(_mm_sub_ps(b, a), f(t)));
    }

    inline __m256 load_m256_2f(const float *x) {
      __m256 x0 = _mm256_loadu_ps(x + 0);
      __m256 x1 = _mm256_loadu_ps(x + 8);
      __m256 x2 = _mm256_permute2f128_ps(x1, x0, 0b11);

      x0 = _mm256_shuffle_ps(x0, x1, _MM_SHUFFLE(2, 0, 2, 0));
      x2 = _mm256_shuffle_ps(x2, x2, _MM_SHUFFLE(2, 0, 2, 0));

      return _mm256_blend_ps(x0, x2, 0b00111100);
    }

    inline __m256 marge_m256_2f(__m128 x, __m128 y) {
      __m256 x_ = _mm256_castps128_ps256(x);
      __m256 y_ = _mm256_castps128_ps256(y);
      __m256 _y = _mm256_permute2f128_ps(y_, x_ /* any __m256 */, 0b11);
      __m256 _x = _mm256_permute2f128_ps(x_, y_ /* any __m256 */, 0b11);

      __m256 xy = _mm256_blend_ps(x_, _y, 0b11110000);
      __m256 yx = _mm256_blend_ps(y_, _x, 0b11110000);

      __m256 oxxoxoox = _mm256_shuffle_ps(yx, xy, _MM_SHUFFLE(3, 1, 0, 2));
      __m256 xooxoxxo = _mm256_shuffle_ps(xy, yx, _MM_SHUFFLE(1, 3, 2, 0));

      return _mm256_blend_ps(oxxoxoox, xooxoxxo, 0b01101001);
    }

    inline __m256 exp_m256(__m256 x) {
      const __m256 _exp_hi = _mm256_set1_ps( 88.3762626647949f);
      const __m256 _exp_lo = _mm256_set1_ps(-88.3762626647949f);

      const __m256 _cephes_ln2     = _mm256_set1_ps(0.693359375); // ln2
      const __m256 _cephes_ln2_inv = _mm256_set1_ps(1.44269504088896341); // 1 / ln2

      const __m256 _cephes_exp_p0 = _mm256_set1_ps(1.9875691500e-4); // 1 / 7!
      const __m256 _cephes_exp_p1 = _mm256_set1_ps(1.3981999507e-3); // 1 / 6!
      const __m256 _cephes_exp_p2 = _mm256_set1_ps(8.3334519073e-3); // 1 / 5!
      const __m256 _cephes_exp_p3 = _mm256_set1_ps(4.1665795894e-2); // 1 / 4!
      const __m256 _cephes_exp_p4 = _mm256_set1_ps(1.6666665459e-1); // 1 / 3!
      const __m256 _cephes_exp_p5 = _mm256_set1_ps(5.0000001201e-1); // 1 / 2!

      const __m256 _one  = _mm256_set1_ps(1);
      const __m256i _0x7f = _mm256_set1_epi32(0x7f);

      // clip x
      x = clip_m256(x, _exp_lo, _exp_hi);

      // exp(x) = exp(g + n * ln2) = exp(g) * 2 ^ n
      // n <- round(g / ln2 + n) , (|z / ln2| < 1)
      __m256 n_f = _mm256_mul_ps(x, _cephes_ln2_inv);
      n_f = _mm256_round_ps(n_f, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
      // g <- x - n * ln2
      __m256 g = _mm256_sub_ps(x, _mm256_mul_ps(n_f, _cephes_ln2));
      __m256 gg = _mm256_mul_ps(g, g);

      // calculate taylor series
      __m256 y = _cephes_exp_p0;
      y = _mm256_add_ps(_mm256_mul_ps(y, g), _cephes_exp_p1);
      y = _mm256_add_ps(_mm256_mul_ps(y, g), _cephes_exp_p2);
      y = _mm256_add_ps(_mm256_mul_ps(y, g), _cephes_exp_p3);
      y = _mm256_add_ps(_mm256_mul_ps(y, g), _cephes_exp_p4);
      y = _mm256_add_ps(_mm256_mul_ps(y, g), _cephes_exp_p5);
      y = _mm256_add_ps(_mm256_mul_ps(y, gg), _mm256_add_ps(g, _one));

      // 2 ^ n
      __m256i n_i = _mm256_cvttps_epi32(n_f);
      n_i = _mm256_add_epi32(n_i, _0x7f);
      n_i = _mm256_slli_epi32(n_i, 23);
      __m256 pow2n = _mm256_castsi256_ps(n_i);

      return _mm256_mul_ps(y, pow2n);
    }
  }
}

#endif //FLUID_MATH_UTIL_H
