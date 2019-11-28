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

  template <typename T>
  void show(T *v, int size, const char* prefix) {
    std::cout << "[debug] (" << prefix << ")";
    for (int i = 0; i < size; i++) std::cout << " " << i << ":" << v[i];
    std::cout << std::endl;
  }

  inline int at(int m, int i, int j) {
    return i + j * m;
  }

  inline int at(int dim, int m, int i, int j, int k) {
    return i + (j + k * m) * dim;
  }

  inline int at2_x(int m, int i, int j) {
    return at(2, m, 0, i, j);
  }

  inline int at2_y(int m, int i, int j) {
    return at(2, m, 1, i, j);
  }

  template <typename T>
  inline T clip(T t, T min, T max) {
    return std::max(min, std::min(t, max));
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

  inline void lerp2d4f(const float *q, float *q_new, float t_x, float t_y, int m) {
    // q0_x = lerp(q[q_i + 0, q_j + 0].x, q[q_i + 1, q_j + 0].x, t_x)
    // q0_y = lerp(q[q_i + 0, q_j + 0].y, q[q_i + 1, q_j + 0].y, t_x)
    // q0_z = lerp(q[q_i + 0, q_j + 0].z, q[q_i + 1, q_j + 0].z, t_x)
    // q0_w = lerp(q[q_i + 0, q_j + 0].w, q[q_i + 1, q_j + 0].w, t_x)
    float q0_x = lerp(*(q + 0), *(q + 4), t_x);
    float q0_y = lerp(*(q + 1), *(q + 5), t_x);
    float q0_z = lerp(*(q + 2), *(q + 6), t_x);
    float q0_w = lerp(*(q + 3), *(q + 7), t_x);

    q += m * 4;

    // q1_x = lerp(q[q_i + 0, q_j + 1].x, q[q_i + 1, q_j + 1].x, t_x)
    // q1_y = lerp(q[q_i + 0, q_j + 1].y, q[q_i + 1, q_j + 1].y, t_x)
    // q1_z = lerp(q[q_i + 0, q_j + 1].z, q[q_i + 1, q_j + 1].z, t_x)
    // q1_w = lerp(q[q_i + 0, q_j + 1].w, q[q_i + 1, q_j + 1].w, t_x)
    float q1_x = lerp(*(q + 0), *(q + 4), t_x);
    float q1_y = lerp(*(q + 1), *(q + 5), t_x);
    float q1_z = lerp(*(q + 2), *(q + 6), t_x);
    float q1_w = lerp(*(q + 3), *(q + 7), t_x);

    // q_new[i, j].x = lerp(q0_x, q1_x, t_y);
    // q_new[i, j].y = lerp(q0_y, q1_y, t_y);
    // q_new[i, j].z = lerp(q0_z, q1_z, t_y);
    // q_new[i, j].w = lerp(q0_w, q1_w, t_y);
    *(q_new + 0) = lerp(q0_x, q1_x, t_y);
    *(q_new + 1) = lerp(q0_y, q1_y, t_y);
    *(q_new + 2) = lerp(q0_z, q1_z, t_y);
    *(q_new + 3) = lerp(q0_w, q1_w, t_y);
  }

  namespace avx {
    void show_m128(__m128 v, const char* prefix) {
      __attribute__((aligned(32))) float t[4] = {0};
      _mm_store_ps(t, v);
      show(t, 4, prefix);
    }

    void show_m256(__m256 v, const char* prefix) {
      __attribute__((aligned(32))) float t[8] = {0};
      _mm256_store_ps(t, v);
      show(t, 8, prefix);
    }

    void show_m256i(__m256i v, const char* prefix) {
      __attribute__((aligned(32))) int t[8] = {0};
      _mm256_store_si256((__m256i*) t, v);
      show(t, 8, prefix);
    }

    inline __m256 load_m256_2f(const float *x) {
      __m256 x0 = _mm256_loadu_ps(x + 0);
      __m256 x1 = _mm256_loadu_ps(x + 8);
      __m256 x2 = _mm256_permute2f128_ps(x1, x0, 0b11);

      x0 = _mm256_shuffle_ps(x0, x1, _MM_SHUFFLE(2, 0, 2, 0));
      x2 = _mm256_shuffle_ps(x2, x2, _MM_SHUFFLE(2, 0, 2, 0));

      return _mm256_blend_ps(x0, x2, 0b00111100);
    }

    inline __m256 marge_m128_2f(__m128 x, __m128 y) {
      __m256 _x = _mm256_castps128_ps256(x);
      __m256 _y = _mm256_castps128_ps256(y);

      __m256 lo = _mm256_unpacklo_ps(_x, _y);
      __m256 hi = _mm256_unpackhi_ps(_x, _y);

      return _mm256_insertf128_ps(lo, _mm256_castps256_ps128(hi), 1);
    }

    inline void store_m256_2f(float *u, __m256 x, __m256 y) {
      __m256 lo = _mm256_unpacklo_ps(x, y);
      __m256 hi = _mm256_unpackhi_ps(x, y);

      __m256 p0 = _mm256_insertf128_ps(lo, _mm256_castps256_ps128(hi), 1);
      __m256 p1 = _mm256_permute2f128_ps(lo, hi, 0x31);

      _mm256_storeu_ps(u + 0, p0);
      _mm256_storeu_ps(u + 8, p1);
    }

    inline __m256 set_m256_2f(float x, float y) {
      return _mm256_set_ps(y, x, y, x, y, x, y, x);
    }

    inline __m256 set_m256_4f(float x, float y, float z, float w) {
      return _mm256_set_ps(w, z, y, x, w, z, y, x);
    }

    inline __m256 clip_m256(__m256 x, __m256 min, __m256 max) {
      return _mm256_min_ps(_mm256_max_ps(x, min), max);
    }

    inline __m256i clip_m256i(__m256i x, __m256i min, __m256i max) {
      return _mm256_min_epi32(_mm256_max_epi32(x, min), max);
    }

    inline __m256 lerp(__m256 a, __m256 b, __m256 t) {
      // y = a + (b - a) * t;
      return _mm256_add_ps(a, _mm256_mul_ps(_mm256_sub_ps(b, a), t));
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

    inline void lerp2d2f(const float *q, float *q_new, __m256i _idx,
                         __m256 t_x, __m256 t_y, int m) {
      // q0_x = lerp(q[q_i + 0, q_j + 0].x, q[q_i + 1, q_j + 0].x, t_x)
      // q0_y = lerp(q[q_i + 0, q_j + 0].y, q[q_i + 1, q_j + 0].y, t_x)
      __m256 q0_x = lerp(
          _mm256_i32gather_ps(q + 0, _idx, sizeof(float)),
          _mm256_i32gather_ps(q + 2, _idx, sizeof(float)), t_x);
      __m256 q0_y = lerp(
          _mm256_i32gather_ps(q + 1, _idx, sizeof(float)),
          _mm256_i32gather_ps(q + 3, _idx, sizeof(float)), t_x);

      q += m * 2;

      // q1_x = lerp(q[q_i + 0, q_j + 1].x, q[q_i + 1, q_j + 1].x, t_x)
      // q1_y = lerp(q[q_i + 0, q_j + 1].y, q[q_i + 1, q_j + 1].y, t_x)
      __m256 q1_x = lerp(
          _mm256_i32gather_ps(q + 0, _idx, sizeof(float)),
          _mm256_i32gather_ps(q + 2, _idx, sizeof(float)), t_x);
      __m256 q1_y = lerp(
          _mm256_i32gather_ps(q + 1, _idx, sizeof(float)),
          _mm256_i32gather_ps(q + 3, _idx, sizeof(float)), t_x);

      // q_new[i, j].x = lerp(q0_x, q1_x, t_y);
      // q_new[i, j].y = lerp(q0_y, q1_y, t_y);
      store_m256_2f(q_new, lerp(q0_x, q1_x, t_y), lerp(q0_y, q1_y, t_y));
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
