/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 23
 */

#ifndef FLUID_MATH_UTIL_H
#define FLUID_MATH_UTIL_H

namespace fluid::util {
  template<class Process>
  void timer(Process p, const char* prefix, int loop = 1, bool enable = true) {
    using namespace std;

    clock_t bgn = clock();
    for (int i = 0; i < loop; i++) p();
    clock_t end = clock();
    double t = (double) (end - bgn) / CLOCKS_PER_SEC * 1000;

    if (enable) cout << "[debug] (" << prefix << ") time: " << t << "ms" << endl;
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

    inline __m256 load_m256_2f(float *x) {
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
      x = _mm256_min_ps(x, _exp_hi);
      x = _mm256_max_ps(x, _exp_lo);

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
