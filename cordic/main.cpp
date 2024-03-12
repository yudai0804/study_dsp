#include <math.h>
#include <stdio.h>

#include <iostream>
#include <vector>

#define CIRCULAR_MODE 0
#define VECTOR_MODE 1

// 単位はrad
std::vector<double> atan_table;
std::vector<double> atanh_table;
// 参考サイトではn=18 実際にやってみた感じ、n=18で有効数字5~6桁くらい出た。
// ちなみにn=100とかにすると2^i乗が計算できなくなって死ぬので注意
int n;
double m = 1;
double mh = 1;

double degToRad(double deg) { return deg / 180.0 * M_PI; }
double radToDeg(double rad) { return rad * 180.0 / M_PI; }

void generateAtanTable() {
  double j = 1;
  for (int i = 0; i < n; i++) {
    atan_table.push_back(atan(j));
    m *= cos(atan_table[i]);
    j /= 2.0;
  }
}

void generateAtanhTable() {
  // atanは初期値はatan(1)だったが、atanhの場合は初期値はatanh(1/2)とする
  double j = 0.5;
  // create table
  for (int i = 0; i < n; i++) {
    atanh_table.push_back(atanh(j));
    j /= 2.0;
  }
  int k = 3;
  for (int i = 0; i < n; i++) {
    for (int j1 = 0; j1 < 2; j1++) {
      mh *= cosh(atan_table[i]);
      if (k) {
        k--;
        break;
      } else {
        k = 3;
      }
    }
  }
}

struct Vector {
  double x;
  double y;
  double theta;
};

Vector calculateCircularCordic(const Vector v0, double a, int mode) {
  double sign, xi, yi;
  Vector v = v0;
  long long t = 1;
  for (int i = 0; i < n; i++) {
    if (mode) {
      sign = (v.y < a) ? 1.0 : -1.0;
      xi = v.x - v.y * sign / (double)t;
      yi = v.y + v.x * sign / (double)t;
      v.theta += -sign * atan_table[i];
    } else {
      sign = (v.theta < a) ? 1.0 : -1.0;
      xi = v.x - v.y * sign / (double)t;
      yi = v.y + v.x * sign / (double)t;
      v.theta += sign * atan_table[i];
    }
    t *= 2;
    v.x = xi;
    v.y = yi;
  }
  return v;
}

Vector calculateHyperbolicCordic(const Vector v0, double a, int mode) {
  double sign, xi, yi;
  Vector v = v0;
  long long t = 2;
  int k = 3;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < 2; j++) {
      if (mode) {
        sign = (v.y < a) ? 1.0 : -1.0;
        xi = v.x + v.y * sign / (double)t;
        yi = v.y + v.x * sign / (double)t;
        v.theta += sign * atanh_table[i];
      } else {
        sign = (v.theta < a) ? 1.0 : -1.0;
        xi = v.x + v.y * sign / (double)t;
        yi = v.y + v.x * sign / (double)t;
        v.theta += sign * atanh_table[i];
      }
      v.x = xi;
      v.y = yi;
      if (k) {
        k--;
        break;
      } else {
        k = 3;
      }
    }
    t *= 2;
  }
  return v;
}

void calculateHyerbolicGain() {
  Vector v{.x = 1, .y = 0, .theta = 0};
  auto ret = calculateHyperbolicCordic(v, 0, CIRCULAR_MODE);
  mh = 1 / ret.x;
}

/**
 * @brief sin
 *
 * @param a -pi/2 <= a <= pi/2
 * @return double
 */
double cordicSin(double a) {
  Vector v{.x = m, .y = 0, .theta = 0};
  auto ret = calculateCircularCordic(v, a, CIRCULAR_MODE);
  return ret.y;
}

/**
 * @brief cos
 *
 * @param a -pi/2 <= a <= pi/2
 * @return double
 */
double cordicCos(double a) {
  Vector v{.x = m, .y = 0, .theta = 0};
  auto ret = calculateCircularCordic(v, a, CIRCULAR_MODE);
  return ret.x;
}

/**
 * @brief tan
 *
 * @param a -pi/2 <= a <= pi/2
 * @noet |pi/2|に値が近づくと精度が落ちる(発散する)ので注意
 * @return double
 */
double cordicTan(double a) {
  Vector v{.x = m, .y = 0, .theta = 0};
  auto ret = calculateCircularCordic(v, a, CIRCULAR_MODE);
  return ret.y / ret.x;
}

double cordicAtan(double a) {
  Vector v{.x = 1.0, .y = a, .theta = 0.0};
  auto ret = calculateCircularCordic(v, 0.0, VECTOR_MODE);
  return ret.theta;
}

/**
 * @brief
 * @note このアルゴリズムだと、|a| < 0.98の範囲でしか正確に計算できなかった
 * @param a
 * @return double
 */
double cordicAsin(double a) {
  Vector v{.x = m, .y = 0, .theta = 0.0};
  auto ret = calculateCircularCordic(v, -a, VECTOR_MODE);
  return ret.theta;
}

/**
 * @brief
 * @note このアルゴリズムだと、|a| < 0.98の範囲でしか正確に計算できなかった
 * @param a
 * @return double
 */
double cordicAcos(double a) { return M_PI / 2 - cordicAsin(a); }
/**
 * @brief |a| < 1.13くらい
 *
 * @param a
 * @return double
 */
double cordicSinh(double a) {
  // cordic
  Vector v{.x = mh, .y = 0, .theta = 0};
  auto ret = calculateHyperbolicCordic(v, a, CIRCULAR_MODE);
  return ret.y;
}

/**
 * @brief |a| < 1.13くらい
 *
 * @param a
 * @return double
 */
double cordicCosh(double a) {
  // cordic
  Vector v{.x = mh, .y = 0, .theta = 0};
  auto ret = calculateHyperbolicCordic(v, a, CIRCULAR_MODE);
  return ret.x;
}

double cordicTanh(double a) {
  Vector v{.x = mh, .y = 0, .theta = 0};
  auto ret = calculateHyperbolicCordic(v, a, CIRCULAR_MODE);
  return ret.y / ret.x;
}

/**
 * @brief |a| < 0.81
 *
 * @param a
 * @return double
 */
double cordicAtanh(double a) {
  Vector v{.x = 1.0, .y = -a, .theta = 0.0};
  auto ret = calculateHyperbolicCordic(v, 0.0, VECTOR_MODE);
  return ret.theta;
}

/**
 * @brief 1<=a<=9.5
 *
 * @param a
 * @return double
 */
double cordicLog(double a) { return 2 * cordicAtanh((a - 1) / (a + 1)); }

/**
 * @brief |a| < 1.13
 *
 * @param a
 * @return double
 */
double cordicExp(double a) {
  Vector v{.x = mh, .y = 0, .theta = 0};
  auto ret = calculateHyperbolicCordic(v, a, CIRCULAR_MODE);
  return ret.x + ret.y;
}
/**
 * @brief a^b
 * @note a^b = exp(b*log(a))
 * @note 範囲はb*log(a) < 1.13 && 1<=a<=9.5 、bは負の数でも可
 *
 * @param a
 * @param b
 * @return double
 */
double cordicPow(double a, double b) { return cordicExp(b * log(a)); }

/**
 * @brief 0.03 < a < 2.3
 *
 * @param a
 * @return double
 */
double cordicSqrt(double a) {
  Vector v{.x = a + 0.25, .y = a - 0.25, .theta = 0};
  auto ret = calculateHyperbolicCordic(v, 0.0, VECTOR_MODE);
  return ret.x * mh;
}

/**
 * @brief
 * @note 精度が悪かった
 * @param a
 * @return double
 */
double cordicAsinh(double a) { return cordicLog(a + cordicSqrt(a * a + 1)); }

/**
 * @brief
 * @note 精度が悪かった
 *
 * @param a
 * @return double
 */
double cordicAcosh(double a) { return cordicLog(a + cordicSqrt(a * a - 1)); }

void init() {
  generateAtanTable();
  generateAtanhTable();
  calculateHyerbolicGain();
}

int main(void) {
  n = 18;
  init();
  printf("n = %d, m = %f, mh = %f\r\n", n, m, mh);
#if 0
  for (int i = -90; i <= 90; i += 10) {
    double cpp_cos = cos(degToRad(i));
    double cordic_cos = cordicCos(degToRad(i));
    printf(
        "=======\r\ntheta = %d\r\ncpp_cos    = %.15f\r\ncordic_cos = "
        "%.15f\r\n",
        i, cpp_cos, cordic_cos);
  }
#endif
#if 0
  for (double i = -1; i <= 1; i += 0.1) {
    printf("a = %f ", i);
    double cpp_asin = asin(i);
    double cordic_asin = cordicAsin(i);
    double cpp_acos = acos(i);
    double cordic_acos = cordicAcos(i);
    printf("%f %f %f %f\r\n", cpp_asin, cordic_asin, cpp_acos, cordic_acos);
  }
#endif
#if 0
  for (double i = -2.0; i <= 1.5; i += 0.01) {
    double cpp_cosh = cosh(i);
    double cordic_cosh = cordicCosh(i);
    printf(
        "theta = %f\r\ncpp_cosh    = %.15f\r\ncordic_cosh = "
        "%.15f\r\n\r\n",
        i, cpp_cosh, cordic_cosh);
  }
#endif
#if 0
  for (double i = -1; i <= 1; i += 0.01) {
    printf("tan = %f\r\n", i);
    double cpp_atanh = atanh(i);
    double cordic_atanh = cordicAtanh(i);
    printf(
        "cpp_atanh    = %.15f\r\ncordic_atanh = "
        "%.15f\r\n\r\n",
        cpp_atanh, cordic_atanh);
  }
#endif
#if 0
  for (double i = -1; i <= 1; i += 0.1) {
    printf("a = %f ", i);
    double cpp_asin = asin(i);
    double cordic_asin = cordicAsin(i);
    double cpp_acos = acos(i);
    double cordic_acos = cordicAcos(i);
    printf("%f %f %f %f\r\n", cpp_asin, cordic_asin, cpp_acos, cordic_acos);
  }
#endif
#if 0
  for (double i = 1.0; i <= 10; i += 0.01) {
    printf("log %f ", i);
    double a = log(i);
    double b = cordicLog(i);
    printf("%f %f\r\n", a, b);
  }
#endif
#if 0
  for (double i = -1.3; i <= 1.3; i += 0.1) {
    printf("exp %f ", i);
    double a = exp(i);
    double b = cordicExp(i);
    printf("%f %f\r\n", a, b);
  }
#endif
#if 0
  for (double i = 0; i < 2; i += 0.1) {
    for (double j = -2.5; j < 0; j += 0.1) {
      printf("pow(%f, %f) %f %f\r\n", i, j, pow(i, j), cordicPow(i, j));
    }
  }
#endif
#if 0
  for (double i = 0; i < 4; i += 0.1) {
    printf("sqrt(%f), %f, %f\r\n", i, sqrt(i), cordicSqrt(i));
  }
#endif
#if 0
  for (double i = 0; i < 4; i += 0.1) {
    printf("i=%f | %f %f | %f %f\r\n", i, acosh(i), cordicAcosh(i), asinh(i),
           cordicAsinh(i));
  }
#endif
}
