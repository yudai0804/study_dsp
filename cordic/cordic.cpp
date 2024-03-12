#include <stdio.h>

#include <cmath>
#include <vector>
#define CIRCULAR_MODE 0
#define VECTOR_MODE 1

namespace cordic {

// 単位はrad
std::vector<double> atan_table;
std::vector<double> atanh_table;
// 参考サイトではn=18 実際にやってみた感じ、n=18で有効数字5~6桁くらい出た。
// ちなみにn=100とかにすると2^i乗が計算できなくなって死ぬので注意
int n;
double m = 1;
double mh = 1;

struct Vector {
  double x;
  double y;
  double theta;
};

static double degToRad(double deg) { return deg / 180.0 * M_PI; }
static double radToDeg(double rad) { return rad * 180.0 / M_PI; }

void generateAtanTable() {
  double j = 1;
  for (int i = 0; i < n; i++) {
    atan_table.push_back(std::atan(j));
    m *= std::cos(atan_table[i]);
    j /= 2.0;
  }
}

void generateAtanhTable() {
  // atanは初期値はatan(1)だったが、atanhの場合は初期値はatanh(1/2)とする
  double j = 0.5;
  // create table
  for (int i = 0; i < n; i++) {
    atanh_table.push_back(std::atanh(j));
    j /= 2.0;
  }
  int k = 3;
  for (int i = 0; i < n; i++) {
    for (int j1 = 0; j1 < 2; j1++) {
      mh *= std::cosh(atan_table[i]);
      if (k) {
        k--;
        break;
      } else {
        k = 3;
      }
    }
  }
}

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
double sin(double a) {
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
double cos(double a) {
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
double tan(double a) {
  Vector v{.x = m, .y = 0, .theta = 0};
  auto ret = calculateCircularCordic(v, a, CIRCULAR_MODE);
  return ret.y / ret.x;
}

double atan(double a) {
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
double asin(double a) {
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
double acos(double a) { return M_PI / 2 - asin(a); }
/**
 * @brief |a| < 1.13くらい
 *
 * @param a
 * @return double
 */
double sinh(double a) {
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
double cosh(double a) {
  // cordic
  Vector v{.x = mh, .y = 0, .theta = 0};
  auto ret = calculateHyperbolicCordic(v, a, CIRCULAR_MODE);
  return ret.x;
}

double tanh(double a) {
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
double atanh(double a) {
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
double log(double a) { return 2 * atanh((a - 1) / (a + 1)); }

/**
 * @brief |a| < 1.13
 *
 * @param a
 * @return double
 */
double exp(double a) {
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
double pow(double a, double b) { return exp(b * log(a)); }

/**
 * @brief 0.03 < a < 2.3
 *
 * @param a
 * @return double
 */
double sqrt(double a) {
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
double asinh(double a) { return log(a + sqrt(a * a + 1)); }

/**
 * @brief
 * @note 精度が悪かった
 *
 * @param a
 * @return double
 */
double acosh(double a) { return log(a + sqrt(a * a - 1)); }

void init(int _n) {
  n = _n;
  generateAtanTable();
  generateAtanhTable();
  calculateHyerbolicGain();
}

}  // namespace cordic