#include <math.h>
#include <stdio.h>

#include <iostream>
#include <vector>

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

Vector calculateCordic1(const Vector v0, double target) {
  double sign, xi, yi;
  Vector v = v0;
  long long t = 1;
  for (int i = 0; i < n; i++) {
    sign = (v.theta < target) ? 1.0 : -1.0;
    xi = v.x - v.y * sign / (double)t;
    yi = v.y + v.x * sign / (double)t;
    v.theta += sign * atan_table[i];
    t *= 2;
    v.x = xi;
    v.y = yi;
    // printf("sign = %f, x = %f, y = %f, theta = %f\r\n", sign, v.x, v.y,
    //  v.theta);
  }
  return v;
}

Vector calculateCordic2(const Vector v0, double target) {
  double sign, xi, yi;
  Vector v = v0;
  long long t = 1;
  for (int i = 0; i < n; i++) {
    sign = (v.y < target) ? 1.0 : -1.0;
    xi = v.x - v.y * sign / (double)t;
    yi = v.y + v.x * sign / (double)t;
    v.theta += -sign * atan_table[i];
    t *= 2;
    v.x = xi;
    v.y = yi;
    // printf("sign = %f, x = %f, y = %f, theta = %f\r\n", sign, v.x, v.y,
    //  v.theta);
  }
  return v;
}

Vector calculateCordic3(const Vector v0, double a) {
  double sign, xi, yi;
  Vector v = v0;
  long long t = 2;
  int k = 3;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < 2; j++) {
      sign = (v.theta < a) ? 1.0 : -1.0;
      xi = v.x + v.y * sign / (double)t;
      yi = v.y + v.x * sign / (double)t;
      v.theta += sign * atanh_table[i];
      v.x = xi;
      v.y = yi;
      // printf("sign = %f, x = %f, y = %f, theta = %f\r\n", sign, v.x, v.y,
      //  v.theta);
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

void calculateHypoblicGain() {
  Vector v{.x = 1, .y = 0, .theta = 0};
  auto ret = calculateCordic3(v, 0);
  printf("gain = %f\r\n", 1.0 / ret.x);
  mh = 1 / ret.x;
}

Vector calculateCordic4(const Vector v0, double a) {
  double sign, xi, yi;
  Vector v = v0;
  long long t = 2;
  int k = 3;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < 2; j++) {
      sign = (v.y < a) ? 1.0 : -1.0;
      xi = v.x + v.y * sign / (double)t;
      yi = v.y + v.x * sign / (double)t;
      v.theta += sign * atanh_table[i];
      v.x = xi;
      v.y = yi;
      // printf("sign = %f, x = %f, y = %f, theta = %f\r\n", sign, v.x, v.y,
      //  v.theta);
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

/**
 * @brief sin
 * @note cordic自体はthetaが0 <= theta <= 2 / piしか対応していない
 * そのため、範囲外の値が来た場合は、三角関数の対称性を利用して計算する
 *
 * @param theta 0 <= theta <= 2 / pi
 * @return double
 */
double cordicSin(double theta) {
  double sign = 1;
  // 0 ~ 2PI以外だった場合は、0~2PIにする
  while (theta > 2.0 * M_PI) {
    theta -= 2.0 * M_PI;
  }
  while (theta < 0.0) {
    theta += 2.0 * M_PI;
  }
  // 三角関数の対称性を利用して、他の象限でも計算可能にする
  if (theta > 0.5 * M_PI && theta <= M_PI) {
    theta = M_PI - theta;
  } else if (theta > M_PI && theta <= 1.5 * M_PI) {
    theta = theta - M_PI;
    sign = -1.0;
  } else if (theta > 1.5 * M_PI && theta <= 2.0 * M_PI) {
    theta = 2.0 * M_PI - theta;
    sign = -1.0;
  }
  // cordic
  Vector v{.x = m, .y = 0, .theta = 0};
  auto ret = calculateCordic1(v, theta);
  return sign * ret.y;
}

/**
 * @brief cos
 *
 * @param theta 0 <= theta <= 2 / pi
 * @return double
 */
double cordicCos(double theta) {
  double sign = 1;
  // 0 ~ 2PI以外だった場合は、0~2PIにする
  while (theta > 2.0 * M_PI) {
    theta -= 2.0 * M_PI;
  }
  while (theta < 0.0) {
    theta += 2.0 * M_PI;
  }
  // 三角関数の対称性を利用して、他の象限でも計算可能にする
  if (theta > 0.5 * M_PI && theta <= M_PI) {
    theta = M_PI - theta;
    sign = -1.0;
  } else if (theta > M_PI && theta <= 1.5 * M_PI) {
    theta = theta - M_PI;
    sign = -1.0;
  } else if (theta > 1.5 * M_PI && theta <= 2.0 * M_PI) {
    theta = 2.0 * M_PI - theta;
  }
  // cordic
  Vector v{.x = m, .y = 0, .theta = 0};
  auto ret = calculateCordic1(v, theta);
  return sign * ret.x;
}

/**
 * @brief tan
 *
 * @param theta 0 <= theta <= 2 / pi
 * @return double
 */
double cordicTan(double theta) { return cordicSin(theta) / cordicCos(theta); }

double cordicAtan(double t) {
  Vector v{.x = 1.0, .y = t, .theta = 0.0};
  auto ret = calculateCordic2(v, 0.0);
  return ret.theta;
}

/**
 * @brief
 * @note このアルゴリズムだと、|target| < 0.98の範囲でしか正確に計算できなかった
 * @param target
 * @return double
 */
double cordicAsin(double target) {
  double _sign = 1.0;
  Vector v{.x = m, .y = 0, .theta = 0.0};
  auto ret = calculateCordic2(v, -target);
  return ret.theta;
}

/**
 * @brief
 * @note このアルゴリズムだと、|target| < 0.98の範囲でしか正確に計算できなかった
 * @param target
 * @return double
 */
double cordicAcos(double target) { return M_PI / 2 - cordicAsin(target); }
/**
 * @brief |a| < 1.13くらい
 *
 * @param a
 * @return double
 */
double cordicSinh(double a) {
  // cordic
  Vector v{.x = mh, .y = 0, .theta = 0};
  auto ret = calculateCordic3(v, a);
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
  auto ret = calculateCordic3(v, a);
  return ret.x;
}

double cordicTanh(double a) { return cordicSinh(a) / cordicCosh(a); }

/**
 * @brief |a| < 0.81
 *
 * @param a
 * @return double
 */
double cordicAtanh(double a) {
  Vector v{.x = 1.0, .y = -a, .theta = 0.0};
  auto ret = calculateCordic4(v, 0.0);
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
double cordicExp(double a) { return cordicCosh(a) + cordicSinh(a); }

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

int main(void) {
  n = 18;
  generateAtanTable();
  generateAtanhTable();
  calculateHypoblicGain();
  printf("n = %d, m = %f, mh = %f\r\n", n, m, mh);
#if 0
  for (int i = 0; i <= 360; i += 10) {
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
#if 1
  for (double i = 0; i < 2; i += 0.1) {
    for (double j = -2.5; j < 0; j += 0.1) {
      printf("pow(%f, %f) %f %f\r\n", i, j, pow(i, j), cordicPow(i, j));
    }
  }

#endif
}
