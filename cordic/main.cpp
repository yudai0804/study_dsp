#include <math.h>
#include <stdio.h>

#include <iostream>
#include <vector>

// 単位はrad
std::vector<double> atan_table;
// 参考サイトではn=18 実際にやってみた感じ、n=18で有効数字5~6桁くらい出た。
// ちなみにn=100とかにすると2^i乗が計算できなくなって死ぬので注意
int n;
double m = 1;

double degToRad(double deg) { return deg / 180.0 * M_PI; }
double radToDeg(double rad) { return rad * 180.0 / M_PI; }

void generateTable() {
  long long j = 1;
  for (int i = 0; i < n; i++) {
    atan_table.push_back(atan(1.0 / (double)j));
    m *= cos(atan_table[i]);
    printf("atan_table[%d] = %f, m = %f\r\n", atan_table[i], m);
    j *= 2;
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
  long long t = 2;
  for (int i = 1; i < n; i++) {
    sign = (v.theta < target) ? 1.0 : -1.0;
    xi = v.x - v.y * sign / (double)t;
    yi = v.y + v.x * sign / (double)t;
    v.theta += sign * atan_table[i];
    t *= 2;
    v.x = xi;
    v.y = yi;
    printf("sign = %f, x = %f, y = %f, theta = %f\r\n", sign, v.x, v.y,
           v.theta);
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
    printf("sign = %f, x = %f, y = %f, theta = %f\r\n", sign, v.x, v.y,
           v.theta);
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
  Vector v{.x = 1.0, .y = 1.0, .theta = M_PI / 4};
  auto ret = calculateCordic1(v, theta);
  return sign * m * ret.y;
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
  Vector v{.x = 1.0, .y = 1.0, .theta = M_PI / 4};
  auto ret = calculateCordic1(v, theta);
  return sign * m * ret.x;
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

double cordicAsin(double t) {
  double sign = -1.0;
  if (t < 0) {
    t = -t;
    sign = 1;
  }
  Vector v{.x = m, .y = 0, .theta = 0.0};
  auto ret = calculateCordic2(v, t);
  return sign * ret.theta;
}
// TODO: acosうまく動かない
double cordicAcos(double t) {}

int main(void) {
  n = 18;
  generateTable();
  printf("n = %d\r\n", n);
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
#if 1
  for (double i = 0; i <= 1; i += 0.1) {
    double cpp_atan = atan(i);
    double cordic_atan = cordicAtan(i);
    printf(
        "=======\r\ntan = %f\r\ncpp_atan    = %.15f\r\ncordic_atan = "
        "%.15f\r\n",
        i, cpp_atan, cordic_atan);
  }
#endif
}