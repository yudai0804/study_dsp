#include <math.h>
#include <stdio.h>

#include <iostream>
#include <vector>

// 単位はrad
std::vector<double> tan_table;
// 参考サイトではn=18 実際にやってみた感じ、n=18で有効数字5~6桁くらい出た。
// ちなみにn=100とかにすると2^i乗が計算できなくなって死ぬので注意
int n;
double m = 1;

double degToRad(double deg) { return deg / 180.0 * M_PI; }
double radToDeg(double rad) { return rad * 180.0 / M_PI; }

void generateTable() {
  long long j = 1;
  for (int i = 0; i < n; i++) {
    tan_table.push_back(atan(1.0 / (double)j));
    m *= cos(tan_table[i]);
    printf("tan_table[%d] = %f, m = %f\r\n", tan_table[i], m);
    j *= 2;
  }
}
/**
 * @brief
 *
 * @param x0
 * @param y0
 * @param theta 0 <= theta <= 2 / pi
 * @return std::pair<double, double> first = x,second = y
 */
std::pair<double, double> calculateCordic1(const double x0, const double y0,
                                           const double theta) {
  double current_theta = atan(y0 / x0);
  double x = x0;
  double y = y0;
  long long t = 2;
  for (int i = 1; i < n; i++) {
    double sign = (current_theta < theta) ? 1.0 : -1.0;
    double xi = x - y * sign / (double)t;
    double yi = y + x * sign / (double)t;
    current_theta += sign * tan_table[i];
    t *= 2;
    x = xi;
    y = yi;
    // printf("sign = %f, x = %f, y = %f\r\n", sign, x, y);
  }
  return std::make_pair(x, y);
}

/**
 * @brief
 *
 * @param x0
 * @param y0 0 <= theta <= 2 / pi
 * @return std::pair<double, double> first = theta, second = scalar
 */
std::pair<double, double> calculateCordic2(const double x0, const double y0) {
  double theta = 0;
  double x = x0;
  double y = y0;
  long long t = 1;
  for (int i = 0; i < n; i++) {
    double sign = (y < 0) ? 1.0 : -1.0;
    double xi = x - y * sign / (double)t;
    double yi = y + x * sign / (double)t;
    theta += -sign * tan_table[i];
    t *= 2;
    x = xi;
    y = yi;
    printf("sign = %f, x = %f, y = %f, theta = %f\r\n", sign, x, y, theta);
  }
  double scalar = x * m;
  return std::make_pair(theta, scalar);
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
  auto ret = calculateCordic1(1.0, 1.0, theta);
  return sign * m * ret.second;
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
  auto ret = calculateCordic1(1.0, 1.0, theta);
  return sign * m * ret.first;
}

/**
 * @brief tan
 *
 * @param theta 0 <= theta <= 2 / pi
 * @return double
 */
double cordicTan(double theta) { return cordicSin(theta) / cordicCos(theta); }

double cordicAtan(double t) {
  auto ret = calculateCordic2(1, t);
  return ret.first;
}

int main(void) {
  n = 30;
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
  for (double i = 0; i <= 2; i += 0.1) {
    double cpp_atan = atan(i);
    double cordic_atan = cordicAtan(i);
    printf(
        "=======\r\ntan = %f\r\ncpp_atan    = %.15f\r\ncordic_atan = "
        "%.15f\r\n",
        i, cpp_atan, cordic_atan);
  }
#endif
}