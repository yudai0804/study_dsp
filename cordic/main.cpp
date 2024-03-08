#include <math.h>
#include <stdio.h>

#include <iostream>
#include <vector>

// 単位はrad
std::vector<double> tan_table;
// 参考サイトではn=18
int n;
double m = 1;

double degToRad(double deg) { return deg / 180.0 * M_PI; }
double radToDeg(double rad) { return rad * 180.0 / M_PI; }

void generateTable() {
  int j = 1;
  for (int i = 0; i < n; i++) {
    tan_table.push_back(atan(1.0 / (double)j));
    m *= cos(tan_table[i]);
    // printf("tan_table[%d] = %f, m = %f\r\n", tan_table[i], m);
    j *= 2;
  }
}
/**
 * @brief
 *
 * @param x0
 * @param y0
 * @param theta 0 <= theta <= 2 / pi
 * @return std::pair<double, double>
 */
std::pair<double, double> calculateCordic(const double x0, const double y0,
                                          const double theta) {
  double current_theta = atan(y0 / x0);
  double x = x0;
  double y = y0;
  double t = 2.0;
  for (int i = 1; i < n; i++) {
    double sign = (current_theta < theta) ? 1.0 : -1.0;
    double xi = x - y * sign / t;
    double yi = y + x * sign / t;
    current_theta += sign * tan_table[i];
    t *= 2.0;
    x = xi;
    y = yi;
    // printf("sign = %f, x = %f, y = %f\r\n", sign, x, y);
  }
  return std::make_pair(x, y);
}

/**
 * @brief sin
 *
 * @param theta 0 <= theta <= 2 / pi
 * @return double
 */
double cordicSin(double theta) {
  auto ret = calculateCordic(1.0, 1.0, theta);
  return m * ret.second;
}

/**
 * @brief cos
 *
 * @param theta 0 <= theta <= 2 / pi
 * @return double
 */
double cordicCos(double theta) {
  auto ret = calculateCordic(1.0, 1.0, theta);
  return m * ret.first;
}

/**
 * @brief tan
 *
 * @param theta 0 <= theta <= 2 / pi
 * @return double
 */
double cordicTan(double theta) { return cordicSin(theta) / cordicCos(theta); }

int main(void) {
  n = 18;
  generateTable();
  printf("n = %d\r\n", n);
  for (int i = 0; i <= 90; i += 10) {
    double cpp_cos = cos(degToRad(i));
    double cordic_cos = cordicCos(degToRad(i));
    printf(
        "=======\r\ntheta = %d\r\ncpp_cos    = %.15f\r\ncordic_cos = "
        "%.15f\r\n",
        i, cpp_cos, cordic_cos);
  }
}