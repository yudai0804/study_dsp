#include <stdio.h>

#include <cmath>

namespace maclaurin {

static double factorial(int x) {
  if (x == 1 || x == 0) {
    return 1;
  } else
    return x * factorial(x - 1);
}

static double power(double x, int n) {
  double y = 1;
  for (int i = 0; i < n; i++) {
    y *= x;
  }
  return y;
}

static double bernoulli() {}

/**
 * @brief sinのn次マクローリン展開
 *
 * @param x
 * @param n n=1,3,5,...
 * @return double
 */
double sin(double x, int n) {
  if (n % 2 == 0) {
    printf("n = %d is error\r\n", n);
    return 0;
  }
  // nを1,3,5,...から0,1,2,...の並びにする
  n = n / 2;
  double y = 0;
  for (int i = 0; i < n; i++) {
    double sign = (i % 2 == 0) ? 1 : -1;
    y += sign / factorial(2 * i + 1) * power(x, 2 * i + 1);
  }
  return y;
}
/**
 * @brief cosのn次マクローリン展開
 *
 * @param x
 * @param n n=0,2,4,...
 * @return double
 */
double cos(double x, int n) {
  if (n % 2 == 1) {
    printf("n = %d is error\r\n", n);
    return 0;
  }
  // nを0,2,4,...から0,1,2,...の並びにする
  n = n / 2;
  double y = 0;
  for (int i = 0; i < n; i++) {
    double sign = (i % 2 == 0) ? 1 : -1;
    y += sign / factorial(2 * i) * power(x, 2 * i);
  }
  return y;
}
double tan(double a);
double asin(double a);
double acos(double a);
double atan(double a);
double sinh(double a);
double cosh(double a);
double tanh(double a);
double asinh(double a);
double acosh(double a);
double atanh(double a);
double log(double a);
double exp(double a);
double sqrt(double a);

}  // namespace maclaurin

int main(void) {
  for (double i = -M_PI / 2; i < M_PI / 2; i += 0.01) {
    printf("%f %f %f\r\n", i, std::cos(i), maclaurin::cos(i, 6));
  }
}