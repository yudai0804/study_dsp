#include <stdio.h>

#include <cmath>
#include <vector>

namespace maclaurin {

/**
 * @brief
 * @param x
 * @return double x!
 */
static double factorial(int x) {
  if (x == 1 || x == 0) {
    return 1;
  } else
    return x * factorial(x - 1);
}

/**
 * @brief
 * @param x
 * @param n
 * @return double x^n
 */
static double power(double x, int n) {
  double y = 1;
  for (int i = 0; i < n; i++) {
    y *= x;
  }
  return y;
}

/**
 * @brief
 *
 * @param n
 * @param k
 * @return double n!/(k!*(n-k)!)
 */
static double combinetion(double n, double k) {
  return factorial(n) / factorial(k) / factorial(n - k);
}

/**
 * @brief
 * @param n
 * @return double
 */
static double bernoulli(int n) {
  if (n == 0) return 1;
  std::vector<double> b(n + 1);
  b[0] = 1;
  for (int i = 1; i <= n; i++) {
    double sum = 0;
    for (int j = 0; j < i; j++) {
      sum += combinetion(i + 1, j) * b[j];
    }
    b[i] = -1 / ((double)i + 1) * sum;
  }
  return b[n];
}

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
  printf("program start\r\n");
  // for (double i = -M_PI / 2; i < M_PI / 2; i += 0.01) {
  //   printf("%f %f %f\r\n", i, std::cos(i), maclaurin::cos(i, 6));
  // }
  for (int i = 0; i < 7; i++) {
    printf("%d %f\r\n", i, maclaurin::bernoulli(i));
  }
  return 0;
}