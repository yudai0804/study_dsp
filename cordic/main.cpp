#include <stdio.h>

#include <cmath>
#include <iostream>
#include <vector>

#include "cordic.h"

std::string file_name;
std::vector<std::string> function_name;
std::vector<double> start;
std::vector<double> end;
std::vector<double> dx;
std::vector<double> output;

int cordic_n = 18;
int main(int argc, char **argv) {
  printf("program start\r\n");
  cordic::init(cordic_n);
  printf("%f\r\n", cordic::cos(0));
  return 0;
}

void test() {
  //  printf("n = %d, m = %f, mh = %f\r\n", n, m, mh);
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