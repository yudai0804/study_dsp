/**
 * @file cordic.h
 * @author Yamaguchi Yudai
 * @brief cordicの関数
 * @version 0.1
 * @date 2024-03-12
 */

#pragma once

namespace cordic {
void init(int _n);
double sin(double a);
double cos(double a);
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
double pow(double a, double b);
double sqrt(double a);
}  // namespace cordic
