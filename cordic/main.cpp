#include <stdio.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

#include "cordic.h"

double __attribute__((noinline)) fsin(double a) {
  double r;
  asm volatile(
      "movsd %%xmm0, %0\n\t"  // 1.
      "fldl %0\n\t"           // 2.
      "fsin\n\t"              // 3.
      "fstpl %1\n\t"          // 4.
      : "=m"(a)
      : "m"(r));
  return r;
}
double __attribute__((noinline)) fcos(double a) {
  double r;
  asm volatile(
      "movsd %%xmm0, %0\n\t"  // 1.
      "fldl %0\n\t"           // 2.
      "fcos\n\t"              // 3.
      "fstpl %1\n\t"          // 4.
      : "=m"(a)
      : "m"(r));
  return r;
}

std::string data_directory = "data";
std::string file_name;
int cordic_n = 18;
std::string function_name;
// start=endのときはstartの値だけ計算して返す
double start;
double end;
double dx;
std::vector<double> x;
std::vector<double> output;

void calc(std::function<double(double)> func) {
  if (start == end) {
    x.push_back(start);
    output.push_back(func(start));
    return;
  }
  int i;
  if (start == end)
    i = 0;
  else
    i = (end - start) / dx;
  for (int j = 0; j < i + 1; j++) {
    x.push_back(start + dx * j);
    output.push_back(func(start + dx * j));
  }
}
// 引数
// ./main file_name cordic_n function_name start end dx ...
int main(int argc, char **argv) {
  if (argc != 7) {
    printf("argument error\r\n");
    return 1;
  }
  // 変数の初期値を代入
  file_name = std::string(argv[1]);
  cordic_n = atoi(argv[2]);
  if (cordic_n == 0) {
    printf("argument error\r\n");
    return 1;
  }
  function_name = std::string(argv[3]);
  start = atof(argv[4]);
  end = atof(argv[5]);
  dx = atof(argv[6]);

  cordic::init(cordic_n);

  std::string &fn = function_name;
  std::function<double(double)> f;
  // clang-format off
  if(fn == "fsin") f = [](double a) { return fsin(a); };
  else if(fn == "fcos") f = [](double a) { return fcos(a); };
  else if (fn == "sin") f = [](double a) { return std::sin(a); };
  else if (fn == "cos") f = [](double a) { return std::cos(a); };
  else if (fn == "tan") f = [](double a) { return std::tan(a); };
  else if (fn == "asin") f = [](double a) { return std::asin(a); };
  else if (fn == "acos") f = [](double a) { return std::acos(a); };
  else if (fn == "atan") f = [](double a) { return std::atan(a); };
  else if (fn == "sinh") f = [](double a) { return std::sinh(a); };
  else if (fn == "cosh") f = [](double a) { return std::cosh(a); };
  else if (fn == "tanh") f = [](double a) { return std::tanh(a); };
  else if (fn == "asinh") f = [](double a) { return std::asinh(a); };
  else if (fn == "acosh") f = [](double a) { return std::acosh(a); };
  else if (fn == "atanh") f = [](double a) { return std::atanh(a); };
  else if (fn == "log") f = [](double a) { return std::log(a); };
  else if (fn == "exp") f = [](double a) { return std::exp(a); };
  else if (fn == "sqrt") f = [](double a) { return std::sqrt(a); };
  else if (fn == "cordic_sin") f = [](double a) { return cordic::sin(a); };
  else if (fn == "cordic_cos") f = [](double a) { return cordic::cos(a); };
  else if (fn == "cordic_tan") f = [](double a) { return cordic::tan(a); };
  else if (fn == "cordic_asin") f = [](double a) { return cordic::asin(a); };
  else if (fn == "cordic_acos") f = [](double a) { return cordic::acos(a); };
  else if (fn == "cordic_atan") f = [](double a) { return cordic::atan(a); };
  else if (fn == "cordic_sinh") f = [](double a) { return cordic::sinh(a); };
  else if (fn == "cordic_cosh") f = [](double a) { return cordic::cosh(a); };
  else if (fn == "cordic_tanh") f = [](double a) { return cordic::tanh(a); };
  else if (fn == "cordic_asinh") f = [](double a) { return cordic::asinh(a); };
  else if (fn == "cordic_acosh") f = [](double a) { return cordic::acosh(a); };
  else if (fn == "cordic_atanh") f = [](double a) { return cordic::atanh(a); };
  else if (fn == "cordic_log") f = [](double a) { return cordic::log(a); };
  else if (fn == "cordic_exp") f = [](double a) { return cordic::exp(a); };
  else if (fn == "cordic_sqrt") f = [](double a) { return cordic::sqrt(a); };
  // clang-format on
  else {
    printf("argument error\r\n");
    return 1;
  }
  // 計算
  calc(f);
  // ファイルにcsv形式で書き出す
  bool is_dir = std::filesystem::is_directory(data_directory);
  if (!is_dir) std::filesystem::create_directory(data_directory);
  std::ofstream ofs(file_name);
  if (!ofs) {
    printf("open file failed\r\n");
    return 1;
  }
  for (int i = 0; i < x.size(); i++) {
    char c[64];
    sprintf(c, "%.16f,%.16f", x[i], output[i]);
    ofs << c << std::endl;
  }
  return 0;
}
