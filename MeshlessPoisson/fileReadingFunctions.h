#ifndef FILE_READING_FUNCTIONS_H
#define FILE_READING_FUNCTIONS_H
#include <vector>
#include <tuple>
#include <stdio.h>
#include "gridclasses.hpp"
std::vector<std::tuple<double, double, double>> pointsFromMshFile(const char* fname);
std::vector<std::tuple<double, double, double>> pointsFromTxts(const char* fname);
std::vector<int> orderFromTxt(const char* fname, int nv);
#endif