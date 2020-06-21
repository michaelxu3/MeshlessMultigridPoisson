#pragma once
#include <vector>
#include <tuple>
#include <stdio.h>
#include "gridclasses.hpp"

std::vector<std::tuple<double, double, double>> pointsFromMshFile(const char* fname) {
	FILE* file;
	file = fopen(fname, "r");
	double xCoord, yCoord, zCoord;
	std::vector<std::tuple<double, double, double>> points;
	while (fscanf(file, "%lf %lf %lf ", &xCoord, &yCoord, &zCoord) != EOF) {
		points.push_back(std::tuple<double, double, double>(xCoord, yCoord, zCoord));
	}
	fclose(file);
	return points;
}