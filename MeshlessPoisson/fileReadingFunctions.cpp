#pragma once
#include <vector>
#include <tuple>
#include <stdio.h>
#include "gridclasses.hpp"

//stolen from Shantanu's code
std::vector<std::tuple<double, double, double>> pointsFromMshFile(const char* fname) {
	double xCoord, yCoord, zCoord;
	std::vector<std::tuple<double, double, double>> points;
	FILE *file;
	int itemp, nv;
	double xtemp, ytemp, ztemp;
	char temp[50];
	file = fopen(fname, "r");
	while (true)
	{
		fscanf(file, "%s ", temp);
		if (strcmp(temp, "$Nodes") == 0)
		{
			break;
		}
	}
	fscanf(file, "%i ", &nv);
	for (int iv = 0; iv < nv; iv++)
	{
		fscanf(file, "%i ", &itemp);  //vertex number
		fscanf(file, "%lf ", &xtemp); //x co-ordinate
		fscanf(file, "%lf ", &ytemp); //y co-ordinate
		fscanf(file, "%lf ", &ztemp); //z co-ordinate
		points.push_back(std::tuple<double, double, double>(xtemp, ytemp, ztemp));
	}
	return points;
}