#pragma once
#include "fileReadingFunctions.h"
//stolen from Shantanu's code
std::vector<std::tuple<double, double, double>> pointsFromMshFile(const char* fname) {
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
		//std::cout << xtemp << " " << ytemp << " " << ztemp << std::endl;
		points.push_back(std::tuple<double, double, double>(xtemp, ytemp, ztemp));
	}
	return points;
}
std::vector<std::tuple<double, double, double>> pointsFromTxts(const char* fname) {
	std::vector<std::tuple<double, double, double>> points;
	FILE *file;
	int nv;
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
		fscanf(file, "%lf ", &xtemp); //x co-ordinate
		fscanf(file, "%lf ", &ytemp); //y co-ordinate
		fscanf(file, "%lf ", &ztemp); //z co-ordinate
		points.push_back(std::tuple<double, double, double>(xtemp, ytemp, ztemp));
	}
	return points;
}
std::vector<int> orderFromTxt(const char* fname, int nv) {
	std::vector<int> order;
	FILE *file;
	int itemp;
	char temp[50];
	file = fopen(fname, "r");
	for (int iv = 0; iv < nv; iv++)
	{
		fscanf(file, "%i ", &itemp); //x co-ordinate
	}
	return order;
}
void writeVectorToTxt(std::vector<double> vec, const char* filename)
{
	std::ofstream file;
	file.open(filename);
	//f << fixed << setprecision(2) << endl;
	for (size_t i = 0; i < vec.size(); i++) {
		file << vec[i] << "\n";
	}
	file.close();
}
