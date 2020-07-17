#ifndef GENERAL_COMPUTATION_H
#define GENERAL_COMPUTATION_H
#include <tuple>
#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>

typedef std::tuple<double, double, double> Point;
using std::vector;
double distance(Point refPoint, Point queryPoint);
std::pair<double, double> minMaxCoord(vector<Point> pointIDs, char coord);
std::vector<std::tuple<double, double, double>> shifting_scaling(vector<Point> points, Point evalPoint);
void cuthill_mckee_ordering(vector<vector<int>> &adjacency, vector<int> &order);
void reverse_cuthill_mckee_ordering(vector<vector<int>> &adjacency, vector<int> &order);
#endif