#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <list>

using namespace std;

int readDomain(char*,double*&);
int readPoints(char*,double*&,double*&,int&);
string convertToString(double);
bool checkIntersection(double, double, double, double, double, double, double, double);
