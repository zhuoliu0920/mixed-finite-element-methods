#ifndef VORMESH
#define VORMESH

#include <math.h>
#include <string>
#include <list>
#include <iterator>
#include "VoronoiDiagramGenerator.h"
#include "ReadFiles.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"

struct VorEdge 
{
	double length;
	Point start;
	Point end;
	int index;
	int indexNode1;
	int indexNode2;
};

struct VorPolygon
{
	Point center;
	list<Point> vertex;
	int index;
	bool onBound;
};

struct Soln
{
	double *u;
	double *p;
};

struct Error
{
	double u;
	double p;
	double uInt;
};

class VorMesh {
public:
	VorMesh ();
	~VorMesh ();

	void initialize(string,string);
	void cleanup();

	void pwcf(double(*) (double,double), double(*) (double,double));
//	void rt0();
	void checkError(double(*) (double,double), double(*) (double,double), double(*) (double,double), double(*) (double,double,double), double(*) (double,double,double));

	int findEdge(Point,Point);

	void outputResults(string);
	void outputCell(string);
	void outputEdge(string);
	void output(string);
	void printCell();
	void printEdge();

//private:
	int numCell;
	VorPolygon *vCell;
	list<VorEdge> vEdge;

	Soln sol;
	Error err;
};

void insert(list<Point>&, double, double, Point);
double cmp(Point,Point,Point);
double distance(Point,Point);
int checkOnBound(double,double,double*);
double numIntLn(Point,Point,double(*) (double,double));
double numIntLn2(Point,Point,double(*) (double,double,double), double);
//double numIntQuad(Point,Point,Point,Point,double(*) (double,double));
double numIntTri(Point,Point,Point,double(*) (double,double));
double numIntTri2(Point,Point,Point,double(*) (double,double,double),double);
double triArea(Point,Point,Point);
bool compareEdges(const VorEdge&, const VorEdge&); 
bool sameEdges(const VorEdge&, const VorEdge&);
Point midPt(Point,Point);
Point ctrPt(Point,Point,Point);
#endif
