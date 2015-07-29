#include <stdio.h>
#include <search.h>
#include <malloc.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "VorMesh.h"
#include "VoronoiDiagramGenerator.h"
#include "CommandLineParser.h"

const double pi = 3.1415926535897;

using namespace std;

double f(double x, double y) // -Laplacian(p) = f in omega
{
	return 0;
}

double g(double x, double y) // p = g on boundary of omega
{
	return x+y;
}

double h1(double x, double y)
{
	return -1;
}

double h2(double x, double y)
{
	return -1;
}

double h1Sq(double x, double y, double c) // (-partial_x(p)(x,y)-c)^2
{
	return (-1-c)*(-1-c);
}

double h2Sq(double x, double y, double c) // (-partial_y(p)(x,y)-c)^2
{
	return (-1-c)*(-1-c);
}

int main(int argc,char **argv) 
{	
	if (argc == 1)
	{
		cerr << "Error: no argument is passed." << endl;
		return -1;
	}

	CommandLineParser clp(argv[1],',');
	char *inputDomainFile = clp.extract("domain");
	char *inputPointsFile = clp.extract("points");
	char *inputResultFile = clp.extract("result");
//	char *outputEdgesFile = clp.extract("edges");
//	char *outputCellsFile = clp.extract("cells");
	string dom(inputDomainFile);
	string pts(inputPointsFile);
	string rst(inputResultFile);
//	string eds(outputEdgesFile);
//	string cls(outputCellsFile);
	if (inputDomainFile == NULL || inputPointsFile == NULL)
	{
		cerr << "Error: input or output files are not specified." << endl;
		return -1;
	}

	VorMesh mymesh;
	mymesh.initialize(pts,dom);
//	mymesh.printEdge();
//	mymesh.printCell();
	mymesh.pwcf(f,g);
	mymesh.checkError(g,h1,h2,h1Sq,h2Sq);
	mymesh.outputResults(rst);
	mymesh.outputCell("cell.txt");
	mymesh.outputEdge("edge.txt");
//	mymesh.cleanup();
		
	return 0;
}
