#include <stdio.h>
#include <search.h>
#include <malloc.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TdMesh.h"
#include "VorMesh.h"
#include "VoronoiDiagramGenerator.h"
#include "ReadFiles.h"
#include "CommandLineParser.h"

const double pi = 3.1415926535897;

using namespace std;

double f(TPoint p) // -Laplacian(p) = f in omega
{
	return -6;
}

double g(TPoint p) // p = g on boundary of omega
{
	return p.x*p.x+p.y*p.y+p.z*p.z;
}

double h1(TPoint p) // (h1,h2,h3) = -grad(g)
{
	return -2*p.x;
}

double h2(TPoint p) 
{
	return -2*p.y;
}

double h3(TPoint p) 
{
	return -2*p.z;
}

double h1Sq(TPoint p, double c) 
{
	return (-2*p.x-c)*(-2*p.x-c);
}

double h2Sq(TPoint p, double c) 
{
	return (-2*p.y-c)*(-2*p.y-c);
}

double h3Sq(TPoint p, double c) 
{
	return (-2*p.z-c)*(-2*p.z-c);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------


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
	char *numberOfLayers = clp.extract("nlayer");
	char *inputResultFile = clp.extract("result");
	char *perturbPercent = clp.extract("percent");
//	char *outputEdgesFile = clp.extract("edges");
//	char *outputCellsFile = clp.extract("cells");
	string dom(inputDomainFile);
	string pts(inputPointsFile);
	int n = atoi(numberOfLayers);
	double p = atof(perturbPercent);
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

	double domain[4];
	domain[0] = -1; 
	domain[1] = 1; 
	domain[2] = -1; 
	domain[3] = 1; 

	TdMesh mytdmesh;
	mytdmesh.initialize(n-1,mymesh.numCell,domain,mymesh.vCell,mymesh.vEdge,p);
//	mytdmesh.printCell();
//	mytdmesh.printFace();
	mytdmesh.split();
//	mytdmesh.printMyCell();
//	mytdmesh.printMyFace();
//	cout << "number of cells: " << mytdmesh.numCell << endl;
//	cout << "number of faces: " << mytdmesh.myFace.size() << endl;
	mytdmesh.pwcf(f,g);
	mytdmesh.checkError(g,h1,h2,h3,h1Sq,h2Sq,h3Sq);
	mytdmesh.outputResults(rst,n);
		
	return 0;
}
