#include <stdio.h>
#include <search.h>
#include <malloc.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "VoronoiDiagramGenerator.h"
#include "ReadFiles.h"
#include "CommandLineParser.h"

using namespace std;

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
	char *outputEdgesFile = clp.extract("edges");
	if (inputDomainFile == NULL || inputPointsFile == NULL || outputEdgesFile == NULL)
	{
		cerr << "Error: input or output files are not specified." << endl;
		return -1;
	}
	double *domain;
	double *pointsX;
	double *pointsY;
	int count;
	readDomain(inputDomainFile, domain);
	readPoints(inputPointsFile, pointsX, pointsY, count);



	VoronoiDiagramGenerator vdg;
	vdg.generateVoronoi(pointsX,pointsY,count,domain[0],domain[1],domain[2],domain[3],0);

	vdg.resetIterator();

	ofstream outEdges;
	outEdges.open(outputEdgesFile);
	if (!outEdges) {
		cerr << "Error: output file \'" << outputEdgesFile << "\' cannot be created." << endl;
		return -1;
	}
	// (x1,y1) and (x2,y2) are the two end points for each edge. (xx1,yy1) and (xx2,yy2) are two sites whose line segement is perpendicular to edge
	double x1,y1,x2,y2,xx1,yy1,xx2,yy2;
	while(vdg.getNext(x1,y1,x2,y2,xx1,yy1,xx2,yy2))
	{
		outEdges << checkIntersection(x1,y1,x2,y2,xx1,yy1,xx2,yy2) << " " << convertToString(x1) << " " << convertToString(y1) << " " << convertToString(x2) << " " << convertToString(y2) << " " << convertToString(xx1) << " " << convertToString(yy1) << " " << convertToString(xx2) << " " << convertToString(yy2) << endl;
	}
	outEdges.close();
	return 0;
}
