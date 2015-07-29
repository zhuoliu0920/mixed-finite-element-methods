#ifndef TDMESH
#define TDMESH

#include "VorMesh.h"

struct TPoint	
{
	double x,y,z;
};

inline bool operator==(const TPoint& lhs, const TPoint& rhs)
{
    return (fabs(lhs.x-rhs.x)<1e-12 && fabs(lhs.y-rhs.y)<1e-12 && fabs(lhs.z-rhs.z)<1e-12);
}

struct Polygon 
{
	list<TPoint> vertex;
	int index;
	int indexNode1;
	int indexNode2;
};

struct Polyhedron
{
	TPoint center;
	list<TPoint> vertex;
	int index;
	bool onBound;
};

struct Triangle
{
	TPoint originalVertex[4];
	TPoint vertex[3]; // vertex 0, 1, 2
	bool onBound;
	int globalIndex;
	int indexNode1;
	int indexNode2;
	int edgeIndex1[3] = {-1,-1,-1}; // vertex 0 to 1, vertex 1 to 2, vertex 2 to 0
	int edgeIndex2[3] = {-1,-1,-1};
	bool update = false; 
};

struct MeshCell
{
	TPoint center;
	list<int> triangle; //indexes for triangular faces (not on boundary)
	int index;
};
	

class TdMesh {
public:
	TdMesh ();
	~TdMesh ();

	void initialize(int,int,double*,VorPolygon*,list<VorEdge>&,double);
	void split();
	void cleanup();

	void pwcf(double(*) (TPoint), double(*) (TPoint));
//	void rt0();
	void checkError(double(*) (TPoint), double(*) (TPoint), double(*) (TPoint), double(*) (TPoint), double(*) (TPoint,double), double(*) (TPoint,double), double(*) (TPoint,double));

	void outputResults(string,int);
	void outputCell(string);
	void outputFace(string);
	void output(string);
	void printCell();
	void printFace();
	void printMyCell();
	void printMyFace();

//private:
	int numLayer;
	int numCell;
	Polyhedron *tdCell;
	list<Polygon> tdFace;
	MeshCell *myCell;
	list<Triangle> myFace;

	Soln sol;
	Error err;
};

bool operator==(const TPoint& lhs, const TPoint& rhs);
Point perturb(Point,double*,double);
TPoint center(list<TPoint>&);
TPoint* getTrigVertex(TPoint*,int);
Triangle assignTrig(TPoint*,int,int,int);
int trigIndexUpdate(list<Triangle>::iterator,TPoint*,int,int,int);
bool tdCheckOnBound(TPoint*);
bool subset(TPoint*,TPoint*);
bool contain(Triangle,int,int);
TPoint* getVertex(Triangle,int,int);
TPoint getOtherVertex(Triangle,int,int);
TPoint* getVertexInSameFace(Triangle,int,int,int);
bool checkVecDirection(Eigen::Vector3d,TPoint,TPoint);
int getEdgeInSameFace(Triangle,int,int,int);
double tdNumIntTrig(TPoint*,double(*) (TPoint));
double tdNumIntTrig(TPoint,TPoint*,double(*) (TPoint));
double tdNumIntTet(TPoint,TPoint,TPoint,TPoint,double(*) (TPoint));
double tdNumIntTet2(TPoint,TPoint,TPoint,TPoint,double(*) (TPoint,double),double); // calculate int_omega (f(x)-c)
double tdArea(TPoint,TPoint,TPoint);
double tdArea(TPoint*);
double tdVolumn(TPoint*,TPoint);
TPoint midPt(TPoint,TPoint);
TPoint ctrPt(TPoint,TPoint,TPoint);
TPoint min(TPoint,TPoint);
TPoint max(TPoint,TPoint);



#endif
