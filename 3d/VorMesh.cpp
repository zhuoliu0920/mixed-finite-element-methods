#include "VorMesh.h"

using namespace std;

VorMesh::VorMesh() {}

VorMesh::~VorMesh() {VorMesh::cleanup();}

void VorMesh::initialize(string voronoiNodeFile, string voronoiDomainFile)
{
	char *nodeFile = (char*) voronoiNodeFile.c_str();
	char *domainFile = (char*) voronoiDomainFile.c_str();
	double *domain;
	double *pointsX;
	double *pointsY;
	readDomain(domainFile, domain);
	readPoints(nodeFile, pointsX, pointsY, numCell);

	vCell = new VorPolygon[numCell];
	for (int i=0; i<numCell; i++) {
		vCell[i].center.x = pointsX[i];
		vCell[i].center.y = pointsY[i];
		vCell[i].index = i;
		vCell[i].onBound = false;
	}
	VoronoiDiagramGenerator vdg;
	vdg.generateVoronoi(pointsX,pointsY,numCell,domain[0],domain[1],domain[2],domain[3],0);

	int j = 0;
	double x1,y1,x2,y2,xx1,yy1,xx2,yy2;
	VorEdge tmpln;
	vdg.resetIterator();
	while(vdg.getNext(x1,y1,x2,y2,xx1,yy1,xx2,yy2))	{
		if ( fabs(x1-x2) < 1e-12 && fabs(y1-y2) < 1e-12) {
//			cout << "There is degenerate edge at " << j << endl;
			continue;
		}
		for (int i=0; i<numCell; i++) {
			if (xx1 == vCell[i].center.x && yy1 == vCell[i].center.y)  {
				insert(vCell[i].vertex,x1,y1,vCell[i].center);
				insert(vCell[i].vertex,x2,y2,vCell[i].center);
				tmpln.indexNode1 = i;
				if (checkOnBound(x1,y1,domain) || checkOnBound(x2,y2,domain))
					vCell[i].onBound=true;
			}
			else if (xx2 == vCell[i].center.x && yy2 == vCell[i].center.y) {
				insert(vCell[i].vertex,x1,y1,vCell[i].center);
				insert(vCell[i].vertex,x2,y2,vCell[i].center);
				tmpln.indexNode2 = i;
				if (checkOnBound(x1,y1,domain) || checkOnBound(x2,y2,domain)) 
					vCell[i].onBound=true;
			}
		}
		tmpln.start.x = x1;
		tmpln.start.y = y1;
		tmpln.end.x = x2;
		tmpln.end.y = y2;
		tmpln.length = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		vEdge.push_back(tmpln);
	}
	// remove duplicate edges
	vEdge.sort(compareEdges);
	vEdge.unique(sameEdges);
	// denote index to edges
	j = 0;
	list<VorEdge>::iterator vit;
	for (vit = vEdge.begin(); vit != vEdge.end(); ++vit) {
		(*vit).index = j++;
	}

	// add on-bound edge
	list<Point>::iterator it;
	list<Point>::iterator it2;
	Point *corner = new Point[4];
	corner[0].x = domain[0]; corner[0].y = domain[2];
	corner[1].x = domain[1]; corner[1].y = domain[2];
	corner[2].x = domain[1]; corner[2].y = domain[3];
	corner[3].x = domain[0]; corner[3].y = domain[3];
	double m = distance(corner[0],corner[3]);
	double minDist[4] = {m,m,m,m};
	int cornerIn[4] = {-1,-1,-1,-1};
	for (int i=0; i<numCell; i++) {
		if (vCell[i].onBound == true) {
			for (int j=0; j<4; j++) {
				cornerIn[j] = (distance(corner[j],vCell[i].center)<minDist[j])?i:cornerIn[j];
				minDist[j] = (distance(corner[j],vCell[i].center)<minDist[j])?distance(corner[j],vCell[i].center):minDist[j];
			}
		}
	}
	for (int j=0; j<4; j++) {
		insert(vCell[cornerIn[j]].vertex,corner[j].x,corner[j].y,vCell[cornerIn[j]].center);
	}
	delete[] domain;
	delete[] pointsX;
	delete[] pointsY;
	delete[] corner;
}

void VorMesh::cleanup()
{
	for (int i=0; i<numCell; i++) {
		vCell[i].vertex.clear();
	}
	vEdge.clear();
	delete[] vCell;
//	delete[] sol.u;
//	delete[] sol.p;
}

void VorMesh::pwcf(double(*f) (double,double), double(*g) (double,double))
{
	int pSize = vEdge.size();
	int uSize = 0;
	for (int i=0; i<numCell; i++) {
		uSize += vCell[i].vertex.size();
	}
	sol.u = new double[uSize];
	sol.p = new double[pSize];
	Eigen::SparseMatrix<double> MInv(uSize,uSize);
	Eigen::SparseMatrix<double> B(pSize,uSize);
	Eigen::MatrixXd MInvTmp(uSize,uSize);
//	Eigen::MatrixXd B(pSize,uSize);
	Eigen::VectorXd G(uSize);
	Eigen::VectorXd F(pSize);
	Eigen::VectorXd uTmp(uSize);
	Eigen::VectorXd pTmp(pSize);
	MInv.setZero();
	B.setZero();
	G.setZero();
	F.setZero();

	Eigen::MatrixXd Mtmp(1,1), Ntmp(2,2);
	Eigen::Vector3d eMid(0,0,0),eLef(0,0,0),eRig(0,0,0);
	Eigen::Vector2d nMid,nLef,nRig,wMid1,wMid2,wLef,wRig,b;
	int indexM, indexL, indexR;
	double areaL, areaR;
	list<Point>::iterator it;
	int j=0;
	int globalIndex=0;
	int edgeIndex;
	for (int i=0; i<numCell; i++) {
		// find MInv and B by block calculation
		Mtmp.resize(vCell[i].vertex.size(),vCell[i].vertex.size());
		Mtmp.setZero();
		j=0;
		for (it=vCell[i].vertex.begin(); it!=vCell[i].vertex.end(); ++it) {
			eMid(0) = (*it).x-vCell[i].center.x; eMid(1) = (*it).y-vCell[i].center.y;
			nMid(0) = -eMid(1); nMid(1) = eMid(0); nMid = (1/(sqrt(nMid.dot(nMid))))*nMid;
			if (j==0) {
				eLef(0) = (*prev(vCell[i].vertex.end())).x-vCell[i].center.x; eLef(1) = (*prev(vCell[i].vertex.end())).y-vCell[i].center.y;
				eRig(0) = (*next(it)).x-vCell[i].center.x; eRig(1) = (*next(it)).y-vCell[i].center.y;
				indexL = vCell[i].vertex.size()-1;
				indexR = j+1;
			
				edgeIndex = findEdge(*it,*prev(vCell[i].vertex.end()));
				if (edgeIndex < 0) {
					G(globalIndex+j) += numIntLn(*it,*prev(vCell[i].vertex.end()),g)*distance(vCell[i].center,(*it))/distance(*it,*prev(vCell[i].vertex.end()));
				}
				else {B.insert(edgeIndex,globalIndex+j) = -distance(vCell[i].center,(*it));}

				edgeIndex = findEdge(*it,*next(it));
				if (edgeIndex < 0) {
					G(globalIndex+j) += -numIntLn(*it,*next(it),g)*distance(vCell[i].center,(*it))/distance(*it,*next(it));
				}
				else {B.insert(edgeIndex,globalIndex+j) = distance(vCell[i].center,(*it));}
			}
			else if (j==vCell[i].vertex.size()-1){
				eLef(0) = (*prev(it)).x-vCell[i].center.x; eLef(1) = (*prev(it)).y-vCell[i].center.y;
				eRig(0) = (*(vCell[i].vertex.begin())).x-vCell[i].center.x; eRig(1) = (*(vCell[i].vertex.begin())).y-vCell[i].center.y;
				indexL = j-1;
				indexR = 0;

				edgeIndex = findEdge(*it,*prev(it));
				if (edgeIndex < 0) {
					G(globalIndex+j) += numIntLn(*it,*prev(it),g)*distance(vCell[i].center,(*it))/distance(*it,*prev(it));
				}
				else {B.insert(edgeIndex,globalIndex+j) = -distance(vCell[i].center,(*it));}

				edgeIndex = findEdge(*it,*(vCell[i].vertex.begin()));
				if (edgeIndex < 0) {
					G(globalIndex+j) += -numIntLn(*it,*(vCell[i].vertex.begin()),g)*distance(vCell[i].center,(*it))/distance(*it,*(vCell[i].vertex.begin()));
				}
				else {B.insert(edgeIndex,globalIndex+j) = distance(vCell[i].center,(*it));}
			}
			else {
				eLef(0) = (*prev(it)).x-vCell[i].center.x; eLef(1) = (*prev(it)).y-vCell[i].center.y;
				eRig(0) = (*next(it)).x-vCell[i].center.x; eRig(1) = (*next(it)).y-vCell[i].center.y;
				indexL = j-1;
				indexR = j+1;

				edgeIndex = findEdge(*it,*prev(it));
				if (edgeIndex < 0) {
					G(globalIndex+j) += numIntLn(*it,*prev(it),g)*distance(vCell[i].center,(*it))/distance(*it,*prev(it));
				}
				else {B.insert(edgeIndex,globalIndex+j) = -distance(vCell[i].center,(*it));}

				edgeIndex = findEdge(*it,*next(it));
				if (edgeIndex < 0) {
					G(globalIndex+j) += -numIntLn(*it,*next(it),g)*distance(vCell[i].center,(*it))/distance(*it,*next(it));
				}
				else {B.insert(edgeIndex,globalIndex+j) = distance(vCell[i].center,(*it));}
			}
			nLef(0) = -eLef(1); nLef(1) = eLef(0);
			nRig(0) = -eRig(1); nRig(1) = eRig(0);
			nLef = (1/(sqrt(nLef.dot(nLef))))*nLef;
			nRig = (1/(sqrt(nRig.dot(nRig))))*nRig;

			areaL = sqrt((eMid.cross(eLef)).dot(eMid.cross(eLef)))/2;
			areaR = sqrt((eMid.cross(eRig)).dot(eMid.cross(eRig)))/2;

			Ntmp(0,0) = nMid(0); Ntmp(1,0) = nMid(1); Ntmp(0,1) = nLef(0); Ntmp(1,1) = nLef(1);
			b(0) = 1; b(1) = 0; wMid1 = (Ntmp.transpose()).fullPivLu().solve(b);
			b(0) = 0; b(1) = 1; wLef = (Ntmp.transpose()).fullPivLu().solve(b);
			Ntmp(0,1) = nRig(0); Ntmp(1,1) = nRig(1);
			b(0) = 1; b(1) = 0; wMid2 = (Ntmp.transpose()).fullPivLu().solve(b);
			b(0) = 0; b(1) = 1; wRig = (Ntmp.transpose()).fullPivLu().solve(b);
			Mtmp(j,j) = areaL*(wMid1.dot(wMid1)) + areaR*(wMid2.dot(wMid2));
			Mtmp(j,indexL) = areaL*(wMid1.dot(wLef));
			Mtmp(j,indexR) = areaR*(wMid2.dot(wRig));
			j++;
		}
		MInvTmp.block(globalIndex,globalIndex,vCell[i].vertex.size(),vCell[i].vertex.size()) = Mtmp.inverse();
		globalIndex += j;
	}

	list<VorEdge>::iterator it2;
	j = 0;
	for (it2 = vEdge.begin(); it2 != vEdge.end(); ++it2) {
		F(j) = -(numIntTri(vCell[(*it2).indexNode1].center, (*it2).start, (*it2).end, f) + numIntTri(vCell[(*it2).indexNode2].center, (*it2).start, (*it2).end, f)) ;
		j++;
	}

//	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor> > solver;

	MInv = MInvTmp.sparseView();
solver.analyzePattern(B*MInv*(B.transpose()));
solver.factorize(B*MInv*(B.transpose()));
pTmp = solver.solve(B*MInv*G-F);
//	pTmp = solver.compute(B*MInv*(B.transpose())).solve(B*MInv*G-F);
	uTmp = MInv*(G-(B.transpose())*(pTmp));
	for (int i=0; i<uSize; i++) {
		sol.u[i] = uTmp(i);
	}
	for (int i=0; i<pSize; i++) {
		sol.p[i] = pTmp(i);
	}


/*cout << MInv << endl;
cout << endl;
cout << B << endl;
cout << endl;
cout << G << endl;
cout << endl;
cout << F << endl;
cout << endl;*/

/*cout << endl;
for (int k=0; k<uSize; k++) cout << sol.u[k] << " ";
cout << endl;
cout << endl;*/
//int k;
//cout << "num soln: " << endl;
//for (k=0; k<pSize; k++) cout << sol.p[k] << " ";
//cout << endl;
//for (int k=0; k<uSize; k++) cout << sol.u[k] << " ";
//cout << endl;

}


void VorMesh::checkError(double(*g) (double,double), double(*h1) (double,double), double(*h2) (double,double), double(*hh1) (double,double,double), double(*hh2) (double,double,double))
{
	err.u = 0;
	err.p = 0;
	err.uInt = 0;

	// find error ||p_h - p*||
	list<VorEdge>::iterator it2;
	int j = 0;
	double pNorm = 0;
	double area = 0;
//cout << "real soln: " << endl;
	for (it2 = vEdge.begin(); it2 != vEdge.end(); ++it2) {
		area = triArea((*it2).start,(*it2).end,vCell[(*it2).indexNode1].center) + triArea((*it2).start,(*it2).end,vCell[(*it2).indexNode2].center);
		err.p += area * pow( (numIntLn((*it2).start, (*it2).end, g)/distance((*it2).start, (*it2).end) - (sol.p)[j]), 2 );
//		pNorm += area * pow( (numIntLn((*it2).start, (*it2).end, g)/distance((*it2).start, (*it2).end)), 2 );
//cout << numIntLn((*it2).start, (*it2).end, g)/distance((*it2).start, (*it2).end) << " ";
		j++;
	}
//cout << endl;
	err.p = sqrt(err.p);
//	err.p = sqrt(err.p/pNorm);

	// find error ||u_h - u*||
	list<Point>::iterator it;
	list<Point>::iterator nextit;
	Eigen::MatrixXd Ntmp(2,2);
	Eigen::Vector2d nMid,nRig,utmp,uInt,b,c;
	j = 0;
	int k = 0;
	for (int i=0; i<numCell; i++) {
		k = j;
		for (it = vCell[i].vertex.begin(); it != vCell[i].vertex.end(); ++it) {
			if (j-k == vCell[i].vertex.size()-1) {
				nextit = vCell[i].vertex.begin();
				b(1) = (sol.u)[k];
			}
			else {
				nextit = next(it);
				b(1) = (sol.u)[j+1];
			}
			nMid(1) = (*it).x-vCell[i].center.x; nMid(0) = -((*it).y-vCell[i].center.y);
			nMid = (1/(sqrt(nMid.dot(nMid))))*nMid;
			nRig(1) = (*nextit).x-vCell[i].center.x; nRig(0) = -((*nextit).y-vCell[i].center.y);
			nRig = (1/(sqrt(nRig.dot(nRig))))*nRig;
			Ntmp(0,0) = nMid(0); Ntmp(0,1) = nMid(1); Ntmp(1,0) = nRig(0); Ntmp(1,1) = nRig(1);
			b(0) = (sol.u)[j];

			c(0) = (nMid(0)*numIntLn((*it), vCell[i].center, h1) + nMid(1)*numIntLn((*it), vCell[i].center, h2))/distance((*it), vCell[i].center);
			c(1) = (nRig(0)*numIntLn((*nextit), vCell[i].center, h1) + nRig(1)*numIntLn((*nextit), vCell[i].center, h2))/distance((*nextit), vCell[i].center);
//cout << b(0) << " " << c(0) << endl;
			utmp = (Ntmp).fullPivLu().solve(b);
			uInt = (Ntmp).fullPivLu().solve(c);
			err.u += numIntTri2(*it, *nextit, vCell[i].center, hh1, utmp(0));
			err.u += numIntTri2(*it, *nextit, vCell[i].center, hh2, utmp(1));
			err.uInt += numIntTri2(*it, *nextit, vCell[i].center, hh1, uInt(0));
			err.uInt += numIntTri2(*it, *nextit, vCell[i].center, hh2, uInt(1));
//			uNorm += numIntTri(*it, *nextit, vCell[i].center, h1, 0);
//			uNorm += numIntTri(*it, *nextit, vCell[i].center, h2, 0);
			j++;
		}
	}
	err.u = sqrt(err.u);
	err.uInt = sqrt(err.uInt);
//	err.u = sqrt(err.u/uNorm);
	cout << "Error on p is: " << err.p << endl;
	cout << "Error on u is: " << err.u << endl;
	cout << "Error on u interpolant is: " << err.uInt << endl;
//	cout << "p norm is: " << pNorm << endl;
//	cout << "u norm is: " << uNorm << endl;
	
}
		


int VorMesh::findEdge(Point pt1, Point pt2)
{
	list<VorEdge>::iterator it;
	for (it = vEdge.begin(); it != vEdge.end(); ++it) {
		if ((((*it).start == pt1) && ((*it).end == pt2)) || (((*it).start == pt2) && ((*it).end == pt1))) {
			return (*it).index;
		}
	}
	return -1;
}



void VorMesh::outputResults(string filename)
{
	char *outfile;
	outfile = (char*) filename.c_str();
	ofstream rst;
	rst.open(outfile, std::ofstream::app);
	
	int uDOF = 0;
	for (int i=0; i<numCell; i++) {
		uDOF += vCell[i].vertex.size();
	}

	double h = 0;
	double d1, d2;
	list<VorEdge>::iterator it;
	for (it = vEdge.begin(); it != vEdge.end(); ++it) {
		d1 = (*it).length;
		d2 = distance( (vCell[(*it).indexNode1]).center, (vCell[(*it).indexNode2]).center );
		h = (max(d1,d2) > h)?max(d1,d2):h;
	}

	rst << numCell << " & " << vEdge.size() << " & " << uDOF << " & " << h << " & " << err.p << " & " << err.u << " & " << err.uInt << " & " << (err.u/err.uInt) << " \\\\ \\hline" << endl;
	rst.close();
}




void VorMesh::outputEdge(string filename)
{
	char *outfile;
	outfile = (char*) filename.c_str();
	ofstream outEdge;
	outEdge.open(outfile);

	list<VorEdge>::iterator it;
	for (it = vEdge.begin(); it != vEdge.end(); it++) { 
		outEdge << (*it).start.x << " " << (*it).start.y << " " << (*it).end.x << " " << (*it).end.y << endl;
	}
	outEdge.close();
}

void VorMesh::outputCell(string filename)
{
	char *outfile;
	outfile = (char*) filename.c_str();
	ofstream outCell;
	outCell.open(outfile);

	list<Point>::iterator it;
	for (int i=0; i<numCell; i++) {
		for (it = vCell[i].vertex.begin(); it != vCell[i].vertex.end(); it++) 
		outCell << vCell[i].center.x << " " << vCell[i].center.y << " " << (*it).x << " " << (*it).y << endl;
	}
	outCell.close();
}

void VorMesh::printCell()
{
	list<Point>::iterator it;
	for (int i=0; i<numCell; i++) {
		cout << "Node with index " << vCell[i].index << ": (" << vCell[i].center.x << "," << vCell[i].center.y << ")" << endl;
		cout << "On the boundary? " << vCell[i].onBound << endl;
		cout << "Vertex:" << endl;
		for (it = vCell[i].vertex.begin(); it != vCell[i].vertex.end(); it++) 
			cout << "(" << (*it).x << "," << (*it).y << ")" << endl;
	}
}

void VorMesh::printEdge()
{
	list<VorEdge>::iterator it;
	for (it = vEdge.begin(); it != vEdge.end(); it++) {
		cout << "Edge with index " << (*it).index << ": (" << (*it).start.x << "," << (*it).start.y << ") to (" << (*it).end.x << "," << (*it).end.y << "), length is: " << (*it).length << endl;
		cout << "With two adjacent nodes with index " << (*it).indexNode1 << " and " << (*it).indexNode2 << endl;
	}
}

void insert(list<Point>& pt, double xcor, double ycor, Point center)
{
	Point tmppt;
	tmppt.x = xcor;
	tmppt.y = ycor;
	
	if (pt.size() == 0)
		pt.push_back(tmppt);
	else {
		list<Point>::iterator it;
		for (it = pt.begin(); it != pt.end(); ++it) {
			if (fabs(tmppt.x-(*it).x)+fabs(tmppt.y-(*it).y) < 1e-12) {
				return;
			}
		}
		it = pt.begin();
		while(it != pt.end()) {
			if (cmp(tmppt,*it,center) < 0) {
				pt.insert(it,tmppt);
				return;
			}
			it++;
		}
		pt.push_back(tmppt);
	}
	return;
}

double cmp(Point pt1, Point pt2, Point ctr)
{
	double angle1;
	double angle2;
	angle1 = atan2(pt1.y-ctr.y, pt1.x-ctr.x);
	angle2 = atan2(pt2.y-ctr.y, pt2.x-ctr.x);
	return(angle1-angle2);
}

double distance(Point pt1, Point pt2) {
	return (sqrt((pt1.x-pt2.x)*(pt1.x-pt2.x)+(pt1.y-pt2.y)*(pt1.y-pt2.y)));
}

int checkOnBound(double x, double y, double* domain)
{
	if (x==domain[0] && y!=domain[3]) return 0;
	else if (y==domain[2] && x!=domain[0]) return 1;
	else if (x==domain[1] && y!=domain[2]) return 2;
	else if (y==domain[3] && x!=domain[1]) return 3;
	return -1;	
}

double numIntLn(Point pt1, Point pt2, double(*g) (double,double)) //simpson's rule, accurate for p_3
{
	return (g(pt1.x,pt1.y)+g(pt2.x,pt2.y)+4*g((midPt(pt1,pt2)).x,(midPt(pt1,pt2)).y))*distance(pt1,pt2)/6;
}

double numIntLn2(Point pt1, Point pt2, double(*g) (double,double,double), double c)
{
	return (g(pt1.x,pt1.y,c)+g(pt2.x,pt2.y,c)+4*g((midPt(pt1,pt2)).x,(midPt(pt1,pt2)).y,c))*distance(pt1,pt2)/6;
}

//double numIntQuad(Point pt1, Point pt2, Point pt3, Point pt4, double(*f) (double,double))
//{
//	return ((f(pt1.x,pt1.y)+f(pt2.x,pt2.y)+f(pt3.x,pt3.y))*triArea(pt1,pt2,pt3)/3 + (f(pt4.x,pt4.y)+f(pt2.x,pt2.y)+f(pt3.x,pt3.y))*triArea(pt4,pt2,pt3)/3);
//}

double numIntTri(Point pt1, Point pt2, Point pt3, double(*f) (double,double)) //exact for p_3
{
	return (3*(f(pt1.x,pt1.y)+f(pt2.x,pt2.y)+f(pt3.x,pt3.y)) + 8*(f((midPt(pt1,pt2).x),(midPt(pt1,pt2).y))+f((midPt(pt2,pt3).x),(midPt(pt2,pt3).y))+f((midPt(pt1,pt3).x),(midPt(pt1,pt3).y))) + 27*f((ctrPt(pt1,pt2,pt3)).x,(ctrPt(pt1,pt2,pt3)).y)) * triArea(pt1,pt2,pt3) / 60 ;
}

double numIntTri2(Point pt1, Point pt2, Point pt3, double(*f) (double,double,double), double c)
{
	return (3*(f(pt1.x,pt1.y,c)+f(pt2.x,pt2.y,c)+f(pt3.x,pt3.y,c)) + 8*(f((midPt(pt1,pt2).x),(midPt(pt1,pt2).y),c)+f((midPt(pt2,pt3).x),(midPt(pt2,pt3).y),c)+f((midPt(pt1,pt3).x),(midPt(pt1,pt3).y),c)) + 27*f((ctrPt(pt1,pt2,pt3)).x,(ctrPt(pt1,pt2,pt3)).y,c)) * triArea(pt1,pt2,pt3) / 60 ;
}

double triArea(Point pt1, Point pt2, Point pt3)
{
	Eigen::Vector3d e1(0,0,0),e2(0,0,0);
	e1(0) = pt2.x - pt1.x; e1(1) = pt2.y - pt1.y;
	e2(0) = pt3.x - pt2.x; e2(1) = pt3.y - pt2.y;
	return (sqrt((e1.cross(e2)).dot(e1.cross(e2)))/2);
}

bool compareEdges(const VorEdge& v1, const VorEdge& v2)
{
	VorEdge tmp1, tmp2;
	if (v1.start.x > v1.end.x) {
		tmp1.start = v1.end;
		tmp1.end = v1.start;
	}
	else if (v1.start.y > v1.end.y) {
		tmp1.start = v1.end;
		tmp1.end = v1.start;
	}
	else {
		tmp1.start = v1.start;
		tmp1.end = v1.end;
	}
	
	if (v2.start.x > v2.end.x) {
		tmp2.start = v2.end;
		tmp2.end = v2.start;
	}
	else if (v2.start.y > v2.end.y) {
		tmp2.start = v2.end;
		tmp2.end = v2.start;
	}
	else {
		tmp2.start = v2.start;
		tmp2.end = v2.end;
	}

	if (tmp1.start.x < tmp2.start.x)
		return true;
	else if (tmp1.start.x > tmp2.start.x)
		return false;
	else {
		if (tmp1.start.y < tmp2.start.y)
			return true;
		else if (tmp1.start.y > tmp2.start.y)
			return false;
		else {
			if (tmp1.end.x < tmp2.end.x)
				return true;
			else if (tmp1.end.x > tmp2.end.x)
				return false;
			else {
				if (tmp1.end.y < tmp2.end.y)
					return true;
				else if (tmp1.end.y > tmp2.end.y)
					return false;
			}
		}
	}
	return false;
}
	
bool sameEdges(const VorEdge& v1, const VorEdge& v2)
{
	if ( v1.start == v2.start && v1.end == v2.end ) return true;		
	else if ( v1.start == v2.end && v1.end == v2.start ) return true;
	else return false;
}

Point midPt(Point pt1, Point pt2)
{
	Point md;
	md.x = (pt1.x+pt2.x)/2;
	md.y = (pt1.y+pt2.y)/2;
	return md;
}

Point ctrPt(Point pt1, Point pt2, Point pt3)
{
	Point ctr;
	ctr.x = (pt1.x+pt2.x+pt3.x)/3;
	ctr.y = (pt1.y+pt2.y+pt3.y)/3;
	return ctr;
}
