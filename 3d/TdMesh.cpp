#include "TdMesh.h"

using namespace std;

TdMesh::TdMesh() {}

TdMesh::~TdMesh() {TdMesh::cleanup();}

void TdMesh::initialize(int layer, int numVCell, double *domain, VorPolygon *initialVCell, list<VorEdge> &initialVEdge, double pertbPerct)
{
	if (layer < 1) {
		cerr << "Error: at least 1 layers." << endl;
		return;
	}
	numLayer = layer;
	numCell = numVCell*layer;
	tdCell = new Polyhedron[numCell];

	list<VorEdge> topVEdge = initialVEdge;
	list<VorEdge> bottomVEdge;
	VorPolygon *topVCell = new VorPolygon[numVCell];
	VorPolygon *bottomVCell = new VorPolygon[numVCell];
	for (int k = 0; k < numVCell; k++) {
		topVCell[k] = initialVCell[k];
	}

	TPoint tmptpt;
	Polygon tmpface;
	int ind = 0;
// one layer
	for (int i = 0; i < numLayer; i++) {
		// initialize bottom face
		bottomVEdge = topVEdge;
		for (int k = 0; k < numVCell; k++) {
			bottomVCell[k] = topVCell[k];
		}
		
		// do perturbation on bottom face
		list<VorEdge>::iterator eit1;
		list<VorEdge>::iterator eit2;
		list<VorEdge>::iterator peit2;
		for (peit2 = bottomVEdge.begin(); peit2 != bottomVEdge.end(); ++peit2) {
			(*peit2).index = 0;
		}
		list<VorEdge>::iterator peit1 = bottomVEdge.begin();
		for (eit1 = topVEdge.begin(); eit1 != topVEdge.end(); ++eit1) {
			if ((*peit1).index == -3) {}
			else if ((*peit1).index == -2) {
				(*peit1).start = perturb((*eit1).start, domain, pertbPerct);
				(*peit1).index += -1;
			}
			else if ((*peit1).index == -1) {
				(*peit1).end = perturb((*eit1).end, domain, pertbPerct);
				(*peit1).index += -2;
			}
			else {
				(*peit1).start = perturb((*eit1).start, domain, pertbPerct);
				(*peit1).end = perturb((*eit1).end, domain, pertbPerct);
				(*peit1).index += -3;
			}
	
			peit2 = next(peit1);
			for (eit2 = next(eit1); eit2 != topVEdge.end(); ++eit2) {
				if ((*eit2).start == (*eit1).start && (*peit2).index != -1 && (*peit2).index != -3) {
					(*peit2).start = (*peit1).start;
					(*peit2).index += -1;
				}
				else if ((*eit2).start == (*eit1).end && (*peit2).index != -1 && (*peit2).index != -3) {
					(*peit2).start = (*peit1).end;
					(*peit2).index += -1;
				}
				else if ((*eit2).end == (*eit1).end && (*peit2).index != -2 && (*peit2).index != -3) {
					(*peit2).end = (*peit1).end;
					(*peit2).index += -2;
				}
				else if ((*eit2).end == (*eit1).start && (*peit2).index != -2 && (*peit2).index != -3) {
					(*peit2).end = (*peit1).start;
					(*peit2).index += -2;
				}
				++peit2;
			}
			++peit1;
		}
		eit2 = topVEdge.begin();
		list<Point>::iterator it;
		for (peit2 = bottomVEdge.begin(); peit2 != bottomVEdge.end(); ++peit2) {
			for (it = (bottomVCell[(*peit2).indexNode1]).vertex.begin(); it != (bottomVCell[(*peit2).indexNode1]).vertex.end(); ++it) {
				if ((*it) == (*eit2).start)
					(*it) = (*peit2).start;
				if ((*it) == (*eit2).end)
					(*it) = (*peit2).end;
			}
			for (it = (bottomVCell[(*peit2).indexNode2]).vertex.begin(); it != (bottomVCell[(*peit2).indexNode2]).vertex.end(); ++it) {
				if ((*it) == (*eit2).start)
					(*it) = (*peit2).start;
				if ((*it) == (*eit2).end)
					(*it) = (*peit2).end;
			}
			++eit2;
		}

		// assign information to class members
		for (int j = 0; j < numVCell; j++) {
			for (it = topVCell[j].vertex.begin(); it != topVCell[j].vertex.end(); ++it) {
				tmptpt.x = (*it).x;
				tmptpt.y = (*it).y;
				tmptpt.z = 1.0-double(i)*(2.0/double(numLayer));
				tdCell[numVCell*i+j].vertex.push_back(tmptpt);
			}
			for (it = bottomVCell[j].vertex.begin(); it != bottomVCell[j].vertex.end(); ++it) {
				tmptpt.x = (*it).x;
				tmptpt.y = (*it).y;
				tmptpt.z = 1.0-(double(i)+1.0)*(2.0/double(numLayer));
				tdCell[numVCell*i+j].vertex.push_back(tmptpt);
			}
			tdCell[numVCell*i+j].center = center(tdCell[numVCell*i+j].vertex);
			tdCell[numVCell*i+j].index = numVCell*i+j;
			tdCell[numVCell*i+j].onBound = (topVCell[j].onBound || i == 0 || i == numLayer-1);
		}

		peit2 = bottomVEdge.begin();
		for (eit2 = topVEdge.begin(); eit2 != topVEdge.end(); ++eit2) {
			tmptpt.x = (*eit2).start.x;
			tmptpt.y = (*eit2).start.y;
			tmptpt.z = 1.0-double(i)*(2.0/double(numLayer));
			tmpface.vertex.push_back(tmptpt);

			tmptpt.x = (*eit2).end.x;
			tmptpt.y = (*eit2).end.y;
			tmptpt.z = 1.0-double(i)*(2.0/double(numLayer));
			tmpface.vertex.push_back(tmptpt);

			tmptpt.x = (*peit2).end.x;
			tmptpt.y = (*peit2).end.y;
			tmptpt.z = 1.0-(double(i)+1)*(2.0/double(numLayer));
			tmpface.vertex.push_back(tmptpt);

			tmptpt.x = (*peit2).start.x;
			tmptpt.y = (*peit2).start.y;
			tmptpt.z = 1.0-(double(i)+1)*(2.0/double(numLayer));
			tmpface.vertex.push_back(tmptpt);

			tmpface.index = ind++;
			tmpface.indexNode1 = (*eit2).indexNode1 + numVCell*i;
			tmpface.indexNode2 = (*eit2).indexNode2 + numVCell*i;
			tdFace.push_back(tmpface);
			tmpface.vertex.clear();
			peit2++;
		}
		if (i != 0) {
			for (int j = 0; j < numVCell; j++) {
				for (it = topVCell[j].vertex.begin(); it != topVCell[j].vertex.end(); ++it) {
					tmptpt.x = (*it).x;
					tmptpt.y = (*it).y;
					tmptpt.z = 1.0-double(i)*(2.0/double(numLayer));
					tmpface.vertex.push_back(tmptpt);
				}
				tmpface.index = ind++;
				tmpface.indexNode1 = (numVCell*(i-1))+j; 
				tmpface.indexNode2 = numVCell*i+j;
				tdFace.push_back(tmpface);
				tmpface.vertex.clear();
			}
		}

		topVEdge = bottomVEdge;
		for (int k = 0; k < numVCell; k++) {
			topVCell[k] = bottomVCell[k];
		}
	}
	delete[] topVCell;
	delete[] bottomVCell;
}

void TdMesh::split()
{
	// Assume domain is the cubic [-1,1]^3
	myCell = new MeshCell[numCell];

	int gInd = 0;
	int count;
	list<Polygon>::iterator pit;
	list<TPoint>::iterator it;
	list<Triangle>::iterator fit;
	Triangle tmptrig;

	for (pit = tdFace.begin(); pit != tdFace.end(); ++pit) {
		tmptrig.onBound = false;
		tmptrig.indexNode1 = (*pit).indexNode1;
		tmptrig.indexNode2 = (*pit).indexNode2;
		count = 0;
		for (it = (*pit).vertex.begin(); it != (*pit).vertex.end(); ++it) {
			tmptrig.originalVertex[count++] = (*it);
		}
		tmptrig.globalIndex = gInd++;
		myFace.push_back(tmptrig);
		tmptrig.globalIndex = gInd++;
		myFace.push_back(tmptrig);
	}

	TPoint *tmppt = new TPoint[8];
	TPoint *trigvertex;
	int trigIndex;
	int n;
	bool done;

	for (int i = 0; i < numCell; i++) {
		myCell[i].center = tdCell[i].center;
		myCell[i].index = tdCell[i].index;
		
		done = false;
		for (int j = 0; j < 12; j++) {
			copy(tdCell[i].vertex.begin(),tdCell[i].vertex.end(),tmppt);
			trigvertex = getTrigVertex(tmppt,j);
			if (tdCheckOnBound(trigvertex) == true) {
				myFace.push_back(assignTrig(tmppt,j,i,gInd));
				myCell[i].triangle.push_back(gInd++);
			}
			else {
				n = 0;
				for (fit = myFace.begin(); fit != myFace.end() && n < 2*tdFace.size(); ++fit) {
					if ( subset( trigvertex, (*fit).originalVertex ) == true && subset( trigvertex, (*next(fit)).originalVertex ) == true && (*fit).update == false && done == false) {
					// update 1
						trigIndex = trigIndexUpdate(fit, trigvertex, myCell[i].index, j++, 1);
						myCell[i].triangle.push_back(trigIndex);
						done = true;
						break;
					}
					else if ( subset( trigvertex, (*fit).originalVertex ) == true && subset( trigvertex, (*prev(fit)).originalVertex ) == true && (*fit).update == false && done == true) {
					// update 1
						trigIndex = trigIndexUpdate(fit, trigvertex, myCell[i].index, j++, 1);
						myCell[i].triangle.push_back(trigIndex);
						done = false;
						break;
					}
					else if ( subset( trigvertex, (*fit).originalVertex ) == true && subset( trigvertex, (*next(fit)).originalVertex ) == true && (*fit).update == true && done == false) {
					// update 2
						trigIndex = trigIndexUpdate(fit, trigvertex, myCell[i].index, j++, 2);
						myCell[i].triangle.push_back(trigIndex);
						done = true;
						break;
					}
					else if ( subset( trigvertex, (*fit).originalVertex ) == true && subset( trigvertex, (*prev(fit)).originalVertex ) == true && (*fit).update == true && done == true) {
					// update 2
						trigIndex = trigIndexUpdate(fit, trigvertex, myCell[i].index, j++, 2);
						myCell[i].triangle.push_back(trigIndex);
						done = false;
						break;
					}
					n++;
				}
				j--;
			}
			delete[] trigvertex;
		}
	}

	delete[] tmppt;
}

void TdMesh::cleanup()
{
	tdFace.clear();
	delete[] tdCell;
	myFace.clear();
	delete[] myCell;
	delete[] sol.u;
	delete[] sol.p;
}


void TdMesh::pwcf(double(*f) (TPoint), double(*g) (TPoint))
{
	Eigen::IOFormat CommaInitFmt(2, 0, ", ", "\n", "[", "]");
	int pSize = 2*tdFace.size();
	int uSize = numCell*18;
cout << pSize << " and " << uSize << endl;
	sol.u = new double[uSize];
	sol.p = new double[pSize];
	Eigen::SparseMatrix<double> MInv(uSize,uSize);
	Eigen::SparseMatrix<double> B(pSize,uSize);
	Eigen::MatrixXd MInvTmp(uSize,uSize);
	Eigen::VectorXd G(uSize);
	Eigen::VectorXd F(pSize);
	Eigen::VectorXd uTmp(uSize);
	Eigen::VectorXd pTmp(pSize);
	MInv.setZero();
	B.setZero();
	G.setZero();
	F.setZero();

	Eigen::MatrixXd Mtmp(1,1), Ntmp(3,3);
	Eigen::Vector3d eMida(0,0,0),eLef(0,0,0),eRig(0,0,0),eLef1a(0,0,0),eRig1a(0,0,0),eLef2a(0,0,0),eRig2a(0,0,0),eMidb(0,0,0),eLef1b(0,0,0),eRig1b(0,0,0),eLef2b(0,0,0),eRig2b(0,0,0);
	Eigen::Vector3d nMid,nLef1,nRig1,nLef2,nRig2,wMid1,wMid2,wLef1,wLef2,wRig1,wRig2,b;

	Triangle *faceArray = new Triangle[myFace.size()];
	copy(myFace.begin(), myFace.end(), faceArray);

	int k;
	bool updateOneSide;
	TPoint otherVertex;
	TPoint *edgeVertex = new TPoint[2];
	TPoint *edgeVertexInSameFace1 = new TPoint[2];
	TPoint *edgeVertexInSameFace2 = new TPoint[2];

	int indexM, indexL1, indexL2, indexR1, indexR2;
	double volumnL, volumnR;
	list<int>::iterator it;
	int globalIndex=0;
	int edgeIndex;
	bool out;
 
//Eigen::Vector3d grad;
//grad(0) = -1; grad(1) = 0; grad(2) = 0;
//Eigen::VectorXd uRealLocal(18);
//uRealLocal.setZero();
//Eigen::VectorXd Mu(uSize);
//Eigen::VectorXd Bp(pSize);
//Eigen::VectorXd MuBp(uSize);

	for (int i=0; i<numCell; i++) {
		// find MInv and B by block calculation
		Mtmp.resize(18,18); // for cube cell, there will be 18 interfaces inside
		Mtmp.setZero();
		for (k=0; k<18; k++) { // for cube cell, there will be 18 interfaces inside
			indexM = k;
			updateOneSide = false;
			for (it = myCell[i].triangle.begin(); it != myCell[i].triangle.end(); ++it) {
				if (contain(faceArray[*it],k,i) == true && updateOneSide == false) {
					edgeVertex = getVertex(faceArray[*it],k,i);
					otherVertex = getOtherVertex(faceArray[*it],k,i);

					eMida(0) = edgeVertex[0].x-myCell[i].center.x;
					eMida(1) = edgeVertex[0].y-myCell[i].center.y;
					eMida(2) = edgeVertex[0].z-myCell[i].center.z;
					eMidb(0) = edgeVertex[1].x-myCell[i].center.x;
					eMidb(1) = edgeVertex[1].y-myCell[i].center.y;
					eMidb(2) = edgeVertex[1].z-myCell[i].center.z;
					nMid = eMida.cross(eMidb);
					nMid = (1/(sqrt(nMid.dot(nMid))))*nMid;

					out = checkVecDirection(nMid,otherVertex,myCell[i].center);

//cout << "normal flux direction out? " << out << endl;
//uRealLocal(k) = nMid.dot(grad);

					indexL1 = getEdgeInSameFace(faceArray[*it],k,i,1);
					edgeVertexInSameFace1 = getVertexInSameFace(faceArray[*it],k,i,1);
					eLef1a(0) = edgeVertexInSameFace1[0].x-myCell[i].center.x;
					eLef1a(1) = edgeVertexInSameFace1[0].y-myCell[i].center.y;
					eLef1a(2) = edgeVertexInSameFace1[0].z-myCell[i].center.z;
					eLef1b(0) = edgeVertexInSameFace1[1].x-myCell[i].center.x;
					eLef1b(1) = edgeVertexInSameFace1[1].y-myCell[i].center.y;
					eLef1b(2) = edgeVertexInSameFace1[1].z-myCell[i].center.z;
					nLef1 = eLef1a.cross(eLef1b);
					nLef1 = (1/(sqrt(nLef1.dot(nLef1))))*nLef1;

					indexL2 = getEdgeInSameFace(faceArray[*it],k,i,2);
					edgeVertexInSameFace2 = getVertexInSameFace(faceArray[*it],k,i,2);
					eLef2a(0) = edgeVertexInSameFace2[0].x-myCell[i].center.x;
					eLef2a(1) = edgeVertexInSameFace2[0].y-myCell[i].center.y;
					eLef2a(2) = edgeVertexInSameFace2[0].z-myCell[i].center.z;
					eLef2b(0) = edgeVertexInSameFace2[1].x-myCell[i].center.x;
					eLef2b(1) = edgeVertexInSameFace2[1].y-myCell[i].center.y;
					eLef2b(2) = edgeVertexInSameFace2[1].z-myCell[i].center.z;
					nLef2 = eLef2a.cross(eLef2b);
					nLef2 = (1/(sqrt(nLef2.dot(nLef2))))*nLef2;

					updateOneSide = true;

					if (faceArray[*it].onBound == true) {
						if (out == false)
							G(globalIndex+k) -= tdNumIntTrig(faceArray[*it].vertex,g)*tdArea(myCell[i].center,edgeVertex[0],edgeVertex[1])/tdArea(faceArray[*it].vertex);
						else
							G(globalIndex+k) += tdNumIntTrig(faceArray[*it].vertex,g)*tdArea(myCell[i].center,edgeVertex[0],edgeVertex[1])/tdArea(faceArray[*it].vertex);
					}
					else {
						if (out == false)
							B.insert(*it,globalIndex+k) = tdArea(myCell[i].center,edgeVertex[0],edgeVertex[1]);
						else 
							B.insert(*it,globalIndex+k) = -tdArea(myCell[i].center,edgeVertex[0],edgeVertex[1]);
					}
				}
				else if (contain(faceArray[*it],k,i) == true && updateOneSide == true) {
					indexR1 = getEdgeInSameFace(faceArray[*it],k,i,1);
					edgeVertexInSameFace1 = getVertexInSameFace(faceArray[*it],k,i,1);
					eRig1a(0) = edgeVertexInSameFace1[0].x-myCell[i].center.x;
					eRig1a(1) = edgeVertexInSameFace1[0].y-myCell[i].center.y;
					eRig1a(2) = edgeVertexInSameFace1[0].z-myCell[i].center.z;
					eRig1b(0) = edgeVertexInSameFace1[1].x-myCell[i].center.x;
					eRig1b(1) = edgeVertexInSameFace1[1].y-myCell[i].center.y;
					eRig1b(2) = edgeVertexInSameFace1[1].z-myCell[i].center.z;
					nRig1 = eRig1a.cross(eRig1b);
					nRig1 = (1/(sqrt(nRig1.dot(nRig1))))*nRig1;

					indexR2 = getEdgeInSameFace(faceArray[*it],k,i,2);
					edgeVertexInSameFace2 = getVertexInSameFace(faceArray[*it],k,i,2);
					eRig2a(0) = edgeVertexInSameFace2[0].x-myCell[i].center.x;
					eRig2a(1) = edgeVertexInSameFace2[0].y-myCell[i].center.y;
					eRig2a(2) = edgeVertexInSameFace2[0].z-myCell[i].center.z;
					eRig2b(0) = edgeVertexInSameFace2[1].x-myCell[i].center.x;
					eRig2b(1) = edgeVertexInSameFace2[1].y-myCell[i].center.y;
					eRig2b(2) = edgeVertexInSameFace2[1].z-myCell[i].center.z;
					nRig2 = eRig2a.cross(eRig2b);
					nRig2 = (1/(sqrt(nRig2.dot(nRig2))))*nRig2;

					if (faceArray[*it].onBound == true) {
						if (out == false)
							G(globalIndex+k) += tdNumIntTrig(faceArray[*it].vertex,g)*tdArea(myCell[i].center,edgeVertex[0],edgeVertex[1])/tdArea(faceArray[*it].vertex);
						else
							G(globalIndex+k) -= tdNumIntTrig(faceArray[*it].vertex,g)*tdArea(myCell[i].center,edgeVertex[0],edgeVertex[1])/tdArea(faceArray[*it].vertex);
					}
					else {
						if (out == false)
							B.insert(*it,globalIndex+k) = -tdArea(myCell[i].center,edgeVertex[0],edgeVertex[1]);
						else 
							B.insert(*it,globalIndex+k) = tdArea(myCell[i].center,edgeVertex[0],edgeVertex[1]);
					}
				}
			}

			if (eLef1a == eMida || eLef1a == eMidb)
				eLef = eLef1b;
			else
				eLef = eLef1a;

			if (eRig1a == eMida || eRig1a == eMidb)
				eRig = eRig1b;
			else
				eRig = eRig1a;

			volumnL = fabs(eLef.dot(eMida.cross(eMidb))/6);
			volumnR = fabs(eRig.dot(eMida.cross(eMidb))/6);

			Ntmp(0,0) = nMid(0); Ntmp(1,0) = nMid(1); Ntmp(2,0) = nMid(2);
			Ntmp(0,1) = nLef1(0); Ntmp(1,1) = nLef1(1); Ntmp(2,1) = nLef1(2);
			Ntmp(0,2) = nLef2(0); Ntmp(1,2) = nLef2(1); Ntmp(2,2) = nLef2(2);
			b(0) = 1; b(1) = 0; b(2) = 0; wMid1 = (Ntmp.transpose()).fullPivLu().solve(b);
			b(0) = 0; b(1) = 1; b(2) = 0; wLef1 = (Ntmp.transpose()).fullPivLu().solve(b);
			b(0) = 0; b(1) = 0; b(2) = 1; wLef2 = (Ntmp.transpose()).fullPivLu().solve(b);

			Ntmp(0,0) = nMid(0); Ntmp(1,0) = nMid(1); Ntmp(2,0) = nMid(2);
			Ntmp(0,1) = nRig1(0); Ntmp(1,1) = nRig1(1); Ntmp(2,1) = nRig1(2);
			Ntmp(0,2) = nRig2(0); Ntmp(1,2) = nRig2(1); Ntmp(2,2) = nRig2(2);
			b(0) = 1; b(1) = 0; b(2) = 0; wMid2 = (Ntmp.transpose()).fullPivLu().solve(b);
			b(0) = 0; b(1) = 1; b(2) = 0; wRig1 = (Ntmp.transpose()).fullPivLu().solve(b);
			b(0) = 0; b(1) = 0; b(2) = 1; wRig2 = (Ntmp.transpose()).fullPivLu().solve(b);

			Mtmp(k,k) = volumnL*(wMid1.dot(wMid1)) + volumnR*(wMid2.dot(wMid2));
			Mtmp(k,indexL1) = volumnL*(wMid1.dot(wLef1));
			Mtmp(k,indexL2) = volumnL*(wMid1.dot(wLef2));
			Mtmp(k,indexR1) = volumnR*(wMid2.dot(wRig1));
			Mtmp(k,indexR2) = volumnR*(wMid2.dot(wRig2));
//----------------------------------------------------------------------------------------------------------------------------------------------------
//cout << k << ":\n" << nMid.transpose() << "\n" << nLef1.transpose() << "\n" << nLef2.transpose() << "\n" << nRig1.transpose() << "\n" << nRig2.transpose() << endl;
		}
		MInvTmp.block(globalIndex,globalIndex,18,18) = Mtmp.inverse();
//cout << "at cell " << i << endl;
//cout << (Mtmp.inverse()).format(CommaInitFmt) << endl;
//cout << Mtmp.format(CommaInitFmt) << endl;
//cout << endl;
//cout << uRealLocal.transpose() << endl;
//cout << endl;
//uTmp = Mtmp*uRealLocal;
//cout << "in cell " << i << endl;
//cout << uRealLocal.transpose() << endl;
//cout << endl;
//for (int l = 0; l < 18; l++) {
//	Mu(l+globalIndex) = uTmp(l);
//cout <<	Mu(l+globalIndex) << " ";
//}
//cout << endl;
//----------------------------------------------------------------------------------------------------------------------------------------------------
		globalIndex += 18;
	}
//cout << Mu.transpose() << endl;

	for (int i = 0; i < myFace.size(); i++) {
		if (faceArray[i].onBound == true) break;
		else {
		F(faceArray[i].globalIndex) = -(tdNumIntTet(faceArray[i].vertex[0], faceArray[i].vertex[1], faceArray[i].vertex[2], myCell[faceArray[i].indexNode1].center, f) + tdNumIntTet(faceArray[i].vertex[0], faceArray[i].vertex[1], faceArray[i].vertex[2], myCell[faceArray[i].indexNode2].center, f));
		}
	}

	MInv = MInvTmp.sparseView();
//----------------------------------------------------------------------------------------------------------------------------------------------------
//	list<Triangle>::iterator itt;
//	int j = 0;
//	for (itt = myFace.begin(); itt != myFace.end(); ++itt) {
//		if ((*itt).onBound == true) break;
//		pTmp(j) = tdNumIntTrig((*itt).vertex, g)/tdArea((*itt).vertex);
//		j++;
//	} 

//	Eigen::MatrixXd BT(uSize,pSize);
//	BT = B.transpose();
/*	double tmpSum;
	for (int i = 0; i < uSize; i++) {
		tmpSum = 0;
		for (int j = 0; j < pSize; j++) {
			tmpSum += BT(i,j);
		}
cout << "sum is " << tmpSum << " and right hand side is " << G(i) << endl;
		if (fabs(tmpSum-G(i)) > 1e-12) {
			cout << "at row " << i << ", not equal: " << tmpSum << " != " << G(i) << endl;
		}
	} */
//	Bp = BT*pTmp;
//	MuBp = Mu + Bp;
//int m=0;
//for (int i = 0; i < numCell; i++) {
//	cout << "in cell " << i << endl;
//	for (int m = 0; m<18;m++) {
//		cout << Bp(m+i*18) << " ";
//	}
//	cout << endl;
//}
//----------------------------------------------------------------------------------------------------------------------------------------------------
	

//cout << MInv << endl;
//cout << endl;
//cout << B.transpose() << endl;
//cout << endl;
//cout << G.transpose() << endl;
//cout << endl;
//cout << F.transpose() << endl;
//cout << endl; 

//cout << endl;
//for (int i = 0; i < uSize; i++) {
//	if (fabs(MuBp(i)-G(i)) > 1e-12) {
//		cout << "at row " << i << ", not equal: " << MuBp(i) << " v.s. " << G(i) << endl;
//		cout << Mu(i) << " " << Bp(i) << endl;
//	}
//}

//cout << uSize << endl;
//cout << pSize << endl;

//	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor> > solver;

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


//for (int k=0; k<pSize; k++) cout << sol.p[k] << " ";
//cout << endl;
//cout << endl;
//cout << "num soln: " << endl;
//for (int k=0; k<pSize; k++) cout << sol.p[k] << " ";
//cout << endl;
//for (int k=0; k<uSize; k++) {
//	cout << sol.u[k] << " ";
//	if ((k+1) % 18 == 0) {
//		cout << endl;
//	}
//}
//cout << endl; 


	delete[] faceArray;
	delete[] edgeVertex;
	delete[] edgeVertexInSameFace1;
	delete[] edgeVertexInSameFace2;
}


void TdMesh::checkError(double(*g) (TPoint), double(*h1) (TPoint), double(*h2) (TPoint), double(*h3) (TPoint), double(*hh1) (TPoint,double), double(*hh2) (TPoint,double), double(*hh3) (TPoint,double))
{
	err.p = 0;
	err.u = 0;
	err.uInt = 0;

	// find error ||p_h - p*||
	list<Triangle>::iterator it;
	int j = 0;
	double volumn = 0;
//cout << "real soln of p: " << endl;
	for (it = myFace.begin(); it != myFace.end(); ++it) {
		if ((*it).onBound == true) break;
		volumn = tdVolumn((*it).vertex, myCell[(*it).indexNode1].center) + tdVolumn((*it).vertex, myCell[(*it).indexNode2].center);
//cout << tdVolumn((*it).vertex, myCell[(*it).indexNode1].center) << ", " <<  tdVolumn((*it).vertex, myCell[(*it).indexNode2].center) << endl;
		err.p += volumn * pow( (tdNumIntTrig((*it).vertex, g)/tdArea((*it).vertex) - (sol.p)[j]), 2 );
//cout << "On " << j << " face, volumn is " << volumn << " , real soln is " << tdNumIntTrig((*it).vertex, g)/tdArea((*it).vertex) << ", num soln is " << (sol.p)[j] << endl;
//cout << tdNumIntTrig((*it).vertex, g)/tdArea((*it).vertex) << " ";
		j++;
	}
	err.p = sqrt(err.p);
cout << "error is " << err.p << endl;


	// find error ||u_h - u*||
	Triangle *faceArray = new Triangle[myFace.size()];
	copy(myFace.begin(), myFace.end(), faceArray);

	int vertexIndexArray[12][3] = { {0,1,2}, {2,3,0}, {4,5,6}, {6,7,4}, {0,4,5}, {5,1,0}, {1,5,6}, {6,2,1}, {3,7,6}, {6,2,3}, {0,4,7}, {7,3,0} };
	int edgeIndexArray[12][3] = { {0,1,16}, {2,3,16}, {4,5,17}, {6,7,17}, {11,4,12}, {8,0,12}, {8,5,13}, {9,1,13}, {10,6,14}, {9,2,14}, {11,7,15}, {10,3,15} };
	Eigen::MatrixXd Ntmp(3,3);
	Eigen::Vector3d nMid,eMida,eMidb,b,c,uH,uInt;
	list<int>::iterator it2;
	TPoint *edgeVertex;
	int globalIndex = 0;

	for (int i=0; i<numCell; i++) {
		for (it2 = myCell[i].triangle.begin(); it2 != myCell[i].triangle.end(); ++it2) { // traverse all the tetrahedrons
			j = 0;
			for (int k=0; k<18; k++) {
				if (contain(faceArray[*it2],k,i) == true) {
					edgeVertex = getVertex(faceArray[*it2],k,i);
					eMida(0) = edgeVertex[0].x-myCell[i].center.x;
					eMida(1) = edgeVertex[0].y-myCell[i].center.y;
					eMida(2) = edgeVertex[0].z-myCell[i].center.z;
					eMidb(0) = edgeVertex[1].x-myCell[i].center.x;
					eMidb(1) = edgeVertex[1].y-myCell[i].center.y;
					eMidb(2) = edgeVertex[1].z-myCell[i].center.z;
					nMid = eMida.cross(eMidb);
					nMid = (1/(sqrt(nMid.dot(nMid))))*nMid;
					if (j > 2) {
						cout << "error in checking error for u_h" << endl;
					}
					Ntmp(j,0) = nMid(0); Ntmp(j,1) = nMid(1); Ntmp(j,2) = nMid(2);
					b(j) = (sol.u)[globalIndex+k];
					c(j) = (nMid(0)*tdNumIntTrig(myCell[i].center,edgeVertex,h1) + nMid(1)*tdNumIntTrig(myCell[i].center,edgeVertex,h2) + nMid(2)*tdNumIntTrig(myCell[i].center,edgeVertex,h3)) / tdArea(myCell[i].center,edgeVertex[0],edgeVertex[1]);
//cout << j << endl;
//cout << b(j) << " and " << c(j) << endl;
					j++;
				}
			}
			uH = (Ntmp).fullPivLu().solve(b);
//cout << uH(0) << ", " << uH(1) << ", " << uH(2) << endl;
			uInt = (Ntmp).fullPivLu().solve(c);
			err.u += tdNumIntTet2(myCell[i].center,faceArray[*it2].vertex[0],faceArray[*it2].vertex[1],faceArray[*it2].vertex[2],hh1,uH(0));
			err.u += tdNumIntTet2(myCell[i].center,faceArray[*it2].vertex[0],faceArray[*it2].vertex[1],faceArray[*it2].vertex[2],hh2,uH(1));
			err.u += tdNumIntTet2(myCell[i].center,faceArray[*it2].vertex[0],faceArray[*it2].vertex[1],faceArray[*it2].vertex[2],hh3,uH(2));
			err.uInt += tdNumIntTet2(myCell[i].center,faceArray[*it2].vertex[0],faceArray[*it2].vertex[1],faceArray[*it2].vertex[2],hh1,uInt(0));
			err.uInt += tdNumIntTet2(myCell[i].center,faceArray[*it2].vertex[0],faceArray[*it2].vertex[1],faceArray[*it2].vertex[2],hh2,uInt(1));
			err.uInt += tdNumIntTet2(myCell[i].center,faceArray[*it2].vertex[0],faceArray[*it2].vertex[1],faceArray[*it2].vertex[2],hh3,uInt(2));
		}
		globalIndex += 18;
	}
	err.u = sqrt(err.u);
	err.uInt = sqrt(err.uInt);
	cout << "error on uh is: " << err.u << endl;
	cout << "error on uInt is: " << err.uInt << endl;

	delete[] faceArray;
	delete[] edgeVertex;

}
		
/*
void TdMesh::outputResults(string filename)
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




void TdMesh::outputEdge(string filename)
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

void TdMesh::outputCell(string filename)
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
}*/

void TdMesh::printCell()
{
	list<TPoint>::iterator it;
	for (int i=0; i<numCell; i++) {
		cout << "Node with index " << tdCell[i].index << ": (" << tdCell[i].center.x << "," << tdCell[i].center.y << "," << tdCell[i].center.z << ")" << endl;
		cout << "On the boundary? " << tdCell[i].onBound << endl;
		cout << "Vertex:" << endl;
		for (it = tdCell[i].vertex.begin(); it != tdCell[i].vertex.end(); it++) 
			cout << "(" << (*it).x << "," << (*it).y << "," <<  (*it).z << ")  ";
		cout << endl;
	}
}

void TdMesh::printFace()
{
	list<Polygon>::iterator it;
	list<TPoint>::iterator tit;
	for (it = tdFace.begin(); it != tdFace.end(); it++) {
		cout << "Face with index " << (*it).index << ": " << endl;
		for (tit = (*it).vertex.begin(); tit != (*it).vertex.end(); ++tit) {
			cout << "(" << (*tit).x << "," << (*tit).y << "," << (*tit).z << ")  ";
		}
		cout << endl;
		cout << "With two adjacent nodes with index " << (*it).indexNode1 << " and " << (*it).indexNode2 << endl;
	}
} 

void TdMesh::printMyCell()
{
	list<int>::iterator it;
	for (int i=0; i<numCell; i++) {
		cout << "Node with index " << myCell[i].index << ": (" << myCell[i].center.x << "," << myCell[i].center.y << "," << myCell[i].center.z << ")" << endl;
		cout << "has triangular faces with index:" << endl;
		for (it = myCell[i].triangle.begin(); it != myCell[i].triangle.end(); ++it) 
			cout << *it << ", ";
		cout << endl;
	}
}

void TdMesh::printMyFace()
{
	list<Triangle>::iterator it;
	for (it = myFace.begin(); it != myFace.end(); ++it) {
		cout << "Triangular face with index " << (*it).globalIndex << ": " << endl;
		for (int i = 0; i < 3; i++) {
			cout << "(" << (*it).vertex[i].x << "," << (*it).vertex[i].y << "," << (*it).vertex[i].z << ")  ";
		}
		cout << endl;
		cout << "With edge index on the first cell " << (*it).edgeIndex1[0] << ", " << (*it).edgeIndex1[1] << ", " << (*it).edgeIndex1[2] << endl;
		if ((*it).onBound == false) {
			cout << "With edge index on the second cell " << (*it).edgeIndex2[0] << ", " << (*it).edgeIndex2[1] << ", " << (*it).edgeIndex2[2] << endl;
			cout << "With two adjacent cells with index " << (*it).indexNode1 << " and " << (*it).indexNode2 << endl;
		}
		cout << "On the boundary? " << (*it).onBound << endl;
	}
} 

void TdMesh::outputResults(string filename, int nlayer)
{
	char *outfile;
	outfile = (char*) filename.c_str();
	ofstream rst;
	rst.open(outfile, std::ofstream::app);
	
	int uDOF = numCell*18;
	int pDOF = 2*tdFace.size();
	double tmp1 = (double) nlayer;
	double tmp2 = (double) numCell;
	double h1 = 2/(tmp1-1);
	double h2 = 2/sqrt(tmp2/(tmp1-1));
	double h = sqrt(h1*h1+2*h2*h2);

	rst << nlayer-1 << " & " << numCell << " & " << h << " & " << pDOF << " & " << uDOF << " & " << err.p << " & " << err.u << " & " << err.uInt << " & " << (err.u/err.uInt) << " \\\\ \\hline" << endl;
	rst.close();
}


//.....................................................................................................................................................
// functions used in TdMesh::initialize()
Point perturb(Point pt, double *domain, double percent)
{
	Point result = pt;
	double randX = 0;
	double randY = 0;
	randX = double (rand() % 100001 - 50000); randX = percent*randX/5000000;
	randY = double (rand() % 100001 - 50000); randY = percent*randY/5000000;
	if (pt.x > domain[0] && pt.x < domain[1]) {
		result.x = (pt.x+randX > domain[0] && pt.x+randX < domain[1])?(pt.x+randX):(pt.x);
	}
	if (pt.y > domain[2] && pt.y < domain[3]) {
		result.y = (pt.y+randY > domain[2] && pt.y+randY < domain[3])?(pt.y+randY):(pt.y);
	}
	return result;
}

TPoint center(list<TPoint> &pts)
{
	TPoint ctr;
	int numPts = 0;
	ctr.x = 0; ctr.y = 0; ctr.z = 0;
	list<TPoint>::iterator it;
	for (it = pts.begin(); it != pts.end(); it++) {
		ctr.x += (*it).x;
		ctr.y += (*it).y;
		ctr.z += (*it).z;
		numPts++;
	}
	ctr.x = ctr.x/numPts; ctr.y = ctr.y/numPts; ctr.z = ctr.z/numPts;
	return ctr;
}

//.....................................................................................................................................................
//functions used in TdMesh::split()
TPoint* getTrigVertex(TPoint *vertex, int j)
{
	int vertexIndexArray[12][3] = { {0,1,2}, {2,3,0}, {4,5,6}, {6,7,4}, {0,4,5}, {5,1,0}, {1,5,6}, {6,2,1}, {3,7,6}, {6,2,3}, {0,4,7}, {7,3,0} };
	TPoint *result = new TPoint[3];
	for (int i = 0; i < 3; i++) {
		result[i] = vertex[vertexIndexArray[j][i]];
	}
	return result;
}

Triangle assignTrig(TPoint *vertex, int j, int i, int gInd)
{
	int vertexIndexArray[12][3] = { {0,1,2}, {2,3,0}, {4,5,6}, {6,7,4}, {0,4,5}, {5,1,0}, {1,5,6}, {6,2,1}, {3,7,6}, {6,2,3}, {0,4,7}, {7,3,0} };
	int edgeIndexArray[12][3] = { {0,1,16}, {2,3,16}, {4,5,17}, {6,7,17}, {11,4,12}, {8,0,12}, {8,5,13}, {9,1,13}, {10,6,14}, {9,2,14}, {11,7,15}, {10,3,15} };
	Triangle result;
	result.globalIndex = gInd;
	result.indexNode1 = i;
	for (int i = 0; i < 3; i++) {
		result.vertex[i] = vertex[vertexIndexArray[j][i]];
		result.edgeIndex1[i] = edgeIndexArray[j][i];
	}
	result.onBound = tdCheckOnBound(result.vertex);
	result.update = true;
	return result;
}

int trigIndexUpdate(list<Triangle>::iterator fit, TPoint *trigvertex, int cellIndex, int j, int c)
{ 
	if (c != 1 && c != 2) {
		cerr << "Error for trigIndexUpdate!" << endl;
		return -1;
	}
	int edgeIndexArray[12][3] = { {0,1,16}, {2,3,16}, {4,5,17}, {6,7,17}, {11,4,12}, {8,0,12}, {8,5,13}, {9,1,13}, {10,6,14}, {9,2,14}, {11,7,15}, {10,3,15} };
	for (int i = 0; i < 3; i++) {
		(*fit).vertex[i] = trigvertex[i];
		if (c == 1) { 
			(*fit).edgeIndex1[i] = edgeIndexArray[j][i];
			if ((*fit).indexNode1 != cellIndex) {
				(*fit).indexNode2 = (*fit).indexNode1;
				(*fit).indexNode1 = cellIndex;
			}
		}
		else if (c == 2)
			(*fit).edgeIndex2[i] = edgeIndexArray[j][i];
	}
	(*fit).update = true;
	return ( (*fit).globalIndex );
}

bool tdCheckOnBound(TPoint *pts)
{
	if ((pts[0].x == -1 && pts[1].x == -1 && pts[2].x == -1) || (pts[0].x == 1 && pts[1].x == 1 && pts[2].x == 1))
		return true;
	else if ((pts[0].y == -1 && pts[1].y == -1 && pts[2].y == -1) || (pts[0].y == 1 && pts[1].y == 1 && pts[2].y == 1))
		return true;
	else if ((pts[0].z == -1 && pts[1].z == -1 && pts[2].z == -1) || (pts[0].z == 1 && pts[1].z == 1 && pts[2].z == 1))
		return true;
	else
		return false;
}

bool subset(TPoint *set1, TPoint *set2)
{
	bool contain = false;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			if (set1[i] == set2[j])
				contain = true;
		}
		if (contain == false)
			return false;
		else
			contain = false;
	}
	return true;
}

//.....................................................................................................................................................
//functions used in TdMesh::pwcf()
bool contain(Triangle face, int edgeIndex, int cellIndex) 
{
	if (cellIndex == face.indexNode1) {
		for (int i = 0; i < 3; i++) {
			if (face.edgeIndex1[i] == edgeIndex)
				return true;
		}
	}
	else if (cellIndex == face.indexNode2) {
		for (int i = 0; i < 3; i++) {
			if (face.edgeIndex2[i] == edgeIndex)
				return true;
		}
	}
	return false;
}
	

TPoint* getVertex(Triangle face, int edgeIndex, int cellIndex) 
{
	int array[3][2] = { {0,1}, {1,2}, {2,0} };
	TPoint *result = new TPoint[2];
	if (cellIndex == face.indexNode1) {
		for (int i = 0; i < 3; i++) {
			if (face.edgeIndex1[i] == edgeIndex) {
				result[0] = min(face.vertex[array[i][0]], face.vertex[array[i][1]]);
				result[1] = max(face.vertex[array[i][0]], face.vertex[array[i][1]]);
			}
		}
	}
	else if (cellIndex == face.indexNode2) {
		for (int i = 0; i < 3; i++) {
			if (face.edgeIndex2[i] == edgeIndex) {
				result[0] = min(face.vertex[array[i][0]], face.vertex[array[i][1]]);
				result[1] = max(face.vertex[array[i][0]], face.vertex[array[i][1]]);
			}
		}
	}
	return result;
}
	
TPoint getOtherVertex(Triangle face, int edgeIndex, int cellIndex)
{
	int array[3] = { 2, 0, 1 };
	TPoint result;
	if (cellIndex == face.indexNode1) {
		for (int i = 0; i < 3; i++) {
			if (face.edgeIndex1[i] == edgeIndex) {
				result = face.vertex[array[i]];
			}
		}
	}
	else if (cellIndex == face.indexNode2) {
		for (int i = 0; i < 3; i++) {
			if (face.edgeIndex2[i] == edgeIndex) {
				result = face.vertex[array[i]];
			}
		}
	}
	return result;
}

TPoint* getVertexInSameFace(Triangle face, int edgeIndex, int cellIndex, int whichOne) 
{
	int array[3][2] = { {0,1}, {1,2}, {2,0} };
	int otherEdgeIndex[2];
	int k = 0;
	TPoint *result = new TPoint[2];
	if (cellIndex == face.indexNode1) {
		for (int i = 0; i < 3; i++) {
			if (face.edgeIndex1[i] != edgeIndex) {
				otherEdgeIndex[k++] = i;
			}
		}
	}
	else if (cellIndex == face.indexNode2) {
		for (int i = 0; i < 3; i++) {
			if (face.edgeIndex2[i] != edgeIndex) {
				otherEdgeIndex[k++] = i;
			}
		}
	}

	if (whichOne == 1) {
		result[0] = min(face.vertex[array[otherEdgeIndex[0]][0]], face.vertex[array[otherEdgeIndex[0]][1]]);
		result[1] = max(face.vertex[array[otherEdgeIndex[0]][0]], face.vertex[array[otherEdgeIndex[0]][1]]);
	}
	else if (whichOne == 2) {
		result[0] = min(face.vertex[array[otherEdgeIndex[1]][0]], face.vertex[array[otherEdgeIndex[1]][1]]);
		result[1] = max(face.vertex[array[otherEdgeIndex[1]][0]], face.vertex[array[otherEdgeIndex[1]][1]]);
	}
	return result;
}
bool checkVecDirection(Eigen::Vector3d n, TPoint top, TPoint base)
{
	Eigen::Vector3d tmp;
	tmp(0) = top.x - base.x;
	tmp(1) = top.y - base.y;
	tmp(2) = top.z - base.z;
	if (n.dot(tmp) > 0) return false; // normal vector goes out of tetrahedron
	else return true;
}

int getEdgeInSameFace(Triangle face, int edgeIndex, int cellIndex, int whichOne) 
{
	int otherEdgeIndex[2];
	int k = 0;
	int result;
	if (cellIndex == face.indexNode1) {
		for (int i = 0; i < 3; i++) {
			if (face.edgeIndex1[i] != edgeIndex)
				otherEdgeIndex[k++] = i;
		}
		if (whichOne == 1) {
			result = face.edgeIndex1[otherEdgeIndex[0]];
		}
		else if (whichOne == 2) {
			result = face.edgeIndex1[otherEdgeIndex[1]];
		}
	}
	else if (cellIndex == face.indexNode2) {
		for (int i = 0; i < 3; i++) {
			if (face.edgeIndex2[i] != edgeIndex)
				otherEdgeIndex[k++] = i;
		}
		if (whichOne == 1) {
			result = face.edgeIndex2[otherEdgeIndex[0]];
		}
		else if (whichOne == 2) {
			result = face.edgeIndex2[otherEdgeIndex[1]];
		}
	}

	return result;
}

//.....................................................................................................................................................
//functions used in TdMesh::checkError()
double tdNumIntTet2(TPoint p1, TPoint p2, TPoint p3, TPoint p4, double(*f) (TPoint,double), double c) // calculate int_omega (f(x)-c)
{
	Eigen::Vector3d e1,e2,e3;
	e1(0) = p2.x - p1.x; e1(1) = p2.y - p1.y; e1(2) = p2.z - p1.z;
	e2(0) = p3.x - p1.x; e2(1) = p3.y - p1.y; e2(2) = p3.z - p1.z;
	e3(0) = p4.x - p1.x; e3(1) = p4.y - p1.y; e3(2) = p4.z - p1.z;
	double volumn = fabs(e1.dot(e2.cross(e3))/6);
	return ( (f(p1,c) + f(p2,c) + f(p3,c) + f(p4,c) + 16*f(midPt(midPt(p1,p2),midPt(p3,p4)),c)) * volumn / 20 );
}

//.....................................................................................................................................................
//functions that are generally used
double tdNumIntTrig(TPoint *vertex, double(*f) (TPoint))
{
	return ( (3*(f(vertex[0])+f(vertex[1])+f(vertex[2])) + 8*(f(midPt(vertex[0],vertex[1]))+f(midPt(vertex[0],vertex[2]))+f(midPt(vertex[1],vertex[2]))) + 27*f(ctrPt(vertex[0],vertex[1],vertex[2]))) * tdArea(vertex[0],vertex[1],vertex[2]) / 60 );
}

double tdNumIntTrig(TPoint top, TPoint *base, double(*f) (TPoint))
{
	return ( (3*(f(base[0])+f(base[1])+f(top)) + 8*(f(midPt(base[0],base[1]))+f(midPt(base[0],top))+f(midPt(base[1],top))) + 27*f(ctrPt(base[0],base[1],top))) * tdArea(base[0],base[1],top) / 60 );
}

double tdNumIntTet(TPoint p1, TPoint p2, TPoint p3, TPoint p4, double(*f) (TPoint))
{
	Eigen::Vector3d e1,e2,e3;
	e1(0) = p2.x - p1.x; e1(1) = p2.y - p1.y; e1(2) = p2.z - p1.z;
	e2(0) = p3.x - p1.x; e2(1) = p3.y - p1.y; e2(2) = p3.z - p1.z;
	e3(0) = p4.x - p1.x; e3(1) = p4.y - p1.y; e3(2) = p4.z - p1.z;
	double volumn = fabs(e1.dot(e2.cross(e3))/6);
	return ( (f(p1) + f(p2) + f(p3) + f(p4) + 16*f(midPt(midPt(p1,p2),midPt(p3,p4)))) * volumn / 20 );
}

double tdArea(TPoint p1, TPoint p2, TPoint p3)
{
	Eigen::Vector3d e1,e2;
	e1(0) = p2.x - p1.x; e1(1) = p2.y - p1.y; e1(2) = p2.z - p1.z;
	e2(0) = p3.x - p1.x; e2(1) = p3.y - p1.y; e2(2) = p3.z - p1.z;
	return ( sqrt((e1.cross(e2)).dot(e1.cross(e2)))/2 );
}

double tdArea(TPoint *pts)
{
	Eigen::Vector3d e1,e2;
	e1(0) = pts[1].x - pts[0].x; e1(1) = pts[1].y - pts[0].y; e1(2) = pts[1].z - pts[0].z;
	e2(0) = pts[2].x - pts[0].x; e2(1) = pts[2].y - pts[0].y; e2(2) = pts[2].z - pts[0].z;
	return ( sqrt((e1.cross(e2)).dot(e1.cross(e2)))/2 );
}

double tdVolumn(TPoint *base, TPoint top)
{
	Eigen::Vector3d e1,e2,e3;
	e1(0) = base[0].x - top.x; e1(1) = base[0].y - top.y; e1(2) = base[0].z - top.z;
	e2(0) = base[1].x - top.x; e2(1) = base[1].y - top.y; e2(2) = base[1].z - top.z;
	e3(0) = base[2].x - top.x; e3(1) = base[2].y - top.y; e3(2) = base[2].z - top.z;
	return ( fabs(e1.dot(e2.cross(e3))/6) );
}

TPoint midPt(TPoint p1, TPoint p2)
{
	TPoint rst;
	rst.x = (p1.x + p2.x) / 2;
	rst.y = (p1.y + p2.y) / 2;
	rst.z = (p1.z + p2.z) / 2;
	return rst;
}

TPoint ctrPt(TPoint p1, TPoint p2, TPoint p3)
{
	TPoint rst;
	rst.x = (p1.x + p2.x + p3.x) / 3;
	rst.y = (p1.y + p2.y + p3.y) / 3;
	rst.z = (p1.z + p2.z + p3.z) / 3;
	return rst;
}

TPoint min(TPoint p1, TPoint p2)
{
	if (p1.x < p2.x)
		return p1;
	else if (p1.x > p2.x)
		return p2;
	else {
		if (p1.y < p2.y)
			return p1;
		else if (p1.y > p2.y)
			return p2;
		else {
			if (p1.z < p2.z)
				return p1;
			else if (p1.z >= p2.z)
				return p2;
		}
	}
}
		 
TPoint max(TPoint p1, TPoint p2)
{
	if (p1.x < p2.x)
		return p2;
	else if (p1.x > p2.x)
		return p1;
	else {
		if (p1.y < p2.y)
			return p2;
		else if (p1.y > p2.y)
			return p1;
		else {
			if (p1.z < p2.z)
				return p2;
			else if (p1.z >= p2.z)
				return p1;
		}
	}
}
		 
