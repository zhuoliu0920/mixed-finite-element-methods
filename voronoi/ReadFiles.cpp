#include "ReadFiles.h"

int readDomain(char *inputDomainFile, double *& domain) {
	domain = new double[4];
	ifstream inDomain(inputDomainFile);
	if (!inDomain) {
		cerr << "Error: input file \'" << inputDomainFile << "\' cannot be opened." << endl;
		return -1;
	}
	string temp;
	stringstream stline;
	while (getline(inDomain,temp)) {
		stline.str(temp);
		for (int i = 0; i < 4; i++) {
			getline(stline,temp,' ');
			domain[i] = atof(temp.c_str());
		}
	}
	inDomain.close();
	return 0;
}

int readPoints(char *inputPointsFile, double*& pointsX, double*& pointsY, int& count) {
	ifstream inPoints(inputPointsFile);
	if (!inPoints) {
		cerr << "Error: input file \'" << inputPointsFile << "\' cannot be opened." << endl;
		return -1;
	}
	count = 0;
	list<double> pointsX_ls;
	list<double> pointsY_ls;
	string line;
	string temp;
	size_t pos;
	while (getline(inPoints,line)) {
		pos = line.find(' ');
		temp = line.substr(0, pos);
		pointsX_ls.push_back(atof(temp.c_str()));
		temp = line.substr(pos+1);
		pointsY_ls.push_back(atof(temp.c_str()));
		count++;
	}
	pointsX = new double[count];
	pointsY = new double[count];
	list<double>::iterator itX = pointsX_ls.begin();
	list<double>::iterator itY = pointsY_ls.begin();
	for (int i = 0; i < count; i++) {
		pointsX[i] = *itX++;
		pointsY[i] = *itY++;
	}
	inPoints.close();
	return 0;
}

string convertToString(double number) {
	ostringstream buff;
	buff << number;
	return buff.str();
}


bool checkIntersection(double x1, double y1, double x2, double y2, double xx1, double yy1, double xx2, double yy2) {
//	(x1,y1) and (x2,y2) are endpoints for one line segment. (xx1,yy1) and (xx2,yy2) are endpoints for anthoer line segment.
	if (x1 == x2 && y1 == y2) return true;
	if (xx1 == xx2 && yy1 == yy2) return true;
	double d, s, t;
	d = (xx2-xx1)*(y2-y1) - (x2-x1)*(yy2-yy1);
	s = (x1-xx1)*(y2-y1) - (y1-yy1)*(x2-x1);
	t = (x1-xx1)*(yy2-yy1) - (y1-yy1)*(xx2-xx1);
	s = s/d;
	t = t/d;
	if (s >= 0 && s <= 1 && t >= 0 && t <= 1) return true;
	else return false;
}
