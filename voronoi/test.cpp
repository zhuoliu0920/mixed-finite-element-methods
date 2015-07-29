#include "ReadFiles.h"
#include "Eigen/Sparse"

using namespace std;

int main() {
	SparseMatrix<double> A(10,10);
	float *domain;
	readDomain("domain.txt", domain);
	for (int i = 0; i < 4; i++) {
		cout << domain[i] << endl;
	}
	float *pointsX;
	float *pointsY;
	int count;
	readPoints("points.txt", pointsX, pointsY, count);
	for (int i = 0; i < count; i++) {
		cout << pointsX[i] << " and " << pointsY[i] << endl;
	} 
	return 0;
}
