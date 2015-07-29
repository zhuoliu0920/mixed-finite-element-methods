#include <iostream>
#include "Eigen/Sparse"
int main()
{
  Eigen::SparseMatrix<double> sm(5,5);
  Eigen::SparseMatrix<double> sm2(2,2);
  sm2.insert(0,0) = 1;
  sm2.insert(1,1) = 5.02;
  sm.block(3,2,2,2) = sm2;
  std::cout << sm << std::endl;
}
