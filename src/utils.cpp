#include "utils.h"

// Functions to determine data X when it comes with a matrix with the same number of columns as
// the background
void importX(Eigen::MatrixXd x, int nb, int nd,
                    Eigen::VectorXi xI, Eigen::VectorXi xO,
                    std::vector<int> &x_data, Eigen::MatrixXd &zx_data,
                    Eigen::MatrixXd &wx_data)
{
  x_data = std::vector<int>(x.rows());
  Eigen::MatrixXd zX(x.rows(), nb), wX(x.rows(), nd);
  int i, j;
  for (i = 0;i<x.rows();i++)
  {
    for (j = 0; j < (nb); j++)
    {
      if (j == 0) zX(i, 0) = 1.; else zX(i, j) = x(i, xI[j - 1]);
    }
    for (j = 0; j < (nd); j++)
    {
      if (j == 0) wX(i, 0) = 1.; else wX(i, j) = x(i, xO[j - 1]);
    }
  }
  zx_data = zX;
  wx_data = wX;
}
