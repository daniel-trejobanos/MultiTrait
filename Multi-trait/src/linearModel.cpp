#include "linearModel.h"
#include <Eigen/Core>

void LinearModel::simulate() {
  scalar_normal_dist_op<double> randN;
  auto now =std::chrono::system_clock::now();
  auto seed = std::chrono::system_clock::to_time_t(now);

  randN.seed(static_cast<int>(seed));

  int GSize = (int)ceil(M / NumGroups);
  // groups are assigned consecutevely
  G.segment(GSize, GSize).array() = 1;
  G.segment(2 * GSize, GSize).array() = 2;
  

  B.segment(0, M) =
      sqrt(h2 / M) *
      Eigen::Matrix<double, Eigen::Dynamic, -1>::NullaryExpr(M, 1, randN);

  E = sqrt(sigmaE) *
      Eigen::Matrix<double, Eigen::Dynamic, -1>::NullaryExpr(N, 1, randN);

  X = Eigen::Matrix<double, Eigen::Dynamic, -1>::NullaryExpr(N, M, randN);

  Y = X * B + E;
  for (int i = 0; i < NumGroups; i++) {
    BSqN(i) =
        B.segment(i * ceil(M / NumGroups), ceil(M / NumGroups)).squaredNorm();
  }
}
