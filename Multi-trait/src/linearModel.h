
#include <Eigen/Core>
#include <cmath>
#include <math.h>
#include <random>
// functor to sample a i.i.d. normal column vector

template <typename Scalar> struct scalar_normal_dist_op {
  static std::mt19937 rng; // The uniform pseudo-random algorithm
  mutable std::normal_distribution<Scalar> norm; // gaussian combinator

  EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

  template <typename Index>
  inline const Scalar operator()(Index, Index = 0) const {
    return norm(rng);
  }
  inline void seed(const uint64_t &s) { rng.seed(s); }
};

template <typename Scalar> std::mt19937 scalar_normal_dist_op<Scalar>::rng;

template <typename Scalar>
struct Eigen::internal::functor_traits<scalar_normal_dist_op<Scalar>> {
  enum {
    Cost = 50 * NumTraits<Scalar>::MulCost,
    PacketAccess = false,
    IsRepeatable = false
  };
};

class LinearModel {
public:
  LinearModel(){};
  LinearModel(int _Mtot, int _M, int _N, int _NumGroups, double _h2)
      : Mtot(_Mtot), M(_M), N(_N), NumGroups(_NumGroups), h2(_h2) {
    Mtot = 4000;   // total markers
    M = 100;       // causal variants or markers in LD with causal variants
    N = 2000;      // Number of individuals
    NumGroups = 3; // Number of groups
    h2 = 0.6;
    sigmaE = 1 - h2;
    Y.resize(N);
    X.resize(N, M);
    E.resize(N);
    B.resize(M);
    G.resize(M);
    BSqN.resize(NumGroups);

    
    Y.setZero();
    X.setZero();
    E.setZero();
    B.setZero();
    G.setZero();
    BSqN.setZero();
  };
  void simulate();

  int Mtot;
  int M;
  int N;
  int NumGroups;
  double h2;
  double sigmaE;

  Eigen::VectorXd Y;
  Eigen::MatrixXd X;
  Eigen::VectorXi G;
  Eigen::VectorXd E;
  Eigen::VectorXd B;
  Eigen::VectorXd BSqN;
  Eigen::VectorXd sigmaG;
};
