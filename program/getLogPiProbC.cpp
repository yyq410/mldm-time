#include<Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]

// Compute lgamma for vector
VectorXd lgammav(VectorXd a) {
  int length = a.rows();
  for(int i = 0; i < length; ++i)
    a(i) = R::lgammafn(a(i));
  
  return a;
}

// [[Rcpp::export]]
MatrixXd getLogPiProbC(NumericMatrix X, NumericMatrix M, NumericMatrix Z, int K, NumericVector Pi, List parameters) {
  Map<MatrixXd> Xc(as<Map<MatrixXd> >(X));
  Map<MatrixXd> Mc(as<Map<MatrixXd> >(M));
  Map<MatrixXd> Zc(as<Map<MatrixXd> >(Z));
  Map<VectorXd> Pic(as<Map<VectorXd> >(Pi));
  int n = Xc.rows();
  //int p = Xc.cols();
  
  MatrixXd PiProb(n, K);
  
  for(int i = 0; i < K; ++i) {
    List parametersPer = as<List>(parameters[i]);
    Map<MatrixXd> B(as<Map<MatrixXd> >(parametersPer[0]));
    Map<VectorXd> B0(as<Map<VectorXd> >(parametersPer[2]));
    Map<MatrixXd> Theta(as<Map<MatrixXd> >(parametersPer[1]));
    Map<VectorXd> muM(as<Map<VectorXd> >(parametersPer[3]));
    Map<MatrixXd> covM(as<Map<MatrixXd> >(parametersPer[4]));
    
    SelfAdjointEigenSolver<MatrixXd> eigensolverTheta(Theta);
    double LogDetTheta = eigensolverTheta.eigenvalues().array().log().sum();
    SelfAdjointEigenSolver<MatrixXd> eigensolverCovM(covM);
    double LogDetCovM = eigensolverCovM.eigenvalues().array().log().sum();
    double part2 = 0.5 * (LogDetTheta - LogDetCovM);
    
    for(int j = 0; j < n; ++j) {
      VectorXd x = Xc.row(j);
      VectorXd m = Mc.row(j);
      VectorXd z = Zc.row(i*n + j);
      
      VectorXd alpha = (B.transpose() * m + z).array().exp();
      double part1 = (lgammav(alpha + x) - lgammav(alpha)).sum() + (R::lgammafn(alpha.sum()) - R::lgammafn((alpha + x).sum()));
      
      double part4 = - ((z - B0).transpose() * Theta * (z - B0) + (m - muM).transpose() * covM.inverse() * (m - muM))(0,0) * 0.5;
      
      //cout << "part1:"  << part1 << endl;
      //cout << "part2:"  << part2 << endl;
      //cout << "part3:"  << log(Pic(i)) << endl;
      //cout << "part4:"  << part4 << endl;
      
      PiProb(j, i) = part1 + part2 + log(Pic(i)) + part4;
    }
  }
  
  return PiProb;
}
