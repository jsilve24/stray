#ifndef MONGREL_MMODEL_H
#define MONGREL_MMODEL_H

#include <RcppNumerical.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Ref;

namespace mongrel {

class MongrelModel : public Numer::MFuncGrad {
public:
  // Interface for optimization
  virtual double f_grad(Numer::Constvec& x, Numer::Refvec grad) = 0;
  
  // hessian vector multiplication  interface for Spectra library
  virtual int getN() = 0; // rows in hessian
  virtual int getD() = 0; // cols in hessian
  virtual VectorXd calcHessVectorProd(const Ref<const VectorXd>& etavec,
                                      VectorXd v, double r) = 0;
  virtual ~MongrelModel(){}
};

}

#endif