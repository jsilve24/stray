#ifndef MONGREL_ADAM_H
#define MONGREL_ADAM_H

#include <RcppNumerical.h>

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::ArrayXd;
using Eigen::VectorXd;
using Eigen::Ref;


namespace adam{

  // Parallel structure to LBFGSFun 
  // https://github.com/yixuan/RcppNumerical/blob/master/inst/include/optimization/wrapper.h
  // This class is a Functor to calculate gradient for ADAM optimizer
  class ADAMFun
  {
  private:
    Numer::MFuncGrad& f;
  public:
    ADAMFun(Numer::MFuncGrad& f_) : f(f_) {}
    inline double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
    {
      return f.f_grad(x, grad);
    }
  };

  // Class for Adam Optimizer
  class ADAMOptim
  {
    private:
      int p;
      VectorXd thetat; // position
      VectorXd gt; // gradient 
      ArrayXd mt; // first moment
      ArrayXd vt; // second moment
      ArrayXd mthat; // first moment bias correction
      ArrayXd vthat; // second moment bias correction
      double b1;
      double b2;
      int t; // timestep
      ADAMFun& fun;
      double val; // value at current location
      double eta; 
      double epsilon;
      
      // stopping criteria
      double eps_f;
      double eps_g;
      double max_iter;
      
      // verbose behavior
      bool verbose;
      int verbose_rate;
      
    public:
      // Main constructor 
      //   fun_ : ADAMFun object (functor)
      //   thetainit : initial parameter estimates
      //   b1 : 1st moment decay parameter (recomend 0.9) "aka momentum"
      //   b2 : 2nd moment decay parameter (recommend 0.99 or 0.999)
      //   eta : step size (recomend 0.001)
      //   epsilon : parameter to avoid divide by zero in adam
      //   eps_f : normalized funtion improvement for stopping
      //   eps_g : normalized gradient magnitute for stopping
      //   max_iter : maximum number of iterations before stopping
      //   verbose : if true will print stats for stopping criteria and iter no.
      //   verbose_rate : rate to print verbose stats to screen
      ADAMOptim(ADAMFun& fun_, Numer::Refvec& thetainit, 
                double b1, double b2, double eta, double epsilon, 
                double eps_f, double eps_g, int max_iter, 
                bool verbose, int verbose_rate) : fun(fun_){
        p = thetainit.size();
        thetat = thetainit;
        mt = ArrayXd::Zero(p);
        vt = ArrayXd::Zero(p);
        gt = VectorXd::Zero(p);
        mthat = mt;
        vthat = vt;
        this -> b1 = b1;
        this -> b2 = b1;
        this -> eta = eta;
        this -> epsilon = epsilon;
        this -> eps_f = eps_f; 
        this -> eps_g = eps_g; 
        this -> max_iter = max_iter;
        t = 1;
        val=0;
        this -> verbose = verbose;
        this -> verbose_rate = verbose_rate;
      }
      
      int step(){
        R_CheckUserInterrupt();
        double val2 = fun(thetat, gt); // update gradient and value based on init
        double gnorm = gt.norm();
        double xnorm = thetat.norm();
        
        // eval stopping criteria
        if (verbose){
          if (t % verbose_rate == 0){
            Rcout << "iter : " << t << std::endl;
            Rcout << "-Log Like: " << val2 << std::endl;
            Rcout << "normalized rel improvement: " << ((val2 - val)/val2) 
                  << std::endl;
            Rcout << "gnorm, gradient threshold " << gnorm << "," 
                  << eps_g*std::max(xnorm, 1.0) << std::endl;
          }
        }
        if (gnorm <= eps_g*std::max(xnorm, 1.0)){
          return 1; // gradient below threshold
        }
        if ((((val2 - val)/val2) > -eps_f) && (((val2 - val)/val2) != 1)){
          return 2; // function value improvement below threshold
        }
        if (t > max_iter){
          return -1;    // max iter warning is -1
        }
        val = val2;
    
        // if no stopping criteria met -- meat of algorithm
        mt *= b1;
        mt += (1-b1)*gt.array();
        vt *= b2;
        vt += (1-b2)*gt.array().square();
        mthat = mt/(1-pow(b1,t));
        vthat = vt/(1-pow(b2,t));
        ArrayXd tmp(eta/(vthat.sqrt()+epsilon) * mthat);
        thetat -= tmp.matrix();
        t++;
        return 0;
      }
      
      double getVal(){return val;} // get optimal value
      VectorXd getTheta(){return thetat;} // get optimal parameter
  };

  // Main Function to Call from other C++ functions 
  //   fun_ : ADAMFun object (functor)
  //   thetainit : initial parameter estimates
  //   b1 : 1st moment decay parameter (recomend 0.9) "aka momentum"
  //   b2 : 2nd moment decay parameter (recommend 0.99 or 0.999)
  //   eta : step size (recomend 0.001)
  //   epsilon : parameter to avoid divide by zero in adam
  //   eps_f : normalized funtion improvement for stopping
  //   eps_g : normalized gradient magnitute for stopping
  //   max_iter : maximum number of iterations before stopping
  //   verbose : if true will print stats for stopping criteria and iter no.
  //   verbose_rate : rate to print verbose stats to screen
  inline int optim_adam(Numer::MFuncGrad& f, 
                        Numer::Refvec theta, // initial value and thing returned 
                        double& fx_opt, 
                        double b1=0.9, 
                        double b2=0.99,
                        double eta=0.003,
                        double epsilon=10e-7, // had been at 10e-6
                        double eps_f= 1e-8, 
                        double eps_g= 1e-5, 
                        int max_iter= 10000, 
                        bool verbose=false, 
                        int verbose_rate=10){
    
    // create functor
    ADAMFun fun(f);
    
    // Solver
    ADAMOptim optim(fun, theta, b1, b2, eta, epsilon, 
                    eps_f, eps_g, max_iter, verbose, verbose_rate);
    
    int status = 0; 
    while (status == 0){
      status = optim.step();
    }
    if (status == -1){
      Rcpp::warning("Max iterations hit, may not be at optima");
    } else if ((status == 1) && verbose){
      Rcout << "Optimization terminated: change in gradient below threshold" 
            << std::endl;
    } else if ((status ==2) && verbose){
      Rcout << "Optimization terminated: change in function value below threshold" 
            << std::endl;
    }
    fx_opt = optim.getVal();
    theta = optim.getTheta();
    return status;
  }
  
}

#endif
