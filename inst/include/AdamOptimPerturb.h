#ifndef MONGREL_ADAMPERT_H
#define MONGREL_ADAMPERT_H

#include <RcppNumerical.h>
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::ArrayXd;
using Eigen::VectorXd;
using Eigen::Ref;

namespace adamperturb{

  // Parallel to LBFGSFun 
  // https://github.com/yixuan/RcppNumerical/blob/master/inst/include/optimization/wrapper.h
  class ADAMFun
  {
  private:
    Numer::MFuncGrad& f;
  public:
    ADAMFun(Numer::MFuncGrad& f_) : f(f_) {}
    inline double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
    {
      //Rcout << "ADAMFun x" << x.head(10).transpose() << std::endl;
      return f.f_grad(x, grad);
    }
  };

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
      VectorXd propsol; 
      
    public:
      ADAMOptim(ADAMFun& fun_, Numer::Refvec& thetainit, 
                double b1, double b2, double eta, double epsilon, 
                double eps_f, double eps_g, int tthresh, 
                int max_iter) : fun(fun_){
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
        val =-1e100;
        propsol = VectorXd::Zero(p);
      }
      
      int step(){
        R_CheckUserInterrupt();
        double val2 = fun(thetat, gt); // update gradient and value based on init
        //Rcout << "step gt" << gt.head(10).transpose() << std::endl;
        //Rcout << "step thetat" << thetat.head(10).transpose() << std::endl;
        double gnorm = gt.norm();
        double xnorm = thetat.norm();
        
        // eval stopping criteria
        if (t % 10 == 0){
          Rcout << "iter : " << t << std::endl;
          Rcout << "-Log Like: " << val2 << std::endl;
          Rcout << "rel improvement: " << ((val2 - val)/val2) << std::endl;
          Rcout << "gnorm, xnorm " << gnorm << "," << eps_g*std::max(xnorm, 1.0)
            << std::endl;
        }
        if (t > max_iter){
          return -1;    // max iter warning is -1
        }
        if ((gnorm < eps_g*std::max(xnorm, 1.0)) ||
            ((val2- val)/val2 > -eps_f)){
          // NEED SOMETHING ELSE HERE TO DETERMINE WHAT IS "CLOSE"
          propsol = thetat;
          thetat += sampleUniformSphere();
        }
        val = val2;
  
        // standard gradient descent
        thetat -= eta*gt;
        // Rcout << thetat.head(10).transpose() << std::endl;
        // Rcout << gt.head(10).transpose() << std::endl;
        
        // adam
        // if no stopping criteria met -- meat of algorithm
        // mt *= b1;
        // mt += (1-b1)*gt.array();
        // vt *= b2;
        // vt += (1-b2)*gt.array().square();
        // //Rcout << "mt step" << mt.head(10).transpose() << std::endl;
        // //Rcout << "vt step" << vt.head(10).transpose() << std::endl;
        // mthat = mt/(1-pow(b1,t));
        // vthat = vt/(1-pow(b2,t));
        // //Rcout << "mthat step" << mthat.head(10).transpose() << std::endl;
        // //Rcout << "vthat step" << vthat.head(10).transpose() << std::endl;
        // ArrayXd tmp(eta/(vthat.sqrt()+epsilon) * mthat);
        // //Rcout << "step tmp : " << tmp.head(10).transpose() << std::endl;
        // thetat -= tmp.matrix();
        t++;
        return 0;
      }
      
      double getVal(){return val;} // get optimal value
      VectorXd getTheta(){return thetat;} // get optimal parameter
      
      VectorXd sampleUniformSphere(){ // from pg 7 Jin et al (overcoming saddles)
        double u = R::runif(0,1);
        u = pow(u, 1/ ((double) p));
        NumericVector r = rnorm(p, 0, 1);
        Map<VectorXd> rvec(as<Map<VectorXd> >(r));
        return (u*rvec.array()/rvec.norm()).matrix();
      }
  };

  inline int optim_adam(Numer::MFuncGrad& f, 
                        Numer::Refvec theta, // initial value and thing returned 
                        double& fx_opt, 
                        double b1=0.9, 
                        double b2=0.99,
                        //double eta=0.003,
                        double eta = 0.001,
                        double epsilon=10e-7, // had been at 10e-6
                        double eps_f= 1e-8, 
                        double eps_g= 1e-6,
                        int tthresh = 4,
                        int max_iter= 100000){
    
    // create functor
    ADAMFun fun(f);
    
    // Solver
    ADAMOptim optim(fun, theta, b1, b2, eta, epsilon, 
                    eps_f, eps_g, tthresh, max_iter);
    
    int status =0; 
    while (status==0){
      status = optim.step();
    }
    if (status == -1){
      Rcpp::warning("Max iterations hit, may not be at optima");
    } else if (status == 1){
      Rcout << "Optimization terminated: change in gradient below threshold" 
            << std::endl;
    } else if (status ==2){
      Rcout << "Optimization terminated: change in function value below threshold" 
            << std::endl;
    }
    fx_opt = optim.getVal();
    theta = optim.getTheta();
    return status;
  }

  
}

#endif