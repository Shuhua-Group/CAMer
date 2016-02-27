// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//overloaded version for x^v where x is a scalar and v is a arma::vector
arma::colvec pow(const double &x, const arma::colvec &v){
    arma::colvec ret(v);
    for(arma::uword i=0;i<ret.size();++i){
        ret.at(i)=pow(x,ret.at(i));
    }
    return ret;
}

//using simple regression formulae instead of QR factorization
//vec=c(theta0,theta1,ssE)
arma::vec fit_theta(const arma::colvec &Ac, const arma::colvec &Z){
    using namespace arma;
    double xbar=mean(Ac);
    double ybar=mean(Z);
    arma::vec ret(3);
    ret.at(1)=sum((Ac-xbar)%(Z-ybar))/sum(square(Ac-xbar)); //theta1
    ret.at(0)=ybar-ret.at(1)*xbar; //theta0
    ret.at(2)=sum(square(Z-ret.at(0)-ret.at(1)*Ac)); //ssE
    return ret;
}

//for any model and time period [end,start], compute the estimated theta0, theta1 and ssE
arma::vec est_fun(const int &model, const arma::mat &A, const arma::colvec &Z, const double &m1, const double &m2, const int &end, const int &start){
    using namespace arma;
    arma::vec Ac;
    int n;
    switch(model){
    case 1: //HI
        Ac=A.col(static_cast<arma::uword>(start)-1)*m1*m2; //offset -1
        break;
    case 2: //CGF1
        n=start-end+1;
        Ac=A.cols(static_cast<arma::uword>(end)-1,static_cast<arma::uword>(start)-1)
            *pow(m1,arma::linspace<arma::vec>(1-1.0/n,0,n))*(1-pow(m1,1.0/n))*m1; //offset -1
        break;
    case 3: //CGF2
        n=start-end+1;
        Ac=A.cols(static_cast<arma::uword>(end)-1,static_cast<arma::uword>(start)-1)
            *pow(m2,arma::linspace<arma::vec>(1-1.0/n,0,n))*(1-pow(m2,1.0/n))*m2; //offset -1
        break;
    case 4: //GA
        n=start-end+1;
        Ac.set_size(static_cast<arma::uword>(n)); //tempararily Ac=c(rep(n,n-1),1)
        Ac.fill(static_cast<double>(n));
        Ac.at(n-1)=1;
        Ac=A.cols(static_cast<arma::uword>(end)-1,static_cast<arma::uword>(start)-1)
            *(pow(1-1.0/n,arma::linspace<arma::vec>(0,n-1,n))/Ac)*m1*m2; //offset -1
        break;
    }
    return fit_theta(Ac,Z);
}

// [[Rcpp::export]]
List search(const int &model, const int &T, const int &max_duration, const arma::mat &A, const arma::colvec &Z, const double &m1, const double &m2, const bool &isolation, const bool &fast_search){
    using namespace arma;
    if(isolation){
        if(fast_search){
            if(model==1){
                double theta0,theta1,ssE=arma::datum::inf;
                arma::vec coef(3);
                int N;
                for(int start=1;start<=T;++start){
                    if(start%50==0) checkUserInterrupt();
                    coef=est_fun(model,A,Z,m1,m2,start,start);
                    if(coef.at(2)<ssE){
                        theta0=coef.at(0);
                        theta1=coef.at(1);
                        ssE=coef.at(2);
                        N=start;
                    }
                }
                return List::create(_["m"]=N,_["n"]=N,_["start"]=N,_["end"]=N,
                                    _["theta0"]=theta0,_["theta1"]=theta1,
                                    _["ssE"]=ssE,_["msE"]=ssE/(Z.size()-1));
            } else {
                struct Generation{
                    int end,start;
                    double theta0,theta1,ssE;
                    bool changed;
                    Generation(const int &end_, const int &start_, const double &theta0_, const double &theta1_, const double &best_ssE, const bool changed_=true):
                        end(end_),start(start_),theta0(theta0_),theta1(theta1_),ssE(best_ssE),changed(changed_) {}
                    void update_generation(const int &model, const arma::mat &A, const arma::colvec &Z, const double &m1, const double &m2){
                        int end_new=end+1;end_new=(end_new<start) ? end_new : start;
                        int start_new=start-1;start_new=(start_new>end) ? start_new : end;
                        
                        arma::vec coef_end_new=est_fun(model,A,Z,m1,m2,end_new,start);
                        arma::vec coef_start_new=est_fun(model,A,Z,m1,m2,end,start_new);
                        
                        double best_ssE=ssE;
                        if(best_ssE>coef_end_new.at(2)) best_ssE=coef_end_new.at(2);
                        if(best_ssE>coef_start_new.at(2)) best_ssE=coef_start_new.at(2);
                        
                        if(best_ssE==ssE){
                            changed=false;
                        } else if(best_ssE==coef_end_new.at(2)){
                            end=end_new;
                            theta0=coef_end_new.at(0);theta1=coef_end_new.at(1);
                            ssE=best_ssE;changed=true;
                        } else{
                            start=start_new;
                            theta0=coef_start_new.at(0);theta1=coef_start_new.at(1);
                            ssE=best_ssE;changed=true;
                        }
                        return;
                    }
                };
                
                arma::vec coef=est_fun(model,A,Z,m1,m2,1,T);
                Generation gen(1,T,coef.at(0),coef.at(1),coef.at(2));
                int check_interrupt_count=0;
                do{
                    gen.update_generation(model,A,Z,m1,m2);
                    ++check_interrupt_count;
                    if(check_interrupt_count==50){
                        checkUserInterrupt();
                        check_interrupt_count=0;
                    }
                } while (gen.changed);
                
                return List::create(_["m"]=gen.end,_["n"]=gen.start-gen.end+1,
                                    _["start"]=gen.start,_["end"]=gen.end,
                                    _["theta0"]=gen.theta0,_["theta1"]=gen.theta1,
                                    _["ssE"]=gen.ssE,_["msE"]=gen.ssE/(Z.size()-1));
            }
        } else {
            if(model==1){
                double theta0,theta1,ssE=arma::datum::inf;
                arma::vec coef(3);
                int N;
                for(int start=1;start<=T;++start){
                    if(start%50==0) checkUserInterrupt();
                    coef=est_fun(model,A,Z,m1,m2,start,start);
                    if(coef.at(2)<ssE){
                        theta0=coef.at(0);
                        theta1=coef.at(1);
                        ssE=coef.at(2);
                        N=start;
                    }
                }
                return List::create(_["m"]=N,_["n"]=N,_["start"]=N,_["end"]=N,
                                    _["theta0"]=theta0,_["theta1"]=theta1,
                                    _["ssE"]=ssE,_["msE"]=ssE/(Z.size()-1));
            } else {
                double theta0,theta1,ssE=arma::datum::inf;
                int N,M;
                arma::vec coef(3);
                for(int n=1;n<=max_duration;++n){
                    for(int m=1;m<=T-n+1;++m){
                        if(m%50==0 || m==T-n+1) checkUserInterrupt();
                        coef=est_fun(model,A,Z,m1,m2,m,m+n-1);
                        if(coef.at(2)<ssE){
                            theta0=coef.at(0);theta1=coef.at(1);ssE=coef.at(2);
                            M=m;N=n;
                        }
                    }
                }
                return List::create(_["m"]=M,_["n"]=N,_["start"]=M+N-1,_["end"]=M,
                                    _["theta0"]=theta0,_["theta1"]=theta1,
                                    _["ssE"]=ssE,_["msE"]=ssE/(Z.size()-1));
            }
        }
    } else {
        double theta0,theta1,ssE=arma::datum::inf;
        int N;
        arma::vec coef(3);
        for(int n=1;n<=T;++n){
            if(n%50==0) checkUserInterrupt();
            coef=est_fun(model,A,Z,m1,m2,1,n);
            if(coef.at(2)<ssE){
                theta0=coef.at(0);theta1=coef.at(1);ssE=coef.at(2);
                N=n;
            }
        }
        return List::create(_["m"]=(model==1)?N:1,_["n"]=N,_["start"]=N,_["end"]=(model==1)?N:1,
                            _["theta0"]=theta0,_["theta1"]=theta1,
                            _["ssE"]=ssE,_["msE"]=ssE/(Z.size()-1));
    }
}


/////////////////////////////////////////
//function for HI
// [[Rcpp::export]]
List search_HI(const int &T, const arma::mat &A, const arma::colvec &Z, const double &m1, const double &m2){
    using namespace arma;
    return search(1,T,150,A,Z,m1,m2,true,true);
}
