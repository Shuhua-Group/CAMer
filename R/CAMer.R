#' Continuous Admixture Modeler (CAMer)
#' 
#' \pkg{CAMer} includes functions to do Continuous Admixture Modeling (CAM), generate summary plots, select the best-fit model(s), generate statistics to test if the results are credible and miscellaneous functionalities.
#' 
#' See \href{https://github.com/david940408/CAMer/blob/master/inst/doc/intro.md}{An Introduction to CAMer package} or intro.html in under inst/doc/ subdirectory of the package for an introduction. This file demonstrates how to use the functions.
#' 
#' @docType package
#' @name CAMer
NULL

#' Simulated .rawld File for CGF1
#'
#' A data frame read from a .rawld file by \code{read.table}. Forward simulation data using 10 chromosomes where the true model is CGF1 (50 generations) with admixture proportion of population 1 being 0.3.
#'
#' @details
#' The useful variables are as follows:
#' \itemize{
#' \item Distance: genetic distance
#' \item Combined_LD: weighted LD
#' \item Fitted: the fitted LD decay curve using the previous method
#' \item Jack?: Jackknives (?=1,2,...)
#' }
#' @format A data frame with 3497 rows and 21 variables
#' @name CGF_50
NULL

#' Simulated .rawld File for GA-I
#'
#' A data frame read from a .rawld file by \code{read.table}. Forward simulation data using 10 chromosomes where the true model is GA-I (start from the 100th generation and end at the 30th generation) with admixture proportion of population 1 being 0.3.
#'
#' @details
#' The useful variables are as follows:
#' \itemize{
#' \item Distance: genetic distance
#' \item Combined_LD: weighted LD
#' \item Fitted: the fitted LD decay curve using the previous method
#' \item Jack?: Jackknives (?=1,2,...)
#' }
#' @format A data frame with 3497 rows and 21 variables
#' @name GA_I
NULL


distance<-function(v1,v2){
    sum((v1-v2)^2)
}


fit.theta<-function(Ac,Z){
    X<-cbind(rep(1,length(Z)),Ac)
    QR<-qr(X)
    Q<-qr.Q(QR);R<-qr.R(QR)
    theta<-solve(R,t(Q)%*%Z)
    list(theta0=theta[1],theta1=theta[2],ssE=distance(theta[1]+Ac*theta[2],Z))
}

#' Continuous Admixture Modeling (CAM) for a Single LD Decay Curve
#'
#' Find the estimated time intervals/point for HI, CGF1(-I), CGF2(-I) and GA(-I) models and corresponding statictis (ssE, msE, etc.) for a single LD decay curve (e.g. Combined_LD or Jack? in a .rawld file).
#'
#' @param d the numeric vector of genetic distance (Morgan) of LD decay curve
#' @param Z the numeric vector of LD decay curve
#' @param m1 the admixture proportion of population 1. If m2 is the admixing proportion of population 2, then m1+m2=1.
#' @param T the most ancient generation to be searched. Defaults to 500.
#' @param isolation \code{TRUE} if the models used for fitting are HI, CGF1-I, CGF2-I and GA-I; \code{FALSE} if the models used for fitting are HI, CGF1, CGF2 and GA. Defaults to \code{TRUE}.
#' @param fast.search Defaults to \code{TRUE}. See "Details".
#' @param max.duration Defaults to 150. See "Details".
#' @param single.parallel a logical expression indicating whether parallel computation should be used. Defaults to \code{TRUE} if \code{isolation=TRUE,fast.search=FALSE} and \code{FALSE} otherwise.
#' @param single.clusternum the number of clusters in parallel computation. Defaults to 4 for the four models. Used if \code{single.parallel=TRUE}.
#' @return an object of S3 class "CAM.single". A list consisting of:
#' \item{call}{the matched call}
#' \item{maxindex}{the index of the maximal value in \code{Z} See "Details".}
#' \item{d,Z}{identical to function inputs up to some truncation. See "Details"}
#' \item{T,isolation}{identical to function inputs}
#' \item{A}{numeric matrix \eqn{A} with the \eqn{(i,j)}-th entry being \eqn{\text{exp}(-j \cdot d_i)}, \eqn{d_i} meaning the \eqn{i}-th entry of \code{d} and \eqn{j} meaning the genertion.}
#' \item{m1,m2}{admixture proportion of population 1 and 2}
#' \item{estimate}{a list of estimates. Each element contains the estimated parameters \eqn{m}, \eqn{n}, \eqn{\theta_0}, \eqn{\theta_1}, starting generation, ending generation and the corresponding ssE and msE. The time point for HI model is stored in \code{start} variable.}
#' \item{summary}{a data frame containing the information in \code{estimate} in a compact form}
#' @details
#'
#' \code{fast.search} is only used when \code{isolation=TRUE}. \code{TRUE} to use the fast searching algorithm, which sometimes gives slightly wider time intervals than the slow searching algorithm. It is recommended to use \code{fast.search=TRUE} (default), not only because it is significantly faster, but also because according to our experience it can partially solve the over-fitting problem of CGF-I and GA-I models so that HI usually does not perform significantly worse than them.
#'
#' \code{max.duration} is only used when \code{isolation=TRUE} and \code{fast.search=FALSE}. The maximal duration of admixture \eqn{n} to be considered as possible. Smaller values can make the slow searching algorithm faster. If \code{max.duration>T}, it will be set to be \code{T}.
#'
#' Given a single LD decay curve, for each model, this function does the following:
#'
#' If \code{isolation=FALSE}, it goes through all possible time intervals/points in [0,\code{T}], each time estimating \eqn{\theta_0} and \eqn{\theta_1} for the corresponding interval/point, and chooses the time interval/point that achieves the smallest ssE as the estimate for the model. Each corresponding \eqn{\theta=(\theta_0,\theta_1)} is the estimted \eqn{\theta} for each model.
#'
#' If \code{isolation=TRUE,fast.search=FALSE}, it also goes through all possible time intervals/points to estimate parameters. This slow algorithm is not recommended as it takes more than 40 minutes if \code{T=500L,max.duration=150L} and \code{Z} has length 3497 without parallel computation.
#'
#' If \code{isolation=TRUE,fast.search=TRUE}, for CGF1-I, CGF2-I, GA-I models, it uses a fast searching algorithm to search for a local minimum of ssE. This local minimum is not guaranteed to be the global minimum as that in the slow algorithm, but usually it is the same or quite close to that. It is recommended to use the fast algorithm because it takes only about 2 minutes if \code{T=500L,max.duration=150L} and \code{Z} has length 3497 without parallel computation.
#'
#' \code{maxindex} is the index of \code{Z} such that \code{Z[maxindex]} is the maximal value of \code{Z}. If the first few values of \code{Z} are not decreasing as theoretically expected, the \code{1:maxindex} of \code{Z} and \code{d} will be removed in calculation and in returned values.
#'
#' If the last entry of distence is greater than 10, a warning of unit will be given.
#' 
#' If the estimated time intervals/points cover \code{T}, a warning of too small \code{T} is given. The user should re-run the function with a larger \code{T} so that optimal time intervals/points can be reached.
#'
#' Require \pkg{parallel} or \pkg{snow} package installed if \code{single.parallel=TRUE}. For newer versions of \code{R (>=2.14.0)}, \pkg{parallel} is in R-core. If only \pkg{snow} is available, it is recommended to library it before using the parallel computing funcationality. When only \pkg{snow} is available, it will be \code{require}-d and hence the search path will be changed; if \pkg{parallel} is available, it will be used but the search path will not be changed. One may go to \url{https://cran.r-project.org/src/contrib/Archive/snow/} to download and install older versions of \pkg{snow} if the version of \code{R} is too old. If neither of the packages is available but \code{single.parallel=TRUE}, the function will compute sequentially with messages.
#'
#' Be aware that when the computational cost is small (e.g. \code{isolation=FALSE} or \code{T=20L,isoaltion=TRUE,fast.search=FALSE,max.duration=10L}), using parallel computation can result in longer computation time.
#'
#' There is a special method of \code{plot} and \code{print} for this class.
#' @examples
#' data(CGF_50)
#' Z<-CGF_50$Combined_LD
#' d<-CGF_50$Distance
#'
#' #fit models with isolation=FALSE
#' fit<-singleCAM(d,Z,m1=0.3,T=10L,isolation=FALSE) #with warning
#' 
#' #re-run with larger T
#' fit<-singleCAM(d,Z,m1=0.3,T=100L,isolation=FALSE)
#' fit
#'
#' #fit models with isolation=TRUE using fast searching algorithm
#' fit<-singleCAM(d,Z,m1=0.3,T=100L)
#' fit
#'
#' #fit models with isolation=TRUE using slow searching algorithm
#' #with parallel computation
#' fit<-singleCAM(d,Z,m1=0.3,T=100L,fast.search=FALSE,
#'                single.parallel=TRUE,single.clusternum=4L)
#' fit
#'
#' #fit models with isolation=TRUE using slow searching algorithm
#' #without parallel computation
#' fit<-singleCAM(d,Z,m1=0.3,T=70L,fast.search=FALSE,single.parallel=FALSE)
#' fit
#'
#' fitted.curves<-reconstruct.fitted(fit)
#' @note
#' When \code{LD.parallel=TRUE} or \code{single.parallel=TRUE}, it is not recommended to terminate the execution of the function. If \pkg{parallel} package is available, it is said that \code{\link[parallel]{setDefaultCluster}} from \pkg{parallel} can be used to remove the registered cluster, but real experiments does not support this; fortunately, these unused clusters will be removed automatically later, but with warnings. If only \pkg{snow} package is available, according to \url{http://homepage.stat.uiowa.edu/~luke/R/cluster/cluster.html}, "don't interrupt a snow computation". The ultimate method to close the unused clusters is probably to quit the R session.
#' @seealso \code{\link{CAM}}, \code{\link{reconstruct.fitted}}, \code{\link{conclude.model}}
#' @export

singleCAM<-function(d,Z,m1,T=500L,isolation=TRUE,
                    fast.search=TRUE,max.duration=150L,
                    single.parallel=isolation && !fast.search,
                    single.clusternum=4L){
    maxindex<-which(Z==max(Z))[1L]
    if(maxindex>1){
        Z<-Z[-seq_len(maxindex-1L)]
        d<-d[-seq_len(maxindex-1L)]
    }
    A<-exp(-d%*%t(seq_len(T)))
    m2<-1-m1
    s<-length(d)

    if(d[length(d)]>10) warning("The unit of Genetic Distance might not be Morgan.")

    result<-list(call=match.call(),maxindex=maxindex,
                 d=d,T=T,A=A,Z=Z,isolation=isolation,m1=m1,m2=m2)
    if(isolation){
        result$fast.search<-fast.search
        if(!fast.search){
            max.duration<-min(max.duration,T)
            result$max.duration<-max.duration
        }
    }
    class(result)<-"CAM.single"

    est.fun.HI<-function(n){
        Ac<-A[,n]*m1*m2
        fit.theta(Ac,Z)
    }
    est.fun.CGF1<-function(end,start){
        n<-start-end+1L
        alpha<-1-m1^(1/n)
        Ac<-A[,end:start]%*%matrix(m1^((n-1):0/n),ncol=1L)*alpha*m1
        fit.theta(Ac,Z)
    }
    est.fun.CGF2<-function(end,start){
        n<-start-end+1L
        alpha<-1-m2^(1/n)
        Ac<-A[,end:start]%*%matrix(m2^((n-1):0/n),ncol=1L)*alpha*m2
        fit.theta(Ac,Z)
    }
    est.fun.GA<-function(end,start){
        n<-start-end+1L
        Ac<-A[,end:start]%*%matrix((1-1/n)^(0:(n-1))/c(rep(n,n-1),1),ncol=1L)*m1*m2
        fit.theta(Ac,Z)
    }
    est.funs<-list(est.fun.HI,est.fun.CGF1,est.fun.CGF2,est.fun.GA)

    if(isolation){
        if(fast.search){
            upgrade.generation<-function(est.fun,gen.old){
                end.new<-min(gen.old$end+1L,gen.old$start)
                start.new<-max(gen.old$start-1L,gen.old$end)

                gen.end.new<-est.fun(end.new,gen.old$start)
                gen.start.new<-est.fun(gen.old$end,start.new)

                ssE<-min(sapply(list(gen.old,gen.end.new,gen.start.new),function(dummy) dummy$ssE))
                if(gen.old$ssE==ssE){
                    end<-gen.old$end
                    start<-gen.old$start
                    theta0<-gen.old$theta0
                    theta1<-gen.old$theta1
                } else if(gen.end.new$ssE==ssE){
                    end<-end.new
                    start<-gen.old$start
                    theta0<-gen.end.new$theta0
                    theta1<-gen.end.new$theta1
                } else {
                    end<-gen.old$end
                    start<-start.new
                    theta0<-gen.start.new$theta0
                    theta1<-gen.start.new$theta1
                }

                list(start=start,end=end,theta0=theta0,theta1=theta1,ssE=ssE,
                     changed=!(end==gen.old$end && start==gen.old$start))
            }

            search<-function(model){
                est.fun<-est.funs[[model]]
                if(model==1L){
                    ssE<-Inf
                    for(n in seq_len(T)){
                        coef<-est.fun(n)
                        if(coef$ssE<ssE){
                            theta0<-coef$theta0
                            theta1<-coef$theta1
                            ssE<-coef$ssE
                            N<-n
                        }
                    }
                    
                    list(m=NA,n=N,start=N,end=NA,theta0=theta0,theta1=theta1,ssE=ssE,msE=ssE/(length(Z)-1))
                } else {
                    est.old<-est.fun(1,T)
                    gen.old<-list(end=1,start=T,theta0=est.old$theta0,theta1=est.old$theta1,ssE=est.old$ssE)
                    repeat{
                        gen.new<-upgrade.generation(est.fun,gen.old)
                        if(gen.new$changed) gen.old<-gen.new
                        else break
                    }
                    
                    list(m=gen.new$end,n=gen.new$start-gen.new$end+1L,
                         start=gen.new$start,end=gen.new$end,
                         theta0=gen.new$theta0,theta1=gen.new$theta1,
                         ssE=gen.new$ssE,msE=gen.new$ssE/(length(Z)-1))
                }
            }
        } else
            search<-function(model){
                est.fun<-est.funs[[model]]
                if(model==1L){
                    ssE<-Inf
                    for(n in seq_len(T)){
                        coef<-est.fun(n)
                        if(coef$ssE<ssE){
                            theta0<-coef$theta0
                            theta1<-coef$theta1
                            ssE<-coef$ssE
                            N<-n
                        }
                    }
                    
                    list(m=NA,n=N,start=N,end=NA,theta0=theta0,theta1=theta1,ssE=ssE,msE=ssE/(length(Z)-1))
                } else {
                    SSE<-Inf

                    for(n in seq_len(max.duration)){
                        for(m in seq_len(T-n+1L)){
                            coef<-est.fun(m,m+n-1L)
                            if(coef$ssE<SSE){
                                theta0<-coef$theta0
                                theta1<-coef$theta1
                                SSE<-coef$ssE
                                M<-m
                                N<-n
                            }
                        }
                    }

                    list(m=M,n=N,start=M+N-1L,end=M,theta0=theta0,theta1=theta1,ssE=SSE,msE=SSE/(length(Z)-1))
                }
            }
    } else
        search<-function(model){
            est.fun<-est.funs[[model]]
            ssE<-Inf
            for(n in seq_len(T)){
                coef<-if(model==1L) est.fun(n) else est.fun(1L,n)
                if(coef$ssE<ssE){
                    theta0<-coef$theta0
                    theta1<-coef$theta1
                    ssE<-coef$ssE
                    N<-n
                }
            }
            
            list(m=if(model==1L) NA else 1L,n=N,start=N,end=if(model==1L) NA else 1L,theta0=theta0,theta1=theta1,ssE=ssE,msE=ssE/(length(Z)-1))
        }

    if(single.parallel && getRversion()<"2.14.0" && !suppressWarnings(require(snow,quietly=TRUE))){
        message("Cannot find 'parallel' or 'snow' package! Computing sequentially...\n")
        single.parallel<-FALSE
    }
    if(single.parallel){
        if(getRversion()>="2.14.0"){
            cl<-parallel::makeCluster(single.clusternum)
            parallel::clusterExport(cl,c("distance","fit.theta","est.funs","A","Z","m1","m2","T"),envir=environment())
            if(isolation && !fast.search)
                parallel::clusterExport(cl,"max.duration",envir=environment())
            if(isolation && fast.search)
                parallel::clusterExport(cl,"upgrade.generation",envir=environment())
            tryCatch(estimate<-parallel::parLapply(cl,seq_len(4L),search),
                     finally=parallel::stopCluster(cl))
        } else {
            require(snow,quietly=TRUE)
            cl<-makeCluster(single.clusternum)
            clusterExport(cl,c("distance","fit.theta","est.funs","A","Z","m1","m2","T"),envir=environment())
            if(isolation && !fast.search)
                clusterExport(cl,"max.duration",envir=environment())
            if(isolation && fast.search)
                clusterExport(cl,"upgrade.generation",envir=environment())
            tryCatch(estimate<-parLapply(cl,seq_len(4L),search),
                     finally=stopCluster(cl))
        }
    } else estimate<-lapply(seq_len(4L),search)

    names(estimate)<-if(isolation) c("HI","CGF1-I","CGF2-I","GA-I") else c("HI","CGF1","CGF2","GA")
    result$estimate<-estimate
    result$summary<-data.frame(Model=names(estimate),
                             Start=sapply(estimate,function(dummy) dummy$start),
                             End=sapply(estimate,function(dummy) dummy$end),
                             theta0=sapply(estimate,function(dummy) dummy$theta0),
                             theta1=sapply(estimate,function(dummy) dummy$theta1),
                             ssE=sapply(estimate,function(dummy) dummy$ssE),
                             Max.index=rep(result$maxindex,4L),
                             msE=sapply(estimate,function(dummy) dummy$msE))
    row.names(result$summary)<-NULL
    if(any(result$summary$Start==T)) warning("Most Ancient Generation T Reached! Consider Re-running with a Larger T.")

    result
}

#' Continuous Admixture Modeling (CAM)
#'
#' Estimate admixture time intervals/points for HI, CGF1(-I), CGF2(-I) and GA(-I) respectively for all Ld decay curves in a .rawld file.
#'
#' @param rawld a string representing the path of the .rawld file or a data frame read from the .rawld file by \code{read.table}.
#' @param m1 the admixture proportion of population 1. If m2 is the admixing proportion of population 2, then m1+m2=1.
#' @param T the most ancient generation to be searched. Defaults to 500.
#' @param isolation \code{TRUE} if the models used for fitting are HI, CGF1-I, CGF2-I and GA-I; \code{FALSE} if the models used for fitting are HI, CGF1, CGF2 and GA. Defaults to \code{TRUE}.
#' @param fast.search only used when \code{isolation=TRUE}. \code{TRUE} to use the fast searching algorithm, which sometimes gives slightly wider time intervals than the slow searching algorithm. Defaults to \code{TRUE}.
#' @param max.duration Defaults to 150. See "Details".
#' @param LD.parallel a logical expression indicating whether each LD decay curve should be computed parallely. Defaults to \code{TRUE}.
#' @param LD.clusternum the number of clusters in parallel computation. See "Details".
#' @param single.parallel a logical expression. See "Details".
#' @param single.clusternum the number of clusters in parallel computation. Defaults to 4 for the four models. Used if \code{single.parallel=TRUE}.
#' @return an object of S3 class "CAM". A list \code{CAM.list} consisting of some basic information of function call, N objects of "CAM.single" class (where N is the number of LD decay curves in the .rawld file), the fitted Ld decay curve fitted by the previous method (up to some truncation according to \code{max.index}) and a summary table containing the parameter estimates for each model and each curve, and diagnostic statistics msE and quasi-F. For details of "CAM.single" class, see \code{\link{singleCAM}}.
#'
#' There is a special method of \code{plot} and \code{print} for this class.
#' @details
#'
#' \code{max.duration} is only used when \code{isolation=TRUE} and \code{fast.search=FALSE}. The maximal duration of admixture \eqn{n} to be considered as possible. Smaller values can make the slow searching algorithm faster. If \code{max.duration>T}, it will be set to be \code{T}.
#'
#' \code{LD.clusternum} is used if \code{LD.parallel=TRUE}. If not specified, it is set to be the number of LD decay curves in the .rawld file.
#'
#' \code{single.parallel} indicates whether parallel computation should be used when computing a single LD decay curve. Defaults to \code{TRUE} if \code{isolation=TRUE,fast.search=FALSE} and \code{FALSE} otherwise.
#'
#' The .rawld file should include exactly one column nameed "Distance" in Morgan, exactly one column named "Combined_LD", several columns named "Jack?" representing Jackknives where ? is a number and exactly one column named "Fitted" representing the fitted LD decay curve using the previous method. This function fits "Combined_LD" and all Jackknives using all models. See \code{\link{singleCAM}} for further details of fitting algorithm for each LD decay curve.
#'
#' If the last entry of Distence in the .rawld file is greater than 10, a warning of unit will be given.
#' 
#' If the estimated time intervals/points cover \code{T}, a warning of too small \code{T} is given. The user should re-run the function with a larger \code{T} so that optimal time intervals/points can be reached.
#' 
#' Require \pkg{parallel} or \pkg{snow} package installed if \code{LD.parallel=TRUE} or \code{single.parallel=TRUE}. For newer versions of \code{R (>=2.14.0)}, \pkg{parallel} is in R-core. If only \pkg{snow} is available, it is recommended to library it before using the parallel computing funcationality. When only \pkg{snow} is available, it will be \code{require}-d and hence the search path will be changed; if \pkg{parallel} is available, it will be used but the search path will not be changed. One may go to \url{https://cran.r-project.org/src/contrib/Archive/snow/} to download and install older versions of \pkg{snow} if the version of \code{R} is too old. If neither of the packages is available but \code{LD.parallel=TRUE} or \code{single.parallel=TRUE}, the function will compute sequentially with messages.
#' 
#' Be aware that when the computational cost is small (e.g. \code{isolation=FALSE} or \code{T=20L,isoaltion=TRUE,fast.search=FALSE,max.duration=10L}), using parallel computation for single LD decay curves can result in longer computation time.
#'
#' @examples
#' data(GA_I)
#'
#' #fit models with isolation=FALSE.
#' fit<-CAM(GA_I,m1=0.3,T=150L,isolation=FALSE)
#' fit
#' \dontrun{
#' plot(fit) #may not be able to display
#' }
#' #Bad fitting indicates isolation=TRUE should be tried.
#' fit<-CAM(GA_I,m1=0.3,T=150L,isolation=TRUE)
#' fit
#' fit$summary
#' \dontrun{
#' plot(fit) #may not be able to display
#' plot(fit,"plot.pdf") #plot to a .pdf file
#' }
#'
#' data(CGF_50)
#' fit<-CAM(CGF_50,0.3,20L,isolation=FALSE,LD.parallel=FALSE) #with warnings
#' fit
#'
#' \dontrun{
#' #passing a file path to the argument `rawld=`
#' fit<-CAM("CGF_50.rawld",0.3)
#' }
#' @note 
#' When \code{LD.parallel=TRUE} or \code{single.parallel=TRUE}, it is not recommended to terminate the execution of the function. If \pkg{parallel} package is available, it is said that \code{\link[parallel]{setDefaultCluster}} from \pkg{parallel} can be used to remove the registered cluster, but real experiments do not support this; fortunately, these unused clusters will be removed automatically later, but with warnings. If only \pkg{snow} package is available, according to \url{http://homepage.stat.uiowa.edu/~luke/R/cluster/cluster.html}, "don't interrupt a snow computation". The ultimate method to close the unused clusters is probably to quit the R session.
#' 
#' Do care about memory allocation, especially when both \code{LD.parallel=TRUE} and \code{single.parallel=TRUE}.
#' 
#' It is possible that this function opens several nodes but non of them is computing, and hence the execution does not stop, especially when both \code{LD.parallel=TRUE} and \code{single.parallel=TRUE}. The cause has not been identified yet. The current solution is to terminate the function by hand and re-run the function with fewer cores (e.g. set \code{single.parallel=FALSE}).
#' @seealso \code{\link{construct.CAM}}, \code{\link{reconstruct.fitted}}, \code{\link{conclude.model}}
#' @import utils
#' @export
#'

CAM<-function(rawld,m1,T=500L,isolation=TRUE,
              fast.search=TRUE,max.duration=150L,
              LD.parallel=TRUE,LD.clusternum,
              single.parallel=isolation && !fast.search,
              single.clusternum=4L){
    if(is.character(rawld)) rawld<-utils::read.table(rawld,header=TRUE)
    Jack.index<-grep("Jack",names(rawld))
    LD.index<-grep("Combined_LD",names(rawld))
    Zs<-rawld[,c(LD.index,Jack.index)]
    d<-rawld$Distance

    if(d[length(d)]>10) warning("The unit of Genetic Distance might not be Morgan.")

    results<-list(call=match.call(),d=d,Zs=Zs,T=T,isolation=isolation,m1=m1,m2=1-m1)
    if(isolation){
        results$fast.search<-fast.search
        if(!fast.search){
            max.duration<-max(max.duration,T)
            results$max.duration<-max.duration
        }
    }
    class(results)<-"CAM"

    if(LD.parallel && getRversion()<"2.14.0" && !suppressWarnings(require(snow,quietly=TRUE))){
        message("Cannot find 'parallel' or 'snow' package! Computing each LD sequentially...\n")
        LD.parallel<-FALSE
    }
    if(LD.parallel){
        if(missing(LD.clusternum)) LD.clusternum<-ncol(Zs)
        if(getRversion()>="2.14.0"){
        cl<-parallel::makeCluster(LD.clusternum)
        parallel::clusterExport(cl,c("distance","fit.theta","d","m1","T","isolation","single.parallel"),envir=environment())
        if(isolation)
            parallel::clusterExport(cl,"fast.search",envir=environment())
        if(isolation && !fast.search)
            parallel::clusterExport(cl,"max.duration",envir=environment())
        if(single.parallel)
            parallel::clusterExport(cl,"single.clusternum",envir=environment())
        tryCatch(results$CAM.list<-parallel::parCapply(cl,Zs,singleCAM,d=d,m1=m1,T=T,isolation=isolation,fast.search=fast.search,max.duration=max.duration,single.parallel=single.parallel,single.clusternum=single.clusternum),
                 finally=parallel::stopCluster(cl))
        } else {
            require(snow,quietly=TRUE)
            cl<-makeCluster(single.clusternum)
            clusterExport(cl,c("distance","fit.theta","d","m1","T","isolation","single.parallel"),envir=environment())
            if(isolation)
                clusterExport(cl,"fast.search",envir=environment())
            if(isolation && !fast.search)
                clusterExport(cl,"max.duration",envir=environment())
            if(single.parallel)
                clusterExport(cl,"single.clusternum",envir=environment())
            tryCatch(results$CAM.list<-parCapply(cl,Zs,singleCAM,d=d,m1=m1,T=T,isolation=isolation,fast.search=fast.search,max.duration=max.duration,single.parallel=single.parallel,single.clusternum=single.clusternum),
                     finally=stopCluster(cl))
        }
    } else results$CAM.list<-lapply(seq_len(ncol(Zs)),function(ld){
        singleCAM(d,Zs[,ld],m1,T,isolation=isolation,single.parallel=single.parallel,fast.search=fast.search,max.duration=max.duration,single.clusternum=single.clusternum)
    })
    names(results$CAM.list)<-c("Combined_LD",paste("Jack",seq_len(length(Jack.index)),sep=""))

    results$fitted<-rawld$Fitted
    if(results$CAM.list[[1L]]$maxindex>1L) results$fitted<-results$fitted[-seq_len(results$CAM.list[[1L]]$maxindex-1L)]
    v<-distance(results$fitted,results$CAM.list[[1L]]$Z)

    data<-NULL
    for(ld in seq_len(ncol(Zs))){
        data.temp<-data.frame(LD=rep(names(results$CAM.list)[ld],4L),
                              Model=names(results$CAM.list[[ld]]$estimate),
                              Start=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$start),
                              End=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$end),
                              theta0=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$theta0),
                              theta1=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$theta1),
                              ssE=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$ssE),
                              Max.index=rep(results$CAM.list[[ld]]$maxindex,4L),
                              msE=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$msE),
                              quasi.F=if(ld==1L) sapply(results$CAM.list[[1L]]$estimate,function(dummy){dummy$ssE/v}) else rep(NA,4L))
        data<-rbind(data,data.temp)
    }
    data$LD<-as.factor(data$LD)
    data$Model<-factor(data$Model,ordered=TRUE)
    row.names(data)<-NULL
    results$summary<-data
    if(any(data$Start==T)) warning("Most Ancient Generation T Reached! Consider Re-running with a Larger T.")
    results
}

#' Reconstruct the Fitted LD Decay Curves of Four Models
#'
#' Given an object of "CAM.single" class, reconstruct the fitted LD decay curves of the four models (from \code{\link{singleCAM}}) or HI model (from \code{\link{singleHI}}).
#'
#' @param CAM.single an object of "CAM.single" class
#' @return a list consisting of the four fitted curves
#' @examples
#' data(CGF_50)
#' Z<-CGF_50$Combined_LD
#' d<-CGF_50$Distance
#'
#' fit<-singleCAM(d,Z,m1=0.3,T=100L,isolation=FALSE)
#' fitted.curves<-reconstruct.fitted(fit)
#' 
#' fit<-singleHI(d,Z,m1=0.3,T=100L)
#' fitted.curves<-reconstruct.fitted(fit)
#' @seealso \code{\link{singleCAM}}, \code{\link{singleHI}} \code{\link{construct.CAM}}
#' @export

reconstruct.fitted<-function(CAM.single){
    Z1<-CAM.single$estimate[[1L]]$theta0+CAM.single$estimate[[1L]]$theta1*
        CAM.single$A[,CAM.single$estimate[[1L]]$n]*CAM.single$m1*CAM.single$m2

    if(length(CAM.single$estimate)>1){
        alpha<-1-CAM.single$m1^(1/CAM.single$estimate[[2L]]$n)
        Z2<-CAM.single$estimate[[2L]]$theta0+CAM.single$estimate[[2L]]$theta1*
            CAM.single$A[,CAM.single$estimate[[2L]]$end:CAM.single$estimate[[2L]]$start]%*%
            matrix(CAM.single$m1^((CAM.single$estimate[[2L]]$n-1L):0/CAM.single$estimate[[2L]]$n),ncol=1L)*alpha*CAM.single$m1
        
        alpha<-1-CAM.single$m2^(1/CAM.single$estimate[[3L]]$n)
        Z3<-CAM.single$estimate[[3L]]$theta0+CAM.single$estimate[[3L]]$theta1*
            CAM.single$A[,CAM.single$estimate[[3L]]$end:CAM.single$estimate[[3L]]$start]%*%
            matrix(CAM.single$m2^((CAM.single$estimate[[3L]]$n-1):0/CAM.single$estimate[[3L]]$n),ncol=1L)*alpha*CAM.single$m2
        
        Z4<-CAM.single$estimate[[4L]]$theta0+CAM.single$estimate[[4L]]$theta1*
            CAM.single$A[,CAM.single$estimate[[4L]]$end:CAM.single$estimate[[4L]]$start]%*%
            matrix((1-1/CAM.single$estimate[[4L]]$n)^(0:(CAM.single$estimate[[4L]]$n-1L))/c(rep(CAM.single$estimate[[4L]]$n,CAM.single$estimate[[4L]]$n-1L),1L),ncol=1L)*CAM.single$m1*CAM.single$m2
    }

    z<-if(length(CAM.single$estimate)>1) list(Z1,Z2,Z3,Z4) else list(Z1)
    names(z)<-paste(names(CAM.single$estimate),".fitted",sep="")
    z
}

#' Summary Plots for "CAM" Class
#'
#' @param x an object of "CAM" class
#' @param filename pdf file path.
#' @param T.max the most ancient generation to be plotted. If an estimated time interval goes beyond T.max, there will be an arrow at the end of the line. Can be missing.
#' @param model.cols a matrix of colors. See "Details"
#' @param box.log a logical expression. If \code{TRUE}, the scale of the axis will be in log scale. Defaults to \code{TRUE}
#' @param box.lim the limit of msE shown in the boxplot. Defaults to the default of \code{\link[graphics]{boxplot}}
#' @param alpha a numeric value in [0,1] representitng alpha of the fitted LD decay curves for Combined LD. Defaults to 0.6. Can be a vector representing different alphas for different models. Ignored if alpha is specified in the third row of \code{model.cols}.
#' @param fit.lwd line width of the fitted LD decay curve for Combined LD of the four models. Defaults to 3. Can be a vector representing different line widths for different models.
#' @param LD.col color of the original Combined LD curve. Defaults to "black".
#' @param LD.lwd line width of the original Combined LD curve. Defaults to 1.
#' @param LD.lim the limit of the Weighted LD-axis in the third plot (Fitting of Models). Defaults to \code{c(-.002,.25)}.
#' @param ... further arguments
#' @details
#'
#' The function generates three plots in a device. The plot on the top left is the estimated time intervals/points for the four models. The color depth of segments/points corresponds to how many intervals/points covers this part in Jackknives. The deeper the color, the more estimates from Jackknives cover this part. The plot on the top right is the boxplot of msE for the four models. The third plot shows the fitting of four models to \code{Combined_LD} in the .rawld file. The numbers after model names in the legend are quasi-F values of the four models for \code{Combined_LD}.
#'
#' If \code{filename} is set, plot to the .pdf file, otherwise plot to the current device. The function is specially designed for a .pdf plot with width being 9.6 and height being 7.2. To add things to the plot and then save it to a file, better to set the size as above. May not be able to plot directly in an R window.
#'
#' The colors in Column 1/2/3/4 correspond to the colors for HI/CGF1(-I)/CGF2(-I)/GA(-I). The first/second row of \code{model.cols} is the lightest/deepest possible color in the "Time Intervals/Points" plot. The third row of \code{model.cols} is the color for "msE Boxplot" and "Fitting of Models" plot. The colors will be converted to RGB colors by \code{\link[grDevices]{col2rgb}}, so the input should be convertable by this function. The input of \code{model.cols} may also be the numeric vector version of the matrix stated above, i.e. \code{as.numeric(model.cols)}.
#' @examples
#' \dontrun{
#' data(CGF_50)
#' fit<-CAM(CGF_50,0.3,10L,isolation=FALSE,LD.parallel=FALSE) #will give warnings
#' plot(fit,"bad_fitting.pdf")
#' 
#' #may not be able to display
#' #This is not a very informative and user-friendly plot
#' plot(fit,T.max=2L,model.cols=matrix(c("pink","red","pink",
#'                              "lightseagreen","green","green",
#'                              "skyblue","blue","blue",
#'                              "yellow","orange","orange"),ncol=4),
#'      box.log=FALSE,box.lim=c(0,1),
#'      alpha=1,fit.lwd=1,LD.col="gray",LD.lwd=3,LD.lim=c(0,1))
#' }
#' @note
#' It is not recommended to pass other arguments to basic functions like \code{plot} in \code{...}.
#' 
#' It is recommended to select the best-fit model(s) according to the output of \code{\link{conclude.model}} rather than the boxplot of msE because \code{\link{conclude.model}} does paired t-tests, whose conclusion may look quite different from what we see from the boxplot.
#' 
#' @seealso \code{\link{CAM}}, \code{\link{conclude.model}}, \code{\link[grDevices]{rgb}}, \code{\link[grDevices]{col2rgb}}
#' @import graphics
#' @import grDevices
#' @export

plot.CAM<-function(x,filename,T.max,
                   model.cols=matrix(c("#ffa1c2B2","#b50d37B2","#da577c",
                                       "#9bff94","#0fbd02","#0fbd02",
                                       "#e9a1ff","#9111b8","#9111b8",
                                       "#7af8ff","#1ea1a8","#1ea1a8"),ncol=4),
                   box.log=TRUE,box.lim,
                   alpha=0.6,fit.lwd=3,LD.col="black",LD.lwd=1,LD.lim=c(-.002,.25),...){
    model.cols2<-grDevices::col2rgb(model.cols,alpha=TRUE)
    alpha<-alpha*rep(1,4L)
    for(model in 1L:4L){
        if(nchar(model.cols[3L,model])==7L)
            model.cols2[4L,3L*model]<-round(alpha[model]*model.cols2[4L,3L*model])
    }

    data<-x$summary
    data.jack<-data[grep("Jack",data$LD),]
    intervals<-lapply(unique(data$Model),function(i){
        data1<-data.jack[data.jack$Model==i,]
        if(i=="HI"){
            point<-data1$Start
            sapply(sort(unique(point)),function(t) as.numeric(c(t,sum(point==t))))
        } else {
            time<-unlist(lapply(seq_len(nrow(data1)),function(ii)
                data1[ii,"End"]:data1[ii,"Start"]))
            sapply(sort(unique(time)),function(t) as.numeric(c(t,sum(time==t))))
        }
    })
    NJack<-length(levels(data$LD))-1L

    if(!missing(filename)) grDevices::pdf(filename,width=9.6,height=7.2)

    graphics::layout(matrix(c(1,2,3,3),ncol=2L,nrow=2L,byrow=TRUE),widths=c(4,4),heights=c(1,2))

    graphics::par(bty="n",las=1)
    if(missing(T.max)){
        if(max(sapply(intervals,function(dummy) max(dummy[1,])))<=50) T.max<-100
        else if(max(sapply(intervals,function(dummy) max(dummy[1,])))<200) T.max<-200
        else if(max(sapply(intervals,function(dummy) max(dummy[1,])))<500) T.max<-500
        else if(max(sapply(intervals,function(dummy) max(dummy[1,])))<1000) T.max<-1000
        else if(max(sapply(intervals,function(dummy) max(dummy[1,])))<2000) T.max<-2000
        else if(max(sapply(intervals,function(dummy) max(dummy[1,])))<3000) T.max<-3000
        else T.max<-x$T
    }
    graphics::plot(x=1:T.max,y=seq(0,4,length.out=T.max),type="n",ann=FALSE,axes=FALSE,...)

    for(model in 1:4){
        colors<-grDevices::colorRampPalette(model.cols[1:2,model],alpha=TRUE)(NJack)
        dummy<-intervals[[model]]
        if(model==1){
            for(t in seq_len(ncol(dummy)))
                graphics::lines(x=dummy[1,t]+c(-.5,.5),y=rep(4.5-model,2),col=colors[dummy[2,t]],lwd=5,...)
        } else {
            if(max(dummy[1,])<=T.max){
                for(t in seq_len(ncol(dummy)))
                    graphics::lines(x=dummy[1,t]+c(-.5,.5),y=rep(4.5-model,2),col=colors[dummy[2,t]],lwd=5,...)
            } else {
                for(t in seq_len(T.max-min(dummy[1,])+1))
                    graphics::lines(x=dummy[1,t]+c(-.5,.5),y=rep(4.5-model,2),col=colors[dummy[2,t]],lwd=5,...)
                graphics::arrows(x0=T.max-5,y0=2.5,x1=T.max+2,lwd=3,col=colors[dummy[2,T.max-min(dummy[1,])+1]],angle=20,length=.2,...)
            }
        }
    }
    graphics::axis(1);graphics::axis(side=4,labels=NA,at=4:1-.5)
    graphics::title(main="Time Intervals/Points",xlab="Generation",ylab="")
    if(x$isolation){
        graphics::mtext(c("      HI"," CGF1-I"," CGF2-I","    GA-I"),side=4,line=1,at=4:1-.5)
    } else graphics::mtext(c("   HI","CGF1","CGF2","  GA"),side=4,line=2,at=4:1-.5)

    if(missing(box.lim)){
        graphics::boxplot(data.jack$msE[data.jack$Model=="GA"|data.jack$Model=="GA-I"]*1e5,
                          data.jack$msE[data.jack$Model=="CGF2"|data.jack$Model=="CGF2-I"]*1e5,
                          data.jack$msE[data.jack$Model=="CGF1"|data.jack$Model=="CGF1-I"]*1e5,
                          data.jack$msE[data.jack$Model=="HI"]*1e5,
                          horizontal=TRUE,boxwex=.8,pch=20,
                          log=if(box.log) "x" else "",medlwd=1,
                          boxfill=rev(model.cols[3,]),
                          outcol=rev(model.cols[3,]),
                          main="msE Boxplot",xlab=if(box.log) expression(paste("msE (",phantom()%*%10^{-5},") on log scale",sep="")) else expression(paste("msE (",phantom()%*%10^{-5},")",sep="")),axes=FALSE,...)
    } else {
        graphics::boxplot(data.jack$msE[data.jack$Model=="GA"|data.jack$Model=="GA-I"]*1e5,
                          data.jack$msE[data.jack$Model=="CGF2"|data.jack$Model=="CGF2-I"]*1e5,
                          data.jack$msE[data.jack$Model=="CGF1"|data.jack$Model=="CGF1-I"]*1e5,
                          data.jack$msE[data.jack$Model=="HI"]*1e5,
                          horizontal=TRUE,boxwex=.8,pch=20,
                          log=if(box.log) "x" else "",medlwd=1,
                          ylim=box.lim,
                          boxfill=rev(model.cols[3,]),
                          outcol=rev(model.cols[3,]),
                          main="msE Boxplot",xlab=if(box.log) expression(paste("msE (",phantom()%*%10^{-5},") on log scale",sep="")) else expression(paste("msE (",phantom()%*%10^{-5},")",sep="")),axes=FALSE,...)
    }
    graphics::axis(2,labels=NA);graphics::axis(1)

    graphics::par(bty="n",las=1,mar=c(5,17,4,10.5)+.1)
    colors<-apply(model.cols2,2,function(col) grDevices::rgb(col[1],col[2],col[3],col[4],maxColorValue=255))[seq(3,12,by=3)]
    fit.lwd<-fit.lwd*rep(1L,4)
    fitted<-reconstruct.fitted(x$CAM.list[[1]])
    graphics::plot(x$CAM.list[[1]]$d,x$CAM.list[[1]]$Z,
                   col=LD.col,type="l",lwd=LD.lwd,ylim=LD.lim,
                   xlab="Distance (Morgan)",ylab="",main="Fitting of Models",axes=FALSE,...)
    for(model in 1:4)
        graphics::lines(x$CAM.list[[1]]$d,fitted[[model]],col=colors[model],lwd=fit.lwd[model],...)
    graphics::legend(x="topright",
                     legend=paste(names(x$CAM.list[[1]]$estimate)," (",round(x$summary$quasi.F[1:4],3),")",sep=""),
                     col=colors,
                     pch=-1,lty=1,lwd=2,
                     ncol=1,bty="n")
    graphics::axis(2);graphics::axis(1);graphics::text(par("usr")[1]-0.15,y=LD.lim[2]*.6+LD.lim[1]*.4,labels="Weighted LD",xpd=TRUE)

    if(!missing(filename)) grDevices::dev.off()
}

#' Quickly Construct a Simple "CAM" Class Object from .rawld File and Summary Table
#'
#' Construct a simple "CAM" class object which can be passed to \code{\link{plot.CAM}} and whose elements in \code{CAM.list} can be passed to \code{\link{reconstruct.fitted}}. Can be used when only the .rawld file, \eqn{m_1} and the summary table is available (e.g., only the summary table was saved after running \code{\link{CAM}} last time).
#'
#' @param rawld original .rawld filepath or its data frame
#' @param m1 the admixture proportion of population 1.
#' @param dataset summary table as in \code{summary} of an object of "CAM" class
#' @return a simple "CAM" class object.
#' @note 
#' The returned onject is not as complate as the "CAM" class object obtained from \code{\link{CAM}}. Particularly, it does not include the information about how the estimates are found. The \code{T} and \code{A} are not the original ones. The \code{T} is the minial possible one, i.e. the smallest one that is sufficient to do the following analysis and construction.
#' @examples
#' data(GA_I)
#' fit<-CAM(GA_I,m1=0.3,T=150L)
#' dataset<-fit$summary
#' fit2<-construct.CAM(GA_I,m1=0.3,dataset)
#' \dontrun{
#' plot(fit2)
#' }
#' @seealso \code{\link{CAM}}, \code{\link{plot.CAM}}
#' @import utils
#' @export

construct.CAM<-function(rawld,m1,dataset){
    if(is.character(rawld)) rawld<-utils::read.table(rawld,header=TRUE)
    Jack.index<-grep("Jack",names(rawld))
    LD.index<-grep("Combined_LD",names(rawld))
    Zs<-rawld[,c(LD.index,Jack.index)]
    d<-rawld$Distance

    if(d[length(d)]>10) warning("The unit of Genetic Distance might not be Morgan.")

    T<-max(dataset$Start)

    dataset$LD<-as.character(dataset$LD)
    dataset$Model<-as.character(dataset$Model)
    dataset$Model.num<-sapply(dataset$Model,switch,
                              HI=1,CGF1=2,`CGF1-I`=2,CGF2=3,`CGF2-I`=3,GA=4,`GA-I`=4)
    dataset$Model.num<-NULL
    dataset$LD<-as.factor(dataset$LD)
    dataset$Model<-as.factor(dataset$Model)

    m2<-1-m1

    results<-list(isolation="CGF1-I" %in% levels(dataset$Model),d=d,Zs=Zs,T=T,CAM.list=NULL,fitted=rawld$Fitted,summary=dataset)
    class(results)<-"CAM"

    for(ld in seq_along(levels(dataset$LD))){
        d.temp<-d
        Z.temp<-Zs[,ld]
        data2<-dataset[(4*ld-3):(4*ld),]
        maxindex<-data2$Max.index[1L]
        if(maxindex>1){
            Z.temp<-Z.temp[-seq_len(maxindex-1L)]
            d.temp<-d.temp[-seq_len(maxindex-1L)]
        }
        A<-exp(-d.temp%*%t(seq_len(T)))
        results$CAM.list[[ld]]<-list(maxindex=maxindex,d=d.temp,T=T,A=A,Z=Z.temp,m1=m1,m2=m2)
        for(model in 1L:4L){
            data.temp<-data2[model,]
            results$CAM.list[[ld]]$estimate[[model]]<-list(m=data.temp$End,n=if(model==1L) data.temp$Start else data.temp$Start-data.temp$End+1,start=data.temp$Start,end=data.temp$End,theta0=data.temp$theta0,theta1=data.temp$theta1,ssE=data.temp$ssE,msE=data.temp$msE)
        }
        names(results$CAM.list[[ld]]$estimate)<-if(results$isolation) c("HI","CGF1-I","CGF2-I","GA-I") else c("HI","CGF1","CGF2","GA")
        results$CAM.list[[ld]]$summary<-data.frame(Model=names(results$CAM.list[[ld]]$estimate),
                                                   Start=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$start),
                                                   End=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$end),
                                                   theta0=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$theta0),
                                                   theta1=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$theta1),
                                                   ssE=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$ssE),
                                                   Max.index=rep(results$CAM.list[[ld]]$maxindex,4L),
                                                   msE=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$msE))
        class(results$CAM.list[[ld]])<-"CAM.single"
    }
    names(results$CAM.list)<-unique(dataset$LD)
    results
}

#' Print Method for "CAM.single" Class
#'
#' @param x a "CAM.single" class object
#' @param ... further arguments
#' @details
#' Print a very brief summary of a "CAM.single" class object. Include:
#' \itemize{
#' \item Function call
#' \item length of used LD (excluding first few values if necessary)
#' \item a data frame containing the estimated time intervals/points and corresponding msE. The time point for HI model is stored in \code{Start} variable.
#' }

#' @seealso \code{\link{singleCAM}}
#' @export

print.CAM.single<-function(x,...){
    cat("Continuous Admixture Inference (CAM) for a Single LD Decay Curve\n\n")
    cat("Function call: ")
    print(x$call)
    cat("\n")
    cat("Length of Used LD:", length(x$Z),"\n\n")
    print(x$summary[,c("Model","Start","End","msE")],row.names=FALSE,...)
}

#' Print Method for "CAM" Class
#'
#' @param x a "CAM" class object
#' @param ... further arguments
#' @details
#' The print method for "CAM" Class. Include:
#' \itemize{
#' \item Function call (if available from the object)
#' \item Total length of LD decay curve
#' \item a data frame containing the estimated intervals/points, msE and quasi-F for each model and each LD decay curve. The time point for HI model is stored in \code{Start} variable.
#' }
#' @seealso \code{\link{CAM}}
#' @export

print.CAM<-function(x,...){
    cat("Continuous Admixture Inference (CAM) for a .rawlf File\n\n")
    if(!is.null(x$call)){
        cat("Function call:")
        print(x$call)
        cat("\n")
    }
    cat("Total Length of LD:",length(x$d),"\n\n")
    print(x$summary[,c("LD","Model","Start","End","msE","quasi.F")],row.names=FALSE,...)
}

#' Draw Conclusions on Models from a "CAM" class object
#'
#' Draw conclusions on which models are best from a "CAM" class object or its summary table.
#'
#' @param x a "CAM" class object or the summary table of a "CAM" class object
#' @param alpha familywise type-I error rate. Defaults to 0.05
#' @param p.adjust.method method for adjusting p-values to adapt for familywise type-I error rate. Defaults to \code{"holm"}
#' @param log a logical expression. Whether log transformation should be applied to msE. Defaults to \code{TRUE}
#' @return an object of S3 class "CAM.conclusion". A list consisting of:
#' \item{call}{function call}
#' \item{group.means}{a named vector of group means of log(msE)/msE with each model being a group}
#' \item{adjusted.p.value}{a matrix of adjusted p-values for pairwise differences with the i,j-th entry being the adjusted p-value for the difference in log(msE)/msE of Model i and Model j. All entries on the diagonal are \code{NA}.}
#' \item{best.models}{a data frame consisting of best models concluded and their estimated time intervals/points}
#' \item{p.adjust.method}{method for adjusting p values used}
#' @details
#' The function uses pairwise paired Student's t-test on msE based on Jackknives to select the best model(s). If HI model is not significantly worse than any other model, it is chosen as the best model; otherwise, the model(s) with significantly smallest msE are chosen as best model(s).
#' 
#' The estimated interval is the one that include all time points covered by more than half of the intervals estimated from Jackknives. The estimated point (fot HI model) is the nearest integer to the mean of the points estimated from Jackknives.
#'
#' There is a special print method for this class (\code{\link{print.CAM.conclusion}}).
#' @note
#' The function only chooses the best model(s). It does NOT do any diagnostic analysis. Particularly, it does NOT check whether the best model(s) are credible.
#' @examples
#' data(CGF_50)
#' fit<-CAM(CGF_50,0.3,70L,isolation=FALSE)
#' conslusion<-conclude.model(fit)
#' 
#' conslusion<-conclude.model(fit,alpha=0.01,p.adjust.method="bonferroni",log=FALSE)
#' @seealso \code{\link{CAM}}, \code{\link[stats]{p.adjust}}, \code{\link[stats]{pairwise.t.test}}
#' @import stats
#' @export

conclude.model<-function(x,alpha=0.05,p.adjust.method="holm",log=TRUE){
    if(class(x)=="CAM") data<-x$summary else data<-x
    data<-data[grep("Jack",data$LD),]
    NJack<-nrow(data)/4
    if(log) data$msE<-log(data$msE)

    means<-tapply(data$msE,data$Model,mean)
    models<-names(means)

    p.value<-stats::pairwise.t.test(data$msE,data$Model,p.adjust.method=p.adjust.method,paired=TRUE)$p.value
    p.value[is.na(p.value)]<-0
    p.value<-cbind(rbind(0,p.value),0)
    p.value<-p.value+t(p.value)
    p.value[matrix(as.logical(diag(4L)),ncol=4L)]<-rep(NA,4L)
    colnames(p.value)<-row.names(p.value)<-models


    best<-which(means==min(means))[1L]
    if(means[4L]!=means[best] && p.value[4L,best]<alpha){
        best<-c(best,which(p.value[best,]>=alpha))
        best<-models[sort(best)]
    } else best<-"HI"
    
    end<-sapply(best,function(model){
        if(model!="HI"){
            data2<-data[data$Model==model,]
            e<-min(data2$End)
            while(sum(e>=data2$End)<=NJack/2) e<-e+1
            e
        } else NA
    })
    
    start<-sapply(best,function(model){
        data2<-data[data$Model==model,]
        if(model!="HI"){
            s<-max(data$Start)
            while(sum(s<=data2$Start)<=NJack/2) s<-s-1
            s
        } else round(mean(data2$Start))
    })

    conclusion<-list(call=match.call(),alpha=alpha,group.means=means,adjusted.p.value=p.value,best.models=data.frame(Best.Models=best,End=end,Start=start),p.adjust.method=p.adjust.method)
    class(conclusion)<-"CAM.conclusion"
    conclusion
}

#' Print Method for "CAM.conclusion" Class
#'
#' @param x an object of "CAM.conclusion" class
#' @param ... further arguments
#' @seealso \code{\link{conclude.model}}
#' @export

print.CAM.conclusion<-function(x,...){
    cat("CAM Best Model(s) Conclusion:\n\n")
    cat("Function call: ")
    print(x$call)
    cat("\n")
    cat("Familiwise Error Rate: ")
    cat(x$alpha,"\n\n",sep="")
    cat("Best Model(s) and Time Estimation:\n")
    print(x$best.models,row.names=FALSE)
    cat("\n")
    cat("Group Means of log(msE)/msE:\n")
    print(x$group.means,...)
    cat("\n")
    cat("Adjusted p-value:\n")
    print(x$adjusted.p.value,...)
}


#' Time Inference of HI model for a Single LD Decay Curve
#' 
#' Find the estimated time point for HI model and corresponding statictis (ssE, msE, etc.) for a single LD decay curve (e.g. Combined_LD or Jack? in a .rawld file).
#' @param d the numeric vector of genetic distance (Morgan) of LD decay curve
#' @param Z the numeric vector of LD decay curve
#' @param m1 the admixture proportion of population 1. If m2 is the admixing proportion of population 2, then m1+m2=1.
#' @param T the most ancient generation to be searched. Defaults to 500.
#' @return an object of S3 class "CAM.single". A list consisting of:
#' \item{call}{the matched call}
#' \item{maxindex}{the index of the maximal value in \code{Z} See "Details".}
#' \item{d,Z}{identical to function inputs up to some truncation. See "Details"}
#' \item{T,isolation}{identical to function inputs}
#' \item{A}{numeric matrix \eqn{A} with the \eqn{(i,j)}-th entry being \eqn{\text{exp}(-j \cdot d_i)}, \eqn{d_i} meaning the \eqn{i}-th entry of \code{d} and \eqn{j} meaning the genertion.}
#' \item{m1,m2}{admixture proportion of population 1 and 2}
#' \item{estimate}{a list of estimates. Each element contains the estimated parameters \eqn{m}, \eqn{n}, \eqn{\theta_0}, \eqn{\theta_1}, starting generation, ending generation and the corresponding ssE and msE. The time point for HI model is stored in \code{start} variable.}
#' \item{summary}{a data frame containing the information in \code{estimate} in a compact form}
#' @details 
#' This function is similar to \code{\link{singleCAM}}, except that it only considers HI model as the core model.
#' @examples 
#' data(CGF_50)
#' d<-CGF_50$Distance
#' Z<-CGF_50$Combined_LD
#' fitHI<-singleHI(d,Z,m1=.3,T=70L)
#' fitHI
#' @seealso \code{\link{singleCAM}}, \code{\link{HI}}, \code{\link{reconstruct.fitted}}
#' @export

singleHI<-function(d,Z,m1,T=500L){
    maxindex<-which(Z==max(Z))[1L]
    if(maxindex>1){
        Z<-Z[-seq_len(maxindex-1L)]
        d<-d[-seq_len(maxindex-1L)]
    }
    A<-exp(-d%*%t(seq_len(T)))
    m2<-1-m1
    s<-length(d)
    
    if(d[length(d)]>10) warning("The unit of Genetic Distance might not be Morgan.")
    
    result<-list(call=match.call(),maxindex=maxindex,
                 d=d,T=T,A=A,Z=Z,m1=m1,m2=m2)
    class(result)<-"CAM.single"
    
    est.fun<-function(n){
        Ac<-A[,n]*m1*m2
        fit.theta(Ac,Z)
    }
    
    ssE<-Inf
    for(n in seq_len(T)){
        coef<-est.fun(n)
        if(coef$ssE<ssE){
            theta0<-coef$theta0
            theta1<-coef$theta1
            ssE<-coef$ssE
            N<-n
        }
    }
    
    estimate<-list(list(m=NA,n=N,start=N,end=NA,theta0=theta0,theta1=theta1,ssE=ssE,msE=ssE/(length(Z)-1)))
    names(estimate)<-"HI"
    result$estimate<-estimate
    
    result$summary<-data.frame(Model=names(estimate),
                               Start=sapply(estimate,function(dummy) dummy$start),
                               End=sapply(estimate,function(dummy) dummy$end),
                               theta0=sapply(estimate,function(dummy) dummy$theta0),
                               theta1=sapply(estimate,function(dummy) dummy$theta1),
                               ssE=sapply(estimate,function(dummy) dummy$ssE),
                               Max.index=result$maxindex,
                               msE=sapply(estimate,function(dummy) dummy$msE))
    row.names(result$summary)<-NULL
    if(any(result$summary$Start==T)) warning("Most Ancient Generation T Reached! Consider Re-running with a Larger T.")
    
    result
}

#' Time Inference of HI model
#' 
#' Estimate admixture time intervals/points for HI, CGF1(-I), CGF2(-I) and GA(-I) respectively for all Ld decay curves in a .rawld file.
#' @param rawld a string representing the path of the .rawld file or a data frame read from the .rawld file by \code{read.table}.
#' @param m1 the admixture proportion of population 1. If m2 is the admixing proportion of population 2, then m1+m2=1.
#' @param T the most ancient generation to be searched. Defaults to 500.
#' @param LD.parallel a logical expression indicating whether each LD decay curve should be computed parallely. Defaults to \code{TRUE}.
#' @param LD.clusternum the number of clusters in parallel computation. See "Details".
#' @return an object of S3 class "CAM". A list \code{CAM.list} consisting of some basic information of function call, N objects of "CAM.single" class (where N is the number of LD decay curves in the .rawld file) calculated by \code{\link{singleHI}}, the fitted Ld decay curve fitted by the previous method (up to some truncation according to \code{max.index}) and a summary table containing the parameter estimates for each model and each curve, and diagnostic statistics msE and quasi-F. See \code{\link{singleCAM}} for more details of this kind of "CAM.songle" class.
#'
#' There is a special method of \code{print} for this class.
#' @details
#' \code{LD.clusternum} is used if \code{LD.parallel=TRUE}. If not specified, it is set to be the number of LD decay curves in the .rawld file.
#'
#' The .rawld file should include exactly one column nameed "Distance" in Morgan, exactly one column named "Combined_LD", several columns named "Jack?" representing Jackknives where ? is a number and exactly one column named "Fitted" representing the fitted LD decay curve using the previous method. This function fits "Combined_LD" and all Jackknives using HI model. See \code{\link{singleHI}} for further details of fitting algorithm for each LD decay curve.
#'
#' If the last entry of Distence in the .rawld file is greater than 10, a warning of unit will be given.
#' 
#' If the estimated time intervals/points cover \code{T}, a warning of too small \code{T} is given. The user should re-run the function with a larger \code{T} so that optimal time intervals/points can be reached.
#' 
#' Require \pkg{parallel} or \pkg{snow} package installed if \code{LD.parallel=TRUE} or \code{single.parallel=TRUE}. For newer versions of \code{R (>=2.14.0)}, \pkg{parallel} is in R-core. If only \pkg{snow} is available, it is recommended to library it before using the parallel computing funcationality. When only \pkg{snow} is available, it will be \code{require}-d and hence the search path will be changed; if \pkg{parallel} is available, it will be used but the search path will not be changed. One may go to \url{https://cran.r-project.org/src/contrib/Archive/snow/} to download and install older versions of \pkg{snow} if the \code{R} version is too old. If neither of the packages is available but \code{LD.parallel=TRUE} or \code{single.parallel=TRUE}, the function will compute sequentially with messages.
#' 
#' Be aware that when the computational cost is small (e.g. \code{T=20L}), using parallel computation for single LD decay curves can result in longer computation time.
#' @note
#' Although the output is a "CAM" class object, it should \emph{NOT} be passed to \code{\link{plot.CAM}}. Its summary table should \emph{NOT} be passed to \code{\link{construct.CAM}} either.
#' 
#' When \code{LD.parallel=TRUE}, it is not recommended to terminate the execution of the function. If \pkg{parallel} package is available, it is said that \code{\link[parallel]{setDefaultCluster}} from \pkg{parallel} can be used to remove the registered cluster, but real experiments do not support this; fortunately, these unused clusters will be removed automatically later, but with warnings. If only \pkg{snow} package is available, according to \url{http://homepage.stat.uiowa.edu/~luke/R/cluster/cluster.html}, "don't interrupt a snow computation". The ultimate method to close the unused clusters is probably to quit the R session.
#' @examples 
#' data(GA_I)
#' fit<-HI(GA_I,m1=0.3,T=150L)
#' fit
#' 
#' \dontrun{
#' #passing a file path to the argument `rawld=`
#' fit<-CAM("CGF_50.rawld",m1=0.3,T=150L)
#' }
#' @import utils
#' @seealso \code{\link{singleHI}}, \code{\link{CAM}}
#' @export

HI<-function(rawld,m1,T=500L,LD.parallel=TRUE,LD.clusternum){
    if(is.character(rawld)) rawld<-utils::read.table(rawld,header=TRUE)
    Jack.index<-grep("Jack",names(rawld))
    LD.index<-grep("Combined_LD",names(rawld))
    Zs<-rawld[,c(LD.index,Jack.index)]
    d<-rawld$Distance
    
    if(d[length(d)]>10) warning("The unit of Genetic Distance might not be Morgan.")
    
    results<-list(call=match.call(),d=d,Zs=Zs,T=T,m1=m1,m2=1-m1)
    class(results)<-"CAM"
    
    if(LD.parallel && getRversion()<"2.14.0" && !suppressWarnings(require(snow,quietly=TRUE))){
        message("Cannot find 'parallel' or 'snow' package! Computing each LD sequentially...\n")
        LD.parallel<-FALSE
    }
    if(LD.parallel){
        if(missing(LD.clusternum)) LD.clusternum<-ncol(Zs)
        if(getRversion()>="2.14.0"){
            cl<-parallel::makeCluster(LD.clusternum)
            parallel::clusterExport(cl,c("distance","fit.theta","d","m1","T"),envir=environment())
            tryCatch(results$CAM.list<-parallel::parCapply(cl,Zs,singleHI,d=d,m1=m1,T=T),
                     finally=parallel::stopCluster(cl))
        } else {
            require(snow,quietly=TRUE)
            cl<-makeCluster(single.clusternum)
            clusterExport(cl,c("distance","fit.theta","d","m1","T"),envir=environment())
            tryCatch(results$CAM.list<-parCapply(cl,Zs,singleHI,d=d,m1=m1,T=T),
                     finally=stopCluster(cl))
        }
    } else results$CAM.list<-lapply(seq_len(ncol(Zs)),function(ld){
        singleHI(d,Zs[,ld],m1,T)
    })
    names(results$CAM.list)<-c("Combined_LD",paste("Jack",seq_len(length(Jack.index)),sep=""))
    
    results$fitted<-rawld$Fitted
    if(results$CAM.list[[1L]]$maxindex>1L) results$fitted<-results$fitted[-seq_len(results$CAM.list[[1L]]$maxindex-1L)]
    v<-distance(results$fitted,results$CAM.list[[1L]]$Z)
    
    data<-NULL
    for(ld in seq_len(ncol(Zs))){
        data.temp<-data.frame(LD=names(results$CAM.list)[ld],
                              Model=names(results$CAM.list[[ld]]$estimate),
                              Start=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$start),
                              End=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$end),
                              theta0=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$theta0),
                              theta1=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$theta1),
                              ssE=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$ssE),
                              Max.index=results$CAM.list[[ld]]$maxindex,
                              msE=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$msE),
                              quasi.F=if(ld==1L) sapply(results$CAM.list[[1L]]$estimate,function(dummy){dummy$ssE/v}) else NA)
        data<-rbind(data,data.temp)
    }
    data$LD<-as.factor(data$LD)
    data$Model<-factor(data$Model,ordered=TRUE)
    row.names(data)<-NULL
    results$summary<-data
    if(any(data$Start==T)) warning("Most Ancient Generation T Reached! Consider Re-running with a Larger T.")
    results
}