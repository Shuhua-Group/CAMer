#' Continuous Admixture Modeling (CAM)
#' 
#' \pkg{CAM} includes functions to do Continuous Admixture Modeling (CAM), generate summary plots, select the best-fit model(s), generate statistics to test if the results are credible and miscellaneous functionalities.
#' 
#' Type \code{browseVignettes(CAM)} in R for the vignettes to see some examples of how to use the functions.
#' 
#' @docType package
#' @name CAM-package
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


fit.theta<-function(Ac,y){
    X<-cbind(rep(1,length(y)),Ac)
    QR<-qr(X)
    Q<-qr.Q(QR);R<-qr.R(QR)
    theta<-solve(R,t(Q)%*%y)
    list(theta0=theta[1],theta1=theta[2],ssE=distance(theta[1]+Ac*theta[2],y))
}

#' Continuous Admixture Modeling (CAM) for a Single LD Decay Curve
#'
#' Find the estimated time intervals/point for HI, CGF1(-I), CGF2(-I) and GA(-I) models and corresponding ssE for a single LD decay curve (e.g. Combined_LD or Jack? in a .rawld file).
#'
#' @param d the numeric vector of genetic distance (Morgan) of LD decay curve
#' @param T the most ancient generation to be searched. Defaults to 500.
#' @param y the numeric vector of LD decay curve
#' @param m1 the admixture proportion of population 1. If m2 is the admixing proportion of population 2, then m1+m2=1.
#' @param isolation \code{TRUE} if the models used for fitting are HI, CGF1-I, CGF2-I and GA-I; \code{FALSE} if the models used for fitting are HI, CGF1, CGF2 and GA. Defaults to \code{TRUE}.
#' @param fast.search Defaults to \code{TRUE}. See "Details".
#' @param max.duration Defaults to 150. See "Details".
#' @param single.parallel a logical expression indicating whether parallel computation should be used. Defaults to \code{TRUE} if \code{isolation=TRUE,fast.search=FALSE} and \code{FALSE} otherwise.
#' @param single.clusternum the number of clusters in parallel computation. Defaults to 4 for the four models. Used if \code{single.parallel=TRUE}.
#' @return an object of S3 class "CAM.single". A list consisting of:
#' \item{call}{the matched call}
#' \item{maxindex}{the index of the maximal value in \code{y} See "Details".}
#' \item{d,y}{identical to function inputs up to some truncation. See "Details"}
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
#' If \code{isolation=TRUE,fast.search=FALSE}, it also goes through all possible time intervals/points to estimate parameters. This slow algorithm is not recommended as it takes more than 40 minutes if \code{T=500L,max.duration=150L} and \code{y} has length 3497 without parallel computation.
#'
#' If \code{isolation=TRUE,fast.search=TRUE}, for CGF1-I, CGF2-I, GA-I models, it uses a fast searching algorithm to search for a local minimum of ssE. This local minimum is not guaranteed to be the global minimum as that in the slow algorithm, but usually it is the same or quite close to that. It is recommended to use the fast algorithm because it takes only about 2 minutes if \code{T=500L,max.duration=150L} and \code{y} has length 3497 without parallel computation.
#'
#' \code{maxindex} is the index of \code{y} such that \code{y[maxindex]} is the maximal value of \code{y}. If the first few values of \code{y} are not decreasing as theoretically expected, the \code{1:maxindex} of \code{y} and \code{d} will be removed in calculation and in returned values.
#'
#' If the last entry of distence is greater than 10, a warning of unit will be given.
#'
#' Require \pkg{doSNOW}, \pkg{foreach} package and their dependencies if \code{single.parallel=TRUE}. It is recommended to library these required packages before using the parallel functionality.
#'
#' Be aware that when the computational cost is small (e.g. \code{isolation=FALSE} or \code{T=20L,isoaltion=TRUE,fast.search=FALSE,max.duration=10L}), using parallel computation can result in longer computatio time.
#'
#' There is a special method of \code{plot} and \code{print} for this class.
#' @examples
#' data(CGF_50)
#' y<-CGF_50$Combined_LD
#' d<-CGF_50$Distance
#'
#' #fit models with isolation=FALSE
#' fit<-singleCAM(d,y,m1=0.3,T=100L,isolation=FALSE)
#' fit
#'
#' #fit models with isolation=TRUE using fast searching algorithm
#' fit<-singleCAM(d,y,m1=0.3,T=100L)
#' fit
#'
#' #fit models with isolation=TRUE using slow searching algorithm
#' #with parallel computation
#' library(foreach);library(doSNOW)
#' fit<-singleCAM(d,y,m1=0.3,T=100L,fast.search=FALSE,
#'                single.parallel=TRUE,single.clusternum=4L)
#' fit
#'
#' #fit models with isolation=TRUE using slow searching algorithm
#' #without parallel computation
#' fit<-singleCAM(d,y,m1=0.3,T=70L,fast.search=FALSE,single.parallel=FALSE)
#' fit
#'
#' fitted.curves<-reconstruct.fitted(fit)
#' @seealso \code{\link{CAM}}, \code{\link{reconstruct.fitted}}, \code{\link{conclude.model}}
#' @export

singleCAM<-function(d,y,m1,T=500L,isolation=TRUE,
                    fast.search=TRUE,max.duration=150L,
                    single.parallel=isolation && !fast.search,
                    single.clusternum=4L){
    maxindex<-which(y==max(y))[1]
    if(maxindex>1){
        y<-y[-seq_len(maxindex-1L)]
        d<-d[-seq_len(maxindex-1L)]
    }
    A<-exp(-d%*%t(seq_len(T)))
    m2<-1-m1
    s<-length(d)

    if(d[length(d)]>10) warning("The unit of Genetic Distance might not be Morgan.")

    result<-list(call=match.call(),maxindex=maxindex,
                 d=d,T=T,A=A,y=y,isolation=isolation,m1=m1,m2=m2)
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
        fit.theta(Ac,y)
    }
    est.fun.CGF1<-function(end,start){
        n<-start-end+1L
        alpha<-1-m1^(1/n)
        Ac<-A[,end:start]%*%matrix(m1^((n-1):0/n),ncol=1L)*alpha*m1
        fit.theta(Ac,y)
    }
    est.fun.CGF2<-function(end,start){
        n<-start-end+1L
        alpha<-1-m2^(1/n)
        Ac<-A[,end:start]%*%matrix(m2^((n-1):0/n),ncol=1L)*alpha*m2
        fit.theta(Ac,y)
    }
    est.fun.GA<-function(end,start){
        n<-start-end+1L
        Ac<-A[,end:start]%*%matrix((1-1/n)^(0:(n-1))/c(rep(n,n-1),1),ncol=1L)*m1*m2
        fit.theta(Ac,y)
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
                    ssE<-theta0<-theta1<-rep(NA,T)
                    for(n in seq_len(T)){
                        coef<-est.fun(n)
                        theta0[n]<-coef$theta0
                        theta1[n]<-coef$theta1
                        ssE[n]<-coef$ssE
                    }
                    n<-which(ssE==min(ssE,na.rm=TRUE))[1L]
                    list(m=NA,n=n,start=n,end=NA,theta0=theta0[n],theta1=theta1[n],ssE=ssE[n],msE=ssE[n]/(length(y)-1))
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
                         ssE=gen.new$ssE,msE=gen.new$ssE/(length(y)-1))
                }
            }
        } else
            search<-function(model){
                est.fun<-est.funs[[model]]
                if(model==1L){
                    ssE<-theta0<-theta1<-rep(NA,T)
                    for(n in seq_len(T)){
                        coef<-est.fun(n)
                        theta0[n]<-coef$theta0
                        theta1[n]<-coef$theta1
                        ssE[n]<-coef$ssE
                    }
                    n<-which(ssE==min(ssE,na.rm=TRUE))[1L]
                    list(m=NA,n=n,start=n,end=NA,theta0=theta0[n],theta1=theta1[n],ssE=ssE[n],msE=ssE[n]/(length(y)-1))
                } else {
                    ssE<-theta0<-theta1<-matrix(nrow=T,ncol=T)

                    for(n in seq_len(max.duration)){
                        for(m in seq_len(T-n+1L)){
                            coef<-est.fun(m,m+n-1L)
                            theta0[n,m]<-coef$theta0
                            theta1[n,m]<-coef$theta1
                            ssE[n,m]<-coef$ssE
                        }
                    }

                    SSE<-Inf
                    for(n in seq_len(max.duration))
                        for(m in seq_len(T-n+1L))
                            if(ssE[n,m]<SSE){
                                M<-m;N<-n;SSE<-ssE[n,m]
                            }
                    list(m=M,n=N,start=M+N-1L,end=M,theta0=theta0[N,M],theta1=theta1[N,M],ssE=SSE,msE=SSE/(length(y)-1))
                }
            }
    } else
        search<-function(model){
            est.fun<-est.funs[[model]]
            ssE<-theta0<-theta1<-rep(NA,T)
            for(n in seq_len(T)){
                coef<-if(model==1L) est.fun(n) else est.fun(1L,n)
                theta0[n]<-coef$theta0
                theta1[n]<-coef$theta1
                ssE[n]<-coef$ssE
            }
            n<-which(ssE==min(ssE,na.rm=TRUE))[1]
            list(m=if(model==1L) NA else 1L,n=n,start=n,end=if(model==1L) NA else 1L,theta0=theta0[n],theta1=theta1[n],ssE=ssE[n],msE=ssE[n]/(length(y)-1))
        }

    if(single.parallel){
        require(doSNOW);require(foreach)
        cl<-makeCluster(single.clusternum)
        registerDoSNOW(cl)
        clusterExport(cl,c("distance","fit.theta"),envir=environment())
        estimate<-foreach(model=seq_len(4L))%dopar% search(model)
        stopCluster(cl)
    } else estimate<-lapply(seq_len(4L),search)

    names(estimate)<-if(isolation) c("HI","CGF1-I","CGF2-I","GA-I") else c("HI","CGF1","CGF2","GA")
    result$estimate<-estimate
    result$summary<-data.frame(Model=names(estimate),
                             Start=sapply(estimate,function(dummy) dummy$start),
                             End=sapply(estimate,function(dummy) dummy$end),
                             theta0=sapply(estimate,function(dummy) dummy$theta0),
                             theta1=sapply(estimate,function(dummy) dummy$theta1),
                             ssE=sapply(estimate,function(dummy) dummy$ssE),
                             Max.index=rep(result$maxindex,4),
                             msE=sapply(estimate,function(dummy) dummy$msE))

    result
}

#' Continuous Admixture Modeling (CAM)
#'
#' Estimate admixture time intervals/points for HI, CGF1(-I), CGF2(-I) and GA(-I) respectively for all Ld decay curves in a .rawld file.
#'
#' @param rawld a string representing the path of the .rawld file or a data frame read from the .rawld file by \code{read.table}.
#' @param T the most ancient generation to be searched. Defaults to 500.
#' @param m1 the admixture proportion of population 1. If m2 is the admixing proportion of population 2, then m1+m2=1.
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
#' Require \pkg{doSNOW}, \pkg{foreach} package and their dependencies if \code{LD.parallel=TRUE} or \code{single.parallel=TRUE}. It is recommended to library these required packages before using the parallel functionality.
#'
#' \code{max.duration} is only used when \code{isolation=TRUE} and \code{fast.search=FALSE}. The maximal duration of admixture \eqn{n} to be considered as possible. Smaller values can make the slow searching algorithm faster. If \code{max.duration>T}, it will be set to be \code{T}.
#'
#' \code{LD.clusternum} is used if \code{LD.parallel=TRUE}. If not specified, it is set to be the number of LD decay curves in the .rawld file.
#'
#' \code{single.parallel} indicates whether parallel computation should be used when computing a single LD decay curve. Defaults to \code{TRUE} if \code{isolation=TRUE,fast.search=FALSE} and \code{FALSE} otherwise.
#'
#' The .rawld file should include exactly one column nameed "Distance" in Morgan, exactly one column named "Combined_LD", several columns named "Jack?" representing Jackknives where ? is a number and exactly one column named "Fitted" representing the fitted LD decay curve using the previous method. This function fits "Combined_LD" and all Jackknives using all models. See \code{\link{singleCAM}} for further details of fitting algorithm for each LD decay curve.
#'
#' Be aware that when the computational cost is small (e.g. \code{isolation=FALSE} or \code{T=20L,isoaltion=TRUE,fast.search=FALSE,max.duration=10L}), using parallel computation for single LD decay curves can result in longer computation time.
#'
#' If the last entry of Distence in the .rawld file is greater than 10, a warning of unit will be given.
#' @examples
#' data(GA_I)
#' library(foreach);library(doSNOW)
#'
#' #fit models with isolation=FALSE.
#' fit<-CAM(GA_I,m1=0.3,T=150L,isolation=FALSE)
#' \dontrun{
#' plot(fit) #may not be able to display
#' }
#' #Bad fitting indicates isolation=TRUE should be tried.
#' fit<-CAM(GA_I,m1=0.3,T=150L,isolation=TRUE)
#' fit
#' fit$summary
#' \dontrun{
#' plot(fit) #may not be able to display
#' plot(fit,"D:/plot.pdf") #plot to a .pdf file
#' }
#'
#' data(CGF_50)
#' fit<-CAM(CGF_50,0.3,20L,isolation=FALSE,LD.parallel=FALSE)
#' fit #Lengths of intervals being 20 indicates larger T should be tried.
#'
#' \dontrun{
#' #passing a file path to the argument `rawld=`
#' fit<-CAM("CGF_50.rawld",0.3)
#' }
#' @seealso \code{\link{singleCAM}}, \code{\link{plot.CAM}}, \code{\link{construct.CAM}}
#' @import utils
#' @export
#'

CAM<-function(rawld,m1,T=500L,isolation=TRUE,
              fast.search=TRUE,max.duration=150L,
              LD.parallel=TRUE,LD.clusternum,
              single.parallel=isolation && !fast.search,
              single.clusternum=4L){
    if(is.character(rawld)) rawld<-read.table(rawld,header=TRUE)
    Jack.index<-grep("Jack",names(rawld))
    LD.index<-grep("Combined_LD",names(rawld))
    Y<-rawld[,c(LD.index,Jack.index)]
    d<-rawld$Distance

    if(d[length(d)]>10) warning("The unit of Genetic Distance might not be Morgan.")

    results<-list(call=match.call(),d=d,Y=Y,T=T,isolation=isolation,m1=m1,m2=1-m1)
    if(isolation){
        results$fast.search<-fast.search
        if(!fast.search){
            max.duration<-max(max.duration,T)
            results$max.duration<-max.duration
        }
    }
    class(results)<-"CAM"

    if(LD.parallel){
        if(missing(LD.clusternum)) LD.clusternum<-ncol(Y)
        require(doSNOW);require(foreach)
        cl<-makeCluster(LD.clusternum)
        registerDoSNOW(cl)
        clusterExport(cl,c("distance","fit.theta","singleCAM"),envir=environment())
        results$CAM.list<-foreach(ld=seq_len(ncol(Y)))%dopar%{
            singleCAM(d,Y[,ld],m1,T,isolation=isolation,fast.search=fast.search,max.duration=max.duration,single.parallel=single.parallel,single.clusternum=single.clusternum)
        }
        stopCluster(cl)
    } else results$CAM.list<-lapply(seq_len(ncol(Y)),function(ld){
        singleCAM(d,Y[,ld],m1,T,isolation=isolation,single.parallel=single.parallel,fast.search=fast.search,max.duration=max.duration,single.clusternum=single.clusternum)
    })
    names(results$CAM.list)<-c("Combined_LD",paste("Jack",seq_len(length(Jack.index)),sep=""))

    results$fitted<-rawld$Fitted
    if(results$CAM.list[[1]]$maxindex>1) results$fitted<-results$fitted[-seq_len(tempresults[[1]]$maxindex-1)]
    v<-distance(results$fitted,results$CAM.list[[1]]$y)

    data<-NULL
    for(ld in seq_len(ncol(Y))){
        data.temp<-data.frame(LD=rep(names(results$CAM.list)[ld],4),
                              Model=names(results$CAM.list[[ld]]$estimate),
                              Start=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$start),
                              End=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$end),
                              theta0=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$theta0),
                              theta1=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$theta1),
                              ssE=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$ssE),
                              Max.index=rep(results$CAM.list[[ld]]$maxindex,4),
                              msE=sapply(results$CAM.list[[ld]]$estimate,function(dummy) dummy$msE),
                              quasi.F=if(ld==1L) sapply(results$CAM.list[[1L]]$estimate,function(dummy){dummy$ssE/v}) else rep(NA,4))
        data<-rbind(data,data.temp)
    }
    data$LD<-as.factor(data$LD)
    data$Model<-factor(data$Model,ordered=TRUE)
    results$summary<-data
    results
}

#' Reconstruct the Fitted LD Decay Curves of Four Models
#'
#' Given an object of "CAM.single" class, reconstruct the fitted LD decay curves of the four models.
#'
#' @param CAM.single an object of "CAM.single" class
#' @return a list consisting of the four fitted curves
#' @examples
#' data(CGF_50)
#' y<-CGF_50$Combined_LD
#' d<-CGF_50$Distance
#'
#' fit<-singleCAM(d,y,m1=0.3,T=100L,isolation=FALSE)
#' fitted.curves<-reconstruct.fitted(fit)
#' @seealso \code{\link{singleCAM}}, \code{\link{construct.CAM}}
#' @export

reconstruct.fitted<-function(CAM.single){
    y1<-CAM.single$estimate[[1]]$theta0+CAM.single$estimate[[1]]$theta1*
        CAM.single$A[,CAM.single$estimate[[1]]$n]*CAM.single$m1*CAM.single$m2

    alpha<-1-CAM.single$m1^(1/CAM.single$estimate[[2]]$n)
    y2<-CAM.single$estimate[[2]]$theta0+CAM.single$estimate[[2]]$theta1*
        CAM.single$A[,CAM.single$estimate[[2]]$end:CAM.single$estimate[[2]]$start]%*%
        matrix(CAM.single$m1^((CAM.single$estimate[[2]]$n-1):0/CAM.single$estimate[[2]]$n),ncol=1)*alpha*CAM.single$m1
    y2<-as.numeric(y2)

    alpha<-1-CAM.single$m2^(1/CAM.single$estimate[[3]]$n)
    y3<-CAM.single$estimate[[3]]$theta0+CAM.single$estimate[[3]]$theta1*
        CAM.single$A[,CAM.single$estimate[[3]]$end:CAM.single$estimate[[3]]$start]%*%
        matrix(CAM.single$m2^((CAM.single$estimate[[3]]$n-1):0/CAM.single$estimate[[3]]$n),ncol=1)*alpha*CAM.single$m2
    y3<-as.numeric(y3)

    y4<-CAM.single$estimate[[4]]$theta0+CAM.single$estimate[[4]]$theta1*
        CAM.single$A[,CAM.single$estimate[[4]]$end:CAM.single$estimate[[4]]$start]%*%
        matrix((1-1/CAM.single$estimate[[4]]$n)^(0:(CAM.single$estimate[[4]]$n-1))/c(rep(CAM.single$estimate[[4]]$n,CAM.single$estimate[[4]]$n-1),1),ncol=1)*CAM.single$m1*CAM.single$m2
    y4<-as.numeric(y4)

    z<-list(y1,y2,y3,y4)
    names(z)<-paste(names(CAM.single$estimate),".fitted",sep="")
    z
}

#' Summary Plots for "CAM" Class
#'
#' @param x an object of "CAM" class
#' @param filename pdf file path.
#' @param model.cols a matrix of colors. See "Details"
#' @param box.log a logical expression. If \code{TRUE}, the scale of the axis will be in log scale. Defaults to \code{TRUE}
#' @param alpha a numeric value in [0,1] representitng alpha of the fitted LD decay curves for Combined LD. Defaults to 0.6. Can be a vector representing different alphas for different models. Ignored if alpha is specified in the third row of \code{model.cols}.
#' @param fit.lwd line width of the fitted LD decay curve for Combined LD of the four models. Defaults to 3. Can be a vector representing different line widths for different models.
#' @param LD.col color of the original Combined LD curve. Defaults to "black".
#' @param LD.lwd line width of the original Combined LD curve. Defaults to 1.
#' @param ... other graphical arguments passed to basic functions like \code{plot}.
#' @details
#'
#' The function generates three plots in a device. The plot on the top left is the estimated time intervals/points for the four models. The color depth of segments/points corresponds to how many intervals/points covers this part in Jackknives. The deeper the color, the more estimates from Jackknives cover this part. The plot on the top right is the boxplot of msE for the four models. The third plot shows the fitting of four models to \code{Combined_LD} in the .rawld file. The numbers after model names in the legend are quasi-F values of the four models for \code{Combined_LD}.
#'
#' If \code{filename} is set, plot to the .pdf file, otherwise plot to the current device. The function is specially designed for a .pdf plot with width being 9.6 and height being 7.2. To add things to the plot and then save it to a file, better to set the size as above. May not be able to plot directly in an R window.
#'
#' The colors in Column 1/2/3/4 correspond to the colors for HI/CGF1(-I)/CGF2(-I)/GA(-I). The first/second row of \code{model.cols} is the lightest/deepest possible color in the "Time Intervals/Points" plot. The third row of \code{model.cols} is the color for "msE Boxplot" and "Fitting of Models" plot. The colors will be converted to RGB colors by \code{col2rgb}.
#' @examples
#' \dontrun{
#' data(CGF_50)
#' fit<-CAM(CGF_50,0.3,10L,isolation=FALSE,LD.parallel=FALSE)
#' plot(fit,"bad_fitting.pdf")
#' 
#' #may not be able to display
#' plot(fit,model.cols=matrix(c("pink","red","pink",
#'                              "lightseagreen","green","green",
#'                              "skyblue","blue","blue",
#'                              "yellow","orange","orange"),ncol=4),
#'      box.log=FALSE,alpha=1,fit.lwd=1,LD.col="gray",LD.lwd=3)
#' }
#' @note
#' It is not recommended to pass other arguments to basic functions like \code{plot} in \code{...}.

#' @seealso \code{\link{CAM}}, \code{\link[grDevices]{rgb}}, \code{\link[grDevices]{col2rgb}}
#' @import graphics
#' @import grDevices
#' @export

plot.CAM<-function(x,filename,
                   model.cols=matrix(c("#ffa1c2B2","#b50d37B2","#da577c",
                                       "#9bff94","#0fbd02","#0fbd02",
                                       "#e9a1ff","#9111b8","#9111b8",
                                       "#7af8ff","#1ea1a8","#1ea1a8"),ncol=4),
                   box.log=TRUE,alpha=0.6,fit.lwd=3,LD.col="black",LD.lwd=1,...){
    model.cols2<-col2rgb(model.cols,alpha=TRUE)
    alpha<-alpha*rep(1,4)
    for(model in 1:4){
        if(nchar(model.cols[3,model])==7)
            model.cols2[4,3*model]<-round(alpha[model]*model.cols2[4,3*model])
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
    NJack<-length(levels(data$LD))-1

    if(!missing(filename)) pdf(filename,width=9.6,height=7.2)

    layout(matrix(c(1,2,3,3),ncol=2,nrow=2,byrow=TRUE),widths=c(4,4),heights=c(1,2))

    par(bty="n",las=1)
    if(max(sapply(intervals,function(dummy) max(dummy[1,])))<=50) xmax<-100
    else if(min(sapply(intervals,function(dummy) min(dummy[1,])))>=200) xmax<-500
    else xmax<-200
    plot(x=1:xmax,y=seq(0,4,length.out=xmax),type="n",ann=FALSE,axes=FALSE,...)

    for(model in 1:4){
        colors<-colorRampPalette(model.cols[1:2,model],alpha=TRUE)(NJack)
        dummy<-intervals[[model]]
        if(model==1){
            points(x=dummy[model,],y=rep(4.5-model,ncol(dummy)),col=colors[dummy[2,]],pch=15,cex=1.2,...)
        } else {
            if(max(dummy[1,])<=xmax){
                for(t in seq_len(ncol(dummy)))
                    lines(x=dummy[1,t]+c(-.5,.5),y=rep(4.5-model,2),col=colors[dummy[2,t]],lwd=5,...)
            } else {
                for(t in seq_len(xmax-min(dummy[1,])+1))
                    lines(x=dummy[1,t]+c(-.5,.5),y=rep(4.5-model,2),col=colors[dummy[2,t]],lwd=5,...)
                arrows(x0=xmax-5,y0=2.5,x1=xmax+2,lwd=3,col=colors[dummy[2,xmax-min(dummy[1,])+1]],angle=20,length=.2,...)
            }
        }
    }
    axis(1);axis(side=4,labels=NA,at=4:1-.5)
    title(main="Time Intervals/Points",xlab="Generation",ylab="")
    if(x$isolation){
        mtext(c("      HI"," CGF1-I"," CGF2-I","    GA-I"),side=4,line=1,at=4:1-.5)
    } else mtext(c("   HI","CGF1","CGF2","  GA"),side=4,line=2,at=4:1-.5)

    boxplot(data.jack$msE[data.jack$Model=="GA"|data.jack$Model=="GA-I"],
            data.jack$msE[data.jack$Model=="CGF2"|data.jack$Model=="CGF2-I"],
            data.jack$msE[data.jack$Model=="CGF1"|data.jack$Model=="CGF1-I"],
            data.jack$msE[data.jack$Model=="HI"],
            horizontal=TRUE,boxwex=.8,pch=20,
            log=if(box.log) "x" else "",medlwd=1,
            boxfill=rev(model.cols[3,]),
            outcol=rev(model.cols[3,]),
            main="msE Boxplot",xlab=if(box.log) "msE on log scale" else "msE",axes=FALSE,...)
    axis(2,labels=NA);axis(1)

    par(bty="n",las=1,mar=c(5,17,4,10.5)+.1)
    colors<-apply(model.cols2,2,function(col) rgb(col[1],col[2],col[3],col[4],maxColorValue=255))[seq(3,12,by=3)]
    fit.lwd<-fit.lwd*rep(1L,4)
    fitted<-reconstruct.fitted(x$CAM.list[[1]])
    plot(x$CAM.list[[1]]$d,x$CAM.list[[1]]$y,
         col=LD.col,type="l",lwd=LD.lwd,ylim=c(-.002,.25),
         xlab="Distance (Morgan)",ylab="",main="Fitting of Models",axes=FALSE,...)
    for(model in 1:4)
        lines(x$CAM.list[[1]]$d,fitted[[model]],col=colors[model],lwd=fit.lwd[model],...)
    legend(x="topright",
           legend=paste(names(x$CAM.list[[1]]$estimate)," (",round(x$summary$quasi.F[1:4],3),")",sep=""),
           col=colors,
           pch=-1,lty=1,lwd=2,
           ncol=1,bty="n")
    axis(2);axis(1);text(par("usr")[1]-0.15,y=.15,labels="Weighted LD",xpd=TRUE)

    if(!missing(filename)) dev.off()
}

#' Quickly Construct a Simple "CAM" Class Object from .rawld File and Summary Table
#'
#' Construct a simple "CAM" class object which can be passed to \code{\link{plot.CAM}} and whose elements in \code{CAM.list} can be passed to \code{\link{reconstruct.fitted}}. Can be used when only the .rawld file, \eqn{m_1} and the summary table is available (e.g., only the summary table was saved after running \code{\link{CAM}} last time).
#'
#' @param rawld original .rawld filepath or its data frame
#' @param m1 the admixture proportion of population 1.
#' @param dataset summary table as in \code{summary} of an object of "CAM" class
#' @return a simple "CAM" class object. It is not as complate as the "CAM" class object obtained from \code{\link{CAM}}. Particularly, it does not include the information about how the estimates are found.
#' @examples
#' library(foreach);library(doSNOW)
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
    if(is.character(rawld)) rawld<-read.table(rawld,header=TRUE)
    Jack.index<-grep("Jack",names(rawld))
    LD.index<-grep("Combined_LD",names(rawld))
    Y<-rawld[,c(LD.index,Jack.index)]
    d<-rawld$Distance

    if(d[length(d)]>10) warning("The unit of Genetic Distance might not be Morgan.")

    T<-max(dataset$Start)

    dataset$LD<-as.character(dataset$LD)
    dataset$Model<-as.character(dataset$Model)
    dataset$Model.num<-sapply(dataset$Model,switch,
                              HI=1,CGF1=2,`CGF1-I`=2,CGF2=3,`CGF2-I`=3,GA=4,`GA-I`=4)
    ordr<-order(dataset$LD,dataset$Model.num)
    dataset<-dataset[ordr,]
    dataset$Model.num<-NULL
    dataset$LD<-as.factor(dataset$LD)
    dataset$Model<-as.factor(dataset$Model)

    m2<-1-m1

    results<-list(d=d,Y=Y,isolation="CGF1-I" %in% levels(dataset$Model),CAM.list=NULL,fitted=rawld$Fitted,summary=dataset)
    class(results)<-"CAM"

    for(ld in seq_along(levels(dataset$LD))){
        d.temp<-d
        y.temp<-Y[,ld]
        data2<-dataset[(4*ld-3):(4*ld),]
        maxindex<-data2$Max.index[1]
        if(maxindex>1){
            y.temp<-y.temp[-seq_len(maxindex-1)]
            d.temp<-d.temp[-seq_len(maxindex-1)]
        }
        A<-exp(-d.temp%*%t(seq_len(T)))
        results$CAM.list[[ld]]<-list(maxindex=maxindex,d=d.temp,A=A,y=y.temp,m1=m1,m2=m2)
        for(model in 1:4){
            data.temp<-data2[model,]
            results$CAM.list[[ld]]$estimate[[model]]<-list(m=data.temp$End,n=if(model=="HI") data.temp$Start else data.temp$Start-data.temp$End+1,start=data.temp$Start,end=data.temp$End,theta0=data.temp$theta0,theta1=data.temp$theta1,ssE=data.temp$ssE,msE=data.temp$msE)
        }
        names(results$CAM.list[[ld]]$estimate)<-if(results$isolation) c("HI","CGF1-I","CGF2-I","GA-I") else c("HI","CGF1","CGF2","GA")
    }
    names(results$CAM.list)<-unique(dataset$LD)
    results
}

#' Print Method for "CAM.single" Class
#'
#' @param x a "CAM.single" class object
#' @details
#' Print a very brief summary of a "CAM.single" class object. Include:
#' \itemize{
#' \item Function call
#' \item length of used LD (excluding first few values if necessary)
#' \item a data frame containing the estimated time intervals/points and corresponding msE. The time point for HI model is stored in \code{Start} variable.
#' }

#' @seealso \code{\link{singleCAM}}
#' @export

print.CAM.single<-function(x){
    cat("Continuous Admixture Inference (CAM) for a Single LD Decay Curve\n\n")
    cat("Function call: ")
    print(x$call)
    cat("\n")
    cat("Length of Used LD:", length(x$y),"\n\n")
    print(x$summary[,c("Model","Start","End","msE")])
}

#' Print Method for "CAM" Class
#'
#' @param x a "CAM" class object
#' @details
#' The print method for "CAM" Class. Include:
#' \itemize{
#' \item Function call (if available from the object)
#' \item Total length of LD decay curve
#' \item a data frame containing the estimated intervals/points, msE and quasi-F for each model and each LD decay curve. The time point for HI model is stored in \code{Start} variable.
#' }
#' @seealso \code{\link{CAM}}
#' @export

print.CAM<-function(x){
    cat("Continuous Admixture Inference (CAM) for a .rawlf File\n\n")
    if(!is.null(x$call)){
        cat("Function call:")
        print(x$call)
        cat("\n")
    }
    cat("Total Length of LD:",length(x$d),"\n\n")
    print(x$summary[,c("LD","Model","Start","End","msE","quasi.F")])
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
#' \item{best.models}{the set of best models concluded}
#' \item{p.adjust.method}{method for adjusting p values used}
#' @details
#' The function uses pairwise Student's t-test on msE to select the best model(s). If HI model is not significantly worse than any other model, it is chosen as the best model; otherwise, the model(s) with significantly smallest msE are chosen as best model(s).
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
    if(log) data$msE<-log(data$msE)

    means<-tapply(data$msE,data$Model,mean)
    models<-names(means)

    p.value<-pairwise.t.test(data$msE,data$Model,p.adjust.method=p.adjust.method,paired=TRUE)$p.value
    p.value[is.na(p.value)]<-0
    p.value<-cbind(rbind(0,p.value),0)
    p.value<-p.value+t(p.value)
    p.value[matrix(as.logical(diag(4)),ncol=4)]<-rep(NA,4)
    colnames(p.value)<-row.names(p.value)<-models


    best<-which(means==min(means))[1]
    if(means[4]!=means[best] && p.value[4,best]<alpha){
        best<-c(best,which(p.value[best,]>=alpha))
        best<-models[sort(best)]
    } else best<-"HI"

    conclusion<-list(call=match.call(),alpha=alpha,group.means=means,adjusted.p.value=p.value,best.models=best,p.adjust.method=p.adjust.method)
    class(conclusion)<-"CAM.conclusion"
    conclusion
}

#' Print Method for "CAM.conclusion" Class
#'
#' @param x an object of "CAM.conclusion" class
#' @seealso \code{\link{conclude.model}}
#' @export

print.CAM.conclusion<-function(x){
    cat("CAM Best Model(s) Conclusion:\n\n")
    cat("Function call: ")
    print(x$call)
    cat("\n")
    cat("Familiwise Error Rate: ")
    cat(x$alpha,"\n\n",sep="")
    cat("Best Model(s): ")
    for(i in seq_along(x$best.models)){
        cat(x$best.models[i])
        if(i==length(x$best.models)) cat("\n\n")
        else cat(", ")
    }
    cat("Group Means of log(msE)/msE:\n")
    print(x$group.means)
    cat("\n")
    cat("Adjusted p-value:\n")
    print(x$adjusted.p.value)
}
