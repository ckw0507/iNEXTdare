#' Main Functions
#'
#' \code{iNEXT_dare}: compute the rarefaction and extrapolation of diveristy.
#' @param x a list consist of N data.frame/matrix describing species-by-assemblage/plot abundance. Note that
#' the species in each element must exactly match including specpes order. Use \code{data(abundata)} to see data example.
#' @param rho a numeric value of the fraction of sample size, which should be between 0 and 1.
#' @param q a numeric value or vector specifying the diversity order of Hill number.
#' @param datatype data type of input data: individual-based abundance data (datatype = "abundance"),
#' sampling-unit-based incidence frequencies data (datatype = "incidence_freq")
#' @param size an integer vector of sample sizes (number of individuals or sampling units) for which diversity estimates will be computed.
#' If NULL, then diversity estimates will be computed for those sample sizes determined by the specified/default endpoint and knots.
#' @param endpoint an integer specifying the sample size that is the endpoint for rarefaction/extrapolation.
#' If NULL, then endpoint = double reference sample size.
#' @param knots an integer specifying the number of equally-spaced knots (say K, default is 40) between size 1 and the endpoint;
#' each knot represents a particular sample size for which diversity estimate will be calculated. If the endpoint is smaller than the reference sample size,
#' @param se a logical variable to calculate the bootstrap standard error and conf confidence interval.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @param nboot an integer specifying the number of replications.
#' @return a list consists of N elements, where N is the number of region. Each element is a list containing 5 tables of gamma. alpha and beta diversities,
#' and 2 types of dissimialrities (Sorensen and Jaccard), respectively.
#' @examples
#' data(abundata)
#' out <- iNEXT_beta(x = abundata, q = c(0,1,2), datatype = "abundance",se = TRUE)
#' @import
#' @export


iNEXTdare <- function(x, rho, q = 0, datatype = "abundance",
                       endpoint = NULL, knots = 40, se = TRUE, conf = 0.95,
                       nboot = 50) {
  datainf.abu <- function(data,rho){
    if (class(data)!="dataframe") data <- as.data.frame(data)
    a1 <- matrix(0,15,1,dimnames=list(1:15, "value"))
    rownames(a1) <- c("n","N","rho", "S.obs", "C.hat","f1","f2","f3","f4","f5","f6","f7","f8","f9","f10")
    n = sum(data)
    N = ceiling(n/rho)
    a1[1,1] <- n
    a1[2,1] <- N
    a1[3,1] <- rho
    a1[4,1] <-  round(c(sum(data!=0)),0)
    a1[6:15,1] <- c(sum(data==1),sum(data==2),sum(data==3),sum(data==4),sum(data==5),sum(data==6),sum(data==7),sum(data==8),sum(data==9),sum(data==10))
    f1 = sum(data==1)
    # f2 = sum(data==2)
    a1[5,1] <- round(1 - (f1/n)*(N-n)/N,4)
    a1 = data.frame(a1)
    return(a1)
  }
  datainf.inc <- function(data,rho){
    if (class(data)!="dataframe") data <- as.data.frame(data)
    a1 <- matrix(0,16,1,dimnames=list(1:16, "value"))
    rownames(a1) <- c("T","T*","rho", "U", "S.obs", "C.hat","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")
    t <- data[1,1]
    ts <- ceiling(t/rho)
    a1[1,1] <- t
    data = data[-1,1]
    U <- sum(data)
    a1[4,1] <- round(U,0)
    a1[2,1] <- ts
    a1[3,1] <- rho
    a1[5,1] <- round(c(sum(data!=0)),0)
    a1[7:16,1] <- c(sum(data==1),sum(data==2),sum(data==3),sum(data==4),sum(data==5),sum(data==6),sum(data==7),sum(data==8),sum(data==9),sum(data==10))
    Q1 = sum(data==1)
    a1[6,1] <- round(1-(ts-t)/ts*Q1/U,4)
    a1 = data.frame(a1)
    return(a1)
  }
  bootstrap.inc <- function(data,rho){
    t <- data[1]
    ts <- ceiling(t/rho)
    data <- data[-1];data <- data[data>0]
    Q1 <- sum(data==1)
    Q2 <- sum(data==2)
    Q0 <- ifelse(Q2 != 0,(Q1^2)/(t/(t-1)*2*Q2+rho/(1-rho)*Q1),Q1*(Q1-1)/(t/(t-1)*2+rho/(1-rho)*Q1))
    chat <- ifelse(Q2 != 0,1-Q1/t*(1-rho)*(t-1)*Q1/((t-1)*Q1+2*Q2),1-Q1/t*(1-rho)*(t-1)*Q1/((t-1)*Q1+2))
    lam <- ifelse(rho == 1,(1-chat)/sum(data/t*(1-rho)^(data/rho)),0)
    Mi <- c(ceiling(data/rho*(1-lam*(1-rho)^(data/rho))),rep(ceiling(ts*(1-chat)/Q0),ceiling(Q0)))
    Mi[Mi > ts] <- ts
    dama <- matrix(0,length(Mi),ts)
    for (j in 1:length(Mi)) {
      samp <- sample(1:ts,Mi[j],replace = FALSE)
      dama[j,samp] <- 1
    }
    return(dama)
  }
  bootstrap.abu <- function(data,rho){
    data <- data[data>0]
    n <- sum(sapply(unique(data), function(x){ x* sum(data == x)}))
    N <- ceiling(n/rho)
    f1 <- sum(data==1)
    f2 <- sum(data==2)
    f0 <- ifelse(f2 != 0,(f1^2)/(n/(n-1)*2*f2+rho/(1-rho)*f1),f1*(f1-1)/(n/(n-1)*2+rho/(1-rho)*f1))
    chat <- ifelse(f2 != 0,1-f1/n*(1-rho)*(n-1)*f1/((n-1)*f1+2*f2),1-f1/n*(1-rho)*(n-1)*f1/((n-1)*f1+2))
    lam <- ifelse(rho == 1,(1-chat)/sum(data/n*(1-rho)^(data/rho)),0)
    Ni <- c(ceiling(data/rho*(1-lam*(1-rho)^(data/rho))),rep(ceiling(N*(1-chat)/f0),ceiling(f0)))
    temp <- unlist(sapply(1:length(Ni), function(x){
      rep(x,Ni[x])
    })
    )
    return(temp)
  }
  S.hat.abu <- function(data,rho) {
    data <- data[data>0]
    n <- sum(sapply( unique(data) , function(x) { x * sum(data == x)}))
    s <- length(data)
    N <- ceiling(n/rho)
    f1 <- sum(data==1)
    f2 <- sum(data==2)
    return(s+(f1^2)/(n/(n-1)*2*f2+rho/(1-rho)*f1))
  }
  S.hat.inc <- function(data,rho) {
    t <- data[1];data <- data[-1]
    data <- data[data>0]
    U <- sum(sapply( unique(data) , function(x) { x * sum(data == x)}))
    s <- length(data)
    ts <- ceiling(t/rho)
    q1 <- sum(data==1)
    q2 <- sum(data==2)
    return(s+(q1^2)/(t/(t-1)*2*q2+rho/(1-rho)*q1))
  }
  entropy.abu <- function(data,rho) {
    data <- data[data>0]
    n <- sum(sapply(unique(data), function(y) {y * sum(data == y)}))
    N <- ceiling(n/rho)
    f1 <- sum(data==1)
    f2 <- sum(data==2)
    r <- 1:(n-1)
    sor <- sort(unique(data))
    tab <- table(data)
    temp <- sapply(sor, function(z){
      k <- 1:(n-z)
      sum(1/k*(N-k)/N*z/n*exp(lchoose(n-z,k)-lchoose(n-1,k)))
    })
    ha <- sum(temp  * tab)
    A <- ifelse(f1 == 0 & f2 ==0, 0,((1-rho)*2*f2+rho*f1)/((n-1)*f1+2*f2))
    hb <- ifelse(A == 0,0,(N-n)/N*f1/n*(1-A)^(1-n)*(-log(A)-sum(1/r*(1-A)^r)))
    return(ha + hb)
  }
  entropy.inc <- function(data,rho) {
    t <- data[1]
    ts <- ceiling(t/rho)
    data <- data[-1];data <- data[data>0]
    U <- sum(sapply(unique(data), function(y){y * sum(data == y)}))
    Q1 <- sum(data==1)
    Q2 <- sum(data==2)
    r <- 1:(t-1)
    sor <- sort(unique(data))
    tab <- table(data)
    temp <- sapply(sor, function(z){
      # k <- ifelse(t==z,1,1:(t-z))
      k <- 1:(t-z);k <- k[k>0]
      sum(1/k*(ts-k)/ts*z/t*exp(lchoose(t-z,k)-lchoose(t-1,k)))
    })
    ha <- sum(temp  * tab)
    B <- ifelse(Q1 == 0 & Q2 ==0, 0,((1-rho)*2*Q2+rho*Q1)/((t-1)*Q1+2*Q2))
    hb <- ifelse(B == 0,0,(ts-t)/ts*Q1/t*(1-B)^(1-t)*(-log(B)-sum(1/r*(1-B)^r)))
    h0 <- ha + hb
    return(t/U*h0 + log(U/t))
  }
  simp.inc <- function(data,rho) {
    t <- data[1]
    ts <- ceiling(t/rho)
    data <- data[-1];data <- data[data>0]
    U <- sum(sapply(unique(data), function(y){y * sum(data == y)}))
    sor <- sort(unique(data))
    tab <- table(data)
    pi <- (sum(tab * sor*(sor-1))+rho*U)/(t*(t-1)+t*rho)
    return((U/t)^2/pi)
  }
  simp.abu <- function(data,rho) {
    data <- data[data>0]
    n <- sum(sapply(unique(data), function(y){y * sum(data == y)}))
    N <- ceiling(n/rho)
    sor <- sort(unique(data))
    tab <- table(data)
    pi <- (sum(tab * sor*(sor-1))+rho*n)/(n*(n-1)+n*rho)
    return(1/pi)
  }
  rarext.inc <- function(x,rho,et,q){
    Qk.hat <- function(y, nT, t) {
      y <- y[y > 0]
      Qj <- table(y)
      js <- as.numeric(names(Qj))
      if (t <= nT) {
        Sub <- function(k) {
          Qj_k <- Qj[js>=k]
          js_k <- js[js>=k]
          ifelse(length(js_k) == 0, 0, sum(exp(lchoose(js_k, k) + lchoose(nT - js_k, t - k) - lchoose(nT, t)) * Qj_k) )
        }
        sapply(1:t, Sub)
      }
      else {
        p.hat <- EstiBootComm.Sam(c(nT, y))
        Sub <- function(k) sum((choose(t, k) * p.hat^k *
                                  (1 - p.hat)^(t - k))/(1 - (1 - p.hat)^T))
        sapply(1:t, Sub)
      }
    }
    entropy.inc <- function(x,rho) {
      t <- x[1]
      yi <- x[-1];yi <- yi[yi>0]
      U <- sum(sapply(unique(yi), function(y){y * sum(yi == y)}))
      Q1 <- sum(yi==1)
      Q2 <- sum(yi==2)
      r <- 1:(t-1)
      sor <- sort(unique(yi))
      tab <- table(yi)
      temp <- sapply(sor, function(z){
        k <- 1:(t-z);k <- k[k>0]
        sum(1/k*(ts-k)/ts*z/t*exp(lchoose(t-z,k)-lchoose(t-1,k)))
      })
      ha <- sum(temp  * tab)
      ts <- ceiling(t/rho)
      B <- ifelse(Q1 == 0 & Q2 ==0, 0,((1-rho)*2*Q2+rho*Q1)/((t-1)*Q1+2*Q2))
      hb <- ifelse(B == 0,0,(ts-t)/ts*Q1/t*(1-B)^(1-t)*(-log(B)-sum(1/r*(1-B)^r)))
      h0 <- ha + hb
      t/U*h0 + log(U/t)
    }
    t <- x[1]
    ts <- ceiling(t/rho)
    x <- x[-1];x <- x[x>0]
    U <- sum(x)
    Q1 <- sum(x == 1)
    Q2 <- sum(x == 2)
    Q0 <- (Q1^2)/(t/(t-1)*2*Q2+rho/(1-rho)*Q1)
    M0 <- (1/rho-1)*Q1/Q0

    d0.hat <- function(y){
      if (y > t){
        y <- min(y,ts) - t
        return(sum(x != 0)+ Q0*(1-(1-y/(ts-t))^M0))
      }
      if (y <= t){
        Fun <- function(z) {
          if (z <= (t - y))
            exp(lgamma(t - z + 1) + lgamma(t - y + 1) - lgamma(t - z - y + 1) - lgamma(t + 1))
          else 0
        }
        return(sum(1 - sapply(x, Fun)))
      }
    }
    d1.hat <- function(y){
      if (y > t){
        y <- y - t
        return((exp(entropy.inc(c(t,x),rho))-exp(sum(-x/sum(x)*log(x/sum(x)))))*(1-(1-y/(ts-t))^M0)+exp(sum(-x/sum(x)*log(x/sum(x)))))
      }
      if (y <= t){
        k <- 1:y
        Ut.hat <- y/t * U
        return(exp(-sum(k/Ut.hat * log(k/Ut.hat) * Qk.hat(x, t, y))))
      }
    }
    d2.hat <- function(y){
      if (y > t){
        y <- y - t
        pi <- (sum(x*(x-1))+U*rho)/(t*(t-1)+t*rho)
        fra <- (t/(t+y)/U)^2*((1-(t+y)/ts)*(t+y)*U/t-(1-(t+y)/ts)*(t+y)*pi+(t+y)^2*pi)
        return(1/fra)
      }
      if (y <= t){
        Sub <- function(z) {
          1/(1/z * t/U + (1 - 1/z) * sum(x * (x - 1)/U^2/(1 - 1/t)))
        }
        return(Sub(y))
      }
    }
    if (q == 0) {
      sapply(et, d0.hat)
    }  else if (q == 1) {
      sapply(et, d1.hat)
    } else if (q == 2) {
      sapply(et, d2.hat)
    }
  }
  rarext.abu <- function(x,rho,m,q){
    fk.hat <- function(x, m) {
      x <- x[x > 0]
      n <- sum(x)
      fj <- table(x)
      js <- as.numeric(names(fj))
      if (m <= n) {
        Sub <- function(k) {
          fj_k <- fj[js>=k]
          js_k <- js[js>=k]
          ifelse(length(js_k) == 0, 0, sum(exp(lchoose(js_k, k) + lchoose(n - js_k, m - k) - lchoose(n, m)) * fj_k) )
        }
        sapply(1:m, Sub)
      }
      else {
        # f1 <- sum(x == 1)
        # f2 <- sum(x == 2)
        # A <- ifelse(f2 > 0, (n - 1) * f1/((n - 1) * f1 +
        #                                     2 * f2), (n - 1) * f1/((n - 1) * f1 + 2))
        # C.hat <- 1 - f1/n * A
        # p.hat <- x/n * C.hat
        p.hat = DetAbu(x)
        Sub <- function(k) sum((choose(m, k) * p.hat^k *
                                  (1 - p.hat)^(m - k))/(1 - (1 - p.hat)^n))
        sapply(1:m, Sub)
      }
    }
    entropy.abu <- function(x,rho) {
      x <- x[x>0]
      n <- sum(sapply(unique(x), function(y) {y * sum(x == y)}))
      f1 <- sum(x==1)
      f2 <- sum(x==2)
      r <- 1:(n-1)
      sor <- sort(unique(x))
      tab <- table(x)
      temp <- sapply(sor, function(z){
        k <- 1:(n-z)
        sum(1/k*(N-k)/N*z/n*exp(lchoose(n-z,k)-lchoose(n-1,k)))
      })
      ha <- sum(temp  * tab)
      N <- ceiling(n/rho)
      A <- ifelse(f1 == 0 & f2 ==0, 0,((1-rho)*2*f2+rho*f1)/((n-1)*f1+2*f2))
      hb <- ifelse(A == 0,0,(N-n)/N*f1/n*(1-A)^(1-n)*(-log(A)-sum(1/r*(1-A)^r)))
      return(ha + hb)
    }
    n <- sum(x)
    N <- ceiling(n/rho)
    x <- x[x>0]
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    f0 <- (f1^2)/(n/(n-1)*2*f2+rho/(1-rho)*f1)
    N0 <- ifelse(f1 != 0,(1/rho-1)*f1/f0,0)
    d0.hat <- function(m){
      if (m > n){
        m <- m - n
        return(sum(x != 0)+ N0*(1-(1-m/(N-n))^N0))
      }
      if (m <= n) {
        Fun <- function(x) {
          if (x <= (n - m))
            exp(lgamma(n - x + 1) + lgamma(n - m + 1) -
                  lgamma(n - x - m + 1) - lgamma(n + 1))
          else 0
        }
        sum(1 - sapply(x, Fun))
      }
    }
    d1.hat <- function(m){
      if (m > n){
        m <- m - n
        return((exp(entropy.abu(x,rho))-exp(sum(-x/n*log(x/n))))*(1-(1-m/(N-n))^N0)+exp(sum(-x/n*log(x/n))))
      }
      if (m <= n){
        k <- 1:m
        return(exp(-sum(k/m * log(k/m) * fk.hat(x, m))))
      }
    }
    d2.hat <- function(m){
      if (m > n) {
        # m <- m - n
        pi <- (sum(x*(x-1))+n*rho)/(n*(n-1)+n*rho)
        fra <- (1-m/N)/m+(m+m/N-1)/m*pi
        return(1/fra)
      }
      if (m <= n) {
        Sub <- function(m) {1/(1/m + (1 - 1/m) * sum(x * (x - 1)/n/(n - 1)))}
        return(Sub(m))
      }
    }
    if (q == 0) {
      sapply(m, d0.hat)
    }  else if (q == 1) {
      sapply(m, d1.hat)
    } else if (q == 2) {
      sapply(m, d2.hat)
    }
  }
  cover.abu <- function(data,rho,em) {
    data <- data[data>0]
    n <- sum(sapply( unique(data) , function(x) { x * sum(data == x)}))
    s <- length(data)
    N <- ceiling(n/rho)
    f1 <- sum(data == 1)
    f2 <- sum(data == 2)
    sor <- sort(unique(data))
    tab <- table(data)
    N1 <- 1+(N-n)*2*f2/f1/(n-1)
    # A <- ifelse(f2 != 0,2*f2/((n-1) * f1 + 2 * f2),ifelse(f1 != 0,2/((n-1)*(f1 - 1)+2),1))
    sub <- function(y){
      if (y >= N) {return(1)}
      if (y >= n & y < N) {
        y <- y - n
        return( 1-(N-n)/N*f1/n*(1-y/(N-n))^N1)
      }
      if (y < n) {
        return(1-(N-y)/N*sum(tab * sapply(sor,function(x) exp(lchoose(n-x,y)-lchoose(n-1,y))*x/n)))
      }
    }
    sapply(em, sub)
  }
  cover.inc <- function(data,rho,et) {
    t <- data[1]
    data <- data[-1];data <- data[data>0]
    ts <- ceiling(t/rho)
    U <- sum(sapply(unique(data), function(x){x * sum(data == x)}))
    Q1 <- sum(data == 1)
    Q2 <- sum(data == 2)
    sor <- sort(unique(data))
    tab <- table(data)
    M1 <- 1+(ts-t)*2*Q2/Q1/(t-1)
    # A <- ifelse(Q2 != 0,2*Q2/((t-1) * Q1 + 2 * Q2),ifelse(Q1 != 0,2/((t-1)*(Q1 - 1)+2),1))
    sub <- function(y){
      if (y >= ts) {return(1)}
      if (y >= t & y < ts) {
        y <- y - t
        return( 1-(ts-t)/ts*Q1/U*(1-y/(ts-t))^M1)
      }
      if (y < t) {
        return(1-(ts-y)/ts*sum(tab * sapply(sor,function(x) exp(lchoose(t-x,y)-lchoose(t-1,y))*x/U)))
      }
    }
    sapply(et, sub)
  }
  CI <- round(qnorm(0.5+conf/2),3)
  if (datatype == "incidence_freq") {
    if (is.null(names(x)) == 1) {
      names(x) <- paste("Site",1:length(x),sep = ".")
    }
    Site <- names(x)
    DataInfo <- c()
    for (i in 1:length(x)) {
      temp <- c()
      temp <- cbind(Site[i],t(datainf.inc(x[[i]],rho)))
      colnames(temp)[1] <- "Site"
      DataInfo <- rbind(DataInfo,temp)
    }
    rownames(DataInfo) <- 1:length(x)
    DataInfo <- as.data.frame(DataInfo)
    if (se == F | nboot == 0) {
      if (is.null(endpoint) == 1) {
        est <- lapply(1:length(x),function(z){
          da <- x[[z]]
          st <- da[1]
          sts <- ceiling(st/rho)
          et <- sort(unique(c(round(seq(1,min(sts,2*st),length.out = knots)),st)))
          sc <- round(cover.inc(da,rho,et),3)
          method <- ifelse(et<st,"interpolated",ifelse(et>st,"extrapolated","observed"))
          temp <- c()
          for (j in 1:length(q)) {
            qD <- round(rarext.inc(da,rho,et,q[j]),3)
            a1 <- data.frame(t=et,method,order = q[j],qD = qD, SC =sc)
            temp <- rbind(temp,a1)
          }
          temp
        })
      }
      if (is.null(endpoint) == 0) {
        est <- lapply(1:length(x),function(z){
          da <- x[[z]]
          st <- da[1]
          sts <- ceiling(st/rho)
          et <- sort(unique(c(round(seq(1,min(sts,endpoint),length.out = knots)),st)))
          sc <- round(cover.inc(da,rho,et),3)
          method <- ifelse(et<st,"interpolated",ifelse(et>st,"extrapolated","observed"))
          temp <- c()
          for (j in 1:length(q)) {
            qD <- round(rarext.inc(da,rho,et,q[z]),3)
            a1 <- data.frame(t=et,method,order = q[z],qD = qD, SC =sc)
            temp <- rbind(temp,a1)
          }
          temp
        })
      }
      names(est) <- Site
      out$iNextEst <- est
    } else if (se == T) {
      if (is.null(endpoint) == 1) {
        est <- lapply(1:length(x),function(z){
          da <- x[[z]]
          st <- da[1]
          sts <- ceiling(st/rho)
          et <- sort(unique(c(round(seq(1,min(sts,2*st),length.out = knots)),st)))
          sc <- cover.inc(da,rho,et)
          method <- ifelse(et<st,"interpolated",ifelse(et>st,"extrapolated","observed"))
          MI <- bootstrap.inc(da,rho)
          sam <- sapply(1:nboot, function(x){
            c(st,rowSums(MI[,sample(1:ncol(MI),st,replace = F)]))
          })
          scse <-  sapply(1:nboot, function(x){
            da <- sam[,x]
            cover.inc(da,rho,et)
          })
          scsd <- sapply(1:length(et), function(z) {sd(scse[z,])})
          temp <- c()
          for (j in 1:length(q)) {
            qD <- rarext.inc(da,rho,et,q[j])
            qDse <- sapply(1:nboot, function(x){
              da <- sam[,x]
              rarext.inc(da,rho,et,q[j])
            })
            qDsd <- sapply(1:length(et), function(z) {sd(qDse[z,])})
            a1 <- data.frame(t=et,method,order = q[j],qD = qD,
                             qD.LCL = qD - CI*qDsd,
                             qD.UCL = qD + CI*qDsd, SC =sc,
                             SC.LCL = sc - CI*scsd,
                             SC.UCL = sc + CI*scsd)
            a1[,4:9] <- round(a1[,4:9],3)
            temp <- rbind(temp,a1)
          }
          temp
        })
      }
      if (is.null(endpoint) == 0) {
        est <- lapply(1:length(x),function(z){
          da <- x[[z]]
          st <- da[1]
          sts <- ceiling(st/rho)
          et <- sort(unique(c(round(seq(1,min(sts,endpoint),length.out = knots)),st)))
          sc <- cover.inc(da,rho,et)
          method <- ifelse(et<st,"interpolated",ifelse(et>st,"extrapolated","observed"))
          MI <- bootstrap.inc(da,rho)
          sam <- sapply(1:nboot, function(x){
            c(st,rowSums(MI[,sample(1:ncol(MI),st,replace = F)]))
          })
          scse <-  sapply(1:nboot, function(x){
            da <- sam[,x]
            cover.inc(da,rho,et)
          })
          scsd <- sapply(1:length(et), function(z) {sd(scse[z,])})
          temp <- c()
          for (j in 1:length(q)) {
            qD <- rarext.inc(da,rho,et,q[j])
            qDse <- sapply(1:nboot, function(x){
              da <- sam[,x]
              rarext.inc(da,rho,et,q[j])
            })
            qDsd <- sapply(1:length(et), function(z) {sd(qDse[z,])})
            a1 <- data.frame(t=et,method,order = q[j],qD = qD,
                             qD.LCL = qD - CI*qDsd,
                             qD.UCL = qD + CI*qDsd, SC =sc,
                             SC.LCL = sc - CI*scsd,
                             SC.UCL = sc + CI*scsd)
            a1[,4:9] <- round(a1[,4:9],3)
            temp <- rbind(temp,a1)
          }
          temp
        })
      }
      names(est) <- Site

    }
    asy <- c()
    for (i in 1:length(x)) {
      da <- x[[i]]
      ob <- da[-1];ob <- ob[ob>0];st <- da[1]
      temp <- matrix(0,3,7)
      temp <- as.data.frame(temp)
      colnames(temp) <- c("Site","Diversity","Observed","Estimator","s.e.","LCL","UCL")
      temp[,1] <- Site[i]
      temp[,2] <- c("Species richness","Shannon diversity","Simpson diversity")
      temp[,3] <- round(c(sum(ob != 0),exp(-sum(ob/sum(ob)*log(ob/sum(ob)))),sum((ob/sum(ob))^2)^(-1)),3)
      temp[,4] <- round(c(S.hat.inc(da,rho),exp(entropy.inc(da,rho)),simp.inc(da,rho)),3)
      MI <- bootstrap.inc(da,rho)
      sam <- sapply(1:nboot, function(x){
        c(st,rowSums(MI[,sample(1:ncol(MI),st,replace = F)]))
      })
      se <-  sapply(1:nboot, function(x){
        da <- sam[,x]
        c(S.hat.inc(da,rho),exp(entropy.inc(da,rho)),simp.inc(da,rho))
      })
      temp[,5] <- round(c(sd(se[1,]),sd(se[2,]),sd(se[3,])),3)
      temp[,6] <- temp[,4] - qnorm(0.5+conf/2)*temp[,5]
      temp[,7] <- temp[,4] + qnorm(0.5+conf/2)*temp[,5]
      asy <- rbind(asy,temp)
    }
    out <- list("DateInfo" = DataInfo, "iNextEst" = est, "AsyEst" = asy)

  }
  if (datatype == "abundance") {
    if (is.null(names(x)) == 1) {
      names(x) <- paste("Site.",1:length(x))
    }
    Site <- names(x)
    DataInfo <- c()
    for (i in 1:length(x)) {
      temp <- c()
      temp <- cbind(Site[i],t(datainf.abu(x[[i]],rho)))
      colnames(temp)[1] <- "Site"
      DataInfo <- rbind(DataInfo,temp)
    }
    rownames(DataInfo) <- 1:length(x)
    DataInfo <- as.data.frame(DataInfo)
    if (se == F | nboot == 0) {
      if (is.null(endpoint) == 1) {
        est <- lapply(1:length(x),function(z){
          da <- x[[z]]
          n <- sum(da)
          N <- ceiling(n/rho)
          em <- sort(unique(c(round(seq(1,min(N,2*n),length.out = knots)),n)))
          sc <- round(cover.abu(da,rho,em),3)
          method <- ifelse(em<n,"interpolated",ifelse(em>n,"extrapolated","observed"))
          temp <- c()
          for (j in 1:length(q)) {
            qD <- round(rarext.abu(da,rho,em,q[j]),3)
            a1 <- data.frame(n=em,method,order = q[j],qD = qD, SC =sc)
            temp <- rbind(temp,a1)
          }
          temp
        })
      }
      if (is.null(endpoint) == 0) {
        est <- lapply(1:length(x),function(z){
          da <- x[[z]]
          n <- sum(da)
          N <- ceiling(n/rho)
          em <- sort(unique(c(round(seq(1,min(N,endpoint),length.out = knots)),n)))
          sc <- round(cover.inc(da,rho,em),3)
          method <- ifelse(em<n,"interpolated",ifelse(em>n,"extrapolated","observed"))
          temp <- c()
          for (j in 1:length(q)) {
            qD <- round(rarext.abu(da,rho,em,q[z]),3)
            a1 <- data.frame(n=et,method,order = q[z],qD = qD, SC =sc)
            temp <- rbind(temp,a1)
          }
          temp
        })
      }
      names(est) <- Site
      out$iNextEst <- est
    } else if (se == T) {

      if (is.null(endpoint) == 1) {
        est <- lapply(1:length(x),function(z){
          da <- x[[z]]
          n <- sum(da)
          N <- ceiling(n/rho)
          em <- sort(unique(c(round(seq(1,min(N,2*n),length.out = knots)),n)))
          sc <- cover.abu(da,rho,em)
          method <- ifelse(em<n,"interpolated",ifelse(em>n,"extrapolated","observed"))
          NI <- bootstrap.abu(da,rho)
          sam <- sapply(1:nboot, function(x){
            table(NI[sample(1:length(NI),n,replace = F)])
          })

          scse <-  sapply(1:nboot, function(x){
            da <- sam[[x]]
            cover.abu(da,rho,em)
          })

          scsd <- sapply(1:length(em), function(z) {sd(scse[z,])})
          temp <- c()
          for (j in 1:length(q)) {
            qD <- rarext.abu(da,rho,em,q[j])
            qDse <- sapply(1:nboot, function(x){
              da <- sam[[x]]
              rarext.inc(da,rho,em,q[j])
            })
            qDsd <- round(sapply(1:length(em), function(z) {sd(qDse[z,])}),3)
            a1 <- data.frame(n=em,method,order = q[j],qD = qD,
                             qD.LCL = qD - CI*qDsd,
                             qD.UCL = qD + CI*qDsd, SC =sc,
                             SC.LCL = sc - CI*scsd,
                             SC.UCL = sc + CI*scsd)
            a1[,4:9] <- round(a1[,4:9],3)
            temp <- rbind(temp,a1)
          }
          temp
        })
      }
      if (is.null(endpoint) == 0) {
        est <- lapply(1:length(x),function(z){
          da <- x[[z]]
          n <- sum(da)
          N <- ceiling(n/rho)
          em <- sort(unique(c(round(seq(1,min(N,endpoint),length.out = knots)),n)))
          sc <- cover.abu(da,rho,em)
          method <- ifelse(em<n,"interpolated",ifelse(em>n,"extrapolated","observed"))
          NI <- bootstrap.abu(da,rho)
          sam <- sapply(1:nboot, function(x){
            table(NI[sample(1:length(NI),n,replace = F)])
          })

          scse <-  sapply(1:nboot, function(x){
            da <- sam[[x]]
            cover.abu(da,rho,em)
          })

          scsd <- sapply(1:length(em), function(z) {sd(scse[z,])})
          temp <- c()
          for (j in 1:length(q)) {
            qD <- rarext.abu(da,rho,em,q[j])
            qDse <- sapply(1:nboot, function(x){
              da <- sam[[x]]
              rarext.inc(da,rho,em,q[j])
            })
            qDsd <- round(sapply(1:length(em), function(z) {sd(qDse[z,])}),3)
            a1 <- data.frame(n=em,method,order = q[j],qD = qD,
                             qD.LCL = qD - CI*qDsd,
                             qD.UCL = qD + CI*qDsd, SC =sc,
                             SC.LCL = sc - CI*scsd,
                             SC.UCL = sc + CI*scsd)
            a1[,4:9] <- round(a1[,4:9],3)
            temp <- rbind(temp,a1)
          }
          temp
        })
      }
      names(est) <- Site

    }
    asy <- c()
    for (i in 1:length(x)) {
      da <- x[[i]]
      da <- da[da>0]
      n <- sum(sapply(unique(da), function(x){ x* sum(da == x)}))
      temp <- matrix(0,3,7)
      temp <- as.data.frame(temp)
      colnames(temp) <- c("Site","Diversity","Observed","Estimator","s.e.","LCL","UCL")
      temp[,1] <- Site[i]
      temp[,2] <- c("Species richness","Shannon diversity","Simpson diversity")
      temp[,3] <- round(c(sum(da != 0),exp(-sum(da/sum(da)*log(da/sum(da)))),sum((da/sum(da))^2)^(-1)),3)
      temp[,4] <- round(c(S.hat.abu(da,rho),exp(entropy.abu(da,rho)),simp.abu(da,rho)),4)
      NI <- bootstrap.abu(da,rho)
      sam <- sapply(1:nboot, function(x){
        table(NI[sample(1:length(NI),n,replace = F)])
      })
      se <-  sapply(1:nboot, function(x){
        da <- sam[[x]]
        c(S.hat.abu(da,rho),exp(entropy.abu(da,rho)),simp.abu(da,rho))
      })
      temp[,5] <- round(c(sd(se[1,]),sd(se[2,]),sd(se[3,])),3)
      temp[,6] <- temp[,4] - CI*temp[,5]
      temp[,7] <- temp[,4] + CI*temp[,5]
      asy <- rbind(asy,temp)
    }
    out <- list("DataInfo" = DataInfo, "iNextEst" = est, "AsyEst" = asy)
  }
  class(out) <- "iNEXT"
  return(out)
}
