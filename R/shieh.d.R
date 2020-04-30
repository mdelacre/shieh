#' shieh.d function to compute Shieh's d and its confidence intervals
#'
#' @param x a (non-empty) numeric vector of data values (first sample)
#' @param y a (non-empty) numeric vector of data values (second sample)
#' @param conf.level confidence level of the interval
#' @param na.rm set whether Missing Values should be excluded (na.rm = TRUE) or not (na.rm = FALSE) - defaults to TRUE
#' @param alternative a character string describing the alternative hypothesis (two.sided, less or greater)
#'
#' @export shieh.d
#' @S3method shieh.d default
#' @S3method print shieh.d
#'
#' @keywords Two means comparison, effect size, Shieh's d, confidence interval
#' @return Returns Call, Shieh's d, confidence interval, sample size and sample estimates
#' @examples
#' #### Run Shieh.d
#' {# The default is to use a 95% confidence interval around point estimate
#' x <- c(1,4,3,5)
#' y <- c(1,9,4,5,0,2)
#' res <- shieh.d(x,y)
#' res
#'
#' # Moreover, a list of elements can be extracted from the function,
#' # such as the position of outliers in the dataset
#' # and the coordinates of outliers
#' res$statistic
#' res$parameter
#' res$estimate
#' res$conf.int
#'
#' @importFrom stats na.omit t.test uniroot sd pt
#'
# Create a generic function
shieh.d <- function(x, y,conf.level,alternative, na.rm) UseMethod("shieh.d")

shieh.dEst <- function(x,
                    y,
                    conf.level=.95,
                    alternative="two.sided",
                    na.rm=TRUE)
{

  if (na.rm == TRUE ) {
    Group.1 <- na.omit(x)
    Group.2 <- na.omit(y)
  } else {
    Group.1 <- x
    Group.2 <- y
  }

    if(inherits(Group.1,c("numeric","integer")) == FALSE)
      stop("x should be either numeric or integer")

    if(inherits(Group.2,c("numeric","integer")) == FALSE)
      stop("y should be either numeric or integer")

  n1 <- length(Group.1)
  n2 <- length(Group.2)
  s1 <- sd(Group.1) # standard deviation of the first group
  s2 <- sd(Group.2) # standard deviation of the second
  m1 <- mean(Group.1)
  m2 <- mean(Group.2)

  # Computing Shieh's d value
  N <- n1+n2
  q1 <- n1/N
  q2 <- n2/N
  shieh_d <- (m1-m2)/sqrt(s1^2/q1+s2^2/q2)

  # Computing confidence limits

  if(alternative=="two.sided"){

    #  perform two-sample Welch t-test (same assumptions of Shieh's d)
    test <- t.test(Group.1, Group.2, alternative = "two.sided", var.equal = FALSE)
    w_obs <- test$statistic
    #  sample estimates for degrees of freedom DF of noncentral t distribution
    DF <- test$parameter

    # lower limit = limit of lambda such as 1-pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level)/2 = alpha/2
    # with DF = (s1^2/n1 + s2^2/n2)^2 / ((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))

    f=function(lambda,rep) 1-pt(q=w_obs, df=DF, ncp = lambda)-rep
    out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")

    lambda.1 <- out$root
    delta.1 <- lambda.1/sqrt(n1+n2)

    # upper limit = limit of lambda such as pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level)/2 = alpha/2
    # with DF = (s1^2/n1 + s2^2/n2)^2 / ((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))

    f=function(lambda,rep) pt(q=w_obs, df=DF, ncp = lambda)-rep
    out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
    lambda.2 <- out$root
    delta.2 <- lambda.2/sqrt(n1+n2)

    result <- c(delta.1, delta.2)

  } else if (alternative == "greater"){

    #  perform two-sample Welch t-test (same assumptions of Shieh's d)
    test <- t.test(Group.1, Group.2, alternative = "greater", var.equal = FALSE)
    w_obs <- test$statistic
    #  sample estimates for degrees of freedom DF of noncentral t distribution
    DF <- test$parameter

    # lower limit = limit of lambda such as 1-pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
    # with DF = (s1^2/n1 + s2^2/n2)^2 / ((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))

    f=function(lambda,rep) 1-pt(q=w_obs, df=DF, ncp = lambda)-rep
    out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")

    lambda.1 <- out$root
    delta.1 <- lambda.1/sqrt(n1+n2)

    # upper limit = limit of lambda such as pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
    # with DF = (s1^2/n1 + s2^2/n2)^2 / ((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))

    delta.2 <- +Inf

    result <- c(delta.1, delta.2)

  } else if (alternative == "less"){

    #  perform two-sample Welch t-test (same assumptions of Shieh's d)
    test <- t.test(Group.1, Group.2, alternative = "less", var.equal = FALSE)
    w_obs <- test$statistic
    #  sample estimates for degrees of freedom DF of noncentral t distribution
    DF <- test$parameter

    # lower limit = limit of lambda such as 1-pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
    # with DF = (s1^2/n1 + s2^2/n2)^2 / ((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))

    delta.1 <- -Inf

    # upper limit = limit of lambda such as pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
    # with DF = (s1^2/n1 + s2^2/n2)^2 / ((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))

    f=function(lambda,rep) pt(q=w_obs, df=DF, ncp = lambda)-rep
    out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
    lambda.2 <- out$root
    delta.2 <- lambda.2/sqrt(n1+n2)

    result <- c(delta.1, delta.2)

  }

  # Return results in list()
  invisible(
    list(statistic = shieh_d,
         parameter = c(n1,n2),
         estimate = c(m1,m2,s1,s2),
         conf.int = result)
  )

}

# Adding a default method in defining a function called shieh.d.default

shieh.d.default <- function(x, y,conf.level=.95,alternative="two.sided", na.rm=TRUE){
  out <- shieh.dEst(x, y,conf.level,alternative, na.rm)
  out$statistic <- out$statistic
  out$parameter <- out$parameter
  out$estimate <- out$estimate
  out$conf.int <- out$conf.int
  out$call <- match.call()

  class(out) <- "shieh.d"
  out
}

print.shieh.d <- function(x,...){
  cat("\n             Shieh's d\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\n Shieh's d :\n")
  print(round(x$statistic,3))

  cat("\n",paste0(x$conf.level*100," percent confidence interval:"),"\n")
  print(round(x$conf.int,3))

  cat("\n Sample size :\n")
  cat("      n1","n2","\n")
  print(x$parameter)

  cat("\n Sample estimates :\n")
  cat(c("    mean of x","mean of y","\n"))
  print(x$estimate[1:2])
  cat(c("     sd of x","sd of y","\n"))
  print(x$estimate[3:4])


}




