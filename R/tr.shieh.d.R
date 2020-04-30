#' tr.shieh.d function to compute transformed Shieh's d and its confidence intervals
#'
#' @param x a (non-empty) numeric vector of data values (first sample)
#' @param y a (non-empty) numeric vector of data values (second sample)
#' @param conf.level confidence level of the interval
#' @param na.rm set whether Missing Values should be excluded (na.rm = TRUE) or not (na.rm = FALSE) - defaults to TRUE
#' @param alternative a character string describing the alternative hypothesis (two.sided, less or greater)
#'
#' @export tr.shieh.d
#' @S3method tr.shieh.d default
#' @S3method print tr.shieh.d
#'
#' @keywords Two means comparison, effect size, Shieh's d, confidence interval
#' @return Returns Call, Shieh's d, confidence interval, sample size and sample estimates
#' @examples
#' #### Run tr.shieh.d
#' {# The default is to use a 95% confidence interval around point estimate
#' x <- c(1,4,3,5)
#' y <- c(1,9,4,5,0,2)
#' res <- tr.shieh.d(x,y)
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
tr.shieh.d <- function(x, y,conf.level,alternative, na.rm) UseMethod("shieh.d")

tr.shieh.dEst <- function(x,
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

  # Computing transformed Shieh's d value
  N <- n1+n2
  q1 <- n1/N
  q2 <- n2/N
  nratio <- n1/n2
  sigma_bal <- sqrt((s1^2+s2^2)/2)
  sigma_unbal <- sqrt((1-q1)*s1^2+(1-q2)*s2^2) # we give more weight to the variance of the smallest group

  shieh_d_corr <- (m1-m2)/sqrt(s1^2/q1+s2^2/q2)*(((nratio+1)*sigma_unbal)/(2*sigma_bal*sqrt(nratio))) # what value of Shieh's delta would be obtain if n1=n2?

  # Computing confidence limits
  orig.ci <- shieh.d(x,y,conf.level=conf.level,alternative=alternative,na.rm=na.rm)
  corr_shieh_lim <- orig.ci$conf.int*(((nratio+1)*sigma_unbal)/(2*sigma_bal*sqrt(nratio)))

  # Return results in list()
  invisible(
    list(statistic = shieh_d_corr,
         parameter = c(n1,n2),
         estimate = c(m1,m2,s1,s2),
         conf.int = corr_shieh_lim)
  )

}

# Adding a default method in defining a function called shieh.d.default

tr.shieh.d.default <- function(x, y,conf.level=.95,alternative="two.sided", na.rm=TRUE){
  out <- tr.shieh.dEst(x, y,conf.level,alternative, na.rm)
  out$statistic <- out$statistic
  out$parameter <- out$parameter
  out$estimate <- out$estimate
  out$conf.int <- out$conf.int
  out$call <- match.call()

  class(out) <- "tr.shieh.d"
  out
}

print.tr.shieh.d <- function(x,...){
  cat("\n             transformed Shieh's d\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\n Transformed Shieh's d :\n")
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
