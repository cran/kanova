\name{kanova}
\alias{kanova}
\title{
    Pseudo analysis of variance of K-functions
}
\description{
    One-way and two-way pseudo analysis of variance for replicated
    point patterns, grouped by one or two classification factors.
    The analysis is based on the values of a summary function.
    This summary function may be specified by the user.  It is
    usually one of the four standard summary functions, most often
    the K-function, which is the default.
}
\usage{
    kanova(fmla,data,expo=1,rsteps=128,r=NULL,sumFnNm=NULL,
           warnSFN=TRUE,test=TRUE,bylevel=FALSE,
           permtype=c("stdres","data"),nperm=99,brief=TRUE,verb=TRUE,
           keepdata=FALSE,divByVar=TRUE)
}
\arguments{
  \item{fmla}{A formula specifying the test to be conducted.
  (See \bold{Details}.)  There can be at most two main effect
  predictors (possibly with interaction between them).  The left
  hand side of \code{fmla} may be omitted.  If so, it is taken to
  be the name of the first column of \code{data}.
}
  \item{data}{A hyperframe (see
  \code{\link[spatstat.geom]{hyperframe}()}) containing the data to
  be analysed.  If the left hand side of \code{fmla} is supplied,
  then \code{data} must have a column with a name which matches
  this left hand side. If no such column exists, an error is thrown.
  The response column may be a list of point patterns, or it may be
  a list of numeric vectors (all of which must have the same length).

  These numeric vectors are notionally values of a diagnostic
  or summary function, applied to point patterns, and then
  evaluated at the attribute \code{"r"} (see below) of \code{data}.
  However that the entries in the list need not actually \emph{be}
  summary functions.  They are simply numeric vectors and \emph{not}
  objects of class \code{"fv"} as returned by a summary function.
  In particular they may be scalars, which allows the \code{kanova()}
  function to be applied to analysis of variance of scalars.

  Argument \code{data} may also have one or two columns with names
  matching those of the main effect predictors as specified by
  \code{fmla}.  If such columns cannot be found, then objects
  with these names are sought in the \dQuote{parent frame}
  (see \code{\link{parent.frame}()}.)  If the appropriate objects
  still cannot be found, then an error is thrown.  These columns
  (or objects) are internally coerced to factors.

  If the response consists of numeric vectors, and if these vectors
  are \emph{not} scalars, then \code{data} must have an attribute
  \code{"r"}  which specifies the \dQuote{argument} of the pseudo
  summary functions in the list.  This attribute should be a
  numeric vector and must be of the same length as the numeric
  vectors in the list.  If the response functions are scalars,
  then no \code{"r"} attribute is needed, and, if present, is ignored.

  Finally, \code{data} may have a columns named \code{"wts"}
  (\dQuote{weights}) whose entries are positive scalars.  Notionally
  they may be thought of as a power of the numbers of
  points in the corresponding point patterns.  If there is no column
  named \code{"wts"} then the weights are all taken to equal 1.

  If the response column consists of a list of point patterns then
  any column named \code{"wts"} is ignored, and the weights are
  indeed a power of the numbers of points in the
  corresponding point patterns.  The impact of using these weights
  (when \code{expo}, see below, is strictly positive) is to diminish
  the influence of quantities corresponding to patterns having few
  points, and conversely to emphasise the influence of quantities
  corresponding to patterns having many points.
}
  \item{expo}{Non-negative numeric scalar.  Ignored unless the
  response is a list of point patterns.  Statistics in \code{kanova}
  are in general calculated using \emph{weighted} means and in this
  context the weights are the counts of points in
  the patterns, raised to the power \code{expo}.  To use unweighted
  means, set \code{expo=0}.  If \code{expo} is equal to \code{1}
  then the weights are simply the pattern counts.
}
  \item{rsteps}{Integer scalar.  Ignored if argument \code{r}
  (see below) is supplied (i.e. is not \code{NULL}), or if the
  response consists of a list of numeric vectors (rather than a
  list of point patterns). The argument \code{rsteps} specifies the
  number of (equal) steps between values of the vector \code{r} at
  which the summary function is evaluated.  The values of \code{r}
  are equispaced on the interval from 0 to \code{rtop}, the latter
  being calculated internally.  The value of \code{rtop} depends
  on the observation windows of the patterns in the response and
  on their intensities.  It also depends on the summary function
  being used.
}
  \item{r}{Numeric vector.  Ignored unless the response consists of
  a list of point patterns.  Note that if the response consists of
  a list of vectors, representing notional diagnostic functions,
  then argument \code{data} must have an attribute \code{r} to
  provide the relevant \code{r} object.  In this setting an error
  is thrown if the required attribute is not present (unless the
  response vectors are in fact scalars).

  The argument \code{r} specifies the values (distances) at
  which the summary/diagnostic function (see below) is evaluated.
  Errors may be thrown if the first entry of \code{r} is not 0,
  or if the values of \code{r} are not sufficiently finely spaced.
  Generally users should not specify \code{r}.  unless they have
  a sound understanding of what they are doing.
}
  \item{sumFnNm}{Character string naming the summary/diagnostic
  function to be used.  If this is not one of the \dQuote{standard
  four}, i.e. \code{"Kest"}, \code{"Fest"},  \code{"Gest"}, or
  \code{"Jest"},and if \code{warnSFN} (see below) is \code{TRUE},
  then a warning is issued.  In this case \bold{there may
  be problems}; the code is not robust in this respect.
  If \code{sumFnNm} is not specified (left \code{NULL}) then it
  defaults to \code{"Kest"}.  Users should apply summary functions
  other than \code{Kest()} only if they have a sound understanding
  of what they are doing.  See \bold{Notes}.
}
  \item{warnSFN}{Logical scalar.  Should a warning be issued if
  \code{sumFnNm} is not one of the \dQuote{standard four}?
}
  \item{test}{Logical scalar.  Should a Monte Carlo test of the
  null hypothesis be carried out?
} 
  \item{bylevel}{Logical scalar.  Should a test of the model
  \code{y ~ A + B} be carried out as \code{b} one-variable tests
  of the significance of A, within each level of B?  That is,
  should we perform \code{b} tests specified by \code{y ~ A,
  data=sdata[[i]]}, \code{i = 1, ..., b} where \code{b} is the
  number of levels of the factor B, and \code{sdata[[i]]} consists
  of those data corresponding to the i-th level of factor B?

  If \code{bylevel} is \code{TRUE}, then the object returned by
  \code{kanova()} is an object of class \code{"multi.kanova"}, which is
  a list of length \code{b}, each entry of which is an object of
  class \code{"kanova"}.
}
  \item{permtype}{Character string specifying what sort of
  permutations should be done to produce the Monte Carlo test
  statistics.  Ignored if \code{test} is \code{FALSE}.

  If \code{permtype} is \code{"stdres"} then the Monte Carlo data
  are formed by permuting the residuals, from the \emph{saturated}
  model, which have been standardised by dividing them by their
  (estimated) standard deviations.  After permutation, the residuals
  are \dQuote{unstandardised} by multiplying them by the appropriate
  standard deviation.  The results are then added back to the fitted
  values from the null model.

  If \code{permtype} is \code{"data"} then the Monte Carlo data
  are formed by permuting the original data sets. 

  In the two-way setting, when the test is for the main effect A,
  then the data (if \code{permtype} is \code{"data"}) are permuted
  \emph{within} the levels of B.

  If \code{fmla} specifies an  interaction between the main effect
  predictors, then \code{permtype} cannot be \code{"data"} and
  \emph{must} be \code{"stdres"}, otherwise an error is thrown.
}
  \item{nperm}{
  The number of permutations to be used to determine the Monte Carlo
  \eqn{p}-value.  Ignored if \code{test} is \code{FALSE}
}
  \item{brief}{Logical scalar.  Should the object returned by this
  function be \dQuote{brief}?  See \bold{Value}.
}
  \item{verb}{
  I.e. \dQuote{verbose}.  Logical scalar.  Should rudimentary
  \dQuote{progress reports} be printed out (in the course of
  conducting the permutation test for \dQuote{significance} of
  the test statistic)?  Such \dQuote{reports} consists simply of
  indications of how many permutations have been effected so far.
  Ignored if \code{test} is \code{FALSE}.
}
  \item{keepdata}{Logical scalar.  Should a copy of the data, to
  which the model has been fitted, be included as a component of
  of the object returned by \code{kanova()}?  See \bold{Value}
  for more detail.
}
  \item{divByVar}{Logical scalar.  Should the components of the
  test statistic be divided by the variances of the underlying
  un-squared values? (See \bold{Details}.)
}
}
\details{
  \describe{
  \item{formulae:}{
  \eqn{~}\cr
  The formulae used in the \code{fmla} argument may take the
  form \code{y ~ A}, \code{y ~ A + B}, or \code{y ~ A * B}.
  These \dQuote{look like} those used in \dQuote{ordinary} analysis
  of variance, but in the second instance the interpretation
  is different.  The formula \code{y ~ A + B} does not actually
  fit the additive \code{A + B} model.  It effects a test of the
  \dQuote{significance} of A, \dQuote{allowing for B} (see below).
  It is thereby obvious that \code{y ~ A + B} is not equivalent
  to \code{y ~ B + A} (whereas in \dQuote{ordinary} analysis of
  variance) they \emph{are} equivalent).

  The formula \code{y ~ A * B} tests the model with interaction
  against the additive model, as in ordinary analysis of variance.
  The additive model is not often meaningful in the current
  context, so such a test for interaction is probably not meaningful
  either.  The test is included in the code basically for the
  sake of completeness, and to allow for the possibility that a
  user may encounter a circumstance in which the additive model
  actually is meaningful.
  }
  \item{allowing for a second effect:}{
  \eqn{~}\cr
  This concept is pertinent only when the formula in question is of
  the form \code{y ~ A + B}.  In this setting, if \code{permtype}
  is \code{"data"}, the factor B is ``allowed for'' by permuting
  the data \emph{within} the levels of B.  If \code{permtype}
  is \code{"stdres"}) the standardised residuals are permuted.
  These are residuals from the saturated model, which includes a
  B effect which \dQuote{allows for} B.
  }
  \item{integration:}{
  \eqn{~}\cr
  The value of the test statistic is obtained as a sum of numerical
  integrals of certain sums of squares.  The integrals are computed
  via a rough trapezoid rule.
  The integration is carried out over the value of \code{r}, the
  argument of the summary functions that are being analysed.  If the
  response consists of numeric vectors of length 1, i.e. of scalars,
  then no integration is in fact performed, and the corresponding
  (downweighted) sum of squares is returned.  You may, if you like,
  think of this as integrating with respect to a measure which has
  a point mass of 1 at a conceptual single value of \code{r}.
  }
  \item{the \code{divByVar} argument:}{
  \eqn{~}\cr
  The \code{divByVar} argument exists essentially to allow the
  package developers to conduct certain simulation experiments.
  In normal circumstances the user would not set the value of
  this argument (i.e. would leave the value of this argument as
  \code{TRUE}).

  If the argument \code{divByVar} is \code{TRUE}, then the sums
  of squares, from which the test statistic is formed, are are
  divided by the estimated variance of the quantity being squared.
  This procedure is analogous to the studentisation procedure
  used by Hahn, 2012.  If \code{divByVar} is \code{FALSE} then
  no such division takes place and the sums of squares involved
  are \dQuote{raw}.

  The quantities that are involved in the sums of squares are
  formed from certain \dQuote{fitted values} which are weighted
  means of the observations (observed values of summary functions).
  The variance referred to is formed as a weighted mean of squares
  of the residuals, which are of course equal to the observations
  minus the fitted values.

  }
  \item{more detail:}{
  \eqn{~}\cr
  More detail about the test statistic, the fitted values, the
  residuals and the estimated variance, can be found in the vignette
  \code{"testStat"}.
  }
}}

\section{Warning}{
  When \code{keepdata} is \code{TRUE}, the way that the data
  is \dQuote{kept} depends on the value of \code{bylevel}.
  See \bold{Value}.
}
\section{Notes}{
  \itemize{
  \item Simulation experiments have given some evidence that
  \code{Fest()} and \code{Gest()} and \code{Jest()} lead to tests
  that have lower power than that obtained than tests obtained by
  using \code{Kest()}.  The power obtained seems to be substantially
  lower in the case of \code{Fest()} and \code{Gest()}, somewhat
  lower in the case of \code{Jest()}.

  Consequently, users are advised to eschew the use of \code{Fest()},
  \code{Gest()} and \code{Jest()} (despite their ready availability)
  and to stick with the default summary function \code{Kest()}.
  Users should ignore this advice \emph{only} if they have a sound
  reason for doing so and a sound understanding of the consequences.

  \item Only one-way and two-way (pseudo) analyses of variance
  are accommodated.  If you feel inclined to ask why there is no
  provision for higher order analysis of variance, just look at
  the code and the answer should be obvious.  It \emph{might} be
  possible to implement higher order pseudo analysis of variance of
  summary functions, but this is unlikely to have any practical use.
  Writing the code would, for me, be a nightmare!
  }
}
\value{If \code{bylevel} is \code{TRUE} then the object returned is
   of class \code{"multi.kanova"} which is a list of length equal
   to then number of levels of the second predictor in the model.
   Each entry in this list is an object of class \code{"kanova"}.
   If \code{keepdata} is \code{TRUE} then the object in question
   (of class \code{"multi.kanova"}) has an
   attribute \code{"data"}.  This is equal to the \code{data}
   argument, possibly having been augmented by the value of the second
   predictor in the model if this predictor was not found in the
   original \code{data} and was located in the parent frame.

   Note that the \code{"kanova"} components of a \code{"multi.kanova"}
   object \emph{never} themselves have a list entry named
   \code{"data"}, irrespective of the value of \code{keepdata}.
   If \code{keepdata} is \code{TRUE}, then the appropriate
   \code{data} object is returned as an \emph{attribute} of the
   \code{"multi.kanova"} object.

   If \code{test} is \code{TRUE} then the \code{"multi.kanova}"}
   object returned has an attribute named \code{"oapv"}
   (\dQuote{overall \eqn{p}-value}).  This is calculated as \code{1 -
   (1 - pvmin)^b} where \code{pvmin} is the minimum of the \code{b}
   \eqn{p}-values obtained when each of the levels of the second
   predictor is tested individually for an A effect.  Note that
   \deqn{F(y) = 1 - (1-y)^b} is the pdf of the minimum of \eqn{b}
   independent observations that are uniformly distributed on [0,1].
   Alternatively one may think of the expression for \code{oapv}
   as arising from the Sidak adjustment to the minimum \eqn{p}-value
   to allow for multiple comparisons.

   If \code{bylevel} is \code{FALSE} then the object returned is
   of class \code{"kanova"}, and is described as follows:

   If \code{brief} is \code{TRUE}, then the object in question is
   a list with components:

  \item{Effectname}{Character string naming the effect being tested for.}

  \item{stat}{Numerical scalar equal to the value of the test statistic
  calculated from the original data.}

  \item{pvalue}{The Monte Carlo \eqn{p}-value of the test calculated
  by comparing \code{stat} with test statistics formed from simulated
  data, generated by permutation, which satisfy the null hypothesis.}

  If \code{brief} is \code{FALSE} then the \code{"kanova"} object in
  question has additional components:

  \item{nperm}{The \code{nperm} argument.}

  \item{permtype}{The \code{permtype} argument.}

  \item{Tstar}{The vector of \code{nperm} values of the test statistic
  calculated from the simulated data sets.}

  \item{fmla}{The \code{fmla} argument.}

  \item{sumFnNm}{The \code{sumFnNm} argument.}

  \item{data}{The \code{data} argument.  This may possibly have
  been augmented by any predictor values which were not found in
  the original \code{data} and were located in the parent frame.}

  Components \code{pvalue}, \code{nperm} and \code{Tstar} are
  present in an object of class \code{"kanova"} only if \code{test}
  is \code{TRUE}.  Component \code{data} is present only if
  \code{keepdata} is \code{TRUE}.
}
\references{
Diggle, Peter J., Mateu, Jorge and Clough, Helen E. (2000) A
comparison between parametric and non-parametric approaches to the
analysis of replicated spatial point patterns, \emph{Advances in
Applied Probability} \bold{32}, pp. 331 -- 343.

Diggle, P. J., Lange, N. and Benes, F. (1991) Analysis of variance
for replicated spatial point patterns in clinical neuroanatomy,
\emph{Journal of the American Statistical Association}, \bold{86},
pp. 618 -- 625.

Hahn, Ute (2012) A studentized permutation test for the
comparison of spatial point patterns, \emph{Journal of the
American Statistical Association}, \bold{107}, pp. 754 -- 764,
DOI: 10.1080/01621459.2012.688463.
}
\author{Rolf Turner
  \email{rolfturner@posteo.net}
}

\seealso{
\code{\link[spatstat.explore]{studpermu.test}()}
}

\examples{
# The following is inappropriate since there is a second
# classification factor.
set.seed(104)
s1 <- kanova(patterns ~ Pos,data=stomata,permtype="d",nperm=9)

# Here we are testing for a Layer effect allowing for a Pos effect.
set.seed(7)
s2a <- kanova(patterns ~ Layer + Pos, data=stomata,permtype="d",nperm=9)
s2b <- kanova(patterns ~ Layer + Pos, data=stomata,permtype="s",nperm=9)

# Here we are testing for a Pos effect allowing for a Layer effect.
set.seed(78)
s3a <- kanova(patterns ~ Pos + Layer, data=stomata,nperm=9)
# permtype defaults to "stdres".
\dontrun{ # Takes too long.
    set.seed(24)
    s3b <- kanova(patterns ~ Pos + Layer, data=stomata,nperm=999)
    # Get a p-value of 0.001
}

# Here we are testing for a Pos effect by testing for such an
# effect within each level of Layer.
set.seed(770)
s3c <- kanova(patterns ~ Pos + Layer, bylevel=TRUE,data=stomata,nperm=9)
# attr(s3c,"oapv") is 0.2172 --- not significant.

# Here, we are testing for a Layer by Pos interaction.  Unlikely to
# be meaningful.
set.seed(2)
s4 <- kanova(patterns ~ Layer * Pos, data=stomata,nperm=9) # permtype must be "s"

# Artificial data.
\dontrun{ # Takes too long.
if(requireNamespace("spatstat.geom")) {
    set.seed(3)
    r    <- seq(0,25,length=129)
    rsp  <- lapply(1:144,function(k){pi*r^2 + runif(129,-0.1,0.1)})
    rsp  <- lapply(rsp,function(x){pmax(0,x)})
    fctr1 <- factor(rep(1:4,12,each=3))
    fctr2 <- factor(rep(1:3,48))
    wts  <- sample(50:100,144,replace=TRUE)
    X    <- spatstat.geom::hyperframe(rsp=rsp,fctr1=fctr1,fctr2=fctr2,wts=wts)
    attr(X,"r") <- r
    set.seed(118)
# Testing for a fctr1 effect, allowing for a fctr2 effect.
    s5a  <- kanova(rsp ~ fctr1 + fctr2, data=X,brief=FALSE,nperm=9)
# Testing for an interaction; meaningful in this artificial data
# context.
    s5b  <- kanova(rsp ~ fctr1*fctr2,data=X,brief=FALSE,nperm=9)
}
}

# Scalar data.
\dontrun{ # Takes too long.
if(requireNamespace("Devore7")) {
    X   <- spatstat.geom::as.hyperframe(Devore7::xmp11.10)
    s6a <- kanova(Tempr ~ Period*Strain,data=X,nperm=999)
    s6b <- kanova(Tempr ~ Period+Strain,data=X,nperm=999)
    s6c <- kanova(Tempr ~ Strain+Period,data=X,nperm=999)
    chk <- lm(Tempr ~ Period*Strain,data=X)
# Executing anova(chk) reveals p-values that are
# at least roughly similar to those in s6a, s6b, and s6c.
}
}
}

\keyword{ htest }
