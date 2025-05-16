\name{kanova-internal}
\alias{buildM1}
\alias{buildV}
\alias{datagenkv}
\alias{datagenpp}
\alias{datagensto}
\alias{estSigsq}
\alias{getWts}
\alias{iEngine}
\alias{initPrep}
\alias{ordinal}
\alias{oEngine}
\alias{ordinalsuffix}
\alias{permSumFns}
\alias{permWithin}
\alias{reenlist}
\alias{resAndFit}
\alias{resAndFitCmpnts}
\alias{stilsqFn}
\alias{testStat}
\alias{tEngine}
\alias{trapint}
\alias{wtdMean}
\alias{wtdSS}
\alias{wtdVar}
\title{Internal kanova functions.}
\description{
  Internal kanova functions.
}
\usage{
buildM1(Khat,Khati,Khatj,Khatij)
buildV(s2ij,wtm)
datagenkv(av,bv=NULL,interac=NULL,nrep=5,rlen=129,sigma=1,
          seed=NULL,pseudonoise=NULL)
datagenpp(kapmin=30,kapmax=45,nkap=6,scale=0.1,lambda=300,
          nrep=10,mdlEff=TRUE,seed=NULL)
datagensto(nrep=10,sigma=1,interac=0,seed=NULL,
           meansonly=FALSE,perturbLayer=0)
estSigsq(sumFns,satMod=FALSE)
getWts(x)
iEngine(sumFns,divByVar)
initPrep(data,rspNm,Anm,Bnm=NULL,sumFnNm,type,expo,rsteps,r)
oEngine(sumFns,divByVar)
ordinal(k)
ordinalsuffix(k)
permSumFns(sumFns,rAndF,permtype)
permWithin(G)
reenlist(x,f)
resAndFit(sumFns)
resAndFitCmpnts(sumFns)
stilsqFn(s2j,wtm)
testStat(sumFns,divByVar)
tEngine(sumFns,divByVar)
trapint(y,r)
wtdMean(x)
wtdSS(x)
wtdVar(x)
}
\value{
    \describe{
        \item{buildM1}{An array.}
        \item{buildV}{An array.}
        \item{datagenkv}{A hyperframe.}
        \item{datagenpp}{A hyperframe.}
        \item{datagensto}{A hyperframe.}
        \item{estSigsq}{A numeric scalar or a list.}
        \item{getWts}{A numeric vector.}
        \item{iEngine}{A numeric scalar.}
        \item{initPrep}{A list.}
        \item{ordinal}{A character vector.}
        \item{oEngine}{A numeric scalar.}
        \item{ordinalsuffix}{A character vector.}
        \item{permSumFns}{A list.}
        \item{permWithin}{A list.}
        \item{reenlist}{A list.}
        \item{resAndFit}{A list.}
        \item{resAndFitCmpnts}{A list.}
        \item{stilsqFn}{A numeric matrix.}
        \item{testStat}{A numeric scalar.}
        \item{tEngine}{A numeric scalar.}
        \item{trapint}{A numeric scalar.}
        \item{wtdMean}{A numeric scalar.}
        \item{wtdSS}{A numeric scalar.}
    }
}
\details{
   These functions are auxiliary and are not intended to be called
   by the user.
}
\keyword{internal}
