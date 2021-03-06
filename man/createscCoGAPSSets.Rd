% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createscCoGAPSSets.R
\name{createscCoGAPSSets}
\alias{createscCoGAPSSets}
\title{Create Gene Sets for scCoGAPS}
\usage{
createscCoGAPSSets(D, nSets, simulationName, samplingRatio = NULL,
  path = "", anotionObj = NULL)
}
\arguments{
\item{D}{data matrix}

\item{nSets}{number of sets to partition the data into}

\item{simulationName}{name used to identify files created by this simulation}

\item{samplingRatio}{vector of relative quantities to use for sampling celltypes}

\item{path}{character string indicating were to save resulting data objects. default is current working dir}

\item{anotionObj}{vector of same length as number of columns of D}
}
\value{
simulationName used to identify saved files
}
\description{
factors whole genome data into randomly generated sets for indexing
}
\examples{
data(SimpSim)
createscCoGAPSSets(SimpSim.D, nSets=2, simulationName="example")
}

