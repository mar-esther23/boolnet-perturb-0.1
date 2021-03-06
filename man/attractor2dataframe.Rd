% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoolNetPerturb.R
\name{attractor2dataframe}
\alias{attractor2dataframe}
\title{Convert attractor to data frame.}
\usage{
attractor2dataframe(attr, sep = "/")
}
\arguments{
\item{attr}{BoolNet attractor object}

\item{sep}{string to join elements with length > 1, default "/"}
}
\value{
Dataframe, each column corresponds to a property of the attractor
}
\description{
Convert an BoolNet attractor object to a data frame. Each property of attr$attractors corresponds to a dataframe column. If the propertie has elements with length > 1 it converts them to a string and joins them with sep.
}
\examples{
> attr <- getAttractors(cellcycle$genes)
> attractor2dataframe(attr)
              involvedStates basinSize
1                        162       512
2 25/785/849/449/389/141/157       512
}

