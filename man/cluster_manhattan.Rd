% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{cluster_manhattan}
\alias{cluster_manhattan}
\title{Title Clustering Using Manhattan Distances}
\usage{
cluster_manhattan(x, k.gap, plot.clustering)
}
\arguments{
\item{x}{data matrix}

\item{k.gap}{positive integer specifying the number of clusters, less than the number of observation. Default value is 10.}

\item{plot.clustering}{if TRUE (default) a 2-dimensional “clusplot” plot will be printed}
}
\value{
vector with the new assigned clusters
}
\description{
Title Clustering Using Manhattan Distances
}
\examples{
# cluster_manhattan(subx, k.gap=4, plot.clustering=TRUE)

}
