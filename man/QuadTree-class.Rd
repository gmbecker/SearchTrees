\name{QuadTree-class}
\Rdversion{1.1}
\docType{class}
\alias{QuadTree-class}
\title{Class \code{"QuadTree"}}
\description{
A class representing a Quad Tree object for storing 2 dimensional points
for efficient rectangular range and knn lookup.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("QuadTree", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{ref}:}{Object of class \code{"externalptr"} Pointer to
      the internal representation of the tree }
    \item{\code{numNodes}:}{Object of class \code{"integer"} Number of
      nodes in the tree }
    \item{\code{dataNodes}:}{Object of class \code{"integer"} Number of
      nodes in the tree which are storing data}
    \item{\code{maxDepth}:}{Object of class \code{"integer"} Maximum
      depth of the tree.}
    \item{\code{maxBucket}:}{Object of class \code{"integer"} Maximum
      number of data points which are stored at a single node}
    \item{\code{totalData}:}{Object of class \code{"integer"} Number of
      objects stored in the tree }
    \item{\code{dataType}:}{Object of class \code{"character"} Indicates
      type of data stored in the tree.
    }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SearchTree}"}, directly.
}
\section{Methods}{
  \describe{
    \item{knnLookup}{\code{signature(tree = "QuadTree")}: ... }
    \item{rectLookup}{\code{signature(tree = "QuadTree")}: ... }
  }
}
\author{
Gabriel Becker
}
\note{
When using createIndex to create a quadTree, only two columns of the
matrix/data.frame passed to the function will be used to create the
tree. See the columns argument in \code{\link{createTree}}
}
\seealso{
  \code{\link{createTree}}
}
\examples{
showClass("QuadTree")
}
\keyword{classes}
