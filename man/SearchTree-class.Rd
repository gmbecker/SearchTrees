\name{SearchTree-class}
\Rdversion{1.1}
\docType{class}
\alias{SearchTree-class}

\title{Class \code{"SearchTree"}}
\description{
A virtual class representing a search tree for storing geometric points
in a manner designed for efficient lookup. 
}
\section{Objects from the Class}{
  This is a virtual class so objects of class SearchTree cannot be created directly.No methods defined with class "SearchTree" in the signature.
}
\section{Slots}{
  \describe{
    \item{\code{ref}:}{Object of class \code{"externalptr"} Pointer to
      the internal representation of the tree. }
    \item{\code{numNodes}:}{Object of class \code{"integer"} Number of
      nodes in the tree }
    \item{\code{dataNodes}:}{Object of class \code{"integer"} Number of
      nodes in the tree which are storing data.}
    \item{\code{maxDepth}:}{Object of class \code{"integer"} Maximum
      depth of the tree}
    \item{\code{maxBucket}:}{Object of class \code{"integer"} Maximum
      number of data points stored in a single node}
    \item{\code{totalData}:}{Object of class \code{"integer"} Number of
      data objects stored in the tree.}
    \item{\code{dataType}:}{Object of class \code{"character"} Indicates
      type of data stored in the tree.
    }
  }
}
\section{Methods}{
  \code{knnLookup}, \code{rectLookup}
}
\author{
  Gabriel Becker
}
\seealso{
  \code{\linkS4class{QuadTree}} \code{\link{createTree}}
}
\keyword{classes}
\keyword{lookup}