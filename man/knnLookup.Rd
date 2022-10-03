\name{knnLookup}
\alias{knnLookup}
\alias{knnLookup, QuadTree-method}
\title{
Perform k-Nearest Neighbors Lookup Using a Search Tree
}
\description{
  This function performs fast k-Nearest Neighbors lookup on a
  SearchTree object
}
\usage{
knnLookup(tree, newx, newy, newdat, columns = 1:2, k = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
    An object which inherits from the \code{SearchTree} S4 class.
  }
  \item{newx}{
    Numeric. Vector of x values for the points to look up neighbors for.
  }
  \item{newy}{
    Numeric. Vector of x values for the points to look up neighbors for.
  }
  \item{newdat}{
    Matrix or data.frame. Data containing x and y values of the points
    to look up neighbors for. Ignored if \code{x} and \code{y} are specified.
  }
  \item{columns}{
    Numeric. Columns x and y values can be found in within \code{newdat}.
  }
  \item{k}{
    Numeric. Number of neighbors to find for each point.
  }
}
\value{
The return value is an integer matrix indicating the indices in
  the original data used to create \code{treE} where the nearest neighbors were found. Row indicates
  the indice of the new point, while column indicates the order of the k neighbors.
}
\author{
  Gabriel Becker
}
\note{
  No defined order is specified for exact ties in distance.
}
\seealso{
  \code{\link{createTree}} \code{\link{rectLookup}}
}
\examples{
x = rnorm(100)
y = rnorm(100)
tree = createTree(cbind(x,y))
newx = c(0, .5)
newy = c(.5, 0)
inds = knnLookup(tree, newx, newy, k=7)

ch = rep(1, times=100)
ch[inds[1:7]] = 3
ch[inds[8:14]] = 5
cls = rep("black", times=100)
cls[inds[1:7]] = "red"
cls[inds[8:14]] ="blue"

plot(x,y, pch=ch, col = cls)
abline(v=newx[1], h = newy[1] , col="red")
abline(v=newx[2], h = newy[2], col = "blue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ knn }
\keyword{ lookup }
\keyword{ neighbors }
