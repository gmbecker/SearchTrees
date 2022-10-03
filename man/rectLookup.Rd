\name{rectLookup}
\alias{rectLookup}
\alias{rectLookup, QuadTree-method}
\title{
  Perform Rectangular Lookup in 2d Space
}
\description{
  Determine which objects, stored in a SearchTrees indexing object, fall
  within a given rectangle in two-dimensional space.
}
\usage{
rectLookup(tree, ptOne, ptTwo, xlims, ylims)
}
\arguments{
  \item{tree}{
    SearchTree. A SearchTree object to perform the lookup on.
  }
  \item{ptOne}{
    Numeric. A numeric of length two indicating x and y values for one
    corner of the rectangle.
  }
  \item{ptTwo}{
    Numeric. A numeric of length two indicating x and y values for the
    corner of the rectangle opposite to \code{ptOne}
  }
  \item{xlims}{
    Numeric. A numeric vector indicating the minimum and maximum x value
    for the rectangle. Overrides \code{ptOne} and \code{ptTwo}
  }
  \item{ylims}{
    Numeric. A numeric vector indicating the minimum and maximum y value for
    the rectangle. Overrides \code{ptOne} and \code{ptTwo}
  }
}
\details{
  In the case of lookup for rectangular objects, any rectangle which
  overlaps the query rectangle will be returned. 
}
\value{
  A numeric vector indicating the indicies of the object (in the order
  they were in when the SearchTree object was created) which fall (at
  least partially) within the rectangular query.
}
\author{
  Gabriel Becker
}
\seealso{
  \code{\linkS4class{QuadTree}} \code{\link{knnLookup}}
}
\examples{
x = rnorm(100)
y = rnorm(100)
x2 = x + runif(100, .5, 2)
y2 = y + runif(100, .5, 2)
dat2 = cbind(x, y, x2, y2)
tree2 = createTree(dat2, dataType="rect", columns= 1:4)
inrect = rectLookup(tree2, xlim = c(0,1), ylim=c(0, 1))
col = rgb(0, 1, 0, alpha=.5)
plot(x, y2, col="white")
rect(x[inrect], y[inrect], x2[inrect], y2[inrect], col=col)
rect(0, 0, 1, 1, col="blue", lwd=3)
}
\keyword{ rectangular }
\keyword{ lookup }
\keyword{ query }
