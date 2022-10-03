\name{createTree}
\alias{createTree}
\title{
Create a Search Tree Index
}
\description{
Create a search tree from the supplied data for use
  in during future lookups.
}
\usage{
createTree(data, treeType = "quad", dataType = "point",
           columns = if (dataType=="point") 1:2 else 1:4, ...)
}
\arguments{
  \item{data}{
    data.frame or matrix. Data to be indexed.
}
  \item{treeType}{
    Character. Indicates type of index tree to be created. Currently
    only "quad" (quad trees) is supported.
}
  \item{dataType}{
    Character. Indicates type of data being indexed. Currently "point",
    and "rect" are supported corresponding to points and rectangles,
    respectively.  Defaults to "point".
}
  \item{columns}{
    Numeric. Indicates columns in \code{data} the information to be
    indexed can be found. Length depends on value of \code{dataType} (2
    for "point" and 4 for "rect"). Defaults to columns 1 and 2. See Details.
}
  \item{\dots}{
    Any additional/type specific parameters to be passed to the tree
    creation function. These include:
    \describe{
      \item{maxDepth:}{
	Numeric. Maximum depth of the tree. Defaults to 7.
      }
      \item{minNodeArea:}{
	Numeric. Minimum (rectangular) area to be represented by a
	single node. When set, this overrides \code{maxDepth}
      }
   }
    
  }
}
\details{
  For a point based tree, the two columns specified in \code{columns}
  represent the x and y values of the points.

  For a rectangle based tree, four columns must be specified. These
  columns represent the x and y coordinates of point 1 and the x and y
  coordinates of point 2, in that order (where point 1 and point 2
  specify the rectangle to be stored).
}
\value{
  The class of the returned object depends on the tree type created,
  though all will inherit from the \code{SearchTree} S4 class and have the
  following slots:
  \item{ref}{An external pointer to the C level data structure.}
  \item{numNodes}{Total number of nodes comprising the tree.}
  \item{dataNodes}{Number of nodes which store at least
    one data point.}
  \item{maxDepth}{Maximum depth of the tree.}
  \item{maxBucket}{Maximum number of data points stored in a single
     node.}
   \item{totalData}{Number of items indexed in the tree.}
   \item{dataType}{Type of objects stored in the tree.}
 }
 \references{
   Finkel, R. A. and Bentley, J. L. "Quad Trees, a Data Structure for Retrieval on Composite Keys." Acta Informatica 4, 1-9, 1974.}
 \author{
   Gabriel Becker
}
\seealso{
  \code{\linkS4class{SearchTree}} \code{linkS4Class{QuadTree}}
}
\examples{
x = rnorm(100)
y = rnorm(100)
dat = cbind(x,y)
tree = createTree(dat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ indexing }
\keyword{ quadtree }% __ONLY ONE__ keyword per line
