#include <stdio.h>
#include <R.h>
#include <Rdefines.h>


typedef struct qtree2 {
  double upper;
  double lower;
  double left;
  double right;
  int numdata;
  char pos; //to indicate which sector in the parent node this node is. 1 = uleft 2 = uright 3=lright 4 = lleft;
  unsigned char depth;
  int **data;
  char datatype;//indicates the type of objects stored in the tree. 1=points 2=rectangle
  struct qtree2 *parent;
  struct qtree2 *uleft;
  struct qtree2 *uright;
  struct qtree2 *lleft;
  struct qtree2 *lright;
} qtree2_t;


typedef struct point {
  double x;
  double y;
  int index;
} point_t;

typedef struct rect {
  //parallel to the axes only, other rectangles will need to wait until we do the general polygon case;
  double left; //lower x bound
  double right; //upper x bound
  double low; //lower y bound
  double high; //lower x bound
  struct rect *parent;
  int index;
} rect_t;


typedef struct kdtree {
  int depth;
  double split;
  int data;
  struct kdtree *parent;
  struct kdtree *left;
  struct kdtree *right;
} kdtree_t;

SEXP R_Build_Quadtree_Pt(SEXP Rx, SEXP Ry, SEXP RxMax, SEXP RxMin, SEXP RyMax, SEXP RyMin, SEXP RmaxDepth);

int Add_To_Bucket(qtree2_t *tree, void *obj, int type, int numdata, unsigned char maxDepth);

int Add_Pts_To_Bucket(qtree2_t *node, point_t **pt, int numdata, unsigned char maxDepth);

void Increase_Depth( qtree2_t *curnode);

qtree2_t *Create_Branch(double up, double low, double left, double right, qtree2_t *parent, char pos, unsigned char depth, char type);

SEXP R_Find_Neighbors_Pts(SEXP Rtree, SEXP Rnewx, SEXP Rnewy, SEXP Rk);

qtree2_t *Descend_To_Bucket(qtree2_t *node, void *obj, int type);

qtree2_t *Descend_To_Bucket_Pts(qtree2_t *node, point_t *pt);

void insert_dist(double *dists, double newdist, point_t **pts, point_t *newpt, int k, int start);
void Harvest_KNN_Pts(qtree2_t *node, int excludepos, double leftbound, double rightbound, double lowbound, double highbound, point_t **pts, double *dists, point_t *newpt, int k, int start);

int check_bounds(qtree2_t *node, double left, double right, double down, double up);

void R_free_quad_tree( SEXP treeptr);

void Free_QuadTree( qtree2_t *tree);

int Find_MaxDepth(qtree2_t *tree, unsigned char curdepth);

SEXP R_Find_MaxDepth(SEXP Rtree);

//point_t *New_Point_Struct(double x, double y, int ind);
double eucl_dist_pts(point_t *pt1, point_t *pt2);

void Grow_ReturnArray(int **found, int *cursize, char type);
void Rectangle_Pt_Lookup(qtree2_t *tree, double left, double right, double down, double up, int **found, int *pos, int *cursize);

int CheckBounds(qtree2_t *node, double left, double right, double down, double up);
int overlap(double left, double right, double down, double up, void *obj, char type);
void Rectangle_Lookup(qtree2_t *tree, double left, double right, double down, double up, int **found, int *pos, int *cursize, char type);
int Add_Rects_To_Bucket(qtree2_t *node, rect_t **rect, int numdata, unsigned char maxDepth);
void get_tree_attributes(qtree2_t *tree, int *curattr);
double get_max_dist(double *dists, int start, int k);
void allocInternalMem(int **arr, int start, int num, char type);
int get_struct_size(char type);
