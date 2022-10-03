#include "Rsearchtrees2.h"

point_t *New_Point_Struct(double x, double y, int ind)
{
  point_t *toret = malloc(sizeof(point_t));
  toret -> x = x;
  toret -> y = y;
  toret -> index = ind;
  return toret;
}

rect_t *New_Rect_Struct(double x1, double x2, double y1, double y2, int ind)
{
  rect_t *toret = malloc(sizeof(rect_t));
  
  toret->left = (x1 < x2) ? x1 : x2;
  toret->right = (x1 < x2) ? x2 : x1;
  toret->low = (y1 < y2) ? y1 : y2;
  toret->high = (y1 < y2) ? y2 : y1;
  toret->index = ind;
  
  return toret;
}
  
int points_equal(point_t *pt1, point_t *pt2)
{
  int ret = 1;
  if(!pt1 || !pt2)
    return 0;
  if(pt1 == pt2)
    return 1;
  if(pt1 -> x != pt2 -> x)
    ret = 0;
  else if (pt1 -> y != pt2 -> y)
    ret = 0;
  else if (pt1 -> index != pt2 -> index)
    ret = 0;
  
  return ret;
}
SEXP
R_Build_Quadtree_Rect(SEXP Rx1, SEXP Ry1, SEXP Rx2, SEXP Ry2, SEXP RxMax, SEXP RxMin, SEXP RyMax, SEXP RyMin, SEXP RmaxDepth)
{
  double *x1 = REAL(Rx1);
  double *x2 = REAL(Rx2);
  double *y1 = REAL(Ry1);
  double *y2 = REAL(Ry2);

  int len = LENGTH(Rx1);

  qtree2_t *tree;
  
  double upper = REAL(RyMax)[0];
  double lower = REAL(RyMin)[0];
  double left = REAL(RxMin)[0];
  double right = REAL(RxMax)[0];
  int maxDepth =  INTEGER(RmaxDepth)[0];

  tree = Create_Branch(upper, lower, left, right, NULL, 0, 0, 2); 
  rect_t **allrects = calloc(len, sizeof(rect_t *));
  for(int i = 0; i < len; i++)
    allrects[i] = New_Rect_Struct(x1[i], x2[i], y1[i], y2[i], i);

  Add_To_Bucket(tree, allrects, 2, len, maxDepth);

  SEXP klass, ans, ptr, ptr2;
  int *attr = calloc(4, sizeof(int));
  get_tree_attributes(tree, attr);
  PROTECT( klass = MAKE_CLASS( "QuadTree" ) );
  PROTECT( ans = NEW( klass ) );
  PROTECT( ptr = R_MakeExternalPtr( tree ,
				    Rf_install( "QuadTree" ),
				    R_NilValue ) );
  R_RegisterCFinalizerEx(ptr, &R_free_quad_tree, 1);
  SET_SLOT( ans, Rf_install( "ref" ), ptr );
  SET_SLOT( ans, Rf_install( "points" ), PROTECT(ScalarInteger( len ) ));
  
  SET_SLOT( ans, Rf_install( "numNodes" ), PROTECT(ScalarInteger(attr[0] )));
  SET_SLOT( ans, Rf_install( "dataNodes" ), PROTECT(ScalarInteger(attr[1] ) ));
  SET_SLOT( ans, Rf_install( "maxDepth" ), PROTECT(ScalarInteger(attr[2] ) ));

  SET_SLOT( ans, Rf_install( "maxBucket" ), PROTECT(ScalarInteger( attr[3] ) ));
  
  UNPROTECT(8);
  return ans;  

}

SEXP
R_Build_Quadtree_Pt(SEXP Rx, SEXP Ry, SEXP RxMax, SEXP RxMin, SEXP RyMax, SEXP RyMin, SEXP RmaxDepth)
{
  double *x = REAL(Rx);
  double *y = REAL(Ry);
  int len = LENGTH(Rx);
  unsigned char maxDepth = (unsigned char) INTEGER( RmaxDepth)[ 0 ];
  
  //qtree2_t tree;
  qtree2_t *tree;
  double upper = REAL(RyMax)[0];
  double lower = REAL(RyMin)[0];
  double left = REAL(RxMin)[0];
  double right = REAL(RxMax)[0];

  tree = Create_Branch(upper, lower, left, right, NULL, 0, 0, 1); //1 for points;
  int res;
  point_t **allpoints = calloc(len, sizeof(point_t*));

  for (int i=0; i < len; i++)
    {
      allpoints[i] = New_Point_Struct(x[i], y[i], i);
      //int Add_To_Bucket(qtree2_t *tree, void *obj, int type, int numdata, unsigned char maxDepth)

    }
  Add_To_Bucket(tree, allpoints, 1, len, maxDepth); //1 for point;
  //get info about the tree to return
  
  //first int is number of nodes, second int is number of nodes with data, third int is maxDepth, fourth int is max bucketsize

  int *attr = calloc(4, sizeof(int));
  get_tree_attributes(tree, attr);
  SEXP klass, ans, ptr, ptr2;
  PROTECT( klass = MAKE_CLASS( "QuadTree" ) );
  PROTECT( ans = NEW( klass ) );
  PROTECT( ptr = R_MakeExternalPtr( tree ,
				    Rf_install( "QuadTree" ),
				    R_NilValue ) );
  PROTECT( ptr2 = R_MakeExternalPtr( allpoints ,
				    Rf_install( "Data" ),
				    R_NilValue ) );
  
  R_RegisterCFinalizerEx(ptr, &R_free_quad_tree, 1);
  SET_SLOT( ans, Rf_install( "ref" ), ptr );
  SET_SLOT( ans, Rf_install( "data" ), ptr2 );
  SET_SLOT( ans, Rf_install( "points" ), PROTECT(ScalarInteger( len ) ));
  SET_SLOT( ans, Rf_install( "numNodes" ), PROTECT(ScalarInteger(attr[0] )));
  SET_SLOT( ans, Rf_install( "dataNodes" ), PROTECT(ScalarInteger(attr[1] ) ));
  SET_SLOT( ans, Rf_install( "maxDepth" ), PROTECT(ScalarInteger(attr[2] ) ));

  SET_SLOT( ans, Rf_install( "maxBucket" ), PROTECT(ScalarInteger( attr[3] ) ));
  UNPROTECT(9);
  free(attr);
  return ans;  
}

int Add_To_Bucket(qtree2_t *tree, void *obj, int type, int numdata, unsigned char maxDepth)
//single point of contact for all types
{
  int ret = 0;
  switch(type) 
    {
    case 1:
      ret = Add_Pts_To_Bucket(tree, (point_t **) obj, numdata, maxDepth);
      break;
    case 2:
      ret = Add_Rects_To_Bucket(tree, (rect_t **) obj, numdata, maxDepth); 
      break;
    }

  return ret;
}

int contains(void *obj, double left, double right, double lower, double upper, char type)
{
  int toret = 0;
  switch(type)
    {
    case 1:
      break;
    case 2:
      {
	rect_t *rec = (rect_t *) obj;
	if(rec->left >= left && rec->right <= right && rec->low >=lower && rec->high <= upper)
	  toret = 1;
	break;
      }
    }

  return toret;
}

int Can_Split(qtree2_t *node, void *rec, int maxDepth)
{
  int toret = 0;

  if(node->depth < maxDepth)
    {
      double left = node->left;
      double right = node->right;
      double down = node -> lower;
      double up = node->upper;
      double midx = (node->left + node->right)/2.0;
      double midy = (node->lower + node->upper)/2.0;
      if(contains(rec, left, midx, down, midy, 2) )
	toret = 1;
      else if ( contains(rec, left, midx, midy, up, 2))
	toret = 1;
      else if ( contains(rec, midx, right, down, midy, 2))
	toret = 1;
      else if ( contains(rec, midx, right, midy, up, 2))
	toret = 1;
    }
  return toret;

}

int Add_Rects_To_Bucket(qtree2_t *node, rect_t **rect, int numdata, unsigned char maxDepth)
{
  int ret = 0;
  qtree2_t *curnode = node;
  rect_t **olddat;
  int olddatsize;
  int found =0;
  
  for( int i =0; i<numdata; i++)
    {
      found=0;
      curnode = Descend_To_Bucket(node, rect[i], 2 ); //2 for rects;
      while(!found)
	{
	  if(contains(rect[i], curnode->left, curnode -> right, curnode->lower, curnode->upper, 2) || curnode->parent == NULL) //2 for rectangles
	    {
	      while(Can_Split(curnode, (void *) rect[i], maxDepth))
		{
		  Increase_Depth(curnode);
		  curnode = Descend_To_Bucket(curnode, rect[i], 2);
		}
	      
	      found = 1;
	      //if there is no data already there no memory has been allocated yet.
	      if (curnode -> numdata == 0)
		{
		  curnode -> data  = calloc(1, sizeof(rect_t));
		  curnode -> data[0] = (int *) rect[i];
		  curnode -> numdata = 1;
		} else  { 
		//we are at max depth so we must add the point to this bucket.
	  //realloc preserves values of previous elements
		curnode -> data = realloc(curnode -> data, sizeof(rect_t) * (curnode -> numdata + 1));
		curnode -> data [ curnode -> numdata ] = (int *) rect[i];
		curnode -> numdata ++;
		
		ret = 1 ;
		
	      } 
	    }else {
	    curnode = curnode->parent;
	  }
	}
    
    }

  return ret;

}
int Add_Pts_To_Bucket(qtree2_t *node, point_t **pt, int numdata, unsigned char maxDepth)
{
  int ret = 0;
  qtree2_t *curnode = node;
  point_t **olddat;
  int olddatsize;
  
  for( int i =0; i<numdata; i++)
    {
      curnode = Descend_To_Bucket(node, pt[i], 1 ); //1 for points;

      //if there is no data already there no memory has been allocated yet.
      if (curnode -> numdata == 0)
	{
	  curnode -> data  = calloc(1, sizeof(point_t));
	  curnode -> data[0] = (int *) pt[i];
	  curnode -> numdata = 1;
	} else if(curnode -> depth == maxDepth) { 
	  //we are at max depth so we must add the point to this bucket.
	  //realloc preserves values of previous elements
	  curnode -> data = realloc(curnode -> data, sizeof(point_t) * (curnode -> numdata + 1));
	  curnode -> data [ curnode -> numdata ] = (int *) pt[i];
	  curnode -> numdata ++;
	  
	  ret = 1 ;
      } else {
	//we are not yet at max depth, so we add depth and try again.
	olddat = (point_t **) curnode -> data;
	olddatsize = curnode -> numdata;
	Increase_Depth(curnode);
	curnode -> numdata = 0;
	curnode -> data = NULL;
	Add_Pts_To_Bucket(curnode, olddat, olddatsize, maxDepth); 
	ret = Add_Pts_To_Bucket(curnode, &pt[i], 1, maxDepth);
      }
    }
  return ret;
}

void Increase_Depth( qtree2_t *curnode)
{

  double up = curnode -> upper;
  double low = curnode -> lower;
  double left = curnode -> left;
  double right = curnode -> right;
  double midvert = (up + low) / 2.0;
  double midhoriz = (left + right) / 2.0;
  unsigned char depth = curnode -> depth + 1;
  char type = curnode -> datatype;
  curnode -> uleft = Create_Branch(up, midvert, left, midhoriz, curnode, 1, depth, type);
  curnode -> uright  = Create_Branch(up, midvert, midhoriz, right, curnode, 2, depth, type);
  curnode -> lleft =  Create_Branch(midvert, low, left, midhoriz, curnode, 4, depth, type) ;
  curnode -> lright =  Create_Branch(midvert, low, midhoriz, right, curnode, 3, depth, type);

  return;
}

qtree2_t *Create_Branch(double up, double low, double left, double right, qtree2_t *parent, char pos, unsigned char depth, char type)
{
  qtree2_t *toret = malloc(sizeof(qtree2_t));
  toret->upper = up;
  toret->left = left;
  toret->right = right;
  toret->lower = low;
  toret->parent = parent;
  //malloc doesn't initialize the allocated memory to 0s.
  toret->uleft = NULL;
  toret->uright = NULL;
  toret->lleft = NULL;
  toret->lright = NULL;
  toret->numdata = 0;
  toret->data = NULL;
  toret->pos = pos;
  toret->depth = depth;
  toret->datatype = type;
  return toret;
}

SEXP 
//Any handling of multiple data types/forms should be handled on the R side
R_Find_Neighbors_Pts(SEXP Rtree, SEXP Rnewx, SEXP Rnewy, SEXP Rk)
{
  qtree2_t *tree = (qtree2_t *) R_ExternalPtrAddr( GET_SLOT( Rtree, Rf_install( "ref" ) ) );

  double *x = REAL( Rnewx );
  double *y = REAL( Rnewy );
  qtree2_t *curnode;
  int k = INTEGER( Rk ) [ 0 ];
  int len = LENGTH( Rnewx );
  //double dists[ len*k ];
  double *dists = calloc(len*k, sizeof(double));
  point_t **chosen = calloc( len * k , sizeof(point_t*));
  point_t *newpt;
  point_t *oldpt;
  int filled;
  double tmpdist;
  int initcnt = 0;
  int tmpind;
  double curmaxdist;
  char pos;
  int inserted = 0;
  //initialized distances
  for(int i = 0; i< len*k; i++)
    dists[i] = DBL_MAX;
  
  for(int l = 0; l < len; l ++)
    {
      newpt = New_Point_Struct( x[ l ] , y[ l ], l );
      curnode = Descend_To_Bucket(tree, newpt, 1); //1 for pts;
      if( curnode -> numdata > 0)
	{
	  for (int i  = 0; i < curnode -> numdata; i++)
	    {
	      tmpdist = eucl_dist_pts(newpt, (point_t *) curnode->data[i]);
	      oldpt = (point_t *) curnode -> data [ i ] ;
	      //insert_dist( (double *) &dists, tmpdist, chosen, oldpt, k, l*k);
	      insert_dist( dists, tmpdist, chosen, oldpt, k, l*k);
	      initcnt ++;
	    }
	}

      //The rest are initialized to NULL pointers, which should be fine for now.

      //position of the maximum found distance relevant to this point
      tmpind = (l + 1 ) * k - 1;
      while (curnode -> parent != NULL)
	{
	  pos = curnode -> pos;
	  curnode = curnode -> parent;
	  curmaxdist= dists[ tmpind ];
	  //Harvest_KNN_Pts(curnode, pos, newpt -> x - curmaxdist, newpt -> x + curmaxdist, newpt -> y - curmaxdist, newpt -> y + curmaxdist, chosen, (double *) &dists, newpt, k, l*k);
	  Harvest_KNN_Pts(curnode, pos, newpt -> x - curmaxdist, newpt -> x + curmaxdist, newpt -> y - curmaxdist, newpt -> y + curmaxdist, chosen, dists, newpt, k, l*k);
	}	  
      free(newpt);
    }
  
  SEXP ans;

  //XXX At some point do we want the option to return the values instead of the indexes?
  PROTECT( ans = NEW_INTEGER( len * k ) );
  for (int i = 0; i < len * k; i++)
    INTEGER(ans) [ i ] = chosen[i]->index + 1; //+1 because it we start he counting at 0?
  free(chosen);
  free(dists);
  UNPROTECT(1);
  return ans;
}

qtree2_t *Descend_To_Bucket(qtree2_t *node, void *obj, int type)
{
  qtree2_t *ret;
  switch(type) 
    {
    case 1:
      //points
      ret = Descend_To_Bucket_Pts(node, (point_t*) obj);
      break;
    case 2:
      {
	rect_t *rec  = (rect_t *)obj;
	point_t *pt = New_Point_Struct((rec->left + rec->right)/2.0, (rec->low + rec->high)/2.0, 0);
	ret = Descend_To_Bucket(node, pt, 1);
	free(pt);
      }
    }
  return ret;
}

qtree2_t *Descend_To_Bucket_Pts(qtree2_t *node, point_t *pt)
{
  qtree2_t *curnode = node;
  while(curnode -> uleft != NULL)
      while(curnode -> uleft != NULL)
    {
      if (pt -> x < ( curnode -> lleft ) -> right)
	{
	  if ( pt -> y < (curnode -> lleft ) -> upper)
	    curnode = curnode -> lleft;
	  else
	    curnode = curnode -> uleft;
	} else 
	{
	  if ( pt -> y < (curnode -> lleft ) -> upper)
	    curnode = curnode -> lright;
	  else
	    curnode = curnode -> uright;
	}
    }
  return curnode;

}


void insert_dist(double *dists, double newdist, point_t **pts, point_t *newpt, int k, int start)
{
  int pos = start;
  int done = 0;
  while (!done && pos < start + k)
    {
      //if the point is already in there we're done!
      if (points_equal(pts[pos], newpt))
	break;

      if(dists[ pos ] >= newdist )
	done = 1;
      else
	pos ++;
    }
  //if we found a spot for it, insert it;
  if(done)
    {
      point_t *tmppt;
      double tmpdist;
      for (int curpos = pos; curpos < start + k; curpos ++)
	{
	  tmppt = pts[curpos]; //old value at curpos;
	  pts [ curpos ] = newpt; //newvalue in pts;
	  tmpdist = dists [ curpos ];
	  dists[ curpos ] = newdist;
	  if (curpos + 1 < start  + k)
	    {
	      newpt = tmppt; //old value going to curpos  + 1
	      newdist = tmpdist; 
	    }
	}
    }
  return;
}
/*
SEXP R_Get_KNN(SEXP Rtree, SEXP Rnewdat, SEXP Rnewcols, SEXP Rk)
//KNN only makes sense for points, so we can assume points throughout.
{
  qtree2_t *tree = (qtree2_t *) R_ExternalPtrAddr( GET_SLOT( Rtree, Rf_install( "ref" ) ) );

  int k = INTEGER( Rk ) [ 0 ];
  //assuming data.frame
  int *cols = INTEGER(Rnewcols);
  int newn = LENGTH( VECTOR_ELT( Rnewdat, cols[0] - 1));
  double dists[newn * k];
  point_t **found = calloc(newn * k, sizeof(point_t*));
  double *x = REAL( VECTOR_ELT( Rnewdat, cols[0] - 1));
  double *y = REAL( VECTOR_ELT( Rnewdat, cols[1] - 1));
  int pos, cnt = 0;
  point_t *newpt;
  qtree2_t * cur;
  double tmpmaxx = tree->right - tree->left;
  double tmpmaxy = tree->upper - tree->lower;
  for (int i =0; i < newn * k; i++)
    dists[i] = -1.0;
  for (int l = 0; l < newn; l++)
    {
      cnt = 0 ;
      newpt = New_Point_Struct(x[l], y[l], -1); //-1 because no meaningful index.
      cur = Descend_To_Bucket(tree, newpt);
      Harvest_KNN_Pts(cur, -1, newpt->x - tmpmaxx, newpt->x + tmpmaxx, newpt->y - tmpmaxy, newpt->y + tmpmaxy, found, &dists, newpt, l*k);
      
      }
*/
//	  Harvest_KNN_Pts(curnode, pos, newpt -> x - curmaxdist, newpt -> x + curmaxdist, newpt -> y - curmaxdist, newpt -> y + curmaxdist, &chosen, &dists, newpt, k, l*k);

void Harvest_KNN_Pts(qtree2_t *node, int excludepos, double leftbound, double rightbound, double lowbound, double highbound, point_t **pts, double *dists, point_t *newpt, int k, int start)
{

  point_t *oldpt;
  double dist;
  if(node -> numdata > 0)
    {
      for(int i = 0; i < node -> numdata ; i++)
	{
	  oldpt = (point_t *) node -> data[i];
	  dist = eucl_dist_pts(newpt, oldpt);
	  insert_dist(dists, dist, pts, oldpt, k, start);
	}
    } else {
    //each valid node has 0 or 4 children
    if(node -> uleft != NULL)
      {
	if(excludepos != 1)
	  {
	    if( check_bounds( node -> uleft, leftbound, rightbound, lowbound, highbound))
	      Harvest_KNN_Pts(node -> uleft, 0, leftbound, rightbound, lowbound, highbound, pts, dists, newpt, k, start);
	  }
	if(excludepos != 2)
	  {
	    if( check_bounds( node -> uright, leftbound, rightbound, lowbound, highbound))
	      Harvest_KNN_Pts(node -> uright, 0, leftbound, rightbound, lowbound, highbound, pts, dists, newpt, k, start);
	  }
	if(excludepos != 3)
	  {
	    if( check_bounds( node -> lright, leftbound, rightbound, lowbound, highbound))
	      Harvest_KNN_Pts(node -> lright, 0, leftbound, rightbound, lowbound, highbound, pts, dists, newpt, k, start);
	  }
	if(excludepos != 4)
	  {
	    if( check_bounds( node -> lleft, leftbound, rightbound, lowbound, highbound))
	      Harvest_KNN_Pts(node -> lleft, 0, leftbound, rightbound, lowbound, highbound, pts, dists, newpt, k, start);
	  }
	
      }
    
  }
  return;
}


int 
check_bounds(qtree2_t *node, double left, double right, double down, double up)
{
  int toret = 0;
  if ( ! (node -> left > right || node -> right < left || node -> upper < down || node -> lower > up) )
    {
      toret = 1;
    }
  return toret;
}



void
//XXX we need to be able to deal with the pts set as well!!
R_free_quad_tree( SEXP treeptr)
{
  qtree2_t *tree  = R_ExternalPtrAddr(treeptr);
  Free_QuadTree(tree);
  return;
}

void Free_QuadTree( qtree2_t *tree)
{

  if (tree -> parent != NULL)
    {
      qtree2_t *par= tree -> parent;
      switch (tree -> pos)
	{
	case 1:
	  par -> uleft = NULL;
	  break;
	case 2:
	  par -> uright = NULL;
	  break;
	case 3:
	  par -> lright = NULL;
	  break;
	case 4:
	  par -> lleft = NULL;
	  break;
	}
    }
  
  if( tree -> uleft != NULL)
    {
      Free_QuadTree( tree -> uleft);
      tree -> uleft = NULL;
      Free_QuadTree( tree -> uright);
      tree -> uright = NULL;
      Free_QuadTree( tree -> lright);
      tree -> lright = NULL;
      Free_QuadTree( tree -> lleft);
      tree -> lleft = NULL;
    }
  for(int i = 0; i < tree->numdata; i++)
    free(tree->data[i]);

  free(tree->data);
  
  free( tree );
  return;
}


int
Find_MaxDepth(qtree2_t *tree, unsigned char curdepth)
{
  if (tree -> uleft != NULL)
    {
      curdepth = Find_MaxDepth( tree -> uleft, curdepth);
      curdepth = Find_MaxDepth( tree -> uright, curdepth);
      curdepth = Find_MaxDepth( tree -> lright, curdepth);
      curdepth = Find_MaxDepth( tree -> lleft, curdepth);
    } else {
    curdepth = tree->depth > curdepth ? tree->depth : curdepth;
  }
  return curdepth;
}

SEXP
R_Find_MaxDepth(SEXP Rtree)
{
  qtree2_t *tree = (qtree2_t *) R_ExternalPtrAddr( GET_SLOT( Rtree, Rf_install( "ref" ) ) );
  int toret = Find_MaxDepth( tree, 0) ;
  return ScalarInteger(toret);
}

double eucl_dist_pts(point_t *pt1, point_t *pt2)
{
  return sqrt( pow(pt1->x - pt2->x, 2.0) + pow(pt1->y - pt2->y, 2));
}

SEXP
R_Rectangle_Lookup(SEXP Rtree, SEXP Rxlims, SEXP Rylims)
{
  double up, down, left, right;
  qtree2_t *tree = (qtree2_t *) R_ExternalPtrAddr( GET_SLOT( Rtree, Rf_install( "ref" ) ) );
  double *xlim, *ylim;
  xlim = REAL(Rxlims);
  ylim = REAL(Rylims);
  left = xlim[0];
  right = xlim[1];
  down = ylim[0];
  up = ylim[1];
  int size = 100;
  int structsize = get_struct_size(tree->datatype);
  int *found = malloc(size*structsize);
  int pos = 0;
  //fprintf(stderr, "\n"); fflush(stderr);
  Rectangle_Lookup(tree, left, right, down, up, &found, &pos, &size, tree->datatype);
  SEXP ans;
  PROTECT( ans = NEW_INTEGER( pos ) );
  for (int i = 0; i < pos; i ++)
    {
      switch(tree->datatype)
	{
	case 1:
	  INTEGER( ans )[ i ] = ( (point_t *)found)[i].index + 1; 
	  break;
	case 2:
	  INTEGER( ans )[ i ] = ( (rect_t *)found)[i].index + 1; 
	  break;
	}
    }
  UNPROTECT(1);
  free(found);
  return ans;
}

int get_struct_size(char type)
{
  int ret = 0;
  switch(type)
    {
    case 1:
      ret = sizeof(point_t);
      break;
    case 2:
      ret = sizeof(rect_t);
      break;
    }
  return ret;
}

void Rectangle_Lookup(qtree2_t *tree, double left, double right, double down, double up, int **found, int *pos, int *cursize, char type)
{
  

  void *cur;
  if(tree -> numdata)
    {
      for(int i = 0; i < tree -> numdata; i++)
	{
	  cur=  tree->data[i];
	  if(overlap(left, right, down, up, cur, type))
	    {
	      if( *pos >= *cursize)
		{
		  //fprintf(stderr, "\nOld Address: %lx", found);fflush(stderr);
		  Grow_ReturnArray(found, cursize, type);
		  //fprintf(stderr, "\nNew Address: %lx", found);fflush(stderr);
		}
	      switch(type)
		{
		case 1:
		  ((point_t*)*found)[*pos] =*((point_t *) cur);
		  break;
		case 2:
		  ((rect_t*)*found)[*pos] = *((rect_t *) cur);
		}
	      *pos += 1;
	      
	      
	    }
	}
    }
  if (tree -> uleft != NULL)
    {
      qtree2_t *tmptree = tree -> uleft;
      if (CheckBounds(tmptree, left, right, down, up))
	{
	  Rectangle_Lookup(tmptree, left, right, down, up, found, pos, cursize, type);
	}
      tmptree = tree -> uright;
      if (CheckBounds(tmptree, left, right, down, up))
	{
	  Rectangle_Lookup(tmptree, left, right, down, up, found, pos, cursize, type);
	}
      tmptree = tree ->lleft;
      if (CheckBounds(tmptree, left, right, down, up))
	{
	  Rectangle_Lookup(tmptree, left, right, down, up, found, pos, cursize, type);
	}
      tmptree = tree -> lright;
      if (CheckBounds(tmptree, left, right, down, up))
	{
	  Rectangle_Lookup(tmptree, left, right, down, up, found, pos, cursize, type);
	}
      

    }
}

int overlap(double left, double right, double down, double up, void *obj, char type)
{
  int toret = 0;
  switch(type)
    {
    case 1:
      {
	point_t *pt = (point_t *) obj;
	if (pt ->x >= left && pt->x <= right && pt->y >= down &&  pt->y <= up)
	  toret = 1;
	break;
      }
    case 2:
      {
	rect_t *cur = (rect_t *) obj;
	if ( (cur->left <= right && cur-> right >= left) && (cur-> low <= up && cur->high >= down))
	  toret = 1; 
	break;
      }
    }

  return toret;
}
      
 void Grow_ReturnArray(int **found, int *cursize, char type)
{
  int *toret = NULL;
  int oldsize = *cursize;
  int newsize;
  if(oldsize < 1000)
    newsize = oldsize *2;
  else
    newsize = 1.1*oldsize;
  
  //toret = realloc( found, newsize * sizeof(int*));
  
  int structsize = get_struct_size(type);
  toret = realloc(*found, newsize*structsize);  
  *cursize = newsize;
  //fprintf(stderr, "Grow_ReturnArray successful. New size: %d\n", newsize);fflush(stderr);
//toret = memcpy(topret, ar, oldsize * sizeof(int));
  *found = toret;
  //return toret;
  return;
}

int 
CheckBounds(qtree2_t *node, double left, double right, double down, double up)
{
  int toret = 0;
  if ( ! (node -> left > right || node -> right < left || node -> upper < down || node -> lower > up) )
    {
      toret = 1;
    }
  return toret;
}


void get_tree_attributes(qtree2_t *tree, int *curattr)
{
  //first int is number of nodes, second int is number of nodes with data, third int is maxDepth, fourth int is max bucketsize
  curattr[0]++;
  
  
  if(tree -> numdata > 0)
    {
      curattr[1]++;
      curattr[3] = (curattr[3] >= tree -> numdata) ? curattr[3] : tree -> numdata;
      curattr[2] = (curattr[2] >= tree -> depth) ? curattr[2] : (int) tree -> depth;
    }

  if(tree->uleft != NULL)
    {
      get_tree_attributes(tree -> uleft, curattr);
      get_tree_attributes(tree -> uright, curattr);
      get_tree_attributes(tree -> lright, curattr);
      get_tree_attributes(tree -> lleft, curattr);
    }
     
  return; 
}

static const R_CallMethodDef callMethods[]  = {
  {"Build_Quadtree_Rect", (DL_FUNC) &R_Build_Quadtree_Rect, 9},
  {"Build_Quadtree_Pt", (DL_FUNC) &R_Build_Quadtree_Pt, 7},
  {"Find_Neighbors_Pts", (DL_FUNC) &R_Find_Neighbors_Pts, 4},
  {"free_quad_tree", (DL_FUNC) &R_free_quad_tree, 1},
  {"Find_MaxDepth", (DL_FUNC) &R_Find_MaxDepth, 1},
  {"Rectangle_Lookup", (DL_FUNC) &R_Rectangle_Lookup, 3},
  {NULL, NULL, 0}
};

void
R_init_SearchTrees(DllInfo *info)
{
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
   R_useDynamicSymbols(info, FALSE);
   R_forceSymbols(info, TRUE);
}
