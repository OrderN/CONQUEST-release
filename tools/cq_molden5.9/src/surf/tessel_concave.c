/* Copyright amitabh varshney@cs.unc.edu. All Rights Reserved. March 16, 1994 */
/*--------------------------------------------------------------------------------
tessel_concave.c

This file contains the procedures to compute the tessellation of the concave 
spherical triangular patches of the molecular surface. 
----------------------------------------------------------------------------------*/

#include "surf.h"

void gen_spherical_recurse(VertexType **points, double *sq_side_len, float *center, double radius, int convex);

void gen_cuspy_concave_tri(VertexType point0, VertexType point1, VertexType point2, int flip, POINT center, double radius, Vector *plane_eq);

void gen_cuspy_spherical_recurse(VertexType **points, double *sq_side_len, float *center, Vector *plane_eq, double radius);
/*--------------------------------------------------------------------------------
gen_concave generates a concave spherical triangular patch of the molecular surface.
At present it can only handle simple cusps -- those which are due to two concave
triangles on the ends of a single torus. It detects things like whether there is 
a cusp or not and accordingly calls other routines.
----------------------------------------------------------------------------------*/
gen_concave(tor_pts, center, flip_verts, cusp_pts, probe, n, cur_torus, torus_cusp)
VertexType	**tor_pts, cusp_pts[2];
VertexType	center[2];
POINT		probe[2];
int		flip_verts, n;
struct Torus	*cur_torus;
int		torus_cusp;
{
  double	dradius;
  Vector	plane_eq;
  float		mid_probe_pos[3];

  dradius = Probe_radius;

  /* find the eq of the cusp generating (clipping) plane */
  VEC3_V_OP_V(mid_probe_pos, probe[0], +, probe[1]);
  VEC3_V_OP_S(mid_probe_pos, mid_probe_pos, *, HALF);
  VEC3_V_OP_V(plane_eq.coeff, probe[0], -, probe[1]);
  NORMALIZE3(plane_eq.coeff);
  plane_eq.coeff[3] = DOTPROD3(plane_eq.coeff, mid_probe_pos);

  if ((cur_torus->end_atoms[0] == cur_torus->prev->end_atoms[1]) ||
      (cur_torus->end_atoms[0] == cur_torus->prev->end_atoms[0]))
  { /* generate a concave triangle intersected with a plane */ 
    gen_cuspy_concave_tri(tor_pts[0][0],tor_pts[1][0],center[0], flip_verts, probe[0], 
		  dradius,&plane_eq);
    if (torus_cusp) /* generate the extra tri for 2-pt torus cusp */
    {
      gen_cuspy_concave_tri(tor_pts[1][0],cusp_pts[0],center[0],flip_verts, probe[0], 
		  dradius,&plane_eq);
    }
  }
  else
  { gen_concave_tri(tor_pts[0][0],tor_pts[1][0],center[0], flip_verts, probe[0], 
     		  dradius);
    if (torus_cusp) /* generate the extra tri for 2-pt torus cusp */
    { 
      gen_concave_tri(tor_pts[1][0],cusp_pts[0],center[0],flip_verts, probe[0], dradius);
    }
  }

  if ((cur_torus->end_atoms[1] == cur_torus->next->end_atoms[0]) ||
      (cur_torus->end_atoms[1] == cur_torus->next->end_atoms[1]))
  { /* generate a concave triangle intersected with a plane */
    gen_cuspy_concave_tri(tor_pts[1][n-1],tor_pts[0][n-1],center[1],flip_verts,probe[1],
		  dradius,&plane_eq);
    if (torus_cusp) /* generate the extra tris for 2-pt torus cusp */
    { 
      gen_cuspy_concave_tri(cusp_pts[1],tor_pts[1][n-1],center[1],flip_verts, probe[1], 
		  dradius,&plane_eq);
    }
  }
  else
  { gen_concave_tri(tor_pts[1][n-1],tor_pts[0][n-1],center[1],flip_verts,probe[1],
		  dradius);
    if (torus_cusp) /* generate the extra tris for 2-pt torus cusp */
    { 
      gen_concave_tri(cusp_pts[1],tor_pts[1][n-1],center[1],flip_verts, probe[1], 
		  dradius);
    }
  }
}

/*--------------------------------------------------------------------------------
gen_concave_tri generates a concave spherical triangular patch of the molecular 
surface.  This routine assumes that there are no cusps to be taken care of.
----------------------------------------------------------------------------------*/
gen_concave_tri(point0, point1, point2, flip, center, radius)
VertexType	point0, point1, point2;
int		flip;
POINT		center;
double		radius;
{
  float		fcenter[3];

  VEC3_ASN_OP(fcenter, =, center);

  if (flip)
  { 
    gen_sphere_tris(&point2, &point1, &point0, fcenter, radius, 0);
  }
  else
    gen_sphere_tris(&point0, &point1, &point2, fcenter, radius, 0);
}


/*--------------------------------------------------------------------------------
gen_cuspy_concave_tri generates a concave spherical triangular patch which might 
have a cusp.
----------------------------------------------------------------------------------*/
void gen_cuspy_concave_tri(point0, point1, point2, flip, center, radius, plane_eq)
VertexType	point0, point1, point2;
int		flip;
POINT		center;
double		radius;
Vector		*plane_eq;
{
  float		fcenter[3];
  double	dtemp;

  VEC3_ASN_OP(fcenter, =, center);

  dtemp = DOTPROD3(plane_eq->coeff, fcenter) - plane_eq->coeff[3];
  if (FP_EQ_EPS(dtemp, 0, EPS)) /* center lies on the plane so perturb plane */
    plane_eq->coeff[3] += 2*EPS;
    

  if (flip)
  { 
    gen_cuspy_sphere_tris(point2, point1, point0, fcenter, radius, plane_eq);
  }
  else
    gen_cuspy_sphere_tris(point0, point1, point2, fcenter, radius, plane_eq);

}

/*--------------------------------------------------------------------------------
gen_cuspy_sphere_tri generates a spherical triangular patch which might have a 
cusp. This packs the points into an array to be used by recursive calls to 
gen_cuspy_spherical_recurse.
----------------------------------------------------------------------------------*/
gen_cuspy_sphere_tris(pt0, pt1, pt2, center, radius, plane_eq)
VertexType		pt0, pt1, pt2;
float			center[3];
double			radius;
Vector			*plane_eq;
{
  double		sq_side_len[3];
  VertexType*		points[3];
  
  points[0] = &pt0; points[1] = &pt1; points[2] = &pt2;

  /* find the squares of the lengths of the sides */
  sq_side_len[0] = SQ_DIST3(points[0]->Coord, points[1]->Coord);
  sq_side_len[1] = SQ_DIST3(points[1]->Coord, points[2]->Coord);
  sq_side_len[2] = SQ_DIST3(points[2]->Coord, points[0]->Coord);

  /* recursively generate the tris */
  gen_cuspy_spherical_recurse(points, sq_side_len, center, plane_eq, radius);
}

/*--------------------------------------------------------------------------------
gen_cuspy_spherical_recurse generates a spherical triangular patch which might have a 
cusp. This recursively subdivides the patch while clipping it with its intersection 
with plane represented by plane_eq.
----------------------------------------------------------------------------------*/
void gen_cuspy_spherical_recurse(points, sq_side_len, center, plane_eq, radius)
VertexType		*points[3];
double			sq_side_len[3];
float			center[3];
double			radius;
Vector			*plane_eq; 
{ 
  double		new_side_len[3];
  double		max_len, dtemp;
  int			i, j, k, max_side;
  int			itempc, itemp[3], num_fine_points;
  VertexType		new_point1, new_point2;
  VertexType		*temp_points[3];

  if ((sq_side_len[0] < EPS) && 
      (sq_side_len[1] < EPS) && 
      (sq_side_len[2] < EPS)) 
  { return;
  }

  /* evaluate which side of the plane the points lie wrt the probe center */
  /* point[i] lies on the same side as probe center, itemp[i] = 1 */
  dtemp = DOTPROD3(plane_eq->coeff, center) - plane_eq->coeff[3];
  itempc = (dtemp > 0)? 1: 0;
  for(i = 0; i < 3; i++)
  { dtemp = DOTPROD3(plane_eq->coeff, points[i]->Coord)-plane_eq->coeff[3];
    itemp[i] = (dtemp > EPS)? 1: ((dtemp < -EPS)? 0: itempc);
    itemp[i] = (itempc == itemp[i]);
  }
  
  num_fine_points = itemp[0] + itemp[1] + itemp[2];
  max_len = FMAX((FMAX(sq_side_len[0], sq_side_len[1])), sq_side_len[2]);
  max_side = (max_len==sq_side_len[0])?0:((max_len==sq_side_len[1])? 1 : 2);

  if (num_fine_points == 0) return;      /* all points are on the wrong side */

  if (max_len <= Max_Tess_Len_Sq)
  { /* sufficiently small so display but first check to see whether it
       intersects the plane or not and if so, clip it accordingly */
    if (num_fine_points == 3)
    { gen_tris(points[0], points[1], points[2]);
      return;
    }
    else if (num_fine_points == 1) /* two points are on the wrong side */
      i = (itemp[0] == 1)? 0 : ((itemp[1] == 1)? 1: 2);
    else /* if (num_fine_points == 2) one point is on the wrong side */
      i = (itemp[0] == 0)? 0 : ((itemp[1] == 0)? 1: 2);
    /* compute the two intersection points with the plane */
    j = (i+1)%3;
    find_line_seg_plane_int(new_point1.Coord,points[i]->Coord,points[j]->Coord,plane_eq->coeff);
    VEC3_V_OP_V(new_point1.Normal, center, -, new_point1.Coord);
    NORMALIZE3(new_point1.Normal);
    for(k = 0; k < 3; k++)
      new_point1.Coord[k] = center[k] - new_point1.Normal[k]*radius;

    j = (i+2)%3;
    find_line_seg_plane_int(new_point2.Coord,points[i]->Coord,points[j]->Coord,plane_eq->coeff);
    VEC3_V_OP_V(new_point2.Normal, center, -, new_point2.Coord);
    NORMALIZE3(new_point2.Normal);
    for(k = 0; k < 3; k++)
      new_point2.Coord[k] = center[k] - new_point2.Normal[k]*radius;
   
    if (num_fine_points == 1)
      gen_tris(points[i], &new_point1, &new_point2);
    else /* num_fine_points == 2 */
    { j = (i+1)%3; k = (i+2)%3;
      gen_tris(points[j], points[k], &new_point2);
      gen_tris(points[j], &new_point2, &new_point1);
    }
  }
  else	/* split into two triangles and recurse */
  { /* compute new point */
    i = max_side;
    j = (i+1)%3;

    VEC3_V_OP_V(new_point1.Normal,points[i]->Normal,+,points[j]->Normal);
    NORMALIZE3(new_point1.Normal);
    VEC3_V_OP_V_OP_S(new_point1.Coord, center, -, new_point1.Normal, *, radius);
    
    /* recurse on the two triangles */
    k = (i+2)%3;
    temp_points[0] = &new_point1;
    new_side_len[0] = SQ_DIST3(new_point1.Coord, points[j]->Coord);
    temp_points[1] = points[j]; 
    new_side_len[1] = sq_side_len[j];
    temp_points[2] = points[k]; 
    new_side_len[2] = SQ_DIST3(new_point1.Coord, points[k]->Coord);
    gen_cuspy_spherical_recurse(temp_points, new_side_len, center, plane_eq, radius);
  
    temp_points[0] = points[k]; new_side_len[0] = sq_side_len[k];
    temp_points[1] = points[i]; 
    new_side_len[1] = SQ_DIST3(points[i]->Coord, new_point1.Coord);
    temp_points[2] = &new_point1; /* new_side_len[2] is same as above */
    gen_cuspy_spherical_recurse(temp_points, new_side_len, center, plane_eq, radius);
  }
}

/*--------------------------------------------------------------------------------
gen_sphere_tris generates a spherical triangular patch which does not have a cusp.
----------------------------------------------------------------------------------*/
/* generate the sphere triangles */
gen_sphere_tris(pt0, pt1, pt2, center, radius, convex)
VertexType		*pt0, *pt1, *pt2;
float			center[3];
double			radius;
int			convex;
{
  double		sq_side_len[3];
  VertexType*		points[3];
  
  points[0] = pt0; points[1] = pt1; points[2] = pt2;
  /* find the squares of the lengths of the sides */
  sq_side_len[0] = SQ_DIST3(points[0]->Coord, points[1]->Coord);
  sq_side_len[1] = SQ_DIST3(points[1]->Coord, points[2]->Coord);
  sq_side_len[2] = SQ_DIST3(points[2]->Coord, points[0]->Coord);

  /* recursively generate the tris */
  gen_spherical_recurse(points, sq_side_len, center, radius, convex);
}

/*--------------------------------------------------------------------------------
gen_spherical_recurse generates a spherical triangular patch by recursively 
subdividing till all its edge lengths are below a user specified threshold.
----------------------------------------------------------------------------------*/
void gen_spherical_recurse(points, sq_side_len, center, radius, convex)
VertexType		*points[3];
float			center[3];
double			radius;
double			sq_side_len[3];
int			convex;
{ 
  double		new_side_len[3];
  double		max_len;
  int			i, j, k, max_side;
  VertexType		new_point;
  VertexType		*temp_points[3];

  max_len = FMAX((FMAX(sq_side_len[0], sq_side_len[1])), sq_side_len[2]);
  max_side = (max_len==sq_side_len[0])?0:((max_len==sq_side_len[1])? 1 : 2);

  if ((sq_side_len[0] < EPS) && 
      (sq_side_len[1] < EPS) && 
      (sq_side_len[2] < EPS)) 
  { 
    return;
  }

  if (max_len > Max_Tess_Len_Sq)
				/* split into two triangles and recurse */
  { /* compute new point */
    i = max_side; j = (i+1)%3;
    VEC3_V_OP_V(new_point.Normal,points[i]->Normal,+,points[j]->Normal);
    NORMALIZE3(new_point.Normal);

    if (convex)
      VEC3_V_OP_V_OP_S(new_point.Coord, center, +, new_point.Normal, *, radius)
    else /* generating a concave triangle */
      VEC3_V_OP_V_OP_S(new_point.Coord, center, -, new_point.Normal, *, radius)
    
    /* recurse on the two triangles */
    k = (i+2)%3;
    temp_points[0] = &new_point;
    new_side_len[0] = SQ_DIST3(new_point.Coord, points[j]->Coord);
    temp_points[1] = points[j]; 
    new_side_len[1] = sq_side_len[j];
    temp_points[2] = points[k]; 
    new_side_len[2] = SQ_DIST3(new_point.Coord, points[k]->Coord);
    gen_spherical_recurse(temp_points, new_side_len, center, radius, convex);

    temp_points[0] = points[k]; new_side_len[0] = sq_side_len[k];
    temp_points[1] = points[i]; 
    new_side_len[1] = SQ_DIST3(points[i]->Coord, new_point.Coord);
    temp_points[2] = &new_point; /* new_side_len[2] is same as above */
    gen_spherical_recurse(temp_points, new_side_len, center, radius, convex);

  }
  else			     /* sufficiently small, so simply display */
    gen_tris(points[0], points[1], points[2]);
}
