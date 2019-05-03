/* Copyright amitabh varshney@cs.unc.edu. All Rights Reserved. March 16, 1994 */
/*--------------------------------------------------------------------------------
utils.c

This file contains some miscellaneous procedures that are useful in the general
three-dimensional geometry. Some of these have been culled from Graphics Gems,
others, I have written. 

A lot of A-intersects-B type of routines are in here.

The reference for Graphics Gems I is:
"Graphics Gems" edited by Andrew Glassner, Academic Press, 1990.
----------------------------------------------------------------------------------*/

#include "surf.h"

/*----------------------------------------------------------------------------------
find_plane_eq finds the equation of the plane given by pt0, pt1, and pt2.
This returns a normalized plane equation that is in the form: 
[0]x + [1]y + [2]z = [3]
where [3] is +ve.
----------------------------------------------------------------------------------*/
find_plane_eq(pt0, pt1, pt2, plane_eq)
float			pt0[3], pt1[3], pt2[3];
float			plane_eq[4];
{

  POINT			diff1, diff2;

  VEC3_V_OP_V(diff1, pt0, -, pt1);
  VEC3_V_OP_V(diff2, pt0, -, pt2);
  CROSSPROD3(plane_eq, diff1, diff2);
  NORMALIZE3(plane_eq);
  plane_eq[3] = DOTPROD3(plane_eq, pt0);
  if (plane_eq[3] < 0)
    VEC_NEG(plane_eq, plane_eq, 4);
}

/*----------------------------------------------------------------------------------
compute_tri_plane_int  computes the intersection of three planes 
(uses method given in Graphics Gems I, pp 305)
----------------------------------------------------------------------------------*/
compute_tri_plane_int(plane0, plane1, plane2, result)
Vector		plane0, plane1, plane2;
double		result[3];
{
  
  double	denom;
  double	tempv[3];

  /* first compute the 3x3 determinant of plane normals */
  denom = DET3x3(plane0.coeff[X], plane0.coeff[Y], plane0.coeff[Z], 
  		 plane1.coeff[X], plane1.coeff[Y], plane1.coeff[Z], 
  		 plane2.coeff[X], plane2.coeff[Y], plane2.coeff[Z]);

  if (fabs(denom) < EPSILON) 
  { printf("compute_tri_plane_int(): denom close to zero %lf\n",denom);
    denom = (denom > 0)? EPSILON: -EPSILON;
  }
  
  CROSSPROD3(tempv, plane1.coeff, plane2.coeff);
  VEC3_V_OP_S(result, tempv, *, plane0.coeff[3]);

  CROSSPROD3(tempv, plane2.coeff, plane0.coeff);
  VEC3_V_OP_V_OP_S(result, result, +, tempv, *, plane1.coeff[3]);

  CROSSPROD3(tempv, plane0.coeff, plane1.coeff);
  VEC3_V_OP_V_OP_S(result, result, +, tempv, *, plane2.coeff[3]);

  VEC3_V_OP_S(result, result, /, denom);
}


/*----------------------------------------------------------------------------------
find_ray_sphere_int  computes the intersection of a ray with a sphere.
uses method given in Graphics Gems I, pp 388-389.
the direction vector of the ray need not be normalized.
----------------------------------------------------------------------------------*/
find_ray_sphere_int(int_point, ray_pt, ray_dir, center, radius)
float			int_point[3], center[3];
float			ray_pt[3], ray_dir[3];
float			radius;
{

  double		EO[3];
  double		t, v, disc;
  double		one_by_sq_dir_mag;

  /* if the ray_pt is almost on the sphere, return it as the intersection pt*/
  VEC3_V_OP_V(EO, center, -, ray_pt);
  disc = radius*radius - DOTPROD3(EO, EO);
  if (FP_EQ_EPS(disc, 0, GP_EPS*GP_EPS))
  { VEC3_ASN_OP(int_point, =, ray_pt);
    return;
  }

  one_by_sq_dir_mag = 1.0/DOTPROD3(ray_dir, ray_dir);
  v = DOTPROD3(EO, ray_dir);
  disc += v*v*one_by_sq_dir_mag;
  if ((disc < 0) && (disc > -GP_EPS)) disc = 0.0;
  if (disc < 0)	/* no intersection */
  { printf("find_ray_sphere_int(): disc = %lf sqrt probs for atom center %lf %lf %lf\n", 
	    disc, center[X], center[Y], center[Z]);
    VEC3_ZERO(int_point);
  }
  else
  { t = v*one_by_sq_dir_mag - sqrt(disc*one_by_sq_dir_mag);
    VEC3_V_OP_V_OP_S(int_point, ray_pt, +, ray_dir, *, t);
  }
}

/*----------------------------------------------------------------------------------
find_line_plane_int  computes the intersection of a line given by a point and a 
direction vector with a plane. 
This assumes:
  (a) plane eq is normalized and is in the form [0]x + [1]y + [2]z = [3] 
  (b) the direction vector of the line is normalized 
----------------------------------------------------------------------------------*/
find_line_plane_int(int_point, line_pt, line_dir, plane_eq)
float		int_point[3], line_pt[3], line_dir[3];
float		plane_eq[4];
{
  double	dtemp;
  int		k;

  dtemp = plane_eq[3] - DOTPROD3(plane_eq, line_pt);
  dtemp /= DOTPROD3(plane_eq, line_dir);
  for(k = 0; k < 3; k++)
    int_point[k] = line_pt[k] + dtemp*line_dir[k];
}

/*----------------------------------------------------------------------------------
find_line_seg_plane_int  computes the intersection of a line given by two points 
with a plane. This assumes plane eq is normalized and is in the form:
[0]x + [1]y + [2]z = [3] 
----------------------------------------------------------------------------------*/
find_line_seg_plane_int(int_point, point0, point1, plane_eq)
float		int_point[3], point0[3], point1[3];
float		plane_eq[4];
{
  double	dtemp;

  VEC3_V_OP_V(int_point, point1, -, point0);
  dtemp = plane_eq[3] - DOTPROD3(plane_eq, point0);
  dtemp /= DOTPROD3(plane_eq, int_point);
  VEC3_V_OP_V_OP_S(int_point, point0, +, int_point, *, dtemp);
}

/*----------------------------------------------------------------------------------
find_angle returns the angle between (p0'-org) and (p1'-org), where p0' 
and p1' are projections of p0 and p1 on the plane given by plane_eq.
This assumes org lies on plane_eq.
----------------------------------------------------------------------------------*/
double
find_angle(p0, p1, plane_eq, org)
float		org[3], p0[3], p1[3];
float		plane_eq[4];
{
  float		p0_bar[3], p1_bar[3];
  float		diff0[3], diff1[3];
  float		line_dir[3];
  double	result;

  VEC3_ASN_OP(line_dir, =, plane_eq);

  /* parallel project p0 and p1 onto the plane */
  find_line_plane_int(p0_bar, p0, line_dir, plane_eq);
  find_line_plane_int(p1_bar, p1, line_dir, plane_eq);
  
  /* find angle */
  VEC3_V_OP_V(diff0, p0_bar, -, org);
  NORMALIZE3(diff0);
  VEC3_V_OP_V(diff1, p1_bar, -, org);
  NORMALIZE3(diff1);
  
  result = 2 - (1 + DOTPROD3(diff0, diff1));
  /* if angle is zero, result = 0, if angle is PI, result = 2 */

  return(result);
}

/*----------------------------------------------------------------------------------
sph_dist 

This assumes that "p0" and "p1" are points on a sphere whose center is "center"
----------------------------------------------------------------------------------*/
double
sph_dist(p0, p1, center)
float		p0[3], p1[3], center[3];
{
  double	result;
  double	diff0[3], diff1[3];

  VEC3_V_OP_V(diff0, p0, -, center);
  VEC3_V_OP_V(diff1, p1, -, center);
  result = -DOTPROD3(diff0, diff1);
  return(result);
}


