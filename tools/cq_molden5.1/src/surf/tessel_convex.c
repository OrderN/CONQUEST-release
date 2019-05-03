/* Copyright amitabh varshney@cs.unc.edu. All Rights Reserved. March 16, 1994 */
/*--------------------------------------------------------------------------------
tessel.c

This file contains the procedures to compute the tessellation of the smooth 
molecular surface. This computation is done on an atom-by-atom basis so that 
in future, if the patches are to be colored based on atom-type, it would be 
easy to do so.
----------------------------------------------------------------------------------*/

#include "surf.h"

/*--------------------------------------------------------------------------------
gen_convex generates convex spherical patches that join torus ends to the 
regular atom surface.
Here, p = sph_pts, q = tor_pts
----------------------------------------------------------------------------------*/
gen_convex(comp_verts, p, probe_centers, q, num_p, num_q, same_order, flip, atom_id, full_torus)
POINT		*comp_verts, *probe_centers; 
VertexType	*p, *q;
int		num_p, num_q;
int		same_order;
int		flip, full_torus;
int		atom_id;
{
  Big_Point	*sph_side, *torus_side, *big_p;
  int		num_sph_side, guess_num_sph_side;
  POINT		diff0, diff1;
  double	angle = 0;
  int		i;

  for(guess_num_sph_side = 0, i = -1; i < num_p; i++)
  { if (i == -1) {VEC3_V_OP_V(diff0, q[0].Coord, -, atoms[atom_id].center);}
    else {VEC3_V_OP_V(diff0, p[i].Coord, -, atoms[atom_id].center);}
    if (i == num_p-1) 
    {VEC3_V_OP_V(diff1, q[num_q-1].Coord, -, atoms[atom_id].center);}
    else {VEC3_V_OP_V(diff1, p[i+1].Coord, -, atoms[atom_id].center);}
    NORMALIZE3(diff0); NORMALIZE3(diff1);
    angle = DOTPROD3(diff0, diff1); 
    if ((angle >= 1.0) && (angle <= 1.0 + GP_EPS)) angle = 0;
    else if ((angle <= -1.0) && (angle >= -1.0-GP_EPS)) angle = PI;
    else angle = acos(angle);
    guess_num_sph_side += 
       1+2*((int)(angle*atoms[atom_id].radius/Max_Tess_Len+FCEIL));
  }
  
  ALLOCN(sph_side, Big_Point, guess_num_sph_side);
  ALLOCN(torus_side, Big_Point, num_q);
  ALLOCN(big_p, Big_Point, num_p);

  /* initialize torus side */
  for(i = 0; i < num_q; i++)
  { memcpy(&(torus_side[i].vt), &q[i], sizeof(VertexType));
    /* memcpy(torus_side[i].pt, probe_centers[i], sizeof(POINT)); */
    VEC3_V_OP_V(torus_side[i].tes_dir,probe_centers[i],-,atoms[atom_id].tes_origin);
    NORMALIZE3(torus_side[i].tes_dir);
  }

  /* initialize sphere side */
  { for(i = 0; i < num_p; i++)
    { memcpy(&(big_p[i].vt), &p[i], sizeof(VertexType));
      VEC3_V_OP_V(big_p[i].tes_dir,comp_verts[i],-,atoms[atom_id].tes_origin);
      NORMALIZE3(big_p[i].tes_dir);
    }
  }

  initialize_r(big_p, sph_side, q, num_p, &num_sph_side, num_q, same_order, atom_id, full_torus);

  if (num_sph_side >= guess_num_sph_side)
  { fprintf(stderr,"gen_convex(): Too many points on the sphere\n");
    fprintf(stderr,"ft %d num_sph_side %d, guess_num_sph_side %d angle %lf\n", 
     full_torus, num_sph_side, guess_num_sph_side, angle);
    fprintf(stderr,"center (%f %f %f) radius %f\n", atoms[atom_id].center[X], 
     atoms[atom_id].center[Y], atoms[atom_id].center[Z], atoms[atom_id].radius);
    fflush(stderr);
    free(sph_side);
    return;
  }

  
  if (full_torus)
    get_full_torus_mesh(sph_side, num_sph_side, torus_side, num_q,
			flip, atom_id);
  else
    get_partial_torus_mesh(sph_side, num_sph_side, torus_side, num_q, 
			flip, atom_id);
  free(sph_side);
  free(torus_side);
  free(big_p);
} 

/*--------------------------------------------------------------------------------
initialize_r generates a set of points on the sphere to help in tessellation with
the torus.
Here, big_p = sph_pts, q = tor_pts, r = rotated sph_pts 
----------------------------------------------------------------------------------*/
initialize_r(big_p, r, q, num_p, num_r, num_q, same_order, atom_id, full_torus)
VertexType	*q;
Big_Point	*big_p, *r;
int		num_p, num_q, *num_r;
int		same_order, full_torus;
int		atom_id;
{ 
  int		i, j, min_angle_id;
  double	angle, min_angle;
  Big_Point	*temp_array;
  float		ave_q[3], plane_eq[4];

  if (same_order)
  { memcpy(&r[0], &big_p[0], sizeof(Big_Point)); 
    for(i =0, j = 1; i < num_p-1; i++)
      gen_arc_recurse(&big_p[i], &big_p[i+1], r, &j, atom_id, 1);
  }
  else
  { memcpy(&r[0], &big_p[num_p-1], sizeof(Big_Point)); 
    for(j = 1, i = num_p-1; i > 0; i--)
      gen_arc_recurse(&big_p[i], &big_p[i-1], r, &j, atom_id, 1);
  }
  *num_r = j;

  if (full_torus)
  { 
    /* compute the average of q's */
    VEC3_ZERO(ave_q);
    for(i = 0; i < num_q; i++)
      VEC3_ASN_OP(ave_q, +=, q[i].Coord);
    VEC3_V_OP_S(ave_q, ave_q, /, num_q);

    /* find the eq of the plane passing through all the tor_pts */
    find_plane_eq(q[0].Coord, q[1].Coord, q[2].Coord, plane_eq);

    /* find the vertex in r that is closest to q[0] in terms of 
       angle wrt ave_q */
    for(i = min_angle_id = 0, min_angle = HUGE_VAL; i < *num_r; i++)
    { angle = find_angle(r[i].vt.Coord, q[0].Coord, plane_eq, ave_q);
      if (min_angle > angle) 
      { min_angle = angle; min_angle_id = i; }
    }
 
    if (min_angle_id > 0)
    { /* re-adjust r so that r[0], matches with r[min_angle_id] */
      ALLOCN(temp_array, Big_Point, *num_r);

      *num_r -= 1;
      for(i = 0; i < *num_r; i++)
      { j = (i + min_angle_id)%(*num_r);
	memcpy(&temp_array[i], &r[j], sizeof(Big_Point));
      }
      
      memcpy(&temp_array[*num_r], &temp_array[0], sizeof(Big_Point));
      *num_r += 1;

      for(i = 0; i < *num_r; i++)
        memcpy(&r[i], &temp_array[i], sizeof(Big_Point));
      free(temp_array);
    }
  }
}

/*--------------------------------------------------------------------------------
gen_linear_recurse generates a set of points on the sphere recursively, such that
they lie along a great circle of the atom/sphere and connect the points p0 and p1 
on the sphere. The criterion for stopping the subdivision is that the length of the 
longest edge in this arc be smaller than a user-specified max edge length.
----------------------------------------------------------------------------------*/
gen_linear_recurse(p0, p1, r, j, center, radius, convex)
VertexType	p0, p1;
VertexType	*r;
int		*j;
float		center[3];
float		radius;
int		convex;
{ 
  double	len;
  VertexType	new_point;

  len = SQ_DIST3(p0.Coord, p1.Coord);
  if (len > Max_Tess_Len_Sq)
  { /* compute new point */
    VEC3_V_OP_V(new_point.Normal, p0.Normal, +, p1.Normal);
    NORMALIZE3(new_point.Normal);

    if (convex)
      VEC3_V_OP_V_OP_S(new_point.Coord, center, +, new_point.Normal, *, radius)
    else
      VEC3_V_OP_V_OP_S(new_point.Coord, center, -, new_point.Normal, *, radius) 

    gen_linear_recurse(p0, new_point, r, j, center, radius, convex);
    gen_linear_recurse(new_point, p1, r, j, center, radius, convex);
  }
  else
  { VEC3_ASN_OP(r[*j].Normal, =, p1.Normal);
    VEC3_ASN_OP(r[*j].Coord, =, p1.Coord);
    *j += 1;
  }
}

/*--------------------------------------------------------------------------------
gen_flipped_sphere_tris generates sphere tris based on whether they are to be 
flipped or not.
----------------------------------------------------------------------------------*/
gen_flipped_sphere_tris(p0, p1, p2, atom_id, flip)
Big_Point		*p0, *p1, *p2;
int			atom_id;
int			flip;
{
  if (flip)
    new_gen_sphere_tris(p0, p1, p2, atom_id);
  else
    new_gen_sphere_tris(p0, p2, p1, atom_id);
}


/*--------------------------------------------------------------------------------
new_gen_sphere_tris generates sphere tris. this calls multi_gen_spherical_recurse
----------------------------------------------------------------------------------*/
new_gen_sphere_tris(bp0, bp1, bp2, atom_id)
Big_Point		*bp0, *bp1, *bp2;
int			atom_id;
{

  /* recursively generate the tris */
  { POINT		tempv;
    double		planarity, ext_rad = atoms[atom_id].radius + Probe_radius;
    Big_Point		mid_pt;

    CROSSPROD3(tempv,bp0->tes_dir, bp1->tes_dir);
    planarity = DOTPROD3(tempv, bp2->tes_dir);
    if (FP_EQ_EPS(planarity, 0, EPS)) /* three pts coplanar wrt tes_origin */
				      /* compute the triangle mid-point */
    { VEC3_V_OP_V_OP_V(mid_pt.vt.Normal, bp0->vt.Normal, +, bp1->vt.Normal, 
						         +, bp2->vt.Normal);
      NORMALIZE3(mid_pt.vt.Normal);
      VEC3_V_OP_V_OP_S(mid_pt.vt.Coord, atoms[atom_id].center, +,
			      mid_pt.vt.Normal, *, atoms[atom_id].radius);
      VEC3_V_OP_V(tempv, atoms[atom_id].center, -, atoms[atom_id].tes_origin);
      VEC3_V_OP_V_OP_S(mid_pt.tes_dir, tempv, +, mid_pt.vt.Normal, *,ext_rad);
      multi_gen_spherical_recurse(bp0, bp1, &mid_pt, atom_id);
      multi_gen_spherical_recurse(bp1, bp2, &mid_pt, atom_id);
      multi_gen_spherical_recurse(bp2, bp0, &mid_pt, atom_id);
    }
    else 
      multi_gen_spherical_recurse(bp0, bp1, bp2, atom_id); 
  }
}

/*--------------------------------------------------------------------------------
multi_gen_spherical_recurse recursively generates sphere tris. this divides into
2, 3, or 4 triangles in one call depending on how many sides are longer than the 
threshold value.
----------------------------------------------------------------------------------*/
multi_gen_spherical_recurse(pt0, pt1, pt2, atom_id)
Big_Point		*pt0, *pt1, *pt2;
int			atom_id;
{ 
  int			i, j, k;
  int			valid_point = 0, flags[3];
  Big_Point		*tri_pt[3], mid_pt[3];
  double		sq_side;

  /* initialization */
  flags[0] = 0x001; flags[1] = 0x010; flags[2] = 0x100;
  tri_pt[0] = pt0; tri_pt[1] = pt1; tri_pt[2] = pt2;

  /* find the squares of the lengths of the sides */
  for(i = 0; i < 3; i++)
  { j = (i+1)%3; 
    sq_side = SQ_DIST3(tri_pt[i]->vt.Coord, tri_pt[j]->vt.Coord);
    if (sq_side > Max_Tess_Len_Sq) /* find the mid point */
      valid_point  |= flags[i];
  }

  if (valid_point == 0x000) /* each side sufficiently small, simply display */
  { 
    gen_tris(&(tri_pt[0]->vt), &(tri_pt[1]->vt), &(tri_pt[2]->vt)); 
  }
  else /* need to compute one or more edge mid-points */
  { if(valid_point&flags[0]) find_mid_pt(&mid_pt[0],tri_pt[0],tri_pt[1],atom_id);
    if(valid_point&flags[1]) find_mid_pt(&mid_pt[1],tri_pt[1],tri_pt[2],atom_id);
    if(valid_point&flags[2]) find_mid_pt(&mid_pt[2],tri_pt[2],tri_pt[0],atom_id); 

    if (valid_point == 0x111) /* each side big so make 4 tris */
    { multi_gen_spherical_recurse(tri_pt[0], &mid_pt[0], &mid_pt[2], atom_id);
      multi_gen_spherical_recurse(&mid_pt[0], tri_pt[1], &mid_pt[1], atom_id);
      multi_gen_spherical_recurse(&mid_pt[1], tri_pt[2], &mid_pt[2], atom_id);
      multi_gen_spherical_recurse(&mid_pt[0], &mid_pt[1],&mid_pt[2], atom_id);
    }
    else if ((valid_point==0x001)||(valid_point==0x010)||(valid_point==0x100))
    { i = (valid_point == 0x001)? 0 : ((valid_point == 0x010)? 1 : 2);
      j = (i + 1)%3; k = (i + 2)%3;
      multi_gen_spherical_recurse(tri_pt[i], &mid_pt[i], tri_pt[k], atom_id);
      multi_gen_spherical_recurse(&mid_pt[i], tri_pt[j], tri_pt[k], atom_id);
    }
    else if ((valid_point==0x110)||(valid_point==0x101)||(valid_point==0x011))
    { i = (valid_point == 0x110)? 0 : ((valid_point == 0x101)? 1 : 2);
      j = (i + 1)%3; k = (i + 2)%3;
      multi_gen_spherical_recurse(&mid_pt[j], tri_pt[k], &mid_pt[k], atom_id);
      if (SQ_DIST3(mid_pt[j].vt.Coord, tri_pt[i]->vt.Coord) >
  	SQ_DIST3(mid_pt[k].vt.Coord, tri_pt[j]->vt.Coord))
      { multi_gen_spherical_recurse(&mid_pt[j],&mid_pt[k], tri_pt[j], atom_id);
        multi_gen_spherical_recurse(&mid_pt[k], tri_pt[i], tri_pt[j], atom_id);
      }
      else
      { multi_gen_spherical_recurse(&mid_pt[j],&mid_pt[k], tri_pt[i], atom_id);
        multi_gen_spherical_recurse(&mid_pt[j], tri_pt[i], tri_pt[j], atom_id);
      }
    }
  }
}

/*--------------------------------------------------------------------------------
find_mid_pt 
----------------------------------------------------------------------------------*/
find_mid_pt(mid_pt, tri_pt0, tri_pt1, atom_id)
Big_Point		*mid_pt, *tri_pt0, *tri_pt1;
int			atom_id;
{ 
   double		b, c, t, t1, t2, ext_rad;
   POINT		diff;

   VEC3_V_OP_V(mid_pt->tes_dir, tri_pt0->tes_dir, +, tri_pt1->tes_dir);
   NORMALIZE3(mid_pt->tes_dir);
   VEC3_V_OP_V(diff, atoms[atom_id].tes_origin, -, atoms[atom_id].center);
   c = DOTPROD3(diff, diff);
   if (FP_EQ_EPS(c, 0, EPS)) /* tes_origin == center => t = ext_rad; */
     VEC3_ASN_OP(mid_pt->vt.Normal, =, mid_pt->tes_dir)
   else
   { ext_rad = atoms[atom_id].radius + Probe_radius;
     c -= ext_rad*ext_rad;
     b = DOTPROD3(mid_pt->tes_dir, diff);
     t = b*b - c;
     if (t < 0) 
     { printf("find_mid_pt(): t = %g < 0\n",t);
       exit(-1);
     }
     t = sqrt(t);
     t1 = -b + t; t2 = -b -t;
     t = FMAX(t1, t2);
     ext_rad = 1/ext_rad;
     mid_pt->vt.Normal[X] = (diff[X] + t*mid_pt->tes_dir[X])*ext_rad;
     mid_pt->vt.Normal[Y] = (diff[Y] + t*mid_pt->tes_dir[Y])*ext_rad;
     mid_pt->vt.Normal[Z] = (diff[Z] + t*mid_pt->tes_dir[Z])*ext_rad;

     VEC3_V_OP_V(diff, mid_pt->vt.Normal, -, tri_pt0->vt.Normal);
     if (DOTPROD3(diff, diff) == 0.0)
     { printf("Stop here\n");
     }
   }
   VEC3_V_OP_V_OP_S(mid_pt->vt.Coord, atoms[atom_id].center, +, 
		    mid_pt->vt.Normal, *, atoms[atom_id].radius);
    
}

/*--------------------------------------------------------------------------------
gen_tris simply adds a tri to the exisiting surface. The only reason this peculiar
routine exists is to act as a single catching place for further enhancements like 
colors based on atom triangles.
----------------------------------------------------------------------------------*/
gen_tris(pt0, pt1, pt2)
VertexType		*pt0, *pt1, *pt2;
{
   add_tri(pt0, pt1, pt2, Current_atom); 
}


/*--------------------------------------------------------------------------------
add_tri adds the triangle to the list of existing triangles.
----------------------------------------------------------------------------------*/
add_tri(pt0, pt1, pt2, atype)
VertexType		*pt0, *pt1, *pt2;
int			atype;
{

  if (Write_Option == 1) /* write it out immediately */
  { 
    write_archive_tri(pt0, pt1, pt2, atype); 
    Num_polys++;
  }
  else if (Write_Option == 2) /* store the triangle */
  { memcpy(&verts[Num_polys][0], pt0, sizeof(VertexType));
    memcpy(&verts[Num_polys][1], pt1, sizeof(VertexType));
    memcpy(&verts[Num_polys][2], pt2, sizeof(VertexType));
    atom_type[Num_polys] = atype; 
   
    Num_polys++; 

    if (Num_polys >= Max_Gp_Polys)
    { fprintf(stderr,"\n add_tri(): Too many triangles -- Max_Tri (%d) exceeded\n", Max_Gp_Polys);
      exit(-1);
    }
  }
  else /* just count it */
    Num_polys++;
}
