/* Copyright amitabh varshney@cs.unc.edu. All Rights Reserved. March 16, 1994 */
/*--------------------------------------------------------------------------------
tessel_torus.c

This file contains the procedures to compute the tessellation of a toroidal patch.
Every atom computes its half of the toroidal patch. So in future, if the toroidal
patches are to be colored based on atom-type, it should be real easy.
----------------------------------------------------------------------------------*/

#include "surf.h"

/*--------------------------------------------------------------------------------
get_full_torus_mesh is called from gen_convex() to smoothly tessellate the convex
spherical patch where it merges with the toroidal patch. The only reason this 
routine exists is to ensure that there are no cracks in the mesh between the 
spherical and the toroidal patches. This routine generates a mesh for merging 
with a complete torus -- i.e. a torus whose first and last points are the same. 
Such toroidal patches arise, for instance between two close-enough atoms. The toroidal 
patches between three close-enough atoms are _not_ "complete" as they are terminated 
in concave triangular patches -- so they will not be handled by this. 
----------------------------------------------------------------------------------*/
get_full_torus_mesh(sph, num_sph, torus, num_torus, flip, atom_id)
Big_Point	*sph, *torus;
int		num_sph, num_torus;
int		flip;
int		atom_id;
{
  int		i;
  float		plane_eq[4]; /* ax + by + cz = d */
  int		s_index, next_s_index;
  int		t_index, next_t_index;
  float		ave_t[3];

  /* compute the average of torus points */
  VEC3_ZERO(ave_t);
  for(i = 0; i < num_torus; i++)
    VEC3_ASN_OP(ave_t, +=, torus[i].vt.Coord);
  VEC3_V_OP_S(ave_t, ave_t, /, num_torus);

  /* find the eq of the plane passing through all the tor_pts */
  find_plane_eq(torus[0].vt.Coord,torus[1].vt.Coord,torus[2].vt.Coord,plane_eq);

  s_index = t_index = 0;

  if ((num_sph > 1)&&(num_torus > 1))
  { do
    { next_s_index = s_index+1;
      next_t_index = t_index+1;

      if (find_angle(sph[next_s_index].vt.Coord, torus[t_index].vt.Coord, plane_eq, ave_t) < 
	  find_angle(sph[s_index].vt.Coord, torus[next_t_index].vt.Coord, plane_eq, ave_t))
      { gen_flipped_sphere_tris(&torus[t_index], &sph[next_s_index], &sph[s_index], atom_id, flip);
        s_index = next_s_index;
      }
      else
      { gen_flipped_sphere_tris(&sph[s_index], &torus[t_index], &torus[next_t_index], atom_id, flip);
        t_index = next_t_index;
      }
    }
    while((s_index != num_sph-1)&&(t_index != num_torus-1));
  }

  if ((s_index == num_sph-1)&&(t_index != num_torus-1))
  { 
    do
    { next_s_index = (s_index + 1)% num_sph;
      next_t_index = t_index + 1;

      if (find_angle(sph[next_s_index].vt.Coord, torus[t_index].vt.Coord, plane_eq, ave_t) < 
	  find_angle(sph[s_index].vt.Coord, torus[next_t_index].vt.Coord, plane_eq, ave_t))
      { s_index = next_s_index;
      }
      else
      { gen_flipped_sphere_tris(&sph[s_index], &torus[t_index], &torus[next_t_index], atom_id, flip);
        t_index = next_t_index;
      }
    } while(t_index != num_torus-1);
  }
  else if ((t_index == num_torus-1)&&(s_index != num_sph-1))
  {
    do
    { next_s_index = s_index + 1;
      next_t_index = (t_index + 1)%num_torus;

      if (find_angle(sph[next_s_index].vt.Coord, torus[t_index].vt.Coord, plane_eq, ave_t) < 
	  find_angle(sph[s_index].vt.Coord, torus[next_t_index].vt.Coord, plane_eq, ave_t))
      { gen_flipped_sphere_tris(&torus[t_index], &sph[next_s_index], &sph[s_index], atom_id, flip);
        s_index = next_s_index;
      }
      else
      { t_index = next_t_index;
      }
    } while(s_index != num_sph-1);
  }
}
/*--------------------------------------------------------------------------------
get_partial_torus_mesh is called from gen_convex() to smoothly tessellate the convex
spherical patch where it merges with the toroidal patch. The only reason this 
routine exists is to ensure that there are no cracks in the mesh between the 
spherical and the toroidal patches. This routine generates a mesh for merging 
with a partial torus -- i.e. a torus whose first and last points are not the same. 
Such toroidal patches arise, for instance between three close-enough atoms. The
toroidal patches between three close-enough atoms are "partial" and not "complete" 
as they are terminated in concave triangular patches.
----------------------------------------------------------------------------------*/
get_partial_torus_mesh(sph, num_sph, torus, num_torus, flip, atom_id)
Big_Point	*sph, *torus;
int		num_sph, num_torus;
int		flip;
int		atom_id;
{ 
   int		i;
   int		s_index, next_s_index;
   int		t_index, next_t_index;
   float	plane_eq[4];

   find_plane_eq(torus[0].vt.Coord,torus[1].vt.Coord,torus[2].vt.Coord,plane_eq);
   s_index = t_index = 0;

   if ((num_sph > 1)&&(num_torus > 1))
   { do
     { next_s_index = s_index+1;
       next_t_index = t_index+1;

       if (sph_dist(torus[t_index].vt.Coord, sph[next_s_index].vt.Coord, atoms[atom_id].center) < 
	   sph_dist(sph[s_index].vt.Coord, torus[next_t_index].vt.Coord, atoms[atom_id].center))
       { gen_flipped_sphere_tris(&torus[t_index], &sph[next_s_index], &sph[s_index], atom_id, flip); 
         s_index = next_s_index;
       }
       else
       { gen_flipped_sphere_tris(&sph[s_index], &torus[t_index], &torus[next_t_index], atom_id, flip); 
         t_index = next_t_index;
       }
     }
   while((s_index != num_sph-1)&&(t_index != num_torus-1)); 
   }

   if (s_index == num_sph-1)
   { for(i = t_index; i < num_torus-1; i++)
     { gen_flipped_sphere_tris(&sph[s_index], &torus[i], &torus[i+1], atom_id, flip);
     }
   }
   else /* t_index == num_torus-1 */
   { for(i = s_index; i < num_sph-1; i++)
     { gen_flipped_sphere_tris(&torus[t_index], &sph[i+1], &sph[i], atom_id, flip);
     }
   }
}

/*--------------------------------------------------------------------------------
gen_torus is used to smoothly tessellate the toroidal patches.
 radius : how "thick" the torus is (Probe Radius for us)
----------------------------------------------------------------------------------*/
gen_torus(tor_pts, probe_centers, num_verts, flip, radius, cusp)
VertexType	**tor_pts;
POINT		*probe_centers;
int		num_verts;
int		flip;
double		radius;
int		cusp;
{
  int		i, j, k, l, m;
  int		old_j;
  VertexType	**r;
  POINT		diff0, diff1;
  float		fcenter[3], fradius;
  double	angle;
  int		n; /* no of torus verts generated from a single center pos*/
  int		*col;

  VEC3_V_OP_V(diff0, tor_pts[0][0].Coord, -, probe_centers[0]);
  VEC3_V_OP_V(diff1, tor_pts[1][0].Coord, -, probe_centers[0]);
  NORMALIZE3(diff0); NORMALIZE3(diff1);
  angle = DOTPROD3(diff0, diff1); 
  if (FP_EQ_EPS(angle, 1.0, GP_EPS)) angle = 0;
  else if (FP_EQ_EPS(angle, -1.0, GP_EPS)) angle = PI;
  else angle = acos(angle);
  
  n = 1 + 2*((int)(radius*angle/Max_Tess_Len + FCEIL));
  if (n < 2) n = 2;

  ALLOCMAT(r, 2, n, VertexType);
  ALLOCN(col, int, n);
  fradius = radius;

  VEC3_ASN_OP(r[0][0].Normal, =, tor_pts[0][0].Normal);
  VEC3_ASN_OP(r[0][0].Coord, =, tor_pts[0][0].Coord);

  j = 1; VEC3_ASN_OP(fcenter, =, probe_centers[0]);
  gen_linear_recurse(tor_pts[0][0],tor_pts[1][0],r[0],&j,fcenter,fradius,0);

  if (j > n)
  { printf("gen_torus(): Too many torus verts %d (limit %d)\n",j, n);
    return ;
  }
  old_j = j;

  for(i = 1; i < num_verts; i++)
  { j = 1; k = i & 0x01; l = !k; VEC3_ASN_OP(fcenter, =, probe_centers[i]);
    VEC3_ASN_OP(r[k][0].Normal, =, tor_pts[0][i].Normal);
    VEC3_ASN_OP(r[k][0].Coord, =, tor_pts[0][i].Coord);
    gen_linear_recurse(tor_pts[0][i],tor_pts[1][i],r[k],&j,fcenter,fradius,0);
    if (j > n)
    { printf("gen_torus(): Too many torus verts %d (limit %d)\n",j, n);
      return ;
    }

    for(m = 0; m < j-1; m++)
    { if (flip) gen_tris(&r[l][m], &r[k][m+1], &r[k][m]);
      else      gen_tris(&r[l][m], &r[k][m], &r[k][m+1]);
      
      if (((!cusp)&&(j <= old_j))||(m != j-2))
      { if (flip) gen_tris(&r[l][m], &r[l][m+1], &r[k][m+1]);
        else      gen_tris(&r[l][m], &r[k][m+1], &r[l][m+1]);
      }
    }
    if (j < old_j)
    { for(m = j-1; m < old_j - 1; m++)
      { if (flip) gen_tris(&r[l][m], &r[l][m+1], &r[k][j-1]);
        else      gen_tris(&r[l][m], &r[k][j-1], &r[l][m+1]);
      }
    }
    /* should take care of the case when j > old_j + 1 here */

    old_j = j;

  }
  free(col);
  FREEMAT(r, 2, n);
}
