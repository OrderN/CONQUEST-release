/* Copyright amitabh varshney@cs.unc.edu. All Rights Reserved. March 16, 1994 */
/*--------------------------------------------------------------------------------
tessel_patches.c

This file acts as a link between tessel_cases.c and the trio of tessel_concave.c,
tessel_convex.c, and tessel_torus.c. 
----------------------------------------------------------------------------------*/
#include "surf.h"

/*--------------------------------------------------------------------------------
gen_patches_I_II generates the patches for the cases I and II. All the cases I-IV
have been explained in tessel_cases.c
----------------------------------------------------------------------------------*/
gen_patches_I_II(face_atoms, comp_verts, num_verts, concave_center, cir_radius, inv_transf, 
		 angle, cur_torus, torus_cusp)
int		face_atoms[2];
POINT		*comp_verts;
int		num_verts;
VertexType	concave_center[2];
double		cir_radius;
MatrixD4	inv_transf;
double		angle;
struct Torus	*cur_torus;
int		torus_cusp;
{ 
  VertexType	*sph_pts; /* points on the sphere */
  VertexType	**tor_pts; /* points that lie on torus boundaries */
  VertexType	cusp_pts[2]; /* midpoint of the torus in case of a 2pt torus cusp */
  Gp_Atom	*atom_ptr[2];
  POINT		*probe_centers;
  int		i, j, num_tor_verts;
  int		flip_verts;
  double	torus_length, max_radius;
  double	cur_angle, seg_angle;
  double	dpoint1[4], dpoint2[4];
  POINT		concave_probe[2];

  /* initialization */
  ALLOCN(sph_pts, VertexType, num_verts);
  for(i = 0; i < 2; i++)
    atom_ptr[i] = &atoms[face_atoms[i]];

  /* get the points on the sphere for face_atoms[0]  - sph_pts */
  for(j = 0; j < num_verts; j++)/* for each extended-rad component vert*/
    find_sphere_points(sph_pts, j, comp_verts[j], atom_ptr[0]);

  /* find the total arc-length along the seam where the torus 
     meets the neighboring atom. 
     ideally we should have taken the radius to be 
     (cir_radius*atom_ptr[0]->radius)/(atom_ptr[0]->radius + Probe_radius)
     However, this will give different torus_lengths and hence different
     num_tor_verts when the radii of the neighboring atoms are different.
     Hence we just take the maximum of the two radii for computing the 
     torus_length.

  */
  max_radius = FMAX(atom_ptr[0]->radius, atom_ptr[1]->radius);
  torus_length = angle*(cir_radius*max_radius)/(max_radius + Probe_radius);
  num_tor_verts = 1 + (int)(torus_length/Max_Tess_Len+ FCEIL); 

  /* place a lower bound on num_tor_verts based on angle */
  if (angle < 2*PI/3) 
  { if (num_tor_verts < 2) num_tor_verts = 2; }
  else if (angle < 4*PI/3) 
  { if (num_tor_verts < 3) num_tor_verts = 3; }
  else if (angle < TWO_PI) 
  { if (num_tor_verts < 4) num_tor_verts = 4; }

  seg_angle = angle/(num_tor_verts-1);

  /* find the probe_centers */
  ALLOCN(probe_centers, POINT, num_tor_verts);
  dpoint1[Z] = ZERO; dpoint1[3] = ONE;
  for(j = 0, cur_angle = 0; j < num_tor_verts; j++)
  { dpoint1[X] = cos(cur_angle);
    dpoint1[Y] = sin(cur_angle);
    MATVECMULT(dpoint2, inv_transf, dpoint1, 4, 4);
    VEC3_ASN_OP(probe_centers[j], =, dpoint2);
    cur_angle += seg_angle;
  }

  /* get rid of accumulation of errors */
  VEC3_ASN_OP(probe_centers[num_tor_verts-1], =, comp_verts[num_verts-1]);
 
  /* find the torus points on the end bordering face_atoms[0] in tor_pts[0],
     and those in the center of the torus, in tor_pts[1] */
  ALLOCMAT(tor_pts, 2, num_tor_verts, VertexType);
  find_torus_points(tor_pts, cusp_pts, num_tor_verts, atom_ptr, 
		    probe_centers, cir_radius, torus_cusp);

  /* make sure that sph_pts[0] == tor_pts[0][0] and 
		    sph_pts[num_verts-1] ==  tor_pts[0][num_tor_verts-1]
  */
  memcpy(&tor_pts[0][0], &sph_pts[0], sizeof(VertexType));
  memcpy(&tor_pts[0][num_tor_verts-1], &sph_pts[num_verts-1], sizeof(VertexType));

  /* generate the patches */
  flip_verts = find_order(comp_verts, num_verts, face_atoms);
  
  gen_torus(tor_pts, probe_centers, num_tor_verts, flip_verts, Probe_radius, torus_cusp);

  /* sph_pts[0] == tor_pts[0][0] and sph_pts[num_verts-1] ==  tor_pts[0][num_tor_verts-1]
     so we start from sph_pts[1] and specify the no of sph verts as only num_verts-2 */

  gen_convex(comp_verts + 1, &sph_pts[1], probe_centers, tor_pts[0], num_verts-2, 
	     num_tor_verts, TRUE, flip_verts, face_atoms[0], FALSE);

  VEC3_ASN_OP(concave_probe[0], =, comp_verts[0]);
  VEC3_ASN_OP(concave_probe[1], =, comp_verts[num_verts-1]);
  gen_concave(tor_pts, concave_center, flip_verts, cusp_pts,
	      concave_probe, num_tor_verts, cur_torus, torus_cusp);

  /* free up allocated stuff */
  free(sph_pts);
  free(probe_centers);
  FREEMAT(tor_pts, 2, num_tor_verts);
}

/*--------------------------------------------------------------------------------
gen_patches_III generates the patches for the case III. All the cases I-IV have 
been explained in tessel_cases.c 
This assumes that the first and last comp_verts are the same 
----------------------------------------------------------------------------------*/
gen_patches_III(face_atoms, comp_verts, num_verts, cir_center, cir_radius, 
		inv_transf, torus_cusp)
int		face_atoms[2];
POINT		*comp_verts;
POINT		cir_center;
double		cir_radius;
int		num_verts;
MatrixD4	inv_transf;
int		torus_cusp;
{ 
  VertexType	*sph_pts; /* points on the sphere */
  VertexType	**tor_pts; /* points that lie on torus boundaries */
  VertexType	cusp_pts[2]; /* redundant for this case */
  Gp_Atom	*atom_ptr[2];
  int		i, j, num_tor_verts;
  POINT		*probe_centers;
  double	dpoint1[4], dpoint2[4];
  POINT		tmp_normal;
  int		flip_verts, same_verts_order;
  double	torus_length, max_radius;
  double	cur_angle, seg_angle;

  /* initialization */
  ALLOCN(sph_pts, VertexType, num_verts);

  for(i = 0; i < 2; i++)
    atom_ptr[i] = &atoms[face_atoms[i]];

  /* get the points on the sphere for face_atoms[0]  - sph_pts */
  for(j = 0; j < num_verts; j++)/* for each extended-rad component vert*/
    find_sphere_points(sph_pts, j, comp_verts[j], atom_ptr[0]);

  /* find the total arc-length along the seam where the torus 
     meets the face_atoms[0]
  */
  /* since it is a type III patch, it is a total torus, so the 
     torus_length is simply 2*pi*(suitably_weighted_radius) (expl of 
     this given in gen_patches_I_II() 
  */
  max_radius = FMAX(atom_ptr[0]->radius, atom_ptr[1]->radius);
  torus_length = TWO_PI*(cir_radius*max_radius)/(max_radius + Probe_radius);
  num_tor_verts = 1 + (int)(torus_length/Max_Tess_Len+ FCEIL);
  /* place a lower bound on num_tor_verts */
  if (num_tor_verts < 5) num_tor_verts = 5;
  seg_angle = TWO_PI/(num_tor_verts-1);

  /* find the probe centers */
  ALLOCN(probe_centers, POINT, num_tor_verts);
  dpoint1[Z] = ZERO; dpoint1[3] = ONE;
  for(j = cur_angle = 0; j < num_tor_verts; j++)
  { dpoint1[X] = cos(cur_angle);
    dpoint1[Y] = sin(cur_angle);
    MATVECMULT(dpoint2, inv_transf, dpoint1, 4, 4);
    VEC3_ASN_OP(probe_centers[j], =, dpoint2);
    cur_angle += seg_angle;
  }

  /* find the torus points on the end bordering face_atoms[0] in tor_pts[0], 
     and those in the center of the torus, in tor_pts[1] */
  ALLOCMAT(tor_pts, 2, num_tor_verts, VertexType);
  find_torus_points(tor_pts, cusp_pts, num_tor_verts, atom_ptr, probe_centers, 
		    cir_radius, torus_cusp);

  /* find the order of comp_verts wrt to the order of the tor_verts */
  VEC3_V_OP_V(dpoint1, comp_verts[0], -, cir_center);
  NORMALIZE3(dpoint1);
  VEC3_V_OP_V(dpoint2, comp_verts[1], -, cir_center);
  NORMALIZE3(dpoint2);
  CROSSPROD3(tmp_normal, dpoint1, dpoint2);
  /* if the Z-axis (in the transformation matrix inv_transf) is in the same 
     direction as tmp_normal, the vertices of the torus [0..n] occur in the 
     same order as the comp_verts[0..n].
  */
  GETCOL(Z, inv_transf, dpoint1, 3);
  same_verts_order = (DOTPROD3(tmp_normal, dpoint1) > 0)? TRUE: FALSE;

  /* generate the patches */
  flip_verts = find_order(comp_verts, num_verts, face_atoms);
  if (!same_verts_order) flip_verts = !flip_verts;

  gen_torus(tor_pts, probe_centers, num_tor_verts, flip_verts, Probe_radius, torus_cusp);

  gen_convex(comp_verts, sph_pts, probe_centers, tor_pts[0], 
	     num_verts, num_tor_verts, same_verts_order, flip_verts, 
	     face_atoms[0], TRUE);

  /* free up allocated stuff */
  free(sph_pts);
  free(probe_centers);
  FREEMAT(tor_pts, 2, num_tor_verts);
}

/*--------------------------------------------------------------------------------
gen_patches_IV generates the patches for the case IV. All the cases I-IV have been
explained in tessel_cases.c 
----------------------------------------------------------------------------------*/
gen_patches_IV(atom_id, comp_verts, num_verts)
int		atom_id;
POINT		*comp_verts;
int		num_verts;
{
  int		i, j, k;
  int		flip_verts;
  int		guess_num;
  Big_Point	*sp_points;
  Big_Point	*comp_sph_pts; /* points on the sphere */
  Big_Point	center;
  POINT		ave_center;

  /* find guess_num */
  guess_num = 50;
  ALLOCN(sp_points, Big_Point, guess_num);
  ALLOCN(comp_sph_pts, Big_Point, num_verts);

  for(i = 0; i < num_verts; i++)
  { find_sphere_points(&(comp_sph_pts[i].vt),0,comp_verts[i], &atoms[atom_id]);
    VEC3_V_OP_V(comp_sph_pts[i].tes_dir, comp_verts[i], -,
		atoms[atom_id].tes_origin);
    NORMALIZE3(comp_sph_pts[i].tes_dir);
  }  

  memcpy(&sp_points[0], &comp_sph_pts[0], sizeof(Big_Point));

  for(k = 1, j = 0; j < num_verts; j++)
  { i = (j+1)%num_verts;
    gen_arc_recurse(&comp_sph_pts[j], &comp_sph_pts[i], sp_points,&k,atom_id, 1);
  }

  if (k >= guess_num)
  { printf("gen_patches_IV(): too small guess_num %d, k = %d\n", 
            guess_num, k);
    exit(-1);
  }

  VEC3_ZERO(ave_center);
  for(j = 0; j < k; j++) VEC3_ASN_OP(ave_center, +=, sp_points[j].tes_dir);
  VEC3_V_OP_S(ave_center, ave_center, /, k);
  VEC3_ASN_OP(ave_center, +=, atoms[atom_id].tes_origin);
  find_sphere_points(&(center.vt), 0, ave_center, &atoms[atom_id]);
  VEC3_V_OP_V(center.tes_dir, ave_center, -, atoms[atom_id].tes_origin);
  NORMALIZE3(center.tes_dir);
  
  flip_verts = find_vertex_order(&(center.vt), sp_points, k);

  /* recopy the sphere points for the comp verts and start off with
     them for better triangulation 
  */
  for(i = 0; i < num_verts; i++)
    memcpy(&sp_points[i], &comp_sph_pts[i], sizeof(Big_Point));
   
  k = num_verts+1;
  memcpy(&sp_points[k-1], &sp_points[0], sizeof(Big_Point));
#if 0
  for(i = 0; i < k-1; i++)
  { if (flip_verts)
      gen_tris(&sp_points[i+1].vt, &sp_points[i].vt, &center.vt);
    else
      gen_tris(&sp_points[i].vt, &sp_points[i+1].vt, &center.vt);
  }
#else
  for(i = 0; i < k-1; i++)
  { if (flip_verts)
      new_gen_sphere_tris(&sp_points[i+1], &sp_points[i], &center, atom_id);
    else
      new_gen_sphere_tris(&sp_points[i], &sp_points[i+1], &center, atom_id);
  }
#endif

  free(comp_sph_pts);
  free(sp_points);
}

/*--------------------------------------------------------------------------------
gen_arc_recurse recursively generates an arc joining bp0 and bp1 till each of the
segments has a length lesser than a user specified maximum edge length. 
----------------------------------------------------------------------------------*/
gen_arc_recurse(bp0, bp1, sp_points, j, atom_id, convex)
Big_Point       *bp0, *bp1, *sp_points;
int             *j;
int		atom_id;
int             convex;
{
  double        len;
  Big_Point	new_bp;

  len = SQ_DIST3(bp0->vt.Coord, bp1->vt.Coord);
  if (len > Max_Tess_Len_Sq)
  { /* compute new point */
    find_mid_pt(&new_bp, bp0, bp1, atom_id);

    gen_arc_recurse(bp0, &new_bp, sp_points, j, atom_id, convex);
    gen_arc_recurse(&new_bp, bp1, sp_points, j, atom_id, convex);
  }
  else
  { memcpy(&sp_points[*j], bp1, sizeof(Big_Point));
    *j += 1;
  }
}

/*--------------------------------------------------------------------------------
find_sphere_points finds the point of contact between an atom and a probe sphere
given the centers and radii of the atom and the probe sphere.
----------------------------------------------------------------------------------*/
find_sphere_points(sph_pts, j, comp_vert, atom_ptr)
VertexType	*sph_pts;
int		j;
POINT		comp_vert;
Gp_Atom		*atom_ptr;
{ 
  double	dtemp;
  float		ext_rad, ray_dir[3];
  float		tempv[3];

  ext_rad = atom_ptr->radius + Probe_radius;
  /* find the coords on the extended radius sphere for atom i */
  VEC3_V_OP_V(ray_dir, atom_ptr->tes_origin, -, comp_vert);
  /* NORMALIZE3(ray_dir); */
  VEC3_ASN_OP(tempv, =, comp_vert);
  find_ray_sphere_int(sph_pts[j].Coord,tempv,ray_dir,atom_ptr->center,ext_rad);
  /* find the normal */
  VEC3_V_OP_V(sph_pts[j].Normal, sph_pts[j].Coord, -, atom_ptr->center);
  /* NORMALIZE3(sph_pts[j].Normal); */
  VEC3_V_OP_S(sph_pts[j].Normal, sph_pts[j].Normal, /, ext_rad);
  
  /* find the coords on surface of sphere */
  dtemp = atom_ptr->radius;
  VEC3_V_OP_V_OP_S(sph_pts[j].Coord, atom_ptr->center,+,sph_pts[j].Normal, *, dtemp);
}

/*--------------------------------------------------------------------------------
find_torus_points computes the points defining a torus, given a sequence of probe 
sphere positions as it rolls defining the torus. 
tor_pts[0] contains the torus points that are on the side of the sphere
tor_pts[1] contains the torus points that are along the center of the torus 
	   (along the plane that will slice the torus into two symmetric torii)
----------------------------------------------------------------------------------*/
find_torus_points(tor_pts, mid_cusp_pts, num_pts, atom_ptr, probe_centers, 
		  cir_radius, cusp)
VertexType	**tor_pts;
VertexType	mid_cusp_pts[2];
int		num_pts;
Gp_Atom		*atom_ptr[2];
POINT		*probe_centers;
double		cir_radius;
int		cusp;
{ 
  double	dtemp[2];
  int		j;
  POINT		cusp_pt;

  if (cusp)
  { dtemp[0] = sqrt(Probe_radius*Probe_radius - cir_radius*cir_radius);
    dtemp[1] = sqrt(SQ(Probe_radius+atom_ptr[0]->radius) - SQ(cir_radius)) - dtemp[0];
    VEC3_V_OP_V(cusp_pt, atom_ptr[1]->center, -, atom_ptr[0]->center);
    NORMALIZE3(cusp_pt);
    VEC3_V_OP_S(cusp_pt, cusp_pt, *, dtemp[1]);
    VEC3_ASN_OP(cusp_pt, +=, atom_ptr[0]->center);
  }
  dtemp[0] = ONE/(atom_ptr[0]->radius + Probe_radius);
  dtemp[1] = ONE/(atom_ptr[1]->radius + Probe_radius);
  for(j = 0; j < num_pts; j++)
  { /* generate the torus points corresponding to the given probe_center */
    VEC3_V_OP_V(tor_pts[0][j].Normal, probe_centers[j], -, atom_ptr[0]->center);
    VEC3_V_OP_S(tor_pts[0][j].Normal, tor_pts[0][j].Normal, *, dtemp[0]);
    if (cusp)
      VEC3_V_OP_V(tor_pts[1][j].Normal, probe_centers[j], -, cusp_pt) 
    else
    { VEC3_V_OP_V(tor_pts[1][j].Normal, probe_centers[j], -, atom_ptr[1]->center);
      VEC3_V_OP_S(tor_pts[1][j].Normal, tor_pts[1][j].Normal, *, dtemp[1]);
      VEC3_ASN_OP(tor_pts[1][j].Normal, +=, tor_pts[0][j].Normal);
    }
    NORMALIZE3(tor_pts[1][j].Normal);

    /* find the coord */
    VEC3_V_OP_V_OP_S(tor_pts[0][j].Coord, probe_centers[j], -,tor_pts[0][j].Normal,*,Probe_radius);
    if (cusp)
      VEC3_ASN_OP(tor_pts[1][j].Coord, =, cusp_pt)
    else
      VEC3_V_OP_V_OP_S(tor_pts[1][j].Coord, probe_centers[j], -, tor_pts[1][j].Normal,*,Probe_radius);
  }
  /* compute the mid cusp pts if required */
  /* if there is a 2-pt torus cusp, then mid_cusp_pts[0] is the vertex that is on the
     side of tor_pts[0/1][0] and lies on the plane between atoms[0] and atoms[1]
     Sly, mid_cusp_pts[1] is on the side of tor_pts[0/1][num_pts-1]
  */
  if (cusp)
  { /* mid_cusp_pt 0 */
    VEC3_V_OP_V(mid_cusp_pts[0].Normal, probe_centers[0], -, atom_ptr[1]->center);
    VEC3_V_OP_V_OP_S(mid_cusp_pts[0].Normal,tor_pts[0][0].Normal,+,mid_cusp_pts[0].Normal,*,dtemp[1]);
    NORMALIZE3(mid_cusp_pts[0].Normal);
    VEC3_V_OP_V_OP_S(mid_cusp_pts[0].Coord,probe_centers[0], -, mid_cusp_pts[0].Normal,*,Probe_radius);
    /* mid_cusp_pt 1 */
    VEC3_V_OP_V(mid_cusp_pts[1].Normal, probe_centers[num_pts-1], -, atom_ptr[1]->center);
    VEC3_V_OP_V_OP_S(mid_cusp_pts[1].Normal, tor_pts[0][num_pts-1].Normal,+,mid_cusp_pts[1].Normal,*,dtemp[1]);
    NORMALIZE3(mid_cusp_pts[1].Normal);
    VEC3_V_OP_V_OP_S(mid_cusp_pts[1].Coord, probe_centers[num_pts-1],-,mid_cusp_pts[1].Normal,*,Probe_radius);
  }
}


/*--------------------------------------------------------------------------------
find_vertex_order accepts Big_Points and determines the orientation of a sequence 
of triangles that have pt0 as a common point. This assumes that 
verts[0] == verts[num_verts-1]
----------------------------------------------------------------------------------*/
int
find_vertex_order(pt0, verts, num_verts)
VertexType	*pt0; 
Big_Point	*verts;
int		num_verts;
{ 

  POINT		normal, diff[2];
  POINT		ave_normal;
  double	tempd = 0;
  int		i, j;


  /* take the average of consecutive tris for robustness */
  for(i = 0; i < num_verts-1; i++)
  { j = i + 1;
    VEC3_V_OP_V(diff[0], verts[i].vt.Coord, -, pt0->Coord);
    VEC3_V_OP_V(diff[1], verts[j].vt.Coord, -, pt0->Coord);
    CROSSPROD3(normal, diff[0], diff[1]);
    NORMALIZE3(normal);
    VEC3_V_OP_V_OP_V(ave_normal,verts[i].vt.Normal,+,verts[j].vt.Normal,+,pt0->Normal);
    NORMALIZE3(ave_normal);
    tempd += DOTPROD3(normal, ave_normal);
  }
  return (tempd < 0); 
}

/*--------------------------------------------------------------------------------
find_order accepts POINTs and returns 0 if in going from verts[i] to verts[i+1],
right-hand-side atom is atom[0], otherwise it returns 1
----------------------------------------------------------------------------------*/
int
find_order(verts, num_verts, atom)
POINT			*verts;
int			num_verts;
int			atom[2];
{ 
  POINT			cg;
  POINT			vec1, vec2, vec3;
  POINT			atom_diff;
  int			i, max_dist_id;
  double		max_dist, tempd;

  /* find the average of the vertices */
  VEC3_ASN_OP(cg, =, verts[0]);
  for(i = 1; i < num_verts; i++)
    VEC3_ASN_OP(cg, +=, verts[i]);
  tempd = ONE/num_verts;
  VEC3_V_OP_S(cg, cg, *, tempd);

  /* find the consecutive verts that are furthest apart for
     robustness
  */
  for(max_dist = max_dist_id = i = 0; i < num_verts-1; i++)
  { tempd = SQ_DIST3(verts[i], verts[i+1]);
    if (tempd > max_dist)
    { max_dist = tempd;
      max_dist_id = i;
    }
  }
  
  VEC3_V_OP_V(vec1, verts[max_dist_id+1], -, verts[max_dist_id]); 
  VEC3_V_OP_V(vec2, cg, -, verts[max_dist_id]); 
  CROSSPROD3(vec3, vec1, vec2);
  
  VEC3_V_OP_V(atom_diff, atoms[atom[1]].center, -, atoms[atom[0]].center);
  
  return (DOTPROD3(vec3, atom_diff) < 0);
}
