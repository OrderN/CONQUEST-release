/* Copyright amitabh varshney@cs.unc.edu. All Rights Reserved. March 16, 1994 */
/*---------------------------------------------------------------------------------
compute.c

This file controls the processing of the surface. It has the highest-level
routines for surface computation.
----------------------------------------------------------------------------------*/
#include "surf.h"
#define  global
#include  "chull.h"			/* include convex hull structs */
#include  "dual.h"			/* include feasible region stuff */

/* define some vars for this file */

int              Max_cons;
int              Tot_cons;
Vector           *Constraints;
Atom_List        *Atom_Grid;
float            Voxel_Side;
int              Num_Voxels[3];
double           Max_Radius;
Vector           Bounding_Tetra[4];
double           Extents[3][2];
short            Neighbor_list[MAX_CONSTRAINT];

/*---------------------------------------------------------------------------------
init_and_compute initializes the data structures and starts off the computation.
----------------------------------------------------------------------------------*/
init_and_compute()
{ int			i;
  
    if (Max_Tess_Len <= 0) Max_Tess_Len = MAX_TESS_LEN;
    printf("Max edge length = %2.3lf\n", Max_Tess_Len);
    Max_Tess_Len_Sq = SQ(Max_Tess_Len);
    Num_polys = 0;

    ALLOCMAT(verts, Max_Gp_Polys, 3, VertexType);
    ALLOCN(atom_type, short, Max_Gp_Polys);

    /* allocate enough memory to accomodate Constraints and 6 extent planes */
    ALLOCN(Constraints, Vector, MAX_CONSTRAINT+6);
    
    /* compute the extents */
    compute_extents();

    /* compute the bounding tetrahedron */
    compute_bounding_tetra();

    Tot_cons = 0;
    Max_cons = 0;

    for(i = 0; i < Num_atoms; i++) 
    {  compute_components(i, Constraints);
       if (i%100==0) { printf("."); fflush(stdout); }
    }
    printf("\n");

    printf("Total Triangles %d Total constraints %d\n",Num_polys, Tot_cons);
    printf("Max Neighbors per atom %d Average Neighbors per atom %g\n", 
     	    Max_cons, Tot_cons*1.0/Num_atoms); 
}

/*---------------------------------------------------------------------------------
compute_extents computes the extents of the molecule, the number of voxels in the 
global grid and the maximum radius of an atom of the molecule. 
----------------------------------------------------------------------------------*/
compute_extents()
{ 
  int		i, j;
  int		total_voxels;
  int		atom_pos[3];
  int		temp_index;
  float		total_radius, temp_min, temp_max;

  for(j = 0; j < 3; j++)
  { Extents[j][MINC] = HUGE_VAL;
    Extents[j][MAXC] = -HUGE_VAL;
  }

  /* compute the extents and the min, max, & average radius*/
  total_radius = 0;
  Max_Radius = -HUGE_VAL;
  for (i = 0; i < Num_atoms; i++)
  { for( j = 0; j < 3; j++)
    { temp_min = atoms[i].center[j] - atoms[i].radius;
      temp_max = atoms[i].center[j] + atoms[i].radius;
      if (Extents[j][MINC] > temp_min) Extents[j][MINC] = temp_min;
      if (Extents[j][MAXC] < temp_max) Extents[j][MAXC] = temp_max;
    }
    total_radius += atoms[i].radius;
    if (Max_Radius < atoms[i].radius) Max_Radius = atoms[i].radius;
  }

  /* define the Voxel_Side to be the average diameter of an atom */
  Voxel_Side = (total_radius + total_radius)/Num_atoms;

  for(j = 0; j < 3; j++)
  { Num_Voxels[j] = (int)((Extents[j][MAXC] - Extents[j][MINC])
                           /Voxel_Side + FCEIL);
  }

  total_voxels = Num_Voxels[X]*Num_Voxels[Y]*Num_Voxels[Z];

  ALLOCN(Atom_Grid, Atom_List, total_voxels);

  for(i = 0; i < total_voxels; i++)
     Atom_Grid[i].num_neighbors = 0;

  /* initialize the grid */
  for(i = 0; i < Num_atoms; i++)
  { /* locate the grid voxel where this atom's center lies */
    for(j = 0; j < 3; j++)
      atom_pos[j] = (int)((atoms[i].center[j]-Extents[j][MINC])/Voxel_Side);
     
    temp_index = AINDEX(atom_pos,Num_Voxels);
    
    j = Atom_Grid[temp_index].num_neighbors;
    Atom_Grid[temp_index].neighbor[j] = i;
    if (++j > MAX_DENSITY)
    { printf("Atom density > MAX_DENSITY\n");
      exit(-1); 
    }
    else 
      Atom_Grid[temp_index].num_neighbors = j;
  }

  /* initialize the tessellation origins */
  for(i = 0; i < Num_atoms; i++)
  { atoms[i].tes_origin[0] = atoms[i].tes_origin[1] = atoms[i].tes_origin[2]
			   = -1e+22;
  }
}

/*---------------------------------------------------------------------------------
compute_bounding_tetra computes the four planes of a regular tetrahedron that 
bounds the molecule.
----------------------------------------------------------------------------------*/
compute_bounding_tetra()
{ 
  int		i;
  float		center[3];
  double	radius;

  /* compute the center of the extents */
  /* and find radius of circumscribing sphere */
  for(i = radius = 0; i < 3; i++)
  { center[i] = (Extents[i][MINC] + Extents[i][MAXC])*HALF;
    radius += SQ(Extents[i][MAXC] + Probe_radius + ONE - center[i]);
  }

  radius = 9*sqrt(radius);

  /* store the planes in the form [0]x + [1]y + [2]z <= [3] */
  /* plane 0 */
  Bounding_Tetra[0].coeff[X] =
  Bounding_Tetra[0].coeff[Y] = ZERO;
  Bounding_Tetra[0].coeff[Z] = -ONE;
  Bounding_Tetra[0].coeff[3] = radius/sqrt(6.0);
  
  /* plane 1 */
  Bounding_Tetra[1].coeff[X] = ZERO;
  Bounding_Tetra[1].coeff[Y] = -2*sqrt(2.0)/3.0;
  Bounding_Tetra[1].coeff[Z] = 1/3.0;
  Bounding_Tetra[1].coeff[3] = radius/sqrt(6.0);

  /* plane 2 */
  Bounding_Tetra[2].coeff[X] = sqrt(2.0/3.0);
  Bounding_Tetra[2].coeff[Y] = sqrt(2.0)/3.0;
  Bounding_Tetra[2].coeff[Z] = 1/3.0;
  Bounding_Tetra[2].coeff[3] = radius/sqrt(6.0);

  /* plane 3 */
  Bounding_Tetra[3].coeff[X] = -sqrt(2.0/3.0);
  Bounding_Tetra[3].coeff[Y] = sqrt(2.0)/3.0;
  Bounding_Tetra[3].coeff[Z] = 1/3.0;
  Bounding_Tetra[3].coeff[3] = radius/sqrt(6.0);

  for(i = 0; i < 4; i++)
  { /* translate the planes so that the new origin is center[] */
    Bounding_Tetra[i].coeff[3] += DOTPROD3(Bounding_Tetra[i].coeff, center);
  }

  if (Checks_On)
  { for(i = 0; i < Num_atoms; i++)
    { if((fabs(DOTPROD3(Bounding_Tetra[0].coeff,atoms[i].center)-
	       Bounding_Tetra[0].coeff[3]) < atoms[i].radius + Probe_radius)||
         (fabs(DOTPROD3(Bounding_Tetra[1].coeff,atoms[i].center)-
	       Bounding_Tetra[1].coeff[3]) < atoms[i].radius + Probe_radius)||
         (fabs(DOTPROD3(Bounding_Tetra[2].coeff,atoms[i].center)-
	       Bounding_Tetra[2].coeff[3]) < atoms[i].radius + Probe_radius)||
         (fabs(DOTPROD3(Bounding_Tetra[3].coeff,atoms[i].center)-
	       Bounding_Tetra[3].coeff[3]) < atoms[i].radius + Probe_radius))
       { printf("Bounding tetra intersects ext radius sph for atom %d\n",i);
         fflush(stdout);
       }
    }
  }
}


/*---------------------------------------------------------------------------------
compute_components is in some sense the central routine. This makes calls to others
for computing the surface.
----------------------------------------------------------------------------------*/
compute_components(atom_id, constraints)
int		atom_id;
Vector		*constraints;
{ 
  int		num_constraints;

  Current_atom = atom_id;

  if (!compute_neighbors(atom_id, &num_constraints)) return;

  compute_planes(atom_id, num_constraints,constraints);

  Tot_cons += num_constraints;
  if (Max_cons < num_constraints) Max_cons = num_constraints;

  if (find_origin(atom_id, constraints, num_constraints))
  { 
    transform_and_add_extents(atom_id, constraints, &num_constraints);

    atoms[atom_id].num_cons = num_constraints;

    /* transform the planes to points */
    dualize_planes_to_points(constraints, num_constraints);

    /* initialize the convex hull tetrahedron */
    init_tet();

    /* form the convex hull */
    complete_hull();

    /* for the dual feasible region */
    dualize_hull(atom_id, constraints);

    /* verify the Euler's formula for planar graphs */
    if (rvertices_count - redges_count + rfaces_count != 2)
    { printf("Some convex hull bug for atom_id %d case\n", atom_id);
      fflush(stdout);
    }
  
    /* find the origin for tesselation */
    find_tes_origin(atom_id, constraints); 
      
    /* find components */
    find_components(atom_id, constraints); 
  
    /* free up all the memory so that it can be used subsequently */
    clean_fc();
    meta_clean();
  }
}

/*---------------------------------------------------------------------------------
compute_neighbors identifies the atoms that are neighbors (in the molecular surface
sense) to the given atom i. If a probe sphere can be placed touching atoms i and j
simultaneously, (without considering steric hinderance from other atoms) then atoms
i and j are neighbors in the molecular surface sense.
----------------------------------------------------------------------------------*/
int
compute_neighbors(i, num_cons)
int		i;
int		*num_cons;
{
  int		j, k, l, m, n, p, return_val;
  int		voxel_extents[3][2], temp_index[3];
  double	tempf, partial_radius, effective_radius;
  double	dist, twice_probe_radius;

  *num_cons = 0;
  /* find the radius within which we have to look for the centers of
     all "neighboring" atoms.
  */
  twice_probe_radius = 2*Probe_radius;
  partial_radius = atoms[i].radius + twice_probe_radius;
  effective_radius = partial_radius + Max_Radius;

  /* find the range of voxels to be searched for along each of the
     three dimensions
  */
  for(j = 0; j < 3; j++)
  { tempf = atoms[i].center[j] - Extents[j][MINC];
    voxel_extents[j][MINC] = (int)((tempf - effective_radius)/ Voxel_Side);
    /* Since the MAXC voxel should be one more than the ceil of the highest
       voxel covered, we add 1+FCEIL here
    */
    voxel_extents[j][MAXC] = (int)((tempf + effective_radius)/Voxel_Side +
                                FCEIL_PLUS_1);
    /* clamp to make the extents fall in [0, Num_Voxels[j]) */
    if (voxel_extents[j][MINC] < 0) voxel_extents[j][MINC] = 0;
    if (voxel_extents[j][MAXC] > Num_Voxels[j]) 
       voxel_extents[j][MAXC] = Num_Voxels[j];
  }
  
  p = 0; return_val = FALSE; 
  temp_index[0] = voxel_extents[X][MINC]*Num_Voxels[1];
  for(j = voxel_extents[X][MINC]; j < voxel_extents[X][MAXC]; j++)
  { temp_index[1] = (temp_index[0] + voxel_extents[Y][MINC])*Num_Voxels[2];
    for(k = voxel_extents[Y][MINC]; k < voxel_extents[Y][MAXC]; k++)
    { temp_index[2] = temp_index[1] + voxel_extents[Z][MINC];
      for(m = voxel_extents[Z][MINC]; m < voxel_extents[Z][MAXC]; m++)
      { 
        /*if (temp_index[2] != INDEX(j, k, m, Num_Voxels)) printf("gochi \n"); */
        for(l = 0; l < Atom_Grid[temp_index[2]].num_neighbors; l++)
        { n = Atom_Grid[temp_index[2]].neighbor[l];
          /* compute the distance */
          effective_radius = partial_radius + atoms[n].radius;

	  dist = SQ_DIST3(atoms[i].center, atoms[n].center);

          if ((dist <=  SQ(effective_radius) + GP_EPS) && (i != n))
          { /* this atom is close enough to interact */

	      *num_cons += 1;
            /* put n in the neighbor list of atom i, if (i != n) */
            Neighbor_list[p] = n;
            if (++p > MAX_CONSTRAINT)
            { printf("Number of neighbors of atom > MAX_CONSTRAINT\n");
              exit(-1); 
            }
          }
        }
	temp_index[2]++;
      }
      temp_index[1] += Num_Voxels[2];
    }
    temp_index[0] += Num_Voxels[1];
  }
 
  return_val = TRUE;

  return(return_val);
}

/*---------------------------------------------------------------------------------
compute_planes computes all the plane constraints for a given atom due to its 
neighbors. These plane constraints are also known as chordales in some circles 
(no pun intended :-)).
----------------------------------------------------------------------------------*/
compute_planes(atom_id, num_cons, cons)
int		atom_id;
int		num_cons;
Vector		*cons;
{ 
  int		i, k, n;
  int		cur_cons;
  int		num_neighbors;
  static int	tot_k = 0;
  
  num_neighbors = num_cons;
  cur_cons = 0;

  for(k = 0; k < num_neighbors; k++)
  { n = Neighbor_list[k];
    /* Find the equation of the bounding hyperplane between these two 
       atoms. This is of the following form:
       [0]x + [1]y + [2]z <= [3]
    */
    i = cur_cons++;
    VEC3_V_OP_V(cons[i].coeff,atoms[n].center,-,atoms[atom_id].center);
      
    cons[i].coeff[3] = 
     	DOTPROD3(atoms[n].center, atoms[n].center) -
      	DOTPROD3(atoms[atom_id].center, atoms[atom_id].center) -
      	SQ(Probe_radius + atoms[n].radius) + 
      	SQ(Probe_radius + atoms[atom_id].radius);
    cons[i].coeff[3] *= HALF;
    NORMALIZE_PLANE3(cons[i].coeff);
    cons[i].atom_id = n;
  }
  tot_k += num_cons;

}

/*---------------------------------------------------------------------------------
find_origin finds a point that lies inside the feasible region as well as the 
extended-radius sphere for atom atom_id. This is considered as the origin for
tessellation (or tes_origin). If the feasible cell is null, then this is not
possible and the routine returns a FALSE in that case -- it basically means that
this atom does not contribute anything to the surface and is a buried atom.
----------------------------------------------------------------------------------*/
int
find_origin(atom_id, cons, num_cons)
int             atom_id;
Vector          *cons;
int             num_cons;
{ int           i, j;
  float         tempf;
  double        lambda;
  float         extents[3][2];
  Vector        obj, sol;
  int           sol_cons[3];
  int           origin_feasible = TRUE;
  POINT         ave_dir, tempv;

  /* initialize the origin to be the center of the atom */
  VEC3_ASN_OP(New_Origin.coeff, =, atoms[atom_id].center);

  /* check to see if the origin is in the feasible region */
  for(i = 0; i < num_cons; i++)
  { lambda=DOTPROD3(cons[i].coeff, New_Origin.coeff);
    if (lambda > cons[i].coeff[3])
    { origin_feasible = FALSE;
      break;
    }
  }
  
  if (origin_feasible)
  { atoms[atom_id].boundary = TRUE;
    VEC3_ASN_OP(atoms[atom_id].tes_origin, =, New_Origin.coeff);
  }
  else /* check to see if the feasible region is non-null */
  { 
    tempf = atoms[atom_id].radius + Probe_radius + GP_EPS;
    for(i = 0; i < 3; i++)
    { extents[i][MINC] = atoms[atom_id].center[i] - tempf;
      extents[i][MAXC] = atoms[atom_id].center[i] + tempf;
    }

    /* initialize the objective function to x+y+z */
    obj.coeff[0] = obj.coeff[1] = obj.coeff[2] = ONE;
    sol_cons[0] = sol_cons[1] = sol_cons[2] = -1;
    atoms[atom_id].boundary =
                linear_prog(3,obj, extents, cons, num_cons, &sol, sol_cons);
    if (atoms[atom_id].boundary)
    {
      /* add the 6 extent planes also to the list of constraints */
      for(i = 0; i < 3; i++)
      { VEC3_ZERO(cons[num_cons].coeff);
        cons[num_cons].coeff[i] = 1;
        cons[num_cons].coeff[3] = extents[i][MAXC];
        cons[num_cons].atom_id = -1;
        num_cons += 1;
        VEC3_ZERO(cons[num_cons].coeff);
        cons[num_cons].coeff[i] = -1;
        cons[num_cons].coeff[3] = -extents[i][MINC];
        cons[num_cons].atom_id = -1;
        num_cons += 1;
      }

      /* set the sol_cons values correctly here as they are not
         being set correctly right now in lp code. should eventually
	 fix the lp code to correctly set the sol_cons values */
      find_sol_cons(cons, num_cons, sol_cons, sol);

      /* find the ray direction  as the average of the cross-prods of normals 
	 of the planes forming the solution vertex */
      VEC3_ZERO(ave_dir);

      for(i = 0; i < 3; i++)
      { j = (i==2)? 0: (i+1);
        CROSSPROD3(tempv, cons[sol_cons[i]].coeff, cons[sol_cons[j]].coeff);
        if (DOTPROD3(tempv, obj.coeff) > 0) VEC3_NEG(tempv, tempv);
	NORMALIZE3(tempv);
        VEC3_ASN_OP(ave_dir, +=, tempv);
      }
      NORMALIZE3(ave_dir);

      /* find the first intersection of this ray in the feasible region */
      lambda = HUGE_VAL;
      for(i = 0; i < num_cons; i++)
      { tempf = DOTPROD3(cons[i].coeff, ave_dir);
        if ((tempf > GP_EPS) || (tempf < -GP_EPS))
          tempf = (cons[i].coeff[3] - DOTPROD3(cons[i].coeff,sol.coeff))/tempf;
        if (tempf > GP_EPS)
         lambda = FMIN(tempf, lambda);
      }
      /* following will yield a point in the interior of feasible region */
      lambda *= HALF;

      VEC3_V_OP_S(New_Origin.coeff, ave_dir, *, lambda);
      VEC3_ASN_OP(New_Origin.coeff, +=, sol.coeff);

      VEC3_ASN_OP(atoms[atom_id].tes_origin, =, New_Origin.coeff);
  
      if (Checks_On)
      { /* check to see if the sol.coeff is in the feasible region */
        for(i = 0; i < num_cons; i++)
        { lambda=DOTPROD3(cons[i].coeff, sol.coeff) - cons[i].coeff[3];
          if (lambda > GP_EPS)
          { printf("Precision warning(): Solution barely satisfied in atom_id %d, constraint %d lambda %lf\n", atom_id, i, lambda);
          }
        }
        /* check to see if the new origin is in the feasible region */
        for(i = 0; i < num_cons; i++)
        { lambda=DOTPROD3(cons[i].coeff, New_Origin.coeff);
          if (lambda > cons[i].coeff[3])
          { origin_feasible = FALSE;
            printf("Origin outside region in atom %d, constraint %d by %lf\n",
            atom_id, i, lambda - cons[i].coeff[3]);
          }
        }
      }
    }
  }
  return(atoms[atom_id].boundary);
}

/*---------------------------------------------------------------------------------
find_sol_cons finds the three planes/constraints whose intersection is the point
at which the linear program returns the optimal point.
----------------------------------------------------------------------------------*/
find_sol_cons(cons, num_cons, sol_cons, sol)
Vector		*cons;
int		num_cons;
int		sol_cons[3];
Vector		sol;
{
  int		i, j, k, m;
  double	lambda, dtemp;
  double	min_lambda[3]; /* [0] has smallest value and [2] the largest */
  int		itemp;

  min_lambda[0] = min_lambda[1] = min_lambda[2] = HUGE_VAL;
  for(i = j = 0; ((i < num_cons)&&(j < 3)); i++)
  { lambda = fabs(DOTPROD3(cons[i].coeff, sol.coeff) - cons[i].coeff[3]);
    if (lambda < min_lambda[2])
    { min_lambda[2] = lambda;
      sol_cons[2] = i;
      /* re-order min_lambda and sol_cons */
      for(k = 0; k < 3; k++)
      { m = (k & 0x01);
	if (min_lambda[m] > min_lambda[m+1])
        { dtemp = min_lambda[m]; min_lambda[m] = min_lambda[m+1]; min_lambda[m+1] = dtemp;
          itemp = sol_cons[m]; sol_cons[m] = sol_cons[m+1]; sol_cons[m+1] = itemp;
        }
      }
    }
  }
  if (min_lambda[2] > GP_EPS)
  { printf("find_sol_cons() warning: sol is somewhat away from vertex (%lf).\n", min_lambda[2]);
    fflush(stdout); 
  }
}     

/*---------------------------------------------------------------------------------
transform_and_add_extents transforms all the constraints of a given atom such
that the origin of the new coordinate system is the tessellation origin of the 
atom. It also adds the four tetrahedral constraints to ensure that the feasible 
cell is closed.
----------------------------------------------------------------------------------*/
transform_and_add_extents(atom_id, cons, num_cons)
int		atom_id;
Vector		*cons;
int		*num_cons;
{ 
  int		i;
  Vector	tempv;

  /* initialize the New_Origin from the atoms[].tes_origin */
  VEC3_ASN_OP(New_Origin.coeff, =, atoms[atom_id].tes_origin);

  /* transform all constraints to the new coordinate system */
  for(i = 0; i < *num_cons; i++)
    cons[i].coeff[3] -= DOTPROD3(cons[i].coeff, New_Origin.coeff);

  /* sort all these planes in decreasing distance from origin */
  cons[*num_cons].coeff[3] = -HUGE_VAL;
  quick_sort_planes(cons, 0, *num_cons-1);
  
  /* add the 4 planes of a tetrahedron bounding the above extents */
  for(i = 0; i < 4; i++)
  { VEC_ASN_OP(cons[*num_cons].coeff, =, Bounding_Tetra[i].coeff, 4);
    cons[*num_cons].coeff[3] -= DOTPROD3(cons[*num_cons].coeff, 
					 New_Origin.coeff);
    cons[*num_cons].atom_id = -1;
    *num_cons += 1;
  }
  /* swap the first and last constraints (to prevent degeneracies in
     convex-hull code) and thereby make the first four verts in convex 
     hull code to be the verts of this tetrahedron 
  */
  memcpy(&tempv, &cons[0], sizeof(Vector));
  memcpy(&cons[0], &cons[*num_cons-4], sizeof(Vector));
  memcpy(&cons[*num_cons-4], &tempv, sizeof(Vector));
}

/*---------------------------------------------------------------------------------
find_components identifies the arcs formed by the intersection of the feasible cell
with the extended-radius sphere of atom atom_id and accordingly makes appropriate 
calls to generates the patches.
----------------------------------------------------------------------------------*/
find_components(atom_id, cons)
int		atom_id;
Vector		*cons;
{ 
  double			sq_radius;
  int				i, j;
  int				edge_int;
  int				all_in;

  sq_radius = SQ(Probe_radius + atoms[atom_id].radius);
  /* check to make sure at least one vertex is out of ext-radius sphere,
     else trivial reject 
  */
  if (all_verts_in_sphere(rvertices, rvertices_count, sq_radius)) 
  { /* printf("atom_id %d all_in = TRUE \n", atom_id); */
    return;
  }

  for(i = 0; i < rfaces_count; i++)	       /* for each face    */
  { 
    for(j = 0, all_in = TRUE; j < rfaces[i].num_verts; j++)
    { if (rfaces[i].vert[j]->sq_dist > sq_radius) 
      { all_in = FALSE; 
        break;
      }
    }
    if (all_in) continue;


    edge_int = gen_component(atom_id, i, cons);

    if (edge_int == FALSE)	
      gen_case_III_IV(atom_id, i, cons);

  }
}

/*---------------------------------------------------------------------------------
find_tes_origin initializes the tessellation origin appropriately, i.e. it 
should be inside the feasible region as well as inside the extended-radius
sphere of atom atom_id. 
This assumes that int_pts, end_pts, etc, are in local coord system 
(the origin of the system is the center of the atom atom_id).
This returns tes_origin in global coord system.
----------------------------------------------------------------------------------*/
find_tes_origin(atom_id, cons)
int			atom_id;
Vector			*cons;
{ 
  int			i;
  POINT			temp_origin;
  POINT			diff, temp_center;
  double		fac, sq_radius;
  int			found_origin = FALSE;
  int			face_atoms[2];
  int			count;
  double		tempd;

  sq_radius = SQ(Probe_radius + atoms[atom_id].radius);
  VEC3_ASN_OP(atoms[atom_id].tes_origin, =, atoms[atom_id].center);

  VEC3_V_OP_V(diff, New_Origin.coeff, -, atoms[atom_id].center);
  if (DOTPROD3(diff, diff) < 1e-10)
  /* atom center lies in the feasible region, and so can be also
     taken to be the tes_origin - no more computations required */
    return;

  /* try computing the tes_origin as the average of all int_pts */
  VEC3_ZERO(temp_origin); count = 0;

  for(i = 0; i < redges_count; i++)
  { if (redges[i].int_pts[0][0] != HUGE_VAL)
    { VEC3_ASN_OP(temp_origin, +=, redges[i].int_pts[0]);
      count++;
      if (redges[i].int_pts[1][0] != HUGE_VAL)
      { VEC3_ASN_OP(temp_origin, +=, redges[i].int_pts[1]);
	count++;
      }
    }
  }
  if (count > 0)
  { found_origin = TRUE;
    VEC3_V_OP_S(temp_origin, temp_origin, /, count);
    /* convert from local to global coord system */
    VEC3_ASN_OP(temp_origin, +=, atoms[atom_id].center);
  }

  if (!found_origin) /* need to look at the faces */
  { for(i = 0; i < rfaces_count; i++)
    { if (cons[rfaces[i].id].atom_id >= 0)
      { face_atoms[0] = atom_id; face_atoms[1] = cons[rfaces[i].id].atom_id; 
	find_circle(face_atoms, cons[rfaces[i].id], temp_center, &fac);
	if ((fac > GP_EPS) && (point_in_region(temp_center, cons)))
	{ VEC3_ASN_OP(temp_origin, +=, temp_center);
	  count++;
        }
      }
    }
    if (count == 1)
    { found_origin = TRUE;
      tempd = sq_radius - fac*fac;      
      tempd = ((Probe_radius + atoms[atom_id].radius) + sqrt(tempd))*HALF;
      VEC3_V_OP_V(diff, temp_origin, -, atoms[atom_id].center);
      NORMALIZE3(diff);
      VEC3_V_OP_V_OP_S(temp_origin, atoms[atom_id].center, +, 
		       diff, *, tempd);
    }
    else if (count > 1)
    { found_origin = TRUE;
      VEC3_V_OP_S(temp_origin, temp_origin, /, count);
    }
  }

  if (found_origin)
  { VEC3_V_OP_V(diff, temp_origin, -, atoms[atom_id].center); 
    if (DOTPROD3(diff, diff) >= sq_radius)
    { fprintf(stdout,"find_tes_origin(): Bad tes_origin:atom%d:(%lf %lf %lf)\n",
	atom_id, temp_origin[0], temp_origin[1], temp_origin[2]);
      fflush(stdout);
      exit(-1);
    }
    VEC3_ASN_OP(atoms[atom_id].tes_origin, =, temp_origin);
  }
}

/*---------------------------------------------------------------------------------
point_in_region  checks to see if a point is inside the feasible cell as determined
by the cons array.
This assumes:
  (a) The elements of cons are of the form [0]x + [1]y + [2]z <= [3]
  (b) The elements of cons are in the New_Origin coord system 
  (c) The test_point is in the global coord system
----------------------------------------------------------------------------------*/
int
point_in_region(test_point, cons)
POINT		test_point;
Vector		*cons;
{ 
  int		i, j;
  int		result = TRUE;
  double	tempd;
  POINT		tempv;

  /* transform the point to the New Origin coord system */
  VEC3_V_OP_V(tempv, test_point, -, New_Origin.coeff);

  /* check for validity */
  for(i = 0; i < rfaces_count; i++)
  { j = rfaces[i].id;
    tempd = DOTPROD3(cons[j].coeff, tempv) - cons[j].coeff[3];
    if ((tempd <= 0) || (FP_EQ_EPS(tempd, 0, GP_EPS))) continue;
    else 
    { result = FALSE;
      break;
    }
  }
  return(result);
}

/*---------------------------------------------------------------------------------
all_verts_in_sphere checks to see if all the vertices are inside the sq_radius 
(returns 1) or not (returns 0)
----------------------------------------------------------------------------------*/
int
all_verts_in_sphere(rverts, num_rverts, sq_radius)
struct rvertex	*rverts;
int		num_rverts;
double		sq_radius;
{
  int		i;
  int		all_in = TRUE;

  for(i = 0; i < num_rverts; i++)
  { if (rverts[i].sq_dist > sq_radius) 
    { all_in = FALSE; 
      break;
    }
  }
  return(all_in);
}


/*---------------------------------------------------------------------------------
gen_component calls the gen_case_I_II for the possibly several arcs contributed
by a single feasible cell face.
This assumes that all the coordinates of the vertices are in the local coordinate 
system of the atom [atom_id]
----------------------------------------------------------------------------------*/
int
gen_component(atom_id, face_id, cons)
int		atom_id, face_id;
Vector		*cons;
{ 
  int		j, k;
  struct Torus*	torus_list;
  struct Torus*	temp_torus;
  double	sq_radius;
  int		torus_found = FALSE;

  sq_radius = SQ(Probe_radius + atoms[atom_id].radius);
  torus_list = NULL;

  j = 0;
  while (j < rfaces[face_id].num_verts)		/* consider all edges */
  { k = (j == rfaces[face_id].num_verts-1)? 0: j+1;
    if ((rfaces[face_id].edge[j]->int_pts[0][0] != HUGE_VAL) &&
	(rfaces[face_id].vert[k]->sq_dist > sq_radius))
    { ALLOCATE(temp_torus, struct Torus); /* edge intersects sphere */
      get_torus_elt(atom_id, face_id, &j, cons, temp_torus);  
      ADD_QUEUE(torus_list, temp_torus);
    }
    else
      j++;
  }

  if (torus_list)
  { /* generate all the torii one by one */
    temp_torus = torus_list;
    do
    {
      gen_case_I_II(temp_torus, cons); 
      temp_torus = temp_torus->next;
    }
    while(temp_torus != torus_list);
    torus_found = TRUE;
  }

  while(torus_list) 
  { temp_torus = torus_list;
    DEL_QUEUE( torus_list, temp_torus);
  }

  return(torus_found);
}

/*---------------------------------------------------------------------------------
get_torus_elt adds a torus element to the list of torii contributed by a single
feasible face.
----------------------------------------------------------------------------------*/
/* adds a torus element to the list of torii */
get_torus_elt(atom_id, face_id, edge_id, cons, torus_elt)
int             atom_id, face_id;
int		*edge_id;
Vector          *cons;
struct Torus	*torus_elt;
{
  int		i, n, cur_edge;
  int		num_int_verts = 1; /* number of interior verts in component */
  struct redge  *edge[2];

  torus_elt->face_atoms[0] = atom_id;
  torus_elt->face_atoms[1] = cons[rfaces[face_id].id].atom_id;
  torus_elt->face_id = face_id;

  n = rfaces[face_id].num_verts - 1;
  cur_edge = (*edge_id == n)? 0: *edge_id+1;

  if (rfaces[face_id].edge[cur_edge]->int_pts[0][0] == HUGE_VAL)
  { do
    { cur_edge = (cur_edge == n)? 0: cur_edge+1;
      num_int_verts++;
    }
    while (rfaces[face_id].edge[cur_edge]->int_pts[0][0] == HUGE_VAL);
  }
  torus_elt->edge_id[0] = *edge_id;
  torus_elt->edge_id[1] = cur_edge;
  torus_elt->num_int_verts = num_int_verts;
  for(i = 0; i < 2; i++)
  { edge[i] = rfaces[face_id].edge[torus_elt->edge_id[i]];
    torus_elt->end_atoms[i] = (edge[i]->adj_faces[0] == &rfaces[face_id])? 
    			      cons[edge[i]->adj_faces[1]->id].atom_id: 
			      cons[edge[i]->adj_faces[0]->id].atom_id;
  }

  if (cur_edge <= *edge_id) 	/* we have traversed the (n-1, 0) boundary */
    *edge_id = n+1;
  else *edge_id = cur_edge;

}

/*---------------------------------------------------------------------------------
quick_sort_planes uses QuickSort to sort the planes in descending order of their
distance from the origin.
----------------------------------------------------------------------------------*/
/* quick sort the array in descending order */
quick_sort_planes(cons, m, n)
Vector                  *cons;
int                     m, n;
{
  int                   i, j;
  float                 key;
  Vector		temp_vec;
  int			vec_size = sizeof(Vector);

  if (m < n)
  { i = m;
    j = n + 1;
  
    key = cons[m].coeff[3];

    do
    {  while (cons[++i].coeff[3] > key);
       while (cons[--j].coeff[3] < key);
       if (i < j) 
       { /* swap i and j-th elements */
	 memcpy(&temp_vec, &cons[i], vec_size);
	 memcpy(&cons[i], &cons[j], vec_size);
	 memcpy(&cons[j], &temp_vec, vec_size);
       }
    }
    while (i < j);

    /* swap m and j-th elements */
    memcpy(&temp_vec, &cons[m], vec_size);
    memcpy(&cons[m], &cons[j], vec_size);
    memcpy(&cons[j], &temp_vec, vec_size);
    quick_sort_planes(cons, m, j-1);
    quick_sort_planes(cons, j+1, n);
  }
}

