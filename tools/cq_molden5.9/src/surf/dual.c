/* Copyright amitabh varshney@cs.unc.edu. All Rights Reserved. March 16, 1994 */
/*--------------------------------------------------------------------------------
dual.c

This file contains the procedures to compute the dual of the convex hull computed
in chull.c. This yield the intersection of n-halfspaces which we have referred
to in our papers as a feasible cell. The concept of duality and this transformation
are explained well in "Computational Geometry - an Introduction" by F. P. 
Preparata and M. I. Shamos, Springer-Verlag, 1985.

----------------------------------------------------------------------------------*/
#include "surf.h"
#include "dual.h"
#include "chull.h"

/*----------------------------------------------------------------------------------
dualize_hull computes the dual of the convex hull for atom atom_id.
----------------------------------------------------------------------------------*/
dualize_hull(atom_id, cons)
int		atom_id;
Vector		*cons;
{ 
  struct tvertex *temp_v;
  struct tface	 *temp_f;
  int		 i;

  /* number the convex hull faces, vertices, edges and allocate the 
     feasible region/cell vertices, faces, edges.
  */
  alloc_fc_verts();
  alloc_fc_faces();
  alloc_fc_edges();

  /* dualize all the convex hull faces to feasible cell vertices 
     (in rvertices array)
  */
  temp_f = faces;
  do
  { dualize_face_to_point(atom_id, temp_f, cons);
    temp_f = temp_f->next;
  } while(temp_f != faces);

  /* dualize all the convex hull vertices to feasible cell faces 
     (in rfaces array) 
  */
  temp_f = faces;
  do
  { for(i = 0; i < 3; i++)
    { temp_v = temp_f->vert[i];
      if (temp_v->mark == !MARKED)/* generate the feasible cell face for this 
				     vertex*/
      { temp_v->mark = MARKED;
        gen_fc_face(temp_v, temp_f);          
      }
    }
    temp_f = temp_f->next;
  } while(temp_f != faces);


  /* set up ptrs to redges from rvertices */
  setup_vert_edge_ptrs();

  /* compute edge-atom intersections */
  for(i = 0; i < redges_count; i++)
    find_edge_sphere_intersection(atom_id, &redges[i]);
}

/*----------------------------------------------------------------------------------
alloc_fc_verts numbers the convex hull faces and allocates the feasible
cell vertices 
----------------------------------------------------------------------------------*/
alloc_fc_verts()
{
  struct tface	 *temp_f;

  temp_f = faces; rvertices_count = 0;
  do
  { temp_f->fnum = rvertices_count++; 	
    temp_f = temp_f->next;
  } while(temp_f != faces);
  ALLOCN(rvertices, struct rvertex, rvertices_count);
}

/* number the ch vertices and allocate the feasible cell faces */
/* this assumes that alloc_fc_verts() has been called before */
alloc_fc_faces()
{
  struct tvertex *temp_v;
  struct tedge	 *temp_e;
  int		  i;

  temp_v = vertices; rfaces_count = 0;
  do
  { temp_v->vnum[0] = rfaces_count++; 	
    temp_v->mark = !MARKED;
    temp_v = temp_v->next;
  } while(temp_v != vertices);
  ALLOCN(rfaces, struct rface, rfaces_count);

  /* store the number of ch edges incident on a ch vertex in the
     corresponding feasible cell face (in rfaces array) */
  for(i = 0; i < rfaces_count; i++) rfaces[i].num_verts = 0;
  temp_e = edges;
  do
  { rfaces[temp_e->endpts[0]->vnum[0]].num_verts += 1;
    rfaces[temp_e->endpts[1]->vnum[0]].num_verts += 1;
    temp_e = temp_e->next;
  } while(temp_e != edges);
  for(i = 0; i < rfaces_count; i++) 
    if(rfaces[i].num_verts == 0) 
      printf("Vert %d has zero edges in convex hull \n", i);

  for(i = 0; i < rfaces_count; i++) 
  { ALLOCN(rfaces[i].vert, struct rvertex*, rfaces[i].num_verts);
    ALLOCN(rfaces[i].edge, struct redge*,   rfaces[i].num_verts);
  }
}

/*----------------------------------------------------------------------------------
alloc_fc_edges numbers the convex hull edges and allocate the feasible
cell edges. this assumes that alloc_fc_verts() and alloc_fc_faces() have been
called before 
----------------------------------------------------------------------------------*/
alloc_fc_edges()
{
  struct tedge	*temp_e;
  int		i;

  temp_e = edges; redges_count = 0;
  do
  { temp_e->ednum = redges_count++; 	
    temp_e = temp_e->next;
  } while(temp_e != edges);
  ALLOCN(redges, struct redge, redges_count);
  
  temp_e = edges; 
  do
  { i = temp_e->ednum;
    redges[i].adj_faces[0] = &rfaces[temp_e->endpts[0]->vnum[0]];
    redges[i].adj_faces[1] = &rfaces[temp_e->endpts[1]->vnum[0]];

    redges[i].end_pts[0] = &rvertices[temp_e->adjface[0]->fnum];
    redges[i].end_pts[1] = &rvertices[temp_e->adjface[1]->fnum];

    temp_e = temp_e->next;
  } while(temp_e != edges);

}

/*----------------------------------------------------------------------------------
setup_vert_edge_ptrs set up ptrs to redges from rvertices 
----------------------------------------------------------------------------------*/
setup_vert_edge_ptrs()
{
  struct tface  	*temp_f;
  int			i;
   
  temp_f = faces;
  do
  { 
    for(i = 0; i < 3; i++)
      rvertices[temp_f->fnum].in_edges[i] = &redges[temp_f->edg[i]->ednum];

    temp_f = temp_f->next;
  } while (temp_f != faces);
}

/*----------------------------------------------------------------------------------
clean_fc frees up space allocated in a feasible cell
----------------------------------------------------------------------------------*/
clean_fc()
{ 
  int		i;

  for(i = 0; i < rfaces_count; i++)
  { free(rfaces[i].vert);
    free(rfaces[i].edge);
  }

  free(rfaces);
  free(redges);
  free(rvertices);
}

/*----------------------------------------------------------------------------------
gen_fc_face generate a feasible cell face that corresponds to the convex
hull vertex ch_vert
----------------------------------------------------------------------------------*/
gen_fc_face(ch_vert, ch_face)
struct	tvertex		*ch_vert;
struct	tface		*ch_face;
{
  int			fc_face_id = ch_vert->vnum[0];
  int			i, j, k, l;
  struct tface		*temp_f;
  struct tedge		*temp_e = NULL;

  temp_f = ch_face;
  rfaces[fc_face_id].id = ch_vert->data_ptr;
  
  if ((temp_f->edg[0]->endpts[0] == ch_vert) ||
      (temp_f->edg[0]->endpts[1] == ch_vert))
    temp_e = temp_f->edg[0];
  else 
  { if ((temp_f->edg[1]->endpts[0] == ch_vert) ||
        (temp_f->edg[1]->endpts[1] == ch_vert))
        temp_e = temp_f->edg[1];
    else
    { printf("gen_fc_face(): ch_vert != either edge endpt\n");
      exit(-1);
    }
  }


  l = 0;

  do
  { rfaces[fc_face_id].vert[l] = &rvertices[temp_f->fnum];
    rfaces[fc_face_id].edge[l++] = &redges[temp_e->ednum];

    if (temp_e->adjface[0]==temp_f)
      temp_f = temp_e->adjface[1];
    else
    { if (temp_e->adjface[1]==temp_f)
        temp_f = temp_e->adjface[0];
      else
      { printf("gen_fc_face(): temp_f != either edge adjface\n");
        exit(-1);
      }
    }

    i = (temp_f->edg[0] == temp_e)? 0: ((temp_f->edg[1] == temp_e)? 1: 2);
    j = (i+1)%3; k = (i+2)%3; 

    if ((temp_f->edg[j]->endpts[0] == ch_vert) ||
        (temp_f->edg[j]->endpts[1] == ch_vert))
      temp_e = temp_f->edg[j];
    else 
    { if ((temp_f->edg[k]->endpts[0] == ch_vert) ||
          (temp_f->edg[k]->endpts[1] == ch_vert))
          temp_e = temp_f->edg[k];
      else
      { printf("gen_fc_face(): ch_vert != either edge endpt\n");
        exit(-1);
      }
    }

  } while (temp_f != ch_face);
}

/*----------------------------------------------------------------------------------
dualize_face_to_point computes a feasible cell vert from the plane passing through
the three vertices of the convex hull face (face_ptr) in the dual space
----------------------------------------------------------------------------------*/
dualize_face_to_point(atom_id, face_ptr, cons)
int			atom_id;
struct tface		*face_ptr;
Vector			*cons;
{
  int                  		i, j, k;
  double                        result[3];
  int				id = face_ptr->fnum;

  i = face_ptr->vert[0]->data_ptr;
  j = face_ptr->vert[1]->data_ptr;
  k = face_ptr->vert[2]->data_ptr;
  compute_tri_plane_int(cons[i], cons[j], cons[k], result);

  /* convert result to world-space coordinates */
  VEC3_ASN_OP(result, +=, New_Origin.coeff);

  /* store the distance of this point from the center of the atom */
  VEC3_V_OP_V(rvertices[id].v, result, -, atoms[atom_id].center);
  rvertices[id].sq_dist = DOTPROD3(rvertices[id].v, rvertices[id].v);

  for(i = 0; i < 3; i++)
    rvertices[id].adj_faces[i] = &rfaces[face_ptr->vert[i]->vnum[0]];
}


/*----------------------------------------------------------------------------------
dualize_planes_to_points -- dualizes planes to points :-)
----------------------------------------------------------------------------------*/
dualize_planes_to_points(plane, num_planes)
Vector          *plane;
int             num_planes;
{
  register int                  i, j;
  register struct tvertex       *new;
  register int			perturbation;


  for(i = 0; i < num_planes; i++)
  { ALLOCATE(new, struct tvertex);

    /* make sure the plane does not pass through the origin */
    if ((plane[i].coeff[3] < GP_EPS) && (plane[i].coeff[3] >= 0))
      plane[i].coeff[3] = GP_EPS;
    else if ((plane[i].coeff[3] > -GP_EPS) && (plane[i].coeff[3] < 0))
      plane[i].coeff[3] = -GP_EPS;
 
    VEC3_V_OP_S(new->v, plane[i].coeff, /, plane[i].coeff[3]);

    /* add Emiris and Canny's perturbation here */
    for(perturbation = j = 1; j < 4; j++)
    { perturbation *= i;
      new->v[j] += EPSILON*(perturbation % PERT_Q);
    }

    new->duplicate = NULL;
    new->active = !ACTIVE;
    new->mark = !MARKED;
    new->data_ptr = i;
    ADD_QUEUE(vertices, new);

  }
}

/*----------------------------------------------------------------------------------
find_edge_sphere_intersection 
   computes the intersections of an edge with the extended radius 
   sphere around an atom. answer is stored in local coord system of 
   the atom (center of atom is (0,0,0)). 
   . in case of no intersection, the int_pts[0..1] contain HUGE_VAL
   . in case of one intersection, the pt is stored in int_pts[0]
     and int_pts[1] contains HUGE_VAL
   . in case of two intersection pts, the pt that is closer to edge
     endpt[0] is stored in int_pts[0], and the other one (closer to
     endpt[1]) is stored in int_pts[1].
----------------------------------------------------------------------------------*/
find_edge_sphere_intersection(atom_id, edge_ptr)
int		atom_id;
struct redge    *edge_ptr;
{

  double	*point1, *point2;
  double 	a, b, c;
  double	diff_pt_pt, diff_pt_ct;
  double	tempd, t[2];
  int		i;

  a = b = 0;
  c = -(SQ(Probe_radius + atoms[atom_id].radius));

  point1 = edge_ptr->end_pts[0]->v;
  point2 = edge_ptr->end_pts[1]->v;

  for(i = 0; i < 3; i++)
  { diff_pt_pt = point1[i] - point2[i];
    /*point1 and point2 are specified in the local coord system 
      of atom, so the center of the atom can be assumed to be the 
      origin (0, 0, 0) for now.
    */
    diff_pt_ct = point2[i]/* - 0 */; 
    a += SQ(diff_pt_pt);
    b += diff_pt_pt*diff_pt_ct;
    c += SQ(diff_pt_ct);
    /* initialize int_pts to impossible values */
    edge_ptr->int_pts[0][i] = HUGE_VAL;
    edge_ptr->int_pts[1][i] = HUGE_VAL;
  }
  b *= 2;

  /* solve the quad eq: a*t^2 + b*t + c = 0 */
  tempd = b*b - 4*a*c;

  if (tempd >= 0) /* the line given by the edge intersects the sphere */
  { 
    tempd = sqrt(tempd);
    t[0] = (-b + tempd)/(2*a);
    t[1] = (-b - tempd)/(2*a);

    if ((t[0] <= 0)||(t[0] >= 1)) t[0] = -1;
    if ((t[1] <= 0)||(t[1] >= 1)) t[1] = -1;

    if ((t[0] > 0) && (t[1] > 0)) 
    { if (t[0] < t[1]) 
      { tempd = t[1]; t[1] = t[0]; t[0] = tempd; }
    }
    else
    { if (t[1] > 0) 
      { t[0] = t[1]; t[1] = -1; }
    }

    if (t[0] > 0)
    { for(i = 0; i < 3; i++)
       edge_ptr->int_pts[0][i] = point2[i] + t[0]*(point1[i]-point2[i]);
      if (t[1] > 0)
        for(i = 0; i < 3; i++)
          edge_ptr->int_pts[1][i] = point2[i] + t[1]*(point1[i]-point2[i]);
    }
  }
}

/*---------------------------------------------------------------------------------
print_cell prints out the vertices of the current feasible cell.
----------------------------------------------------------------------------------*/
print_cell()
{ 
  int		i;

  for(i = 0; i < rvertices_count; i++)
    printf("Cell vert[%d]: %lf %lf %lf\n", i,
	    rvertices[i].v[X], rvertices[i].v[Y], rvertices[i].v[Z]);
}

/*---------------------------------------------------------------------------------
display_face prints out the triangles for a given face of the current feasible cell.
----------------------------------------------------------------------------------*/
display_face(atom_id, region_face, cons)
int		atom_id;
struct rface	region_face;
Vector		*cons;
{ 
  VertexType	center;
  VertexType	pt0, pt1;
  int		i, j;

  VEC3_ZERO(center.Coord);
  for(i = 0; i < region_face.num_verts; i++)
    VEC3_ASN_OP(center.Coord, +=, region_face.vert[i]->v);
  VEC3_V_OP_V_OP_S(center.Coord, atoms[atom_id].center, +, center.Coord, /, region_face.num_verts);

  for(i = 0; i < region_face.num_verts; i++)
  { j = (i == region_face.num_verts-1)? 0: i+1;
    VEC3_V_OP_V(pt0.Coord, region_face.vert[i]->v, +, atoms[atom_id].center);
    VEC3_V_OP_V(pt1.Coord, region_face.vert[j]->v, +, atoms[atom_id].center);
    VEC3_ASN_OP(pt0.Normal, =, cons[region_face.id].coeff);
    VEC3_ASN_OP(pt1.Normal, =, cons[region_face.id].coeff);
    VEC3_ASN_OP(center.Normal, =, cons[region_face.id].coeff);
    gen_tris(&pt0, &pt1, &center);
  }
}
