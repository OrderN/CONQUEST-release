/*====================================================================
  __file: Make.c
  This file contains the functions that build up the convex hull
  ===================================================================*/

#include "surf.h"
#include "chull.h"

/*----------------------------------------------------------------------
       Make_structs makes a new face and two new edges between the 
  edge and the point that are passed to it. It returns a pointer to
  the new face.
 ----------------------------------------------------------------------*/
struct tface *make_structs( e, p )
struct tedge *e;
struct tvertex *p;
{
	struct tedge *new_edge[2];
        struct tedge *make_edge();
        struct tface *new_face;
        struct tface *make_face();
        int i, j;

        for ( i=0; i < 2; ++i ) 
              if ( !( new_edge[i] = e->endpts[i]->duplicate) ) {
                 /* if the edge does not already exist, make it */
                 new_edge[i] = make_edge();
                 new_edge[i]->endpts[0] = e->endpts[i];
                 new_edge[i]->endpts[1] = p;
                 e->endpts[i]->duplicate = new_edge[i];
                 }

        /* make the new face */
        new_face = make_face();   
        new_face->edg[0] = e;
        new_face->edg[1] = new_edge[0];
        new_face->edg[2] = new_edge[1];
        make_ccw( new_face, e, p ); 
        pre_vol( new_face );
        
        /* set the adjacent face pointers */
        for ( i=0; i < 2; ++i )
              for ( j=0; j < 3; ++j )  
                    if ( !new_edge[i]->adjface[j] ) {
                          new_edge[i]->adjface[j] = new_face;
                          break;
                          }
        
        return new_face;
}

/*------------------------------------------------------------------------
      Make_ccw puts the vetices in the face structure in counter
  clockwise order.  If there is no adjacent face[1] then we know that
  we are working with the first face of the initial tetrahedron.  In this
  case we want to store the vertices in the opposite order from the 
  initial face.  Otherwise, we want to store the vertices in the same order
  as in the visible face.  The third vertex is always p.
  ------------------------------------------------------------------------*/
make_ccw( f, e, p )
struct tface *f;
struct tedge *e;
struct tvertex *p;
{
        register int i;
        register struct tface *fv;
        
        if ( !e->adjface[1] ) {
             /* if this is the initial triangle */
             fv = e->adjface[0];
        
             /* find the index of endpoint[1] */
             for ( i=0; fv->vert[i] != e->endpts[1]; ++i )
                   ;
             /* put the vertices in the opposite order of fv */
             if ( fv->vert[ (i+1) % 3 ] != e->endpts[0] ) {
                  f->vert[0] = e->endpts[1];  
                  f->vert[1] = e->endpts[0];    
                  }
             else {                               
                  f->vert[0] = e->endpts[0];   
                  f->vert[1] = e->endpts[1];      
                  }
             }
        
        else {
             /* otherwise,  set the visible face */
             if  ( e->adjface[0]->visible )      
                   fv = e->adjface[0];
             else fv = e->adjface[1];
        
             /* find the index of endpoint[1] */
             for ( i = 0; fv->vert[i] != e->endpts[1] ; ++i )  
                   ;    				
        
             /* put the vertices in the same order as fv */
             if ( fv->vert[ (i+1) % 3 ] != e->endpts[0] ) {
                  f->vert[0] = e->endpts[0];     
                  f->vert[1] = e->endpts[1];     
                  }
             else {     
                  f->vert[0] = e->endpts[1];     
                  f->vert[1] = e->endpts[0]; 
                  } 
              }
 
        f->vert[2] = p;
}
 
/*---------------------------------------------------------------------
      Make_edge creates a new cell and initializes all pointers to NULL
   and sets all flags to off.  It returns a pointer to the empty cell.
  ---------------------------------------------------------------------*/
struct tedge *make_edge()
{
	register struct tedge *new;

        ALLOCATE( new, struct tedge );
        new->adjface[0] = new->adjface[1] = new->adjface[2] = NULL;
        new->endpts[0] = new->endpts[1] = NULL;
        new->deleted = !DELETED;
        ADD_QUEUE( edges, new );
        return new;
}

/*---------------------------------------------------------------------
     Make_face creates a new face structure and initializes all of its
  flags to NULL and sets all the flags to off.  It returns a pointer
  to the empty cell.
 ----------------------------------------------------------------------*/
struct tface *make_face()
{
	register struct tface *new;
	register int i;

	ALLOCATE( new, struct tface);
	for ( i=0; i < 3; ++i ) {
    	      new->edg[i] = NULL;
    	      new->vert[i] = NULL;
    	      }
	new->visible = 0; /* FALSE */
	ADD_QUEUE( faces, new );
	return new;

}
/*====================================================================
  __file: Initial.c   This file contains the functions that build up the  
  the initial tetrahedron.
  ===================================================================*/

/*-----------------------------------------------------------------------
   Init_tet builds the initial tetrahedron.  It first finds 3 non-colinear
 points and makes a face out of them.   It next finds a fourth point that
 is not co-planar with the face.  The vertices are stored in the face
 structure in counter clockwise order so that the volume between the face
 and the point is negative. Lastly, the 3 newfaces to the fourth point
 are made and the data structures are cleaned up. 
 -----------------------------------------------------------------------*/
init_tet()
{
	register struct tvertex *v1, *v4, *t;
	register struct tface *f;
         	struct tface *make_face(), *make_structs();
	register struct tedge *e1, *e2, *e3, *s;
         	struct tedge *make_edge();
	double vol;
	
	v1 = vertices;

	/* find 3 non-colinear points */
	while ( co_linear( v1, v1->next, v1->next->next ) ) 
      	      if ( ( v1 = v1->next ) == vertices ) {
         	     printf("There are no 3 non-colinear points.\n");
         	     exit(0);
         	     } 
	
	/* mark the original vertices as used */
	v1->mark = MARKED;  v1->next->mark = MARKED; 
	v1->next->next->mark = MARKED;              
                                               
	/* make the edges of the original triangle */
	e1 = make_edge(); 
	e2 = make_edge();
	e3 = make_edge();
	e1->endpts[0] = v1;              e1->endpts[1] = v1->next;
	e2->endpts[0] = v1->next;        e2->endpts[1] = v1->next->next;
	e3->endpts[0] = v1->next->next;  e3->endpts[1] = v1;

	/* make the face of the original triangle */
	f = make_face();   
	f->edg[0] = e1;   f->edg[1] = e2;   f->edg[2] = e3;
	f->vert[0] = v1;  f->vert[1] = v1->next;
	f->vert[2] = v1->next->next;

	/* set the adjacent faces */
	e1->adjface[0] = e2->adjface[0] = e3->adjface[0] = f; 
	
	v4 = v1->next->next->next;
	pre_vol (f);
	
	/* find a fourth non-coplanar point */
	volume( vol, f, v4 );
	while ( vol <= 0.0 )   {
      		if ( ( v4 = v4->next ) == v1 ) {
         	       printf("There are no 4 non-coplanar points.\n");
         	       exit(0);
                       } 
       		volume( vol, f, v4 );
       	        }
	v4->mark = MARKED;

	/* store points in counter clockwise order */
	if( vol < 0 ) {
            SWAP( t, f->vert[0], f->vert[1] );
            SWAP( s, f->edg[1], f->edg[2] );
  	    }

	pre_vol( f );
	
	/* make the faces and edges between the original
   	   triangle and the fourth point. */
	e1->adjface[1] = make_structs( e1, v4 );
	e2->adjface[1] = make_structs( e2, v4 );
	e3->adjface[1] = make_structs( e3, v4 );

	clean_up();
	
	if ( debug )  print_out( v4 );
	
}

/*===================================================================
  __file: Degen.c  This file contains the functions that check for 
  co-linear and co-planar points.
  ===================================================================*/


/*------------------------------------------------------------------
     Co_linear checks to see if the three points given are colinear.
  It checks to see if the cross product  ( rather, each element of
  the cross product  i, j, and k ) between the two vectors formed by
  the points is 0 or not.  If each element is 0 then the area of the
  triangle is 0 and the points are colinear.
---------------------------------------------------------------------*/
co_linear( a, b, c )
struct tvertex *a, *b, *c;
{
	return ( c->v[Z] - a->v[Z] ) * ( b->v[Y] - a->v[Y] ) -
               ( b->v[Z] - a->v[Z] ) * ( c->v[Y] - a->v[Y] ) == 0
            && ( b->v[Z] - a->v[Z] ) * ( c->v[X] - a->v[X] ) -
               ( b->v[X] - a->v[X] ) * ( c->v[Z] - a->v[Z] ) == 0
            && ( b->v[X] - a->v[X] ) * ( c->v[Y] - a->v[Y] ) -
               ( b->v[Y] - a->v[Y] ) * ( c->v[X] - a->v[X] ) == 0  ;
}

      
/*----------------------------------------------------------------------
       Pre_vol calculates the volume information that corresponds to 
   each face and stores it in the face structure.
  ----------------------------------------------------------------------*/
pre_vol( f )
struct tface *f;
{
        double x1, x2, x3, y1, y2, y3, z1, z2, z3;

        x1 = f->vert[0]->v[0];
        y1 = f->vert[0]->v[1];
        z1 = f->vert[0]->v[2];

        x2 = f->vert[1]->v[0] - x1;
        y2 = f->vert[1]->v[1] - y1;
        z2 = f->vert[1]->v[2] - z1;

        x3 = f->vert[2]->v[0] - x1;
        y3 = f->vert[2]->v[1] - y1;
        z3 = f->vert[2]->v[2] - z1;
        
        f->p[X] =  y3*z2 - y2*z3;
        f->p[Y] =  x2*z3 - x3*z2;
        f->p[Z] =  x3*y2 - x2*y3; 

}

/*=========================================================================
  __file: Hull.c    This file contains the functions that make the rest of
  the hull after the initial tetrahedron is built.
 =========================================================================*/

/*-------------------------------------------------------------------------
        Complete_hull goes through the vertex list and adds them to the hull
   if they are not already used.  It will mark the vertex once it is 
   looked at.
  -------------------------------------------------------------------------*/
complete_hull()
{
	struct tvertex *v, *tmp_vert;

	v = vertices;
	do {
           if ( !v->mark ) {
                 v->mark = MARKED;
                 if (add_on( v ))
                 { if ( debug ) print_out( v );
                   clean_up();
                   v = v->next;
		 }
		 else
		 { tmp_vert = v;
		   v = v->next;
		   DEL_QUEUE(vertices, tmp_vert);
                 }
	       }
	       else
		 v = v->next;
            } while ( v != vertices );

	 clean_up();
}

/*---------------------------------------------------------------------------
      Add_on is passed a vertex.  It will first determine all faces that
  are visible from that point.  If none are visible then the point is 
  marked not active and 0 is returned (else a 1 is returned).  Then it will 
  go through the edge list.  If both of the edges adjacent faces are  
  visible then the edge is marked deleted.  If one of the adjacent faces is 
  visible then a new face is made and the edges adjacent face[2] pointer is
  set to point to it. A temporary pointer is used so that any new edge just 
  added to the list will not be examined.
  --------------------------------------------------------------------------*/
int
add_on( p )
struct tvertex *p;
{
	register struct tface *f; 
	         struct tface *make_structs();
	register struct tedge *e;
        struct tedge	*temp;
        int vis = 0; /* FALSE */

	f = faces; 
        do {
             /* mark the visible faces and set a flag */
	     vis |= f->visible = is_visible(f, p);
   	     f = f->next;
   	} while ( f != faces );

        if ( !vis ) {
              /* if no faces are visible then */
              /* mark the vertex inactive */
              p->active = !ACTIVE;  
              return (0); 
              }

        e = edges;
	do {
           temp = e->next;
           if ( e->adjface[0]->visible && e->adjface[1]->visible )
                /* if e is an interior edge, delete it */
                e->deleted = DELETED;

           else if ( e->adjface[0]->visible || e->adjface[1]->visible ) 
                /* if e is a border edge, make a new face */
                e->adjface[2] = make_structs( e, p );
           e = temp;
           } while ( e != edges );

       return(1);
}


/* ================================================================== 
   __file: Debug.c   This file is used whenever the debug flag is set.
    It contains the functions to print out the entire contents of
    each data structure.  It prints into standard error.  
   ==================================================================*/

/*-------------------------------------------------------------------*/
print_out( v )
struct tvertex *v;
{
	fprintf( stderr, "\nAdding vertex %6x :\n", v );
	print_verts();
	print_edges();
	print_fs();
}

/*------------------------------------------------------------------
       Print_verts prints out the contents of each tvertex structure
   that is in the list.  The address of the vertex is printed,
   then the x, y, and z coordinates is printed. Either a 0 or a 1 is
   printed for the active flag and the value of the pointer in
   hexidecimal is printed for the duplicate pointer. 
  -----------------------------------------------------------------*/
print_verts()
{
	struct tvertex *temp;

	temp = vertices;
	fprintf (stderr, "Vertex List\n");
	if (vertices) do {
            fprintf(stderr,"  addr %6x\t", vertices );
            fprintf(stderr,"(%g,%g,%g)",vertices->v[X],
                    vertices->v[Y], vertices->v[Z] );
            fprintf(stderr,"   active:%3d", vertices->active );
            fprintf(stderr,"   duplicate:%5x", vertices->duplicate );
            fprintf(stderr,"   mark:%2d\n", vertices->mark );
            vertices = vertices->next;
            } while ( vertices != temp );

}

/*-----------------------------------------------------------------------
      Print_edges prints out the information of each tedge structure that
  is in the circular list.  For the pointers adjface and endpts, the
  value of the pointer in hexidecimal is printed ( 0 if it is NULL ).
  The flag deleted is either 0 or 1. The address of the edge is also
  printed.
  ------------------------------------------------------------------------*/
print_edges()
{
	register struct tedge *temp;
	register int i;
	
	temp = edges;
	fprintf (stderr, "Edge List\n");
	if (edges) do {
            fprintf( stderr, "  addr: %6x\t", edges );
            fprintf( stderr, "adj: ");
            for (i=0; i<3; ++i) 
                 fprintf( stderr, "%6x", edges->adjface[i] );
            fprintf( stderr, "  endpts:");
            for (i=0; i<2; ++i) 
                 fprintf( stderr, "%8x", edges->endpts[i]);
            fprintf( stderr, "  del:%3d\n", edges->deleted );
            edges = edges->next; 
            } while (edges != temp );

}

/*----------------------------------------------------------------------
       Print_fs prints out the information of each tface structure
   that is in the circular list.  For the pointers edg and vert the 
   value of the pointer in hexidecimal is printed.  For the flag
   visible either a 0 or a 1 is printed. The address of the face is 
   also printed.
  ---------------------------------------------------------------------*/
print_fs()
{
	register int i;
	register struct tface *temp;

	temp = faces;
	fprintf (stderr, "Face List\n");
	if (faces) do {
            fprintf(stderr, "  addr: %6x\t", faces );
            fprintf(stderr, "  edges:");
            for( i=0; i<3; ++i )
                 fprintf(stderr, "%6x", faces->edg[i] );
            fprintf(stderr, "  vert:");
            for ( i=0; i<3; ++i)
                  fprintf(stderr, "%6x", faces->vert[i] );
            fprintf(stderr, "  vis: %d\n", faces->visible );
            faces= faces->next;
            } while ( faces != temp );

}
/* ================================================================== 
    __file: Data.c   This file contains the functions that take care of
    the data structures-  reading the information in, storing it,
    and printing it out.         
   ==================================================================*/

/*--------------------------------------------------------------------
     Read_vertices reads the vertices in and puts them in a circular,
  doubly linked list of tvertex structures.  The duplicate pointer is
  set to NULL and the active flag is set off.  
  --------------------------------------------------------------------*/ 
read_vertices()
{

	register struct tvertex *new;
	double x, y, z;
        int    i = 0;

	while ( scanf ("%lf%lf%lf", &x, &y, &z ) != EOF )  {
                ALLOCATE( new, struct tvertex );
                new->v[X] = x;
                new->v[Y] = y;
                new->v[Z] = z;
                new->vnum[0] = new->vnum[1] = new->vnum[2] = ++i; 
		new->vnum[1] *= new->vnum[0];
		new->vnum[2] *= new->vnum[1];
                new->duplicate = NULL;
                new->active = !ACTIVE;
                new->mark = !MARKED;
                ADD_QUEUE( perm_vert, new );
                }

}

/*--------------------------------------------------------------------
  Copies the vertices from the perm_vert queue into the vertices queue.
  --------------------------------------------------------------------*/ 
copy_vertices()
{

	register struct tvertex *new, *v1;

	v1 = perm_vert;
	do {
          ALLOCATE( new, struct tvertex );
	  memcpy(new, v1, sizeof (struct tvertex));
          ADD_QUEUE( vertices, new );
        } while((v1 = v1->prev)!= perm_vert);

}

/*--------------------------------------------------------------------
       Print prints the output in correct format.  It first determines
  the total number of vertices and prints that number out.  Then the
  x, y, and z coordinates of each vertex are printed.  Next it finds
  the total number of faces and then for each face it prints out the
  cooresponding vertices.  It also finds the total number of edges
  and calls the routines that check the number of edges and faces.
  --------------------------------------------------------------------*/
print()
{
	register struct tvertex *temp_v;
	register struct tface   *temp_f;
	register int i = 0, j = 0 ;

	printf("\n");
        /* total the vertices, then total the faces and print header line */
	temp_v = vertices; 
        temp_f = faces   ;
             do {                                 
                temp_v->vnum[0] = i++;           
                temp_v = temp_v->next;
                } while ( temp_v != vertices );
             do {
                ++j;                              
                temp_f  = temp_f ->next;
                } while ( temp_f  != faces );
       	     printf("%d %d %d\n", i,j, (i+j-2)  );
        
	/* print the vertices and print out each face's vertices */
	temp_v = vertices; 
        temp_f = faces   ;
     	     do {                                 
                printf("%g %g %g\n", temp_v->v[X],
                         temp_v->v[Y], temp_v->v[Z] );
                temp_v = temp_v->next;
                } while ( temp_v != vertices );
             do {                           
                printf("3%5d%6d%6d\n", temp_f->vert[0]->vnum,
                        temp_f->vert[1]->vnum, temp_f->vert[2]->vnum );
                temp_f = temp_f->next;
                } while ( temp_f != faces );
}
/* ================================================================== 
    __file: Clean.c   This file contains the functions that take care of
    cleaning up the data structures. 
   ==================================================================*/

/*-------------------------------------------------------------------------
       Meta_clean goes through each data structure list and clears all
   flags and NULLs out all pointers and frees them.  
  ------------------------------------------------------------------------*/
meta_clean()
{
  register struct tedge *e;
  register struct tvertex *v;
  register struct tface *f;

  while ( edges ) {
    e = edges;
    DEL_QUEUE ( edges, e );
  }  

  while ( vertices ) {
    v = vertices;
    DEL_QUEUE ( vertices, v );
  }  

  while ( faces ) {
    f = faces;
    DEL_QUEUE ( faces, f );
  }  
}
/*-------------------------------------------------------------------------
       Clean_up goes through each data structure list and clears all
   flags and NULLs out some pointers.  In the vertex list,  the active
  flag is set to 0 and the duplicate pointer is set to NULL. In the
  edge list any struture that is marked deleted gets deleted from the 
  edge list, and in the vertex list any face that is marked visible gets
  deleted from the face list.  
  ------------------------------------------------------------------------*/
clean_up()
{

	clean_edges();
	clean_faces();
	clean_vertices();

}

/*------------------------------------------------------------------------
      Clean_edges runs through the edge list and cleans up the structure.
   If there is a face in the edges adjacent face[2] then it will put
   that face in place of the visible face and set adjacent face[2] to 
   NULL.  It will also delete any edge marked deleted.
  -----------------------------------------------------------------------*/
clean_edges()
{
	register struct tedge *e, *tmp_edge;
	
        /* replace the visible face with the new face */
	if ( e = edges )
    	     do {
                if ( e->adjface[2] ) { 
                     if ( e->adjface[0]->visible )
                          e->adjface[0] = e->adjface[2]; 
                     else e->adjface[1] = e->adjface[2];
                     e->adjface[2] = NULL;
                     }
                  e = e->next;
                  } while ( e != edges );

        /* delete any edges marked deleted */
	while ( edges && edges->deleted ) {
                e = edges;
                DEL_QUEUE ( edges, e );
                }
        if ( edges )  {
             e = edges->next;
             do {
                if ( e->deleted ) {
                     tmp_edge = e;
                     e = e->next;
                     DEL_QUEUE( edges, tmp_edge );
                   }
                   else e = e->next;
                } while ( e != edges );

        }

}

/*------------------------------------------------------------------------
         Clean_faces runs through the face list and deletes any face
   marked visible.
  -----------------------------------------------------------------------*/
clean_faces()
{
	register struct tface *f, *tmp_face;

        while ( faces && faces->visible ) { 
                f = faces;
                DEL_QUEUE( faces, f );
                }
	if ( faces ) {
             f = faces->next;
             do {
                if ( f->visible ) {
                     tmp_face = f;
                     f = f->next;
                     DEL_QUEUE( faces, tmp_face );
                   }
                else f = f->next;
                } while ( f != faces );
        }

}

/*-------------------------------------------------------------------------
      Clean_vertices runs through the vertex list and deletes the 
   vertices that are marked as used but are not incident to any undeleted
   edges. It will also set the duplicate pointer to NULL and the active
   flag to not active.
  -------------------------------------------------------------------------*/
clean_vertices()
{
	struct tedge *e;
	struct tvertex *v, *tmp_vert;

        /* mark all vertices that are incident */
        /* to an edge as active */
	if ( e = edges )
             do {
	         e->endpts[0]->active = e->endpts[1]->active = ACTIVE;
	         e = e->next;
             } while (e != edges);

        /* delete all vertices that are */
        /* marked but not active */

        while ( vertices && vertices->mark && !vertices->active ) { 
                v = vertices;
                DEL_QUEUE( vertices, v );
                }

        if ( v = vertices ) { 
               do { 
                 if ( v->mark && !v->active ) {    
                     tmp_vert = v; 
                     v = v->next;
                     DEL_QUEUE( vertices, tmp_vert );
                   }
                 else v = v->next;
               } while ( v != vertices);
        }

        /* reset active flag and duplicate pointer */
	if ( v = vertices )
       	     do {
                v->duplicate = NULL; 
                v->active = !ACTIVE; 
                v = v->next;
                } while ( v != vertices );
   
}
 /*=====================================================================
    __file: Checks.c   This file is used whenever the test flag is set.
    These checks are done after the hull has been made.
   =====================================================================*/

/*------------------------------------------------------------------------
      Consistency runs through the edge list and checks that all
  adjacent faces have their endpoints in opposite order.  This verifies
  that the vertices are in counter clockwise order.
  -----------------------------------------------------------------------*/
consistency()
{
	register struct tedge *e;
	register int i, j;

	e = edges;

	do {
    	   /* find index of endpoint[0] in adjacent face[0] */
           for ( i = 0; e->adjface[0]->vert[i] != e->endpts[0]; ++i )
                 ;
   
   	   /* find index of endpoint[0] in adjacent face[1] */
   	   for ( j = 0; e->adjface[1]->vert[j] != e->endpts[0]; ++j )
         	 ;

           /* check if the endpoints occur in opposite order */
           if ( !( e->adjface[0]->vert[ (i+1) % 3 ] ==
                   e->adjface[1]->vert[ (j+2) % 3 ] ||
                   e->adjface[0]->vert[ (i+2) % 3 ] ==
                   e->adjface[1]->vert[ (j+1) % 3 ] )  )
               break;
           e = e->next;

           } while ( e != edges );

        if ( e != edges )
             fprintf( stderr, "Edges are NOT consistent.\n");

}

/*----------------------------------------------------------------------
       Convexity checks that the volume between every face and every
  point is negative.  This shows that each point is inside every face
  and therefore the hull is convex.
  ---------------------------------------------------------------------*/
convexity()
{
	register struct tface *f;
        register struct tvertex *v;
        int vol;

        f = faces;

        do {
              v = vertices;
              do {
                   volume( vol, f, v );
                   if ( vol < 0 )
                        break;
                   v = v->next;
                   } while ( v != vertices );
              if ( v != vertices )
                   break;
            } while ( f != faces );

        if ( f != faces )
           fprintf( stderr, "It is NOT convex.\n");
        
}

/*----------------------------------------------------------------------
     Check_faces checks that the number of faces is correct according
 to the number of vertices that are in the hull.
  ---------------------------------------------------------------------*/
check_faces(cvertices, cfaces )
int cvertices, cfaces;
{

	if ( cfaces != 2 * cvertices - 4 ) 
             fprintf( stderr, "The number of faces is NOT correct\n");

}

/*----------------------------------------------------------------------
      Check_edges checks that the number of edges is correct according
  to the number of  vertices that are in the hull.
  ---------------------------------------------------------------------*/
check_edges( cvertices, cedges )
int cvertices, cedges;
{

	if ( 2 * cedges != 3 * cvertices ) 
	     fprintf( stderr, "The number of edges is NOT correct.\n");
}
/*-----------------------------------------------------------------------*/
checks()
{
	consistency();
	convexity();
}
