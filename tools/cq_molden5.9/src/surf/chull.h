/*------------------------------------------------------------------    
  chull.h

  defines variables, macros, and structures used in convex hull code
--------------------------------------------------------------------*/


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* Define Boolean type */
typedef int FLAG;

/* general-purpose macros */
#define SWAP(t, x, y)	{ (t) = (x); (x) = (y); (y) = (t); }

#define ALLOCATE(p,type) { p=(type *) malloc (sizeof(type));	\
			   if (p == NULL) {			\
				printf ("Out of Memory!\n");	\
				exit(0);		        \
			   }					\
			 }

#define FREE(p)		if (p) { free ((char *) p); p = NULL; }


#define ADD_QUEUE( head, p )  { if ( head )  { \
                                    p->next = head->next; \
                                    p->prev = head; \
                                    head->next = p; \
                                    p->next->prev = p; \
                                    } \
                                else  { \
				    head = p; \
                                    head->next = head->prev = head; \
				    } \
                              }

#define DEL_QUEUE( head, p )  { if ( head )  { \
                                     if ( head == head->next ) \
                                          head = NULL;  \
                                     else if ( p == head ) \
                                          head = head->next; \
                                     p->next->prev = p->prev;  \
                                     p->prev->next = p->next;  \
                                     FREE( p ); \
                                     } \
                              }

#define is_visible( f, vt )   (f->p[X]*(vt->v[X] - f->vert[0]->v[X]) + \
			       f->p[Y]*(vt->v[Y] - f->vert[0]->v[Y]) + \
			       f->p[Z]*(vt->v[Z] - f->vert[0]->v[Z]) < 0)
#define volume( vol, f, vt )  { vol = f->p[X]*(vt->v[X] - f->vert[0]->v[X]) + \
				      f->p[Y]*(vt->v[Y] - f->vert[0]->v[Y]) + \
				      f->p[Z]*(vt->v[Z] - f->vert[0]->v[Z]);  \
                              }


/* Define structures for vertices, edges and faces */
struct tvertex {
	  double v[3];
          int  vnum[3]; /* store i, i^2 i^3 for use in perturbation*/
	  int  data_ptr;
          struct tedge *duplicate;
          FLAG active, mark;
          struct tvertex *next, *prev;
          };

struct tedge {
	  int	ednum;
          struct tface *adjface[3];
          struct tvertex *endpts[2];
          FLAG deleted;
          struct tedge *next, *prev;
          };

struct tface {
          struct tedge *edg[3];
          struct tvertex *vert[3];
	  int  fnum;
          FLAG visible;
          double p[3];
	  double sq_dist; /* distance of the face from the origin */
          struct tface *next, *prev;
          };

/* Define flags */
#define ACTIVE   1
#define DELETED  1
#define MARKED   1

/* Global variable definitions */
#ifdef global
struct tvertex *vertices = NULL;
struct tvertex *perm_vert= NULL;
struct tedge   *edges    = NULL;
struct tface   *faces    = NULL;
FLAG debug = 0; /* FALSE */
FLAG test  = 0; /* FALSE */
#else
extern struct tvertex *vertices;
extern struct tvertex *perm_vert;
extern struct tedge   *edges;
extern struct tface   *faces;
extern FLAG debug;
extern FLAG test;
#endif
