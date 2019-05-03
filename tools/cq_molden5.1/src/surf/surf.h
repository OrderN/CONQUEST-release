#include	<stdio.h>
#include	<math.h>
#include	<stdlib.h>
#include	<string.h>
#include        <sys/time.h>
#include 	"linalg.h"

#ifdef EXTERN
#undef EXTERN
#define EXTERN
#else
#define EXTERN extern
#endif


/*------------------ various consts ----------------------------------*/

#define 	X   		0
#define 	Y   		1
#define 	Z   		2

#define		TRUE		1
#define		FALSE		0

/* constants defined to prevent unnecessary double->float ops on i860's */
#define         ZERO            (float)(0.0)
#define		EPS		(float)(1e-5)
#define		LP_EPS		(float)(1e-5)
#define		GP_EPS	        (float)(1e-3)
#define		EPSILON		(float)(1e-7)
#define		HALF	        (float)(0.5)
#define         FCEIL           (float)(0.999999)
#define         ONE             (float)(1.0)
#define         FCEIL_PLUS_1    (float)(1.999999)
#define         PI      	(float)(3.14159265358979323846)
#define         TWO_PI      	(float)(6.28318530717958647692)
#define 	DEG_TO_RAD      (float)(.01745329251994329576) /* pi/180 */

#define         MINC            0
#define         MAXC            1

#define		MAX_TESS_LEN	1.2

#define		MAX_TYPES	5

#define         MAX_DENSITY	15 /* max no of atoms in a voxel(odd for algn) */

#define		MAX_CONSTRAINT 	4000 /* 700 */
#define		PERT_Q		4001 /* smallest prime exceeding MAX_CONSTRAINT*/



/*-------------------- general-purpose macros ------------------------*/

#define START gettimeofday(&tm,&tz);\
                et = (tm.tv_sec)+ (0.000001* (tm.tv_usec));

#define STOP gettimeofday(&tm,&tz);\
                et = (tm.tv_sec)+(0.000001*(tm.tv_usec)) - et;

#define  ALLOCN(PTR,TYPE,N) 					\
	{ PTR = (TYPE *) malloc(((unsigned)(N))*sizeof(TYPE));	\
	  if (PTR == NULL) {    				\
	  printf("malloc failed");				\
	  exit(-1);                                             \
	  }							\
	}

#define FMAX(x,y) ((x)>(y) ? (x) : (y))
#define FMIN(x,y) ((x)<(y) ? (x) : (y))
#define FP_EQ_EPS( a, b, c )  ((((a) - (b)) < (c)) && (((a) - (b)) > -(c)))
#define INDEX(x, y, z, num)  (((x)*(num[1]) + (y))*(num[2]) + (z))
#define AINDEX(a, num)       (((a[0])*(num[1]) + (a[1]))*(num[2]) + (a[2]))
#define SQ(a)                ((a)*(a))

/*---------------------------structs-------------------------------------*/
typedef unsigned char   byte;

typedef double		POINT[3];

typedef struct 
{ float         	coeff[4];
  int           	atom_id;
} Vector;


typedef float VectorType[3];	/* vector in three-space */
typedef float PointType[3];     /* point in homogenous three-space */
typedef struct {
    PointType Coord;
    VectorType Normal;
} VertexType;

typedef struct
{ float			radius;
  float			center[3];
  float			tes_origin[3];
  short			num_cons;
  short			type; 
  short			boundary;
} Gp_Atom;

typedef struct 		
{ short			num_neighbors;
  short		        neighbor[MAX_DENSITY];
} Atom_List;

struct Torus
{ int                   face_atoms[2];
  int                   end_atoms[2];
  int                   face_id;
  int                   edge_id[2];
  int                   num_int_verts;
  struct Torus		*next;
  struct Torus		*prev;
};

typedef struct
{ VertexType		vt;
  POINT			tes_dir; 
} Big_Point;


/*---------------------------variables-------------------------------------*/
extern double           find_area();
extern double           find_angle();
extern double		sph_dist();
extern int              linear_prog();
extern int		find_circle();
extern int		find_vertex_order();

EXTERN float            Probe_radius;
EXTERN Gp_Atom          *atoms;
EXTERN int              Num_atoms;
EXTERN int              Checks_On;
EXTERN int              Max_Gp_Polys;
EXTERN double           Max_Tess_Len, Max_Tess_Len_Sq;
EXTERN VertexType       **verts;
EXTERN short            *atom_type;
EXTERN int              Num_polys;
EXTERN Vector           New_Origin;  
EXTERN int              Current_atom;
EXTERN FILE            	*Opf;
EXTERN Vector 		temp_cons[2][MAX_CONSTRAINT+6];
EXTERN int              Write_Option;
EXTERN struct timeval   tm;
EXTERN struct timezone  tz;
EXTERN double           et;
EXTERN_LINALG
