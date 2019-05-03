/* Define structures for vertices, edges and faces of the feasible region*/

struct rvertex {
          double 	v[3];
	  double 	sq_dist;
	  struct rface  *adj_faces[3];
	  struct redge  *in_edges[3];
          };

struct rface {
          struct rvertex **vert;
	  struct redge   **edge; /* edge[i] = (vert[i], vert[(i+1)%num_v]) */
	  int	 	 num_verts;
	  int	         id;
          };

struct redge {
	  struct rface   *adj_faces[2];
	  struct rvertex *end_pts[2];
	  double 	 int_pts[2][3];
	  };

#ifdef global
struct rvertex *rvertices = NULL;
struct rface   *rfaces    = NULL;
struct redge   *redges    = NULL;
int    	        redges_count;
int		rfaces_count;
int		rvertices_count;
#else
extern struct rvertex *rvertices;
extern struct rface   *rfaces;
extern struct redge   *redges;
extern int	       redges_count;
extern int	       rfaces_count;
extern int	       rvertices_count;
#endif
