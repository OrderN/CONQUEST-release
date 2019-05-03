
/* Macros and stuff for linear algebra operations with vectors and matrices. */

/* Initial declarations necessary for linear algebra macros */
#define EXTERN_LINALG	EXTERN double LAmag,LAsum;EXTERN int LAi,LAj,LAk;

/*==================== Matrix Declarations ================*/
typedef double  MatrixD3[3][3];
typedef double  MatrixD4[4][4];
typedef double  PointD3[3];
typedef double  VectorD3[3];

/*==================== 3D vector macros ===================*/
#define VEC3_ZERO(a)	       { a[0]=a[1]=a[2]=0; }
#define VEC3_NEG(a,b)           { a[0]= -b[0]; a[1]= -b[1];a[2]= -b[2];}
#define VEC3_EQ(a,b)           ((a[0]==b[0]) && (a[1]==b[1]) && (a[2]==b[2]))
#define ZERO3(a)               { a[0] = ((a[0]<1e-5)&&(a[0]>-1e-5))?0.0:a[0];\
				 a[1] = ((a[1]<1e-5)&&(a[1]>-1e-5))?0.0:a[1];\
			         a[2] = ((a[2]<1e-5)&&(a[2]>-1e-5))?0.0:a[2];\
			       }

#define VEC3_V_OP_S(a,b,op,c)  {  a[0] = b[0] op c;  \
				  a[1] = b[1] op c;  \
				  a[2] = b[2] op c;  }

#define VEC3_V_OP_V(a,b,op,c)  { a[0] = b[0] op c[0]; \
				 a[1] = b[1] op c[1]; \
				 a[2] = b[2] op c[2]; \
				}
#define VEC3_V_OP_V_OP_S(a,b,op1,c,op2,d)  \
				{ a[0] = b[0] op1 c[0] op2 d; \
				  a[1] = b[1] op1 c[1] op2 d; \
				  a[2] = b[2] op1 c[2] op2 d; }

#define VEC3_V_OP_V_OP_V(a,b,op1,c,op2,d)  \
				{ a[0] = b[0] op1 c[0] op2 d[0]; \
				  a[1] = b[1] op1 c[1] op2 d[1]; \
				  a[2] = b[2] op1 c[2] op2 d[2]; }

#define VEC3_ASN_OP(a,op,b)      {a[0] op b[0]; a[1] op b[1]; a[2] op b[2];}

#define DOTPROD3(a, b)		 (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#define CROSSPROD3(a,b,c)       {a[0]=b[1]*c[2]-b[2]*c[1]; \
                                 a[1]=b[2]*c[0]-b[0]*c[2]; \
                                 a[2]=b[0]*c[1]-b[1]*c[0];}

#define NORMALIZE3(a)		{LAmag=1./sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);\
				 a[0] *= LAmag; a[1] *= LAmag; a[2] *= LAmag;}

#define NORMALIZE_PLANE3(a)	{LAmag=1./sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);\
				 a[0] *= LAmag; a[1] *= LAmag; \
				 a[2] *= LAmag; a[3] *= LAmag;}

#define SQ_DIST3(a, b)          ((a[0]-b[0])*(a[0]-b[0]) +      \
                                 (a[1]-b[1])*(a[1]-b[1]) +      \
                                 (a[2]-b[2])*(a[2]-b[2]))

#define PRINT_VEC3(a,string)    { printf("%s: (%g %g %g)\n",string, a[0], \
				  			   a[1], a[2]);}

/*================= General vector macros ===================*/
#define VEC_ZERO(a,m)	         {for(LAi=0; LAi<m; a[LAi++] = 0); }
#define VEC_NEG(a,b,m)        	 {for(LAi=0; LAi<m; LAi++) a[LAi] = -b[LAi]; }

#define VEC_V_OP_S(a,b,op,c,m)  {for(LAi=0; LAi<m; LAi++)   \
				    a[LAi] = b[LAi] op c; }
#define VEC_V_OP_V(a,b,op,c,m)  {for(LAi=0; LAi<m; LAi++)   \
				    a[LAi] = b[LAi] op c[LAi]; }
#define VEC_V_OP_V_OP_S(a,b,op1,c,op2,d) \
				{for(LAi=0; LAi<m; LAi++)   \
				    a[LAi] = b[LAi] op1 c[LAi] op2 d; }
#define VEC_V_OP_V_OP_V(a,b,op1,c,op2,d) \
				{for(LAi=0; LAi<m; LAi++)   \
				    a[LAi] = b[LAi] op1 c[LAi] op2 d[LAi]; }
#define VEC_ASN_OP(a,op,b,m){for(LAi=0; LAi<m; LAi++) a[LAi] op b[LAi]; }

#define NORMALIZE(a,m)		{ for(LAi=0,LAmag=0.;LAi<m;LAi++) \
				    LAmag += a[LAi]*a[LAi];	  \
				  LAmag=1./sqrt(LAmag);		  \
				  for(LAi=0; LAi<m; a[LAi++] *= LAmag); }
#define PRINT_VEC(a, m, string)	{ printf("%s : (", string);	  \
				  for(LAi=0; LAi<m; LAi++)	  \
					printf("%g ",a[LAi]);	  \
				  printf(")\n"); }
/*================= Matrix macros ===================*/

#define ALLOCMAT(a, m, n, type) {ALLOCN(a, type*, m);	  \
				 for(LAi=0; LAi<m; LAi++) 	  \
				   ALLOCN(a[LAi], type, n);}
#define FREEMAT(a, m, n)        {for(LAi=0; LAi<m; LAi++) free(a[LAi]);	  \
				 free(a); }

/* a = identity */
#define IDENTMAT(a,m)		{for(LAi=0; LAi<m; LAi++) \
                                    for(LAj=0; LAj<m; LAj++) \
                                       a[LAi][LAj] = (LAi==LAj) ? 1. : 0.;}
/* vec a = mat b * vec c; mat b is m x n matrix */
#define MATVECMULT(a,b,c,m,n)	{for(LAi=0; LAi<m; LAi++) {\
				       for(LAj=0,LAsum=0.; LAj<n; LAj++) \
				          LAsum+=b[LAi][LAj]*c[LAj]; \
                                       a[LAi] = LAsum; \
                                    }}
/* a = b * c  (b is m x n and c is n x p) */
#define MATMULT(a,b,c,m,n,p)	{for(LAi=0; LAi<m; LAi++) \
				    for(LAj=0; LAj<p; LAj++) { \
				       for(LAk=0,LAsum=0.; LAk<n; LAk++) \
				          LAsum+=b[LAi][LAk]*c[LAk][LAj]; \
                                       a[LAi][LAj] = LAsum; \
                                    }}
/* a = b (m x n matrices)*/
#define MATASSIGN(a,b,m,n)	{for(LAi=0; LAi<m; LAi++) \
				    for(LAj=0; LAj<n; LAj++) { \
                                       a[LAi][LAj] = b[LAi][LAj]; \
                                    }}
/*      t		 */
/* a = b ( b is m x n )*/
#define MATTRANSPOSE(a,b,m,n)	{for(LAi=0; LAi<m; LAi++) \
				    for(LAj=0; LAj<n; LAj++) { \
                                       a[LAj][LAi] = b[LAi][LAj]; \
                                    }}

/* b[][a] = c, length of column vector is m */
#define SETCOL(a,b,c,m)         {for(LAi=0; LAi<m; LAi++) \
				    b[LAi][a] = c[LAi]; \
				}
/* c = b[][a], length of column vector is m */
#define GETCOL(a,b,c,m)         {for(LAi=0; LAi<m; LAi++) \
				    c[LAi] = b[LAi][a]; \
				}
/* b[a][] = c, length of row vector is n */
#define SETROW(a,b,c,n)         {for(LAi=0; LAi<n; LAi++) \
				    b[a][LAi] = c[LAi]; \
				}
/* c = b[a][], length of row vector is n */
#define GETROW(a,b,c,n)         {for(LAi=0; LAi<n; LAi++) \
				    c[LAi] = b[a][LAi]; \
				}
#define PRINT_MAT(a, p, q, string){ printf("%s :\n", string);	  \
				  for(LAi=0; LAi<p; LAi++)	  \
				  { for(LAj=0; LAj<q; LAj++)	  \
					printf("%g ",a[LAi][LAj]);\
				    printf("\n");		  \
				  }}

#define DET3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3)                \
                ( (a1)*((b2)*(c3) - (c2)*(b3))                    \
                 -(b1)*((a2)*(c3) - (a3)*(c2))                    \
                 +(c1)*((a2)*(b3) - (a3)*(b2)))                   
