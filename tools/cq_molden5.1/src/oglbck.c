
/*

     ##  COPYRIGHT (C) 2000 by G. Schaftenaar  ##

     MOLecular DENsity OpenGL Helper Program
     By G. Schaftenaar
     CMBI, University of Nijmegen, The Netherlands
     (formerly CAOS/CAMM Center)
     2000
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static FILE *out;
static int dowrt = 0;
static double z0[] = {0.0,0.0,0.0};
static double zx[] = {1.0,0.0,0.0};
static double zy[] = {0.0,1.0,0.0};

#define MAXNAM 512
static char ogfil[MAXNAM];

#if defined(VMS) || defined(UNDERSC)
parogf(str,nlen)
#else
#ifdef CRAY
PAROGF(str,nlen)
#else
parogf_(str,nlen)
#endif
#endif

#ifdef VMS
struct dsc$descriptor_s *str;
#else
#ifdef CRAY
_fcd str;
#else
char *str;
#endif
#endif
int *nlen;

{
  int i;

  if (*nlen >= MAXNAM) { 
	fprintf(stderr,"filename too long !\n");
	*nlen = MAXNAM - 1;
  }
#ifdef VMS
  for (i=0; i<*nlen; i++)
         ogfil[i] = str->dsc$a_pointer[i];
#else
#ifdef CRAY
   strncpy(ogfil,_fcdtocp(str),*nlen);
#else
   for (i=0; i<*nlen+1; i++) ogfil[i] = '\0';
   strncpy(ogfil,str,*nlen);
#endif
#endif

}

#if defined(VMS) || defined(UNDERSC)
void ognorm(double *v1, double *v2, double *v3)
#else
#ifdef CRAY
void OGNORM(double *v1, double *v2, double *v3)
#else
void ognorm_(double *v1, double *v2, double *v3)
#endif
#endif
{
      if (dowrt) fprintf(out,"%f %f %f\n",*v1,*v2,*v3);
}

#if defined(VMS) || defined(UNDERSC)
void ogvert(double *v1, double *v2, double *v3)
#else
#ifdef CRAY
void OGVERT(double *v1, double *v2, double *v3)
#else
void ogvert_(double *v1, double *v2, double *v3)
#endif
#endif
{
      if (dowrt) fprintf(out,"%f %f %f\n",*v1,*v2,*v3);
}

#if defined(VMS) || defined(UNDERSC)
void ogcol(double *v1, double *v2, double *v3)
#else
#ifdef CRAY
void OGCOL(double *v1, double *v2, double *v3)
#else
void ogcol_(double *v1, double *v2, double *v3)
#endif
#endif
{
      if (dowrt) fprintf(out,"%f %f %f\n",*v1,*v2,*v3);
}

ogwrmol(double *r, double *adjus, int *natoms, int *nat, 
        double *xsym, double *ysym, double *zsym, double *vdwr)
{
      int i,j,k,ia,ja,nconn;
      double dmaxsq, dijsq, st, tmp1[3], tmp2[3];
      int iconn[30];

      if (!dowrt) return;

      if (*adjus != 1.0) {
         fprintf(out,"[MOLECULE] AU CONN\n");
      } else {
         fprintf(out,"[MOLECULE] CONN\n");
      }

      fprintf(out,"%f %f %f\n",r[0],r[1],r[2]);
      fprintf(out,"%d\n",*natoms);
      for (i=0; i<*natoms; i++) {

         ia = nat[i];
	 tmp1[0] = -1.0 * ysym[i] / r[0];
	 tmp1[1] = -1.0 * xsym[i] / r[0];
	 tmp1[2] = -1.0 * zsym[i] / r[0];


         nconn = 0;
         for (j=0; j<*natoms; j++) {

            ja = nat[j];
	    tmp2[0] = -1.0 * ysym[j] / r[0];
	    tmp2[1] = -1.0 * xsym[j] / r[0];
	    tmp2[2] = -1.0 * zsym[j] / r[0];

            dmaxsq = (vdwr[ia] + vdwr[ja]);
            dmaxsq = dmaxsq * dmaxsq;

            dijsq = 0.0;

            st = (xsym[i] - xsym[j])*(*adjus);
	    st *= st;
	    dijsq = dijsq + st;

	    st = (ysym[i] - ysym[j])*(*adjus);
	    st *= st;
	    dijsq = dijsq + st;

	    st = (zsym[i] - zsym[j])*(*adjus);
	    st *= st;
	    dijsq = dijsq + st;

            if (i != j && dijsq < dmaxsq) {
		iconn[nconn] = j+1;
		nconn++;
            }

         }
         fprintf(out,"%d %f %f %f %d",nat[i],xsym[i],ysym[i],zsym[i],nconn);
	 for (k=0; k<nconn; k++)
	   fprintf(out," %d",iconn[k]);
         fprintf(out,"\n");
      }

}


#if defined(VMS) || defined(UNDERSC)
void ogbeg(isurf)
#else
#ifdef CRAY
void OGBEG(isurf)
#else
void ogbeg_(isurf)
#endif
#endif
int *isurf;
{
      int surf,surfabs;

      if (!dowrt) return;

      surf = *isurf;
      surfabs = abs(surf);

      if (surfabs == 3) {
	   fprintf(out,"[SURFACE] MAPPED\n");
      } else {
	if (surf<0) {
	   fprintf(out,"[SURFACE] TRANS COLOR 1.0 1.0 0.0\n");
	} else {
	   fprintf(out,"[SURFACE]\n");
	}
      }

}

#if defined(VMS) || defined(UNDERSC)
void ogend()
#else
#ifdef CRAY
void OGEND()
#else
void ogend_()
#endif
#endif
{
      int idum;

      idum = 0;
}

double vln(double *a)
{
      double vl;
      double tot;

      tot = a[0]*a[0]+a[1]*a[1]+a[2]*a[2];

      vl = 0.0;
      if (tot > 0.0) vl = sqrt(tot);

      return(vl);

}

crpsin(a,b,c,d)
double *a;
double *b;
double *c;
double *d;
{
/*

       calculates cross product:  (b-a) x (c-a) = d
                                   ---     ---    -
*/
      int i;
      double v1[3],v2[3],dlen;

      for (i=0; i<3; i++) {
         v1[i] = b[i]-a[i];
         v2[i] = c[i]-a[i];
      }

      d[0] = v2[1]*v1[2]-v2[2]*v1[1];
      d[1] = v2[2]*v1[0]-v2[0]*v1[2];
      d[2] = v2[0]*v1[1]-v2[1]*v1[0];
      dlen = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);

      if (dlen > 0 ) for (i=0; i<3; i++) d[i] = d[i] / dlen;
}

void znorm(double rpts, double cnst, double *dens, double *vn,
	   int npts1, int npts2, int i, int j)
{
     double vl,z1[3],z2[3],z3[3],nx[3],ny[3];
     int k;

     if (j+1 < npts2) {
         z1[0] = 1.0; z1[1] = 0.0;
         z1[2] = cnst*rpts*(dens[j+1+i*npts2]-dens[j+i*npts2]);
         crpsin(z0,z1,zy,z2);
     } else {
         for (k=0; k<3; k++) z2[k] = z0[k];
     }

     if (j-1 >= 0) {
         z1[0] = -1.0; z1[1] = 0.0;
         z1[2] = cnst*rpts*(dens[j-1+i*npts2]-dens[j+i*npts2]);
         crpsin(z0,zy,z1,z3);
     } else {
         for (k=0; k<3; k++) z3[k] = z0[k];
     }

     for (k=0; k<3; k++) nx[k] = z3[k] + z2[k];
     vl = vln(nx);
     for (k=0; k<3; k++) nx[k] = nx[k] / vl;

     if (i+1 < npts1) {
         z1[0] = 0.0; z1[1] = 1.0;
         z1[2] = cnst*rpts*(dens[j+(i+1)*npts2]-dens[j+i*npts2]);
         crpsin(z0,z1,zx,z2);
     } else {
         for (k=0; k<3; k++) z2[k] = z0[k];
     }

     if (i-1 >= 0) {
         z1[0] = 0.0; z1[1] = -1.0;
         z1[2] = cnst*rpts*(dens[j+(i-1)*npts2]-dens[j+i*npts2]);
         crpsin(z0,zx,z1,z3);
     } else {
         for (k=0; k<3; k++) z3[k] = z0[k];
     }

     for (k=0; k<3; k++) ny[k] = z3[k] + z2[k];
     vl = vln(ny);
     for (k=0; k<3; k++) ny[k] = ny[k] / vl;

     for (k=0; k<3; k++) vn[k] = ny[k] - nx[k];
     vl = vln(vn);
     for (k=0; k<3; k++) vn[k] = vn[k] / vl;

}

BuildList(r,cnst,nnpts1,nnpts2,dens)
double *r;
double *cnst;
int *nnpts1;
int *nnpts2;
double *dens;
{

     int i,j,k,noff1,noff2,npts1,npts2;
     double v[3], rpts,vl;
     double vec1[3], vn1[3];
     float hinv1,hinv2;

/*
     a grid of n*n points has n-1*n-1 squares
     and twice as much triangular polygons
     so 2*(npts-1)**2
     npts is 80 at maximum so lets make it 12800
*/

      npts1 = *nnpts1;
      npts2 = *nnpts2;
      noff1 = npts1/2;
      noff2 = npts2/2;
      rpts = (double) (npts1-1);
      hinv1 = r[1]/ (npts2*r[0]);
      hinv2 = 1.0/ npts1;

      if (dowrt) fprintf(out,"[ELEVATIONGRID]\n");

      for (i=0; i<npts1-1; i++) {
         for (j=0; j<npts2-1; j++) {
/*
        first triangle
*/
            vec1[0] = (double) (j-noff2);
            vec1[1] = (double) (i-noff1);
            vec1[2] = dens[j+i*npts2]*(*cnst)*rpts;

	    znorm(rpts,*cnst,dens,vn1,npts1,npts2,i,j);

#if defined(VMS) || defined(UNDERSC)
            ognorm(&vn1[0],&vn1[1],&vn1[2]);
#else
#ifdef CRAY
            OGNORM(&vn1[0],&vn1[1],&vn1[2]);
#else
            ognorm_(&vn1[0],&vn1[1],&vn1[2]);
#endif
#endif

	    v[0] = vec1[0]*hinv1;
	    v[1] = vec1[1]*hinv2;
	    v[2] = vec1[2]*hinv1;
#if defined(VMS) || defined(UNDERSC)
	    ogvert(&v[0],&v[1],&v[2]);
#else
#ifdef CRAY
	    OGVERT(&v[0],&v[1],&v[2]);
#else
	    ogvert_(&v[0],&v[1],&v[2]);
#endif
#endif

            vec1[0] = (double) (j+1-noff2);
            vec1[1] = (double) (i-noff1);
            vec1[2] = dens[j+1+i*npts2]*(*cnst)*rpts;

	    znorm(rpts,*cnst,dens,vn1,npts1,npts2,i,j+1);

#if defined(VMS) || defined(UNDERSC)
            ognorm(&vn1[0],&vn1[1],&vn1[2]);
#else
#ifdef CRAY
            OGNORM(&vn1[0],&vn1[1],&vn1[2]);
#else
            ognorm_(&vn1[0],&vn1[1],&vn1[2]);
#endif
#endif

	    v[0] = vec1[0]*hinv1;
	    v[1] = vec1[1]*hinv2;
	    v[2] = vec1[2]*hinv1;
#if defined(VMS) || defined(UNDERSC)
	    ogvert(&v[0],&v[1],&v[2]);
#else
#ifdef CRAY
	    OGVERT(&v[0],&v[1],&v[2]);
#else
	    ogvert_(&v[0],&v[1],&v[2]);
#endif
#endif

            vec1[0] = (double) (j+1-noff2);
            vec1[1] = (double) (i+1-noff1);
            vec1[2] = dens[j+1+(i+1)*npts2]*(*cnst)*rpts;

	    znorm(rpts,*cnst,dens,vn1,npts1,npts2,i+1,j+1);

#if defined(VMS) || defined(UNDERSC)
            ognorm(&vn1[0],&vn1[1],&vn1[2]);
#else
#ifdef CRAY
            OGNORM(&vn1[0],&vn1[1],&vn1[2]);
#else
            ognorm_(&vn1[0],&vn1[1],&vn1[2]);
#endif
#endif

	    v[0] = vec1[0]*hinv1;
	    v[1] = vec1[1]*hinv2;
	    v[2] = vec1[2]*hinv1;
#if defined(VMS) || defined(UNDERSC)
	    ogvert(&v[0],&v[1],&v[2]);
#else
#ifdef CRAY
	    OGVERT(&v[0],&v[1],&v[2]);
#else
	    ogvert_(&v[0],&v[1],&v[2]);
#endif
#endif

            vec1[0] = (double) (j-noff2);
            vec1[1] = (double) (i+1-noff1);
            vec1[2] = dens[j+(i+1)*npts2]*(*cnst)*rpts;

	    znorm(rpts,*cnst,dens,vn1,npts1,npts2,i+1,j);

#if defined(VMS) || defined(UNDERSC)
            ognorm(&vn1[0],&vn1[1],&vn1[2]);
#else
#ifdef CRAY
            OGNORM(&vn1[0],&vn1[1],&vn1[2]);
#else
            ognorm_(&vn1[0],&vn1[1],&vn1[2]);
#endif
#endif

	    v[0] = vec1[0]*hinv1;
	    v[1] = vec1[1]*hinv2;
	    v[2] = vec1[2]*hinv1;
#if defined(VMS) || defined(UNDERSC)
	    ogvert(&v[0],&v[1],&v[2]);
#else
#ifdef CRAY
	    OGVERT(&v[0],&v[1],&v[2]);
#else
	    ogvert_(&v[0],&v[1],&v[2]);
#endif
#endif

         }
      }


}

void launchViewer(opt)
int opt;
{
    int pid;

    pid = fork();
    switch(pid) {
    case -1:
          fprintf(stderr,"Couldnt Fork\n");
          break;
    case 0:       /*child */
#ifdef __CYGWIN__
	  if (opt) {
             execlp("./mogl","mogl","-r","-b","0",ogfil,NULL);
	  } else {
             execlp("./mogl","mogl",ogfil,NULL);
	  }
#else
	  if (opt) {
             execlp("moldenogl","moldenogl","-r","-b","0",ogfil,NULL);
	  } else {
             execlp("moldenogl","moldenogl",ogfil,NULL);
	  }
#endif
	  exit(0);
          break;
    default:      /*parent*/
          break;
    }
}

#if defined(VMS) || defined(UNDERSC)
void oginid(double *r, double *adjus, int *natoms, int *nat, int *icol,
#else
#ifdef CRAY
void OGINID(double *r, double *adjus, int *natoms, int *nat, int *icol,
#else
void oginid_(double *r, double *adjus, int *natoms, int *nat, int *icol,
#endif
#endif
       double *xsym, double *ysym, double *zsym, double *vdwr,
       double *cnst, int *nnpts1, int *nnpts2, double *dens)
{

    dowrt = 1;
    out = fopen(ogfil,"w");

    if (out == NULL) {
          fprintf(stderr,"Unable to open file %s\n",ogfil);
	  dowrt = 0;
	  return;
    }

    fprintf(out,"[MOLDENOGL]\n");
    ogwrmol(r,adjus,natoms,nat,xsym,ysym,zsym,vdwr);
    BuildList(r,cnst,nnpts1,nnpts2,dens);

    if (dowrt) fclose(out);

    launchViewer(1);

}

#if defined(VMS) || defined(UNDERSC)
void oginsp(double *r, double *adjus, int *natoms, int *nat, int *iatclr, int *icol,
#else
#ifdef CRAY
void OGINSP(double *r, double *adjus, int *natoms, int *nat, int *iatclr, int *icol,
#else
void oginsp_(double *r, double *adjus, int *natoms, int *nat, int *iatclr, int *icol,
#endif
#endif

       double *xsym, double *ysym, double *zsym, double *vdwr, 
       int *mopt, int *conn, int *nconn, int *iconn, int *ofil)
{

    dowrt = 1;
    out = fopen(ogfil,"w");

    if (out == NULL) {
          fprintf(stderr,"Unable to open file %s\n",ogfil);
	  dowrt = 0;
	  return;
    }

    fprintf(out,"[MOLDENOGL]\n");
    ogwrmol(r,adjus,natoms,nat,xsym,ysym,zsym,vdwr);

}

#if defined(VMS) || defined(UNDERSC)
void initog(int *iwrt)
#else
#ifdef CRAY
void INITOG(int *iwrt)
#else
void initog_(int *iwrt)
#endif
#endif
{
 int i;

/* just dummy */

 i = *iwrt;

}

#if defined(VMS) || defined(UNDERSC)
void ogspst()
#else
#ifdef CRAY
void OGSPST()
#else
void ogspst_()
#endif
#endif
{
    int pid;
    int stat;

    if (dowrt) fclose(out);

    launchViewer(1);
}

#if defined(VMS) || defined(UNDERSC)
void viewer()
#else
#ifdef CRAY
void VIEWER()
#else
void viewer_()
#endif
#endif
{
    launchViewer(0);
}
