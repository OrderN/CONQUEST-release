/*

     ##  COPYRIGHT (C) 2000 by G. Schaftenaar  ##

     MOLecular DENsity OpenGL Helper Program
     By G. Schaftenaar
     CMBI, University of Nijmegen, The Netherlands
     (formerly CAOS/CAMM Center)
     2000
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef WIN32
#include <GL/gl.h>
#endif
#ifdef DARWIN
/* MacOS X "Panther" OpenGL implementation -- Uses Xcode Tools */
#include <OpenGL/glext.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define NUMATM 20000
#define MXEL 100
#define MXCON 12

static double xsym[NUMATM], ysym[NUMATM], zsym[NUMATM];
static int nat[NUMATM], iatclr[NUMATM], nconn[NUMATM], iconn[NUMATM][MXCON];

#define TRDEF 0.6
extern int ecol = 2;
extern int bgcol = 6;
extern int perspon = 1;
extern double tr_val = TRDEF;
extern int mapped = 0;
extern float diffuseColor[15][4];
extern int hires = 0;
extern GLfloat poszz = 0.0;
extern int AtomColors[15][3];
extern int lines = 0;
extern int spacefill = 0;
extern int atcol = 1;
extern int width = 600;
extern int height = 600;
extern int RotType = 0;

static float ribcol[6][4] =
{
  {0.0,0.0,1.0,0.8},
  {1.0,0.0,1.0,0.8},
  {0.0,1.0,0.0,0.8},
  {0.6,0.6,0.6,0.8},
  {0.5,1.0,0.5,0.8},
  {1.0,1.0,1.0,0.8}
};
 
static FILE *out;

static int icol[] = {
15, 4, 8, 8, 9, 14, 7, 1, 2, 12, 8, 8, 7, 11, 2, 6, 3, 12, 2, 1, 2, 2,
2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 12, 2, 1, 2, 2, 2, 2, 2, 2, 2,
2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2,
2, 2, 2, 2, 2, 1, 13 };

/* vdwr in a.u. */

static double vdwr[] = {
0.430, 0.741, 0.880, 0.550, 1.030, 0.900, 0.880, 0.880, 0.840, 0.815,
1.170, 1.300, 1.550, 1.400, 1.250, 1.220, 1.190, 0.995, 1.530, 1.190,
1.640, 1.670, 1.530, 1.550, 1.550, 1.540, 1.530, 1.700, 1.720, 1.650,
1.420, 1.370, 1.410, 1.420, 1.410, 1.069, 1.670, 1.320, 1.980, 1.760,
1.680, 1.670, 1.550, 1.600, 1.650, 1.700, 1.790, 1.890, 1.830, 1.660,
1.660, 1.670, 1.600, 1.750, 1.870, 1.540, 2.070, 2.030, 2.020, 2.010,
2.000, 2.000, 2.190, 1.990, 1.960, 1.950, 1.940, 1.930, 1.920, 2.140,
1.920, 1.770, 1.630, 1.570, 1.550, 1.570, 1.520, 1.700, 1.700, 1.900,
1.750, 1.740, 1.740, 1.880, 0.200, 0.200, 0.200, 2.100, 2.080, 1.990,
1.810, 1.780, 1.750, 0.200, 1.710, 0.200, 0.200, 1.730, 0.100, 0.200 };

ogrdmol(char *top, double *r, double *adjus, int *natoms, int *nat, int *iatclr,
        double *xsym, double *ysym, double *zsym, int *iopt, int *conn, int *nconn, int iconn[NUMATM][MXCON])
{
      int i;

      *adjus = 1.0;
      if (strstr(top,"AU")) *adjus = 0.529177;
      *conn = 0;
      if (strstr(top,"CONN")) *conn = 1;
      *iopt = 0;
      if (strstr(top,"UNSCALED")) *iopt = 1;
      spacefill = 0;
      if (strstr(top,"SPACEFILL")) spacefill = 1;
      atcol = 1;
      if (strstr(top,"GRPCOL")) atcol = 0;
      if (strstr(top,"HIGH")) hires = 1;
      
      if (*iopt) {
        r[0] = 1.0;
      } else {
	fgets(top,100,out);
	sscanf(top,"%lf %lf %lf",&r[0],&r[1],&r[2]);
      }
      fgets(top,100,out);
      sscanf(top,"%d",natoms);
      for (i=0; i<*natoms; i++) {
	  fgets(top,100,out);
	  if (atcol) {
	    sscanf(top,"%d %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %d",
		&nat[i],&xsym[i],&ysym[i],&zsym[i],&nconn[i],
		&iconn[i][0],&iconn[i][1],&iconn[i][2],&iconn[i][3],
		&iconn[i][4],&iconn[i][5],&iconn[i][6],&iconn[i][7],
		&iconn[i][8],&iconn[i][9], &iconn[i][10],&iconn[i][11]);
	  } else {
	    sscanf(top,"%d %d %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %d",
		&nat[i],&iatclr[i],&xsym[i],&ysym[i],&zsym[i],&nconn[i],
		&iconn[i][0],&iconn[i][1],&iconn[i][2],&iconn[i][3],
		&iconn[i][4],&iconn[i][5],&iconn[i][6],&iconn[i][7],
		&iconn[i][8],&iconn[i][9], &iconn[i][10],&iconn[i][11]);
	  }
      }
}

dist(double *coo, double *dsq)
{
  double d;

  d = coo[0]*coo[0] + coo[1]*coo[1] + coo[2]*coo[2]; 
  if (d > 0.0) d = sqrt(d);
  if (d > *dsq) *dsq = d;
}

ogsize(int *natoms, double *xsym, double *ysym, double *zsym, double *r)
{
   int i;
   double dijsq;

   for (i=0; i<*natoms; i++) {
	dijsq = xsym[i]*xsym[i] + ysym[i]*ysym[i] + zsym[i]*zsym[i];
	dijsq = dijsq / (r[0]*r[0]);
	if (dijsq > 0.0) dijsq = sqrt(dijsq);
	if (dijsq > poszz) poszz = dijsq;

   }
}

int main(int argc, char **argv)
{
   double adjus, r[3];
   int natoms, conn, mopt;
   int iopt,fend,i,ic,elev,ribb,surf,loop,iribc;
   char rdstr[100];
   char infile[100];
   char *colstr;
   double v[3],vt[3],vc[3],dis;

   strcpy(infile,"molden.ogl");

   hires = 0;
   perspon = 1;
   dis = 0.0;

   i = 1;
   while (i < argc) {
      if (strcmp(argv[i], "-h") == 0) {
#ifdef WIN32
         fprintf(stderr,
	"Usage: mogl [-b 0-7] [-c 0-7] [-t 0.75] [-o] [-r] [-h] [-W width] [-H heigth] filename\n");
#else
         fprintf(stderr,
	"Usage: moldenogl [-b 0-7] [-c 0-7] [-t 0.75] [-o] [-r] [-h] [-W width] [-H heigth] filename\n");
#endif
      } else if (strcmp(argv[i], "-c") == 0) {
	 i++;
         if (i<argc) ecol = atoi(argv[i]);
      } else if (strcmp(argv[i], "-t") == 0) {
	 i++;
         if (i<argc) tr_val = atof(argv[i]);
      } else if (strcmp(argv[i], "-l") == 0) {
	 hires = 0;
      } else if (strcmp(argv[i], "-r") == 0) {
	 hires = 1;
      } else if (strcmp(argv[i], "-b") == 0) {
	 i++;
         if (i<argc) bgcol = atoi(argv[i]);
      } else if (strcmp(argv[i], "-p") == 0) {
	 perspon = 1;
      } else if (strcmp(argv[i], "-o") == 0) {
	 perspon = 0;
      } else if (strcmp(argv[i], "-W") == 0) {
	 i++;
	 width = atoi(argv[i]);
      } else if (strcmp(argv[i], "-H") == 0) {
	 i++;
	 height = atoi(argv[i]);
      } else if (strcmp(argv[i], "-R") == 0) {
	 RotType = 1;
      } else {
         strcpy(infile,argv[i]);
      }
      i++;
   }

   out = fopen(infile,"r");
   if (out == NULL) {
      fprintf(stderr,"Couldnt open file %s\n",infile);
      exit(1);
   }

   fgets(rdstr,100,out);
   if (strncmp(rdstr,"[MOLDENOGL]",11) != 0) {
      fprintf(stderr,"This is not an [MOLDENOGL] file !\n");
      exit(1);
   }

   iopt = 0;
#if defined(VMS) || defined(UNDERSC)
      initog(&iopt);
#else
#ifdef CRAY
      INITOG(&iopt);
#else
      initog_(&iopt);
#endif
#endif

   while(!(fend = feof(out))) {
      fgets(rdstr,100,out);
      if (strstr(rdstr,"[MOLECULE]") == NULL) goto Surf;
      ogrdmol(rdstr,r,&adjus,&natoms,nat,iatclr,xsym,ysym,zsym,&mopt,&conn,nconn,iconn);
      ogsize(&natoms,xsym,ysym,zsym,r);

#if defined(VMS) || defined(UNDERSC)
      oginsp(r,&adjus,&natoms,nat,iatclr,icol,xsym,ysym,zsym,vdwr,&mopt,&conn,nconn,iconn,&iopt);
#else
#ifdef CRAY
      OGINSP(r,&adjus,&natoms,nat,iatclr,icol,xsym,ysym,zsym,vdwr,&mopt,&conn,nconn,iconn,&iopt);
#else
      oginsp_(r,&adjus,&natoms,nat,iatclr,icol,xsym,ysym,zsym,vdwr,&mopt,&conn,nconn,iconn,&iopt);
#endif
#endif
   }


Surf:

   loop = 0;

   elev = 0;
   ribb = 0;
   lines = 0;
   surf = 0;
   mapped = 0;

   while(!(fend = feof(out))) {

      if (loop) fgets(rdstr,100,out);
      loop = 1;
      
      if (strncmp(rdstr,"[LINES",6) == 0) {
#if defined(VMS) || defined(UNDERSC)
	if (surf) ogend();
	oglin();
#else
#ifdef CRAY
	if (surf) OGEND();
	OGLIN();
#else
	if (surf) ogend_();
	oglin_();
#endif
#endif
	surf++;
	lines = 1;
	continue;
      }

      iribc = -1;
      if (strstr(rdstr,"[COL STRANDTOP]") != NULL) iribc = 0;
      if (strstr(rdstr,"[COL STRANDBOTTOM]") != NULL) iribc = 1;
      if (strstr(rdstr,"[COL HELIXOUT]") != NULL) iribc = 2;
      if (strstr(rdstr,"[COL HELIXIN]") != NULL) iribc = 3;
      if (strstr(rdstr,"[COL RNA]") != NULL) iribc = 4;
      if (strstr(rdstr,"[COL COIL]") != NULL) iribc = 5;
      if (iribc != -1) {
	colstr = strstr(rdstr,"]");
	if (colstr != NULL) {
	   colstr = colstr + 1;
	   sscanf(colstr,"%f %f %f",&ribcol[iribc][0],
		&ribcol[iribc][1],&ribcol[iribc][2]);
	}
	continue;
      }

      if (strncmp(rdstr,"[SURFACE",8) == 0 || 
	  strncmp(rdstr,"[ELEVATION",10) == 0 || 
	  strncmp(rdstr,"[RIBBON",7) == 0) {
	lines = 0;
#if defined(VMS) || defined(UNDERSC)
	if (surf) ogend();
#else
#ifdef CRAY
	if (surf) OGEND();
#else
	if (surf) ogend_();
#endif
#endif
	colstr = strstr(rdstr,"COLOR");
	if (colstr != NULL) {
	   colstr = colstr + 5;
	   sscanf(colstr,"%f %f %f",&diffuseColor[surf][0],
		&diffuseColor[surf][1],&diffuseColor[surf][2]);
	}

	if (strncmp(rdstr,"[RIBBON",7) == 0) {
           iribc = -1;
	   if (strstr(rdstr,"STRANDTOP") != NULL) iribc = 0;
	   if (strstr(rdstr,"STRANDBOTTOM") != NULL) iribc = 1;
	   if (strstr(rdstr,"HELIXOUT") != NULL) iribc = 2;
	   if (strstr(rdstr,"HELIXIN") != NULL) iribc = 3;
	   if (strstr(rdstr,"RNA") != NULL) iribc = 4;
	   if (strstr(rdstr,"COIL") != NULL) iribc = 5;

	   if (iribc != -1)
	      for (i=0; i<3; i++) diffuseColor[surf][i] = ribcol[iribc][i];
	}


	if (strstr(rdstr,"MAPPED") != NULL) mapped = 1;

	iopt = 1;
	colstr = strstr(rdstr,"TRANS");
	if (colstr != NULL) {
	   colstr = colstr + 5;
	   iopt = -iopt;
	   if (sscanf(colstr,"%lf",&tr_val) <= 0) {
		tr_val = TRDEF;
	   }
	}

	if (strncmp(rdstr,"[ELEVATION",10) == 0 ||
	    strncmp(rdstr,"[RIBBON",7) == 0) {
	   if (strncmp(rdstr,"[ELEVATION",10) == 0) {
		elev = 1;
		ogelev();
	   } else {
		ribb = 1;
		ogribb();
	   }
	} else {
#if defined(VMS) || defined(UNDERSC)
	   ogbeg(&iopt);
#else
#ifdef CRAY
	   OGBEG(&iopt);
#else
	   ogbeg_(&iopt);
#endif
#endif
	}
	surf++;
	fgets(rdstr,100,out);
      }

      if (lines) {

/*
Yes you are right the x and y coordinates are swapped in molden
Have to fix it sometime
*/
	sscanf(rdstr,"%d %lf %lf %lf %lf %lf %lf",
		&ic,&v[0],&v[1],&v[2],&vt[0],&vt[1],&vt[2]);

	for (i=0; i<3; i++) {
	   vc[i] = ((GLdouble) AtomColors[ic-1][i]) / 255.0 ;
	}
#if defined(VMS) || defined(UNDERSC)
	ogcol(&vc[0],&vc[1],&vc[2]);
	ogvert(&v[0],&v[1],&v[2]);
	ogvert(&vt[0],&vt[1],&vt[2]);
#else
#ifdef CRAY
	OGCOL(&vc[0],&vc[1],&vc[2]);
	OGVERT(&v[0],&v[1],&v[2]);
	OGVERT(&vt[0],&vt[1],&vt[2]);
#else
	ogcol_(&vc[0],&vc[1],&vc[2]);
	ogvert_(&v[0],&v[1],&v[2]);
	ogvert_(&vt[0],&vt[1],&vt[2]);
#endif
#endif
	continue;

      }

      if (mapped) {
	sscanf(rdstr,"%lf %lf %lf",&v[0],&v[1],&v[2]);
	v[0] = fabs(v[0]);
#if defined(VMS) || defined(UNDERSC)
	ogcol(&v[0],&v[1],&v[2]);
#else
#ifdef CRAY
	OGCOL(&v[0],&v[1],&v[2]);
#else
	ogcol_(&v[0],&v[1],&v[2]);
#endif
#endif
        fgets(rdstr,100,out);
      }

      sscanf(rdstr,"%lf %lf %lf",&v[0],&v[1],&v[2]);
#if defined(VMS) || defined(UNDERSC)
      ognorm(&v[0],&v[1],&v[2]);
#else
#ifdef CRAY
      OGNORM(&v[0],&v[1],&v[2]);
#else
      ognorm_(&v[0],&v[1],&v[2]);
#endif
#endif

      fgets(rdstr,100,out);
      sscanf(rdstr,"%lf %lf %lf",&v[0],&v[1],&v[2]);
#if defined(VMS) || defined(UNDERSC)
      ogvert(&v[0],&v[1],&v[2]);
#else
#ifdef CRAY
      OGVERT(&v[0],&v[1],&v[2]);
#else
      ogvert_(&v[0],&v[1],&v[2]);
#endif
#endif
      dist(v,&dis);

      if (elev || ribb) {

       for (i=0; i<3; i++) {

         fgets(rdstr,100,out);
         sscanf(rdstr,"%lf %lf %lf",&v[0],&v[1],&v[2]);
#if defined(VMS) || defined(UNDERSC)
         ognorm(&v[0],&v[1],&v[2]);
#else
#ifdef CRAY
         OGNORM(&v[0],&v[1],&v[2]);
#else
         ognorm_(&v[0],&v[1],&v[2]);
#endif
#endif
         fgets(rdstr,100,out);
         sscanf(rdstr,"%lf %lf %lf",&v[0],&v[1],&v[2]);
#if defined(VMS) || defined(UNDERSC)
         ogvert(&v[0],&v[1],&v[2]);
#else
#ifdef CRAY
         OGVERT(&v[0],&v[1],&v[2]);
#else
         ogvert_(&v[0],&v[1],&v[2]);
#endif
#endif
         dist(v,&dis);
        }

      }
   }

   if (dis > poszz) poszz = dis;

#if defined(VMS) || defined(UNDERSC)
   ogend();
#else
#ifdef CRAY
   OGEND();
#else
   ogend_();
#endif
#endif

   fclose(out);

#if defined(VMS) || defined(UNDERSC)
   ogspst();
#else
#ifdef CRAY
   OGSPST();
#else
   ogspst_();
#endif
#endif


}
