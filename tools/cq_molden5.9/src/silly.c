
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

#ifndef WIN32
#include <GL/gl.h>
#endif
#ifdef DARWIN
/*
      MacOS X "Panther" OpenGL implemantation -- Uses Xcode Tools
*/
#include <OpenGL/glext.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

extern int ecol;
extern int bgcol;
extern double tr_val;
extern int mapped;
extern int hires;
extern int perspon;
extern GLfloat poszz;
extern int lines;
extern int lines;
extern int spacefill;
extern int atcol;
extern int width,height;
extern int RotType;
static int TRANS;
static int PERSP = 0;
static int l1on = 0;
static int l2on = 0;
static int l3on = 1;
static int stacks = 20;
static int DoLines = 0;
static int DoCap = 0;
static int DoFog = 0;
static GLfloat fdens = 0.03;
static int DoLights = 0;
static int ColBG = -1;

/* Some <math.h> files do not define M_PI... */
#ifndef M_PI
#define M_PI 3.141592654
#endif
#define TORAD M_PI/180.0

static GLdouble RR[4][4] = { 
{1.0,0.0,0.0,0.0}, {0.0,1.0,0.0,0.0}, {0.0,0.0,1.0,0.0}, {0.0,0.0,0.0,1.0} 
};

#define MAXFNAME 512
static int pixtyp = 0;
static char *pixext[] = {"rgb","ppm"};
#ifdef WIN32
static char *rgbfile = "mogl.rgb";
static char *ppmfile = "mogl.ppm";
#else
static char *rgbfile = "moldenogl.rgb";
static char *ppmfile = "moldenogl.ppm";
#endif
static char pixtmp[MAXFNAME];
#define PICMAX 180
static int picnum = 0;

static double z0[] = {0.0,0.0,0.0};
static double zx[] = {1.0,0.0,0.0};
static double zy[] = {0.0,1.0,0.0};

#define MAXSURF 500
extern float diffuseColor[MAXSURF][4] =
{ 
  {0.15,0.15,1.0,0.8}, 
  {1.0 ,0.1 ,0.1,0.8}, 
  {1.0,1.0,0.0,0.8},
  {1.0,0.5,0.5,0.8},
  {1.0,0.5,0.5,0.8},
  {1.0,0.5,0.5,0.8},
  {1.0,0.5,0.5,0.8},
  {1.0,0.5,0.5,0.8},
  {1.0,0.5,0.5,0.8},
  {1.0,0.5,0.5,0.8},
  {1.0,0.5,0.5,0.8},
  {1.0,0.5,0.5,0.8},
  {1.0,0.5,0.5,0.8},
  {1.0,0.5,0.5,0.8},
  {1.0,0.5,0.5,0.8}
};

static float ambientColor[4] = {0.0,0.0,0.0,0.0};
static float ambientFColor[4] = {1.0,1.0,1.0,1.0};
static float specularColor[MAXSURF][4] =
         { {0.8,0.8,0.8,1.0} , {0.8,0.8,0.8,1.0}, {0.9,0.8,0.8,1.0},
           {0.8,0.8,0.8,1.0} , {0.8,0.8,0.8,1.0}, {0.9,0.8,0.8,1.0},
           {0.8,0.8,0.8,1.0} , {0.8,0.8,0.8,1.0}
};

static float materialColor[8][4] =
{
  {0.8, 0.8, 0.8, 1.0},
  {0.8, 0.0, 0.0, 1.0},
  {0.0, 0.8, 0.0, 1.0},
  {0.0, 0.0, 0.8, 1.0},
  {0.0, 0.8, 0.8, 1.0},
  {0.8, 0.0, 0.8, 1.0},
  {0.8, 0.8, 0.0, 1.0},
  {1.0, 0.5, 0.5, 1.0},
};

static float clrColor[7][4] =
{
  {0.0, 0.0, 0.0, 0.0},
  {1.0, 1.0, 1.0, 1.0},
  {1.0, 0.0, 0.0, 1.0},
  {0.0, 1.0, 0.0, 1.0},
  {0.0, 1.0, 1.0, 1.0},
  {1.0, 0.5, 0.5, 0.0},
  {0.5, 0.5, 0.5, 0.0},
};

static char *clrNames[7] =
{ "Black", "White", "Red", "Green", "Blue", "Peach", "Grey"};

static Nclr = 7;

static float AmbientNul[4] = {0.0,0.0,0.0,0.0};

static float light0_position[] = {-2.0, -2.0, 0.0, 1};
static float light1_position[] = {4.0, 4.0, 4.0, 1};
static float light2_position[] = {0.0, 0.0, 4.0, 1};

#define MXEL 100
#define MXMOL 1000
static FILE *out;
static int dowrt = 0;
static int molon = 1;
static int AmbAndDiff = 1;
static int moving = 0;
static int startx, starty;
static int AutoRot = 0;
GLfloat ang1 = 0.0;
GLfloat ang2 = 0.0;
double angincr1 = 0.0;
double angincr2 = 0.0;
#define HISTMAX 3
#define MULTH 10.0
static double AngIncrHist1[HISTMAX];
static double AngIncrHist2[HISTMAX];
static int mmoving = 0;
static int mstartx, mstarty;
GLfloat posx = 0.0;
GLfloat posy = 0.0;
GLfloat posz = 2.0;
static int smoving = 0;
static int sstarty;

static int NSurf = -1;
GLuint theSurf[MAXSURF];
static int NMols = -1;
static int CMols = -1;
static int Anim = 1;
static int AnimSlow = 0;
static int AnimCount = 1;
GLuint theMol[MXMOL];
static GLUquadricObj *cyl;
static GLUquadricObj *sphere;

#define ROTINCR 5.0
static GLfloat xrot = -45.0;
static GLfloat yrot = 240.0;
static int Vsize = 600;

extern int AtomColors[15][3] =
{
      {255,0,0},
      {255,159,9},
      {0,255,0},
      {78,255,187},
      {0,255,255},
      {255,191,0},
      {132,193,214},
      {115,115,115},
      {255,0,255},
      {16,176,16},
      {239,202,140},
      {255,122,0},
      {230,214,92},
      {184,56,6},
      {255,255,255},
};

/* How many feedback buffer GLfloats each of the three objects need. */
int objectComplexity = 3800000;  /* Teapot requires ~1.5 megabytes for
                           its feedback results! */

void idle(void)
{

  if (AutoRot) {
     ang1 = ang1 + 2.0;
     angincr1 = 2.0;
     glutPostRedisplay();
  }

  AnimCount--;
  if (AnimCount < 0) AnimCount = AnimSlow;
  if (!AnimCount) {
     if (molon && Anim) {
	CMols++;
	if (CMols >= NMols) CMols = 0;
     }
  }
  if (NMols > 1) glutPostRedisplay();

  if (DoCap) {
     picnum++;
     if (picnum < PICMAX) {
	sprintf(pixtmp, "mol%03d.%s",picnum,pixext[pixtyp]);
	switch(pixtyp) {
	case 0:
	  save_rgb(pixtmp);
	  break;
	case 1:
	  save_pixmap(pixtmp);
	  break;
	}

     }
  }

}

static void
setColor(int c)
{
      GLfloat ambientColor[4] = {0.0, 0.0, 0.0, 0.0};
      GLfloat mat_specular[4] = {0.8, 0.8, 0.8, 1.0 };

      glMaterialfv(GL_FRONT_AND_BACK,
        GL_AMBIENT, ambientColor);
      glMaterialfv(GL_FRONT_AND_BACK,
        GL_DIFFUSE, &materialColor[c][0]);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
      glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 130);
}


void cross(float *a,float *b, float *c)
{
/*

       calculates cross product:   a x b = c
                                   -   -   -
*/
 
      c[0] = a[1]*b[2] - a[2]*b[1];
      c[1] = a[2]*b[0] - a[0]*b[2];
      c[2] = a[0]*b[1] - a[1]*b[0];

}

float veclen(float *a)
{
      float vl;
      double tot;

      tot = a[0]*a[0]+a[1]*a[1]+a[2]*a[2];

      vl = 0.0;
      if (tot > 0.0) vl = (float) sqrt(tot);

      return(vl);

}

void improd(float *a, float *b, double *c)
{
      int i;
      float rimp, al, bl;

      rimp = 0.0;

      for (i=0; i<3; i++) rimp = rimp + a[i]*b[i];

      al = veclen(a);
      bl = veclen(b);

      if (al > 0.0 && bl > 0.0) {
         *c = (double) rimp / (veclen(a)*veclen(b));
      } else {
         *c = 0.0;
      }

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
      glNormal3d(*v1,*v2,*v3);
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
      glVertex3d(*v1,*v2,*v3);
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
      glColor4d(*v1,*v2,*v3,tr_val);
      if (dowrt) fprintf(out,"%f %f %f\n",*v1,*v2,*v3);
}

extern void
setAtomColor(int c)
{
      GLfloat tmp_Color[4];
      int i;
      GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };

      for (i=0; i<3; i++) 
	 tmp_Color[i] = ((GLfloat) AtomColors[c][i]) / 255.0 ; 
      tmp_Color[3] = 0.5;

      if (AmbAndDiff) {
         glMaterialfv(GL_FRONT_AND_BACK,
           GL_AMBIENT_AND_DIFFUSE, tmp_Color);
      } else {
         glMaterialfv(GL_FRONT_AND_BACK,
           GL_AMBIENT, AmbientNul);
         glMaterialfv(GL_FRONT_AND_BACK,
           GL_DIFFUSE, tmp_Color);
      }
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
      glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100);

}


void ogrod(int ic, float *p1, float *p2, double rad)
{
      int i;
      double cosa;
      float angle, vl;
      float v1[3], v2[3], v3[3], p3[3];
      float todeg;

      
      todeg = 45.0 / atan(1.0);

      for (i=0; i<3; i++) {
         v1[i] = p2[i] - p1[i];
         v2[i] = 0.0;
         p3[i] = v1[i]/2.0 + p1[i];
      }
      v2[2] = 1.0;

      improd(v1,v2,&cosa);

      if (abs(cosa) == 1.0) {

         for (i=0; i<3; i++) v3[i] = 0.0;
         v3[1] = 1.0;
         angle = 0.0;
	 if (cosa < 0) angle = 180.0;

      } else {

         angle = (float) acos(cosa)*todeg;
         cross(v2,v1,v3);
         vl = veclen(v3);
         for (i=0; i<3; i++) v3[i] = v3[i] / vl;

      }


      glPushMatrix();
      glTranslatef(p1[0],p1[1],p1[2]);
      glRotatef(angle,v3[0],v3[1],v3[2]);
      setAtomColor(ic);
      vl = veclen(v1);
      if (hires) {
	 gluCylinder( cyl, rad, rad, vl, stacks, 10 );
      } else {
	 gluCylinder( cyl, rad, rad, vl, 15, 1 );
      }
      glPopMatrix();

}

void ogsphere(int ic, float *p1, double rad)
{

      glPushMatrix();
      glTranslatef(p1[0],p1[1],p1[2]);
      setAtomColor(ic);
      if (hires) {
	 gluSphere(sphere, rad, stacks, stacks);
      } else {
	 gluSphere(sphere, rad, 10, 10);
      }
      glPopMatrix();

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

#define MXCON 12

ogmol(double *r, double *adjus, int *natoms, int *nat, int* iatclr, int *icol,
      double *xsym, double *ysym, double *zsym, double *vdwr,
      int *mopt, int *conn, int *nconn, int *iconn)
{
    double roddef, dmaxsq, dijsq, st;
    float tmp1[3],tmp2[3], tmp3[3];
    int i,j,k,l, ia, ic, ja, ido, doclr ;

    doclr = 0;
    if (iatclr != NULL && !atcol) doclr = 1;

    if (dowrt) ogwrmol(r,adjus,natoms,nat,xsym,ysym,zsym,vdwr);

    roddef = (0.13/0.52917706)/r[0];

    if (*adjus == 1.0) roddef = roddef*0.52917706;

    glPopMatrix();
    glPushMatrix();

    if (NMols >= MXMOL - 1) return;
    if (*natoms == 0) return;

    NMols++;
    theMol[NMols] = glGenLists(1);
    glNewList(theMol[NMols], GL_COMPILE);

    if (*conn) {

      for (i=0; i<*natoms; i++) {

         ia = nat[i];
	 if (doclr) {
            ic = iatclr[i]-1;
	 } else {
            ic = icol[ia-1]-1;
	 }
	 if (*mopt) {
            tmp1[0] = (float) xsym[i];
            tmp1[1] = (float) ysym[i];
            tmp1[2] = (float) zsym[i];
	 } else {
            tmp1[0] = (float) -1.0 * ysym[i] / r[0];
            tmp1[1] = (float) -1.0 * xsym[i] / r[0];
            tmp1[2] = (float) -1.0 * zsym[i] / r[0];
	 }


	 if (spacefill) {

	   roddef = vdwr[ia]*1.4/(r[0]*0.52917706);
	   if (*adjus == 1.0) roddef = roddef*0.52917706;
           ogsphere(ic,tmp1,roddef);

	 } else {

           ogsphere(ic,tmp1,roddef);

           for (k=0; k<nconn[i]; k++) {

            j = iconn[k+i*MXCON] - 1;
            ja = nat[j];
	    if (*mopt) {
		tmp2[0] = (float) xsym[j];
      		tmp2[1] = (float) ysym[j];
     		tmp2[2] = (float) zsym[j];
	    } else {
		tmp2[0] = (float) -1.0 * ysym[j] / r[0];
		tmp2[1] = (float) -1.0 * xsym[j] / r[0];
		tmp2[2] = (float) -1.0 * zsym[j] / r[0];
	    }

            if (1) {

                 ido = 1;
                 if (ja == ia) {
                    if (j > i) {
                       for (l=0; l<3; l++) tmp3[l] = tmp2[l];
                    } else {
                       ido = 0;
                    }
                 } else {
                    for (l=0; l<3; l++) 
                       tmp3[l] = (tmp2[l] - tmp1[l])/2.0 + tmp1[l];
                 }
  
                 if (ido) ogrod(ic,tmp1,tmp3,roddef);

            }
	   }
	 }
      }

    } else {

      for (i=0; i<*natoms; i++) {

         ia = nat[i];
	 if (doclr) {
            ic = iatclr[i]-1;
	 } else {
            ic = icol[ia-1]-1;
	 }
	 if (*mopt) {
            tmp1[0] = (float) xsym[i];
            tmp1[1] = (float) ysym[i];
            tmp1[2] = (float) zsym[i];
	 } else {
	    tmp1[0] = (float) -1.0 * ysym[i] / r[0];
	    tmp1[1] = (float) -1.0 * xsym[i] / r[0];
	    tmp1[2] = (float) -1.0 * zsym[i] / r[0];
	 }


	 if (spacefill) {

	   roddef = vdwr[ia]*1.4/(r[0]*0.52917706);
	   if (*adjus == 1.0) roddef = roddef*0.52917706;
           ogsphere(ic,tmp1,roddef);

	 } else {

           ogsphere(ic,tmp1,roddef);

           for (j=0; j<*natoms; j++) {

            ja = nat[j];
	    if (*mopt) {
		tmp2[0] = (float) xsym[j];
		tmp2[1] = (float) ysym[j];
		tmp2[2] = (float) zsym[j];
	    } else {
		tmp2[0] = (float) -1.0 * ysym[j] / r[0];
		tmp2[1] = (float) -1.0 * xsym[j] / r[0];
		tmp2[2] = (float) -1.0 * zsym[j] / r[0];
	    }

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

            if (dijsq < dmaxsq) {

                 ido = 1;
                 if (ja == ia) {
                    if (j > i) {
                       for (l=0; l<3; l++) tmp3[l] = tmp2[l];
                    } else {
                       ido = 0;
                    }
                 } else {
                    for (l=0; l<3; l++) 
                       tmp3[l] = (tmp2[l] - tmp1[l])/2.0 + tmp1[l];
                 }
  
                 if (ido) ogrod(ic,tmp1,tmp3,roddef);

            }

           }
         }
      }
    }
    glEndList();

}

#if defined(VMS) || defined(UNDERSC)
void oglin()
#else
#ifdef CRAY
void OGLIN()
#else
void oglin_()
#endif
#endif
{

      glPopMatrix();
      glPushMatrix();

      NSurf++;
      theSurf[NSurf] = glGenLists(1);

      glNewList(theSurf[NSurf], GL_COMPILE);

      glEnable(GL_COLOR_MATERIAL);
      glBegin(GL_LINES);
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
      int itrns;

      itrns = 0;
      if (*isurf < 0) {
	  itrns = 1;
      }
      TRANS = itrns;

      if (dowrt) {
	  if (itrns) {
	     fprintf(out,"[SURFACE] TRANS COLOR 1.0 1.0 0.0\n");
          } else {
             fprintf(out,"[SURFACE]\n");
          }
      }

      NSurf++;
      theSurf[NSurf] = glGenLists(1);

      glNewList(theSurf[NSurf], GL_COMPILE);

      if (mapped) {
        glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
      }

      glBegin(GL_TRIANGLES);
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
      glEnd();
      glDisable(GL_BLEND);
      glDisable(GL_COLOR_MATERIAL);
      glDisable(GL_CULL_FACE);
      glEndList();
}

void ogelev()
{
      int i;

      glPopMatrix();
      glPushMatrix();

      NSurf++;
      theSurf[NSurf] = glGenLists(1); 
      glNewList(theSurf[NSurf], GL_COMPILE);

      for (i=0; i<4; i++)
	diffuseColor[NSurf][i] = materialColor[ecol][i];

      glBegin(GL_QUADS);
}

void ogribb()
{

      glPopMatrix();
      glPushMatrix();

      NSurf++;
      theSurf[NSurf] = glGenLists(1); 
      glNewList(theSurf[NSurf], GL_COMPILE);

      glBegin(GL_QUADS);
}

void crpsin(a,b,c,d)
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

double vln( double *a)
{
      double vl;
      double tot;

      tot = a[0]*a[0]+a[1]*a[1]+a[2]*a[2];

      vl = 0.0;
      if (tot > 0.0) vl = sqrt(tot);

      return(vl);

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

void ps_header(file)
FILE *file;
{
   fprintf(file,"%%true {\n");
   fprintf(file,"systemdict /colorimage known not {\n");
   fprintf(file,"%%\n");
   fprintf(file,"/colorImageDict 50 dict def\n");
   fprintf(file,"/colorimage {\n");
   fprintf(file,"    colorImageDict begin\n");
   fprintf(file,"    /Ncomp exch def\n");
   fprintf(file,"    {\n");
   fprintf(file,"        (Multi-source not implemented\\n) print flush\n");
   fprintf(file,"        limitcheck\n");
   fprintf(file,"    } {\n");
   fprintf(file,"        /Dsrc exch def\n");
   fprintf(file,"        /Matrix exch def\n");
   fprintf(file,"        /Bcomp exch def\n");
   fprintf(file,"        /Height exch def\n");
   fprintf(file,"        /Width exch def\n");
   fprintf(file,"        /Bycomp Bcomp 7 add 8 idiv def\n");
   fprintf(file,"        Bcomp 8 gt { (Only 8 bit per sample images \\n)\n");
   fprintf(file,"                     print flush limitcheck\n");
   fprintf(file,"                   } if\n");
   fprintf(file,"        Width Height Bcomp Matrix\n");
   fprintf(file,"        Ncomp 1 eq {\n");
   fprintf(file,"            { Dsrc exec }\n");
   fprintf(file,"        } if\n");
   fprintf(file,"        Ncomp 3 eq {\n");
   fprintf(file,"          /Gstr Bycomp Width mul string def\n");
   fprintf(file,"          { Dsrc exec\n");
   fprintf(file,"             /Cstr exch def\n");
   fprintf(file,"             0 1 Width 1 sub {\n");
   fprintf(file,"               /I exch def\n");
   fprintf(file,"               /X I 3 mul def\n");
   fprintf(file,"               Gstr I\n");
   fprintf(file,"                 Cstr X       get 0.3  mul\n");
   fprintf(file,"                 Cstr X 1 add get 0.59 mul\n");
   fprintf(file,"                 Cstr X 2 add get 0.11 mul\n");
   fprintf(file,"                 add add cvi\n");
   fprintf(file,"               put\n");
   fprintf(file,"              } for\n");
   fprintf(file,"             Gstr\n");
   fprintf(file,"          }\n");
   fprintf(file,"        } if\n");
   fprintf(file,"        Ncomp 4 eq {\n");
   fprintf(file,"          /Gstr Bycomp Width mul string def\n");
   fprintf(file,"          { Dsrc exec\n");
   fprintf(file,"             /Cstr exch def\n");
   fprintf(file,"             0 1 Width 1 sub {\n");
   fprintf(file,"               /I exch def\n");
   fprintf(file,"               /X I 4 mul def\n");
   fprintf(file,"               Gstr I\n");
   fprintf(file,"                 2 Bcomp exp 1 sub\n");
   fprintf(file,"                 Cstr X       get 0.3  mul\n");
   fprintf(file,"                 Cstr X 1 add get 0.59 mul\n");
   fprintf(file,"                 Cstr X 2 add get 0.11 mul\n");
   fprintf(file,"                 Cstr X 3 add get\n");
   fprintf(file,"                 add add add dup 2 index gt {pop dup} if\n");
   fprintf(file,"                 sub cvi\n");
   fprintf(file,"               put\n");
   fprintf(file,"              } for\n");
   fprintf(file,"             Gstr\n");
   fprintf(file,"          }\n");
   fprintf(file,"        } if\n");
   fprintf(file,"        image\n");
   fprintf(file,"    } ifelse\n");
   fprintf(file,"    end\n");
   fprintf(file,"} bind def\n");
   fprintf(file,"} if\n");
}

#define Max(a,b) (((a)>(b))?(a):(b))
#define Min(a,b) (((a)<(b))?(a):(b))

/* OpenGL's GL_3D_COLOR feedback vertex format. */
typedef struct _Feedback3Dcolor {
  GLfloat x;
  GLfloat y;
  GLfloat z;
  GLfloat red;
  GLfloat green;
  GLfloat blue;
  GLfloat alpha;
} Feedback3Dcolor;




/* Write contents of one vertex to stdout. */
void
print3DcolorVertex(GLint size, GLint * count,
  GLfloat * buffer)
{
  int i;

  printf("  ");
  for (i = 0; i < 7; i++) {
    printf("%4.2f ", buffer[size - (*count)]);
    *count = *count - 1;
  }
  printf("\n");
}

void
printBuffer(GLint size, GLfloat * buffer)
{
  GLint count;
  int token, nvertices;

  count = size;
  while (count) {
    token = (int) buffer[size - count];
    count--;
    switch (token) {
    case GL_PASS_THROUGH_TOKEN:
      printf("GL_PASS_THROUGH_TOKEN\n");
      printf("  %4.2f\n", buffer[size - count]);
      count--;
      break;
    case GL_POINT_TOKEN:
      printf("GL_POINT_TOKEN\n");
      print3DcolorVertex(size, &count, buffer);
      break;
    case GL_LINE_TOKEN:
      printf("GL_LINE_TOKEN\n");
      print3DcolorVertex(size, &count, buffer);
      print3DcolorVertex(size, &count, buffer);
      break;
    case GL_LINE_RESET_TOKEN:
      printf("GL_LINE_RESET_TOKEN\n");
      print3DcolorVertex(size, &count, buffer);
      print3DcolorVertex(size, &count, buffer);
      break;
    case GL_POLYGON_TOKEN:
      printf("GL_POLYGON_TOKEN\n");
      nvertices = (int) buffer[size - count];
      count--;
      for (; nvertices > 0; nvertices--) {
        print3DcolorVertex(size, &count, buffer);
      }
    }
  }
}

GLfloat pointSize;

static char *gouraudtriangleEPS[] =
{
  "/bd{bind def}bind def /triangle { aload pop   setrgbcolor  aload pop 5 3",
  "roll 4 2 roll 3 2 roll exch moveto lineto lineto closepath fill } bd",
  "/computediff1 { 2 copy sub abs threshold ge {pop pop pop true} { exch 2",
  "index sub abs threshold ge { pop pop true} { sub abs threshold ge } ifelse",
  "} ifelse } bd /computediff3 { 3 copy 0 get 3 1 roll 0 get 3 1 roll 0 get",
  "computediff1 {true} { 3 copy 1 get 3 1 roll 1 get 3 1 roll 1 get",
  "computediff1 {true} { 3 copy 2 get 3 1 roll  2 get 3 1 roll 2 get",
  "computediff1 } ifelse } ifelse } bd /middlecolor { aload pop 4 -1 roll",
  "aload pop 4 -1 roll add 2 div 5 1 roll 3 -1 roll add 2 div 3 1 roll add 2",
  "div 3 1 roll exch 3 array astore } bd /gouraudtriangle { computediff3 { 4",
  "-1 roll aload 7 1 roll 6 -1 roll pop 3 -1 roll pop add 2 div 3 1 roll add",
  "2 div exch 3 -1 roll aload 7 1 roll exch pop 4 -1 roll pop add 2 div 3 1",
  "roll add 2 div exch 3 -1 roll aload 7 1 roll pop 3 -1 roll pop add 2 div 3",
  "1 roll add 2 div exch 7 3 roll 10 -3 roll dup 3 index middlecolor 4 1 roll",
  "2 copy middlecolor 4 1 roll 3 copy pop middlecolor 4 1 roll 13 -1 roll",
  "aload pop 17 index 6 index 15 index 19 index 6 index 17 index 6 array",
  "astore 10 index 10 index 14 index gouraudtriangle 17 index 5 index 17",
  "index 19 index 5 index 19 index 6 array astore 10 index 9 index 13 index",
  "gouraudtriangle 13 index 16 index 5 index 15 index 18 index 5 index 6",
  "array astore 12 index 12 index 9 index gouraudtriangle 17 index 16 index",
  "15 index 19 index 18 index 17 index 6 array astore 10 index 12 index 14",
  "index gouraudtriangle 18 {pop} repeat } { aload pop 5 3 roll aload pop 7 3",
  "roll aload pop 9 3 roll 4 index 6 index 4 index add add 3 div 10 1 roll 7",
  "index 5 index 3 index add add 3 div 10 1 roll 6 index 4 index 2 index add",
  "add 3 div 10 1 roll 9 {pop} repeat 3 array astore triangle } ifelse } bd",
  NULL
};

GLfloat *
spewPrimitiveEPS(FILE * file, GLfloat * loc)
{
  int token;
  int nvertices, i;
  GLfloat red, green, blue;
  int smooth;
  GLfloat dx, dy, dr, dg, db, absR, absG, absB, colormax;
  int steps;
  Feedback3Dcolor *vertex;
  GLfloat xstep, ystep, rstep, gstep, bstep;
  GLfloat xnext, ynext, rnext, gnext, bnext, distance;

  token = *loc;
  loc++;
  switch (token) {
  case GL_LINE_RESET_TOKEN:
  case GL_LINE_TOKEN:
    vertex = (Feedback3Dcolor *) loc;

    dr = vertex[1].red - vertex[0].red;
    dg = vertex[1].green - vertex[0].green;
    db = vertex[1].blue - vertex[0].blue;

    if (dr != 0 || dg != 0 || db != 0) {
      /* Smooth shaded line. */
      dx = vertex[1].x - vertex[0].x;
      dy = vertex[1].y - vertex[0].y;
#ifdef DARWIN
      distance = sqrt(dx * dx + dy * dy);
#else
      distance = sqrtf(dx * dx + dy * dy);
#endif

      absR = fabsf(dr);
      absG = fabsf(dg);
      absB = fabsf(db);


#define EPS_SMOOTH_LINE_FACTOR 0.06  /* Lower for better smooth 

                                        lines. */

      colormax = Max(absR, Max(absG, absB));
      steps = Max(1.0, colormax * distance * EPS_SMOOTH_LINE_FACTOR);

      xstep = dx / steps;
      ystep = dy / steps;

      rstep = dr / steps;
      gstep = dg / steps;
      bstep = db / steps;

      xnext = vertex[0].x;
      ynext = vertex[0].y;
      rnext = vertex[0].red;
      gnext = vertex[0].green;
      bnext = vertex[0].blue;

      /* Back up half a step; we want the end points to be
         exactly the their endpoint colors. */
      xnext -= xstep / 2.0;
      ynext -= ystep / 2.0;
      rnext -= rstep / 2.0;
      gnext -= gstep / 2.0;
      bnext -= bstep / 2.0;
    } else {
      /* Single color line. */
      steps = 0;
    }

    fprintf(file, "%g %g %g setrgbcolor\n",
      vertex[0].red, vertex[0].green, vertex[0].blue);
    fprintf(file, "%g %g moveto\n", vertex[0].x, vertex[0].y);

    for (i = 0; i < steps; i++) {
      xnext += xstep;
      ynext += ystep;
      rnext += rstep;
      gnext += gstep;
      bnext += bstep;
      fprintf(file, "%g %g lineto stroke\n", xnext, ynext);
      fprintf(file, "%g %g %g setrgbcolor\n", rnext, gnext, bnext);
      fprintf(file, "%g %g moveto\n", xnext, ynext);
    }
    fprintf(file, "%g %g lineto stroke\n", vertex[1].x, vertex[1].y);

    loc += 14;          /* Each vertex element in the feedback
                           buffer is 7 GLfloats. */

    break;
  case GL_POLYGON_TOKEN:
    nvertices = *loc;
    loc++;

    vertex = (Feedback3Dcolor *) loc;

    if (nvertices > 0) {
      red = vertex[0].red;
      green = vertex[0].green;
      blue = vertex[0].blue;
      smooth = 0;
      for (i = 1; i < nvertices; i++) {
        if (red != vertex[i].red || green != vertex[i].green || blue != vertex[i].blue) {
          smooth = 1;
          break;
        }
      }
      if (smooth) {
        /* Smooth shaded polygon; varying colors at vetices. */

        /* Break polygon into "nvertices-2" triangle fans. */
        for (i = 0; i < nvertices - 2; i++) {
          fprintf(file, "[%g %g %g %g %g %g]",
            vertex[0].x, vertex[i + 1].x, vertex[i + 2].x,
            vertex[0].y, vertex[i + 1].y, vertex[i + 2].y);
          fprintf(file, " [%g %g %g] [%g %g %g] [%g %g %g] gouraudtriangle\n",
            vertex[0].red, vertex[0].green, vertex[0].blue,
            vertex[i + 1].red, vertex[i + 1].green, vertex[i + 1].blue,
            vertex[i + 2].red, vertex[i + 2].green, vertex[i + 2].blue);
        }
      } else {
        /* Flat shaded polygon; all vertex colors the same. */
        fprintf(file, "newpath\n");
        fprintf(file, "%g %g %g setrgbcolor\n", red, green, blue);

        /* Draw a filled triangle. */
        fprintf(file, "%g %g moveto\n", vertex[0].x, vertex[0].y);
        for (i = 1; i < nvertices; i++) {
          fprintf(file, "%g %g lineto\n", vertex[i].x, vertex[i].y);
        }
        fprintf(file, "closepath fill\n\n");
      }
    }
    loc += nvertices * 7;  /* Each vertex element in the
                              feedback buffer is 7 GLfloats. */
    break;
  case GL_POINT_TOKEN:
    vertex = (Feedback3Dcolor *) loc;
    fprintf(file, "%g %g %g setrgbcolor\n", vertex[0].red, vertex[0].green, vertex[0].blue);
    fprintf(file, "%g %g %g 0 360 arc fill\n\n", vertex[0].x, vertex[0].y, pointSize / 2.0);
    loc += 7;           /* Each vertex element in the feedback
                           buffer is 7 GLfloats. */
    break;
  default:
    /* XXX Left as an excersie to the reader. */
    printf("Incomplete implementation.  Unexpected token (%d).\n", token);
    exit(1);
  }
  return loc;
}

void
spewUnsortedFeedback(FILE * file, GLint size, GLfloat * buffer)
{
  GLfloat *loc, *end;

  loc = buffer;
  end = buffer + size;
  while (loc < end) {
    loc = spewPrimitiveEPS(file, loc);
  }
}

typedef struct _DepthIndex {
  GLfloat *ptr;
  GLfloat depth;
} DepthIndex;

static int
compare(const void *a, const void *b)
{
  DepthIndex *p1 = (DepthIndex *) a;
  DepthIndex *p2 = (DepthIndex *) b;
  GLfloat diff = p2->depth - p1->depth;

  if (diff > 0.0) {
    return 1;
  } else if (diff < 0.0) {
    return -1;
  } else {
    return 0;
  }
}

void
spewSortedFeedback(FILE * file, GLint size, GLfloat * buffer)
{
  int token;
  GLfloat *loc, *end;
  Feedback3Dcolor *vertex;
  GLfloat depthSum;
  int nprimitives, item;
  DepthIndex *prims;
  int nvertices, i;

  end = buffer + size;

  /* Count how many primitives there are. */
  nprimitives = 0;
  loc = buffer;
  while (loc < end) {
    token = *loc;
    loc++;
    switch (token) {
    case GL_LINE_TOKEN:
    case GL_LINE_RESET_TOKEN:
      loc += 14;
      nprimitives++;
      break;
    case GL_POLYGON_TOKEN:
      nvertices = *loc;
      loc++;
      loc += (7 * nvertices);
      nprimitives++;
      break;
    case GL_POINT_TOKEN:
      loc += 7;
      nprimitives++;
      break;
    default:
      /* XXX Left as an excersie to the reader. */
      printf("Incomplete implementation.  Unexpected token (%d).\n",
        token);
      exit(1);
    }
  }

  /* Allocate an array of pointers that will point back at
     primitives in the feedback buffer.  There will be one
     entry per primitive.  This array is also where we keep the
     primitive's average depth.  There is one entry per
     primitive  in the feedback buffer. */
  prims = (DepthIndex *) malloc(sizeof(DepthIndex) * nprimitives);

  item = 0;
  loc = buffer;
  while (loc < end) {
    prims[item].ptr = loc;  /* Save this primitive's location. */
    token = *loc;
    loc++;
    switch (token) {
    case GL_LINE_TOKEN:
    case GL_LINE_RESET_TOKEN:
      vertex = (Feedback3Dcolor *) loc;
      depthSum = vertex[0].z + vertex[1].z;
      prims[item].depth = depthSum / 2.0;
      loc += 14;
      break;
    case GL_POLYGON_TOKEN:
      nvertices = *loc;
      loc++;
      vertex = (Feedback3Dcolor *) loc;
      depthSum = vertex[0].z;
      for (i = 1; i < nvertices; i++) {
/*        depthSum += vertex[i].z;*/
        depthSum = Min(depthSum, vertex[i].z);
      }
/*      prims[item].depth = depthSum / nvertices;*/
      prims[item].depth = depthSum;
      loc += (7 * nvertices);
      break;
    case GL_POINT_TOKEN:
      vertex = (Feedback3Dcolor *) loc;
      prims[item].depth = vertex[0].z;
      loc += 7;
      break;
/*
    default:
*/
      /* XXX Left as an excersie to the reader. */
/* assert(1); */
    }
    item++;
  }
/*
  assert(item == nprimitives);
*/

  /* Sort the primitives back to front. */
  qsort(prims, nprimitives, sizeof(DepthIndex), compare);

  /* XXX Understand that sorting by a primitives average depth
     doesn't allow us to disambiguate some cases like self
     intersecting polygons.  Handling these cases would require
     breaking up the primitives.  That's too involved for this
     example.  Sorting by depth is good enough for lots of
     applications. */

  /* Emit the Encapsulated PostScript for the primitives in
     back to front order. */
  for (item = 0; item < nprimitives; item++) {
    (void) spewPrimitiveEPS(file, prims[item].ptr);
  }

  free(prims);
}

#define EPS_GOURAUD_THRESHOLD 0.1  /* Lower for better (slower) 

                                      smooth shading. */

void
spewWireFrameEPS(FILE * file, int doSort, GLint size, GLfloat * buffer, char *creator)
{
  GLfloat clearColor[4], viewport[4];
  GLfloat lineWidth;
  int i;

  /* Read back a bunch of OpenGL state to help make the EPS
     consistent with the OpenGL clear color, line width, point
     size, and viewport. */
  glGetFloatv(GL_VIEWPORT, viewport);
  glGetFloatv(GL_COLOR_CLEAR_VALUE, clearColor);
  glGetFloatv(GL_LINE_WIDTH, &lineWidth);
  glGetFloatv(GL_POINT_SIZE, &pointSize);

  /* Emit EPS header. */
  fputs("%!PS-Adobe-2.0 EPSF-2.0\n", file);
  /* Notice %% for a single % in the fprintf calls. */
  fprintf(file, "%%%%Creator: %s (using OpenGL feedback)\n", file, creator);
  fprintf(file, "%%%%BoundingBox: %g %g %g %g\n",
    viewport[0], viewport[1], viewport[2], viewport[3]);
  fputs("%%EndComments\n", file);
  fputs("\n", file);
  fputs("gsave\n", file);
  fputs("\n", file);

  /* Output Frederic Delhoume's "gouraudtriangle" PostScript
     fragment. */
  fputs("% the gouraudtriangle PostScript fragement below is free\n", file);
  fputs("% written by Frederic Delhoume (delhoume@ilog.fr)\n", file);
  fprintf(file, "/threshold %g def\n", EPS_GOURAUD_THRESHOLD);
  for (i = 0; gouraudtriangleEPS[i]; i++) {
    fprintf(file, "%s\n", gouraudtriangleEPS[i]);
  }

  fprintf(file, "\n%g setlinewidth\n", lineWidth);

  /* Clear the background like OpenGL had it. */
  fprintf(file, "%g %g %g setrgbcolor\n",
    clearColor[0], clearColor[1], clearColor[2]);
  fprintf(file, "%g %g %g %g rectfill\n\n",
    viewport[0], viewport[1], viewport[2], viewport[3]);

  if (doSort) {
    spewSortedFeedback(file, size, buffer);
  } else {
    spewUnsortedFeedback(file, size, buffer);
  }

  /* Emit EPS trailer. */
  fputs("grestore\n\n", file);
  fputs("showpage\n", file);

  fclose(file);
}

/**************************************************************************
*
*	Save the Frame Buffer in a ps/eps format file
*
*	With thanks to Pedro Vazquez vazquez@penelope.iqm.unicamp.br
*
**************************************************************************/
void
save_ps(mode)
int mode;
{
   int x,y;
   int i,k,l,rowlen;
   FILE *file;
   GLubyte *rgbbuf;

        x = glutGet(GLUT_WINDOW_WIDTH);
        y = glutGet(GLUT_WINDOW_HEIGHT);

	rowlen = Vsize;
	if (x < Vsize) rowlen = x;

	glPixelStorei(GL_PACK_ROW_LENGTH,rowlen);
	glPixelStorei(GL_PACK_ALIGNMENT,1);

        rgbbuf = (GLubyte *)malloc(3*rowlen*y*sizeof(GLubyte));
	if (!rgbbuf) {
	   fprintf(stderr,"moldenogl: couldn't allocate memory\n");
	   return;
	}

        glReadBuffer(GL_FRONT);
        glReadPixels(0,0,x,y,GL_RGB,GL_UNSIGNED_BYTE,rgbbuf);

	if(mode)
#ifdef WIN32
        file=fopen("mogl.ps","w");
#else
        file=fopen("moldenogl.ps","w");
#endif
	else
#ifdef WIN32
        file=fopen("mogl.eps","w");
#else
        file=fopen("moldenogl.eps","w");
#endif

	if (!file) {
	   fprintf(stderr,"moldenogl: can't open output file\n");
	   return;
	}

        fprintf(file,"%%!PS-Adobe-2.0 EPSF-2.0\n");
        fprintf(file,"%%%%BoundingBox: 16 16 %d %d\n",x+16,y+16);
        fprintf(file,"%%%%Creator: Moldenogl\n");
        fprintf(file,"%%%%Title: Moldenogl output file\n");
        fprintf(file,"%%%%EndComments\n");

	ps_header(file);

	fprintf(file,"/picstr %d string def\n",y*3);

	fprintf(file,"16 16 translate\n");
        fprintf(file,"%d %d scale\n",x,y);
        fprintf(file,"%d %d 8 [ %d 0 0 %d 0 0] \n",x,y,x,y);
        fprintf(file,"{ currentfile picstr readhexstring  pop }\n");
        fprintf(file,"false 3 colorimage\n");

        k = l = 0;
	for (i=0;i<x*y+1;i++){

	   fprintf(file,"%02x%02x%02x",rgbbuf[l],rgbbuf[l+1],rgbbuf[l+2]);
	   l += 3;
           k += 6;

           if (k>70){ 
                fprintf(file,"\n");
                k=0;
           }
        }
	if (mode)
        fprintf(file,"\nshowpage\n");   
        fprintf(file,"%%%%Trailer\n");
        fclose(file);
        free(rgbbuf);

}

/**************************************************************************
*
*
*       Save the Frame Buffer in a ppm format file
*
*	With thanks to Pedro Vazquez vazquez@penelope.iqm.unicamp.br
*
**************************************************************************/

save_pixmap(ppmfile)
char *ppmfile;
{       
FILE *file;
int i,j,k,x,y,rowlen;
GLubyte *rgbbuf;

	x = glutGet(GLUT_WINDOW_WIDTH);
	y = glutGet(GLUT_WINDOW_HEIGHT);

	rowlen = Vsize;
	if (x < Vsize) rowlen = x;

	glPixelStorei(GL_PACK_ROW_LENGTH,rowlen);
	glPixelStorei(GL_PACK_ALIGNMENT,1);

	rgbbuf = (GLubyte *) malloc(3*rowlen*y*sizeof(GLubyte));
	if (!rgbbuf) {
            fprintf(stderr,"moldenogl: couldn't allocate memory\n");
            return;
	}

	glReadBuffer(GL_FRONT);
	glReadPixels(0,0,x,y,GL_RGB,GL_UNSIGNED_BYTE,rgbbuf);

        file = fopen(ppmfile,"w");

        if (!file) {
            fprintf(stderr,"moldenogl: can't open output file\n");
            return;
        }
        fprintf(file,"P6\n");
        fprintf(file,"#Image rendered with moldenogl\n");
        fprintf(file,"%d\n%d\n255\n", x,y);

	for(i=y-1; i>= 0; i--){
	   for(j=0; j< x; j++){
		k = 3*(j + i*x);
		fwrite( &rgbbuf[k] ,sizeof(*rgbbuf), 1, file);
		fwrite( &rgbbuf[k+1] ,sizeof(*rgbbuf), 1, file);
		fwrite( &rgbbuf[k+2] ,sizeof(*rgbbuf), 1, file);
	   }
	}

	fclose(file);
	free(rgbbuf);
} 



putbyte(outf,val)
FILE *outf;
unsigned char val;
{
        unsigned char buf[1];

        buf[0] = val;
        fwrite(buf,1,1,outf);
}

putshort(outf,val)
FILE *outf;
unsigned short val;
{
        unsigned char buf[2];

        buf[0] = (val>>8);
        buf[1] = (val>>0);
        fwrite(buf,2,1,outf);
    }

static int putlong(outf,val)
FILE *outf;
unsigned long val;
{
	unsigned char buf[4];

	buf[0] = (val>>24);
	buf[1] = (val>>16);
	buf[2] = (val>>8);
	buf[3] = (val>>0);
	return fwrite(buf,4,1,outf);
}

save_rgb(rgbfile)
char *rgbfile;
{
	FILE *of;
        char iname[80];
        int i, k, rowlen;
	int Xsize, Ysize;
	GLubyte *rgbbuf;

        of = fopen(rgbfile,"w");

        if (!of) {
            fprintf(stderr,"moldenogl: can't open output file\n");
            return;
        }

	Xsize = glutGet(GLUT_WINDOW_WIDTH);
	Ysize = glutGet(GLUT_WINDOW_HEIGHT);

	rowlen = Vsize;
	if (Xsize < Vsize) rowlen = Xsize;

	glPixelStorei(GL_PACK_ROW_LENGTH,rowlen);
	glPixelStorei(GL_PACK_ALIGNMENT,1);

	rgbbuf = (GLubyte *)malloc(3*rowlen*Ysize*sizeof(GLubyte));
	if (!rgbbuf) {
	   fprintf(stderr,"moldenogl: couldn't allocate memory\n");
	   close(of);
	   return;
	}
	glReadBuffer(GL_FRONT);
	glReadPixels(0,0,Xsize,Ysize,GL_RGB,GL_UNSIGNED_BYTE,rgbbuf);

        putshort(of,474);	/* MAGIC		*/
        putbyte(of,0);		/* STORAGE is VERBATIM	*/
        putbyte(of,257);	/* BPC is 257          	*/
        putshort(of,3);		/* DIMENSION is 3	*/
        putshort(of,Xsize);	/* XSIZE               	*/
        putshort(of,Ysize);	/* YSIZE               	*/
        putshort(of,3);		/* ZSIZE               	*/
        putlong(of,0);		/* PIXMIN is 0         	*/
        putlong(of,255);	/* PIXMAX is 255       	*/
        for(i=0; i<4; i++)	/* DUMMY 4 bytes 	*/
            putbyte(of,0);
        strcpy(iname,"Moldenogl");
        fwrite(iname,80,1,of);	/* IMAGENAME  		*/
        putlong(of,0);		/* COLORMAP is 0 	*/

        for(i=0; i<404; i++)
            putbyte(of,0);

/* red */
        k=0;
        for(i=0;i< Xsize*Ysize; i++){
         fwrite( &rgbbuf[k] ,sizeof(*rgbbuf), 1, of);
         k=k+3;
        }

/* green */
        k=0;
        for(i=0;i< Xsize*Ysize; i++){
         fwrite( &rgbbuf[k+1] ,sizeof(*rgbbuf), 1, of);
         k=k+3;
        }

/* blue */
        k=0;
        for(i=0;i< Xsize*Ysize; i++){
         fwrite( &rgbbuf[k+2] ,sizeof(*rgbbuf), 1, of);
         k=k+3;
        }

        fclose(of);
	free(rgbbuf);
}

unsigned char bmp_header[]=
{ 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0,
  40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0, 0,0,0,0, 0,0,0,0,
  0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 };

static void WLSBL(val,arr)
    int val;
    char* arr;
{
    arr[0] = (char) (val&0xff);
    arr[1] = (char) ((val>>8) &0xff);
    arr[2] = (char) ((val>>16)&0xff);
    arr[3] = (char) ((val>>24)&0xff);
}

save_bmp()
{ int i,j;
  FILE *fp;
  GLubyte *rgbbuf;
  GLubyte rgbtmp[3];
  int Xsize, Ysize, rowlen;
  int pad;

#ifdef WIN32
  if ((fp=fopen("mogl.bmp","wb"))==NULL) return(1);
#else
  if ((fp=fopen("moldenogl.bmp","wb"))==NULL) return(1);
#endif

  Xsize = glutGet(GLUT_WINDOW_WIDTH);
  Ysize = glutGet(GLUT_WINDOW_HEIGHT);

  rowlen = Vsize;
  if (Xsize < Vsize) rowlen = Xsize;

  glPixelStorei(GL_PACK_ROW_LENGTH,rowlen);
  glPixelStorei(GL_PACK_ALIGNMENT,1);

  rgbbuf = (GLubyte *)malloc(3*rowlen*Ysize*sizeof(GLubyte));
  if (!rgbbuf) {
	fprintf(stderr,"moldenogl: couldn't allocate memory\n");
	close(fp);
	return;
  }
  glReadBuffer(GL_FRONT);
  glReadPixels(0,0,Xsize,Ysize,GL_RGB,GL_UNSIGNED_BYTE,rgbbuf);

/* The number of bytes on a screenline should be wholly devisible by 4 */

  pad = (Xsize*3)%4;
  if (pad) pad = 4 - pad;

  WLSBL((int) (3*Xsize+pad)*Ysize+54,bmp_header+2);
  WLSBL((int) Xsize,bmp_header+18);
  WLSBL((int) Ysize,bmp_header+22);
  WLSBL((int) 3*Xsize*Ysize,bmp_header+34);

  fwrite(bmp_header,1,54,fp);

  for (i=0;i<Ysize;i++) {
    for (j=0;j<Xsize;j++) {
	rgbtmp[0] = rgbbuf[(j+Xsize*i)*3+2];
	rgbtmp[1] = rgbbuf[(j+Xsize*i)*3+1];
	rgbtmp[2] = rgbbuf[(j+Xsize*i)*3+0];
	fwrite(rgbtmp,3,1,fp);
    }
    rgbtmp[0] = (char) 0;
    for (j=0;j<pad;j++) fwrite(rgbtmp,1,1,fp);
  }

  fclose(fp);
  free(rgbbuf);
}

BuildList(r,cnst,nnpts1,nnpts2,dens)
double *r;
double *cnst;
int *nnpts1;
int *nnpts2;
double *dens;
{

     int i,j,noff1,noff2,npts1,npts2;
     double v[3], rpts;
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

      NSurf++;
      theSurf[NSurf] = glGenLists(1); 

      glPopMatrix();
      glPushMatrix();

      glNewList(theSurf[0], GL_COMPILE);

      setColor(2);

      glBegin(GL_QUADS);

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

      glEnd();

      glEndList();

}

void dispsf(void)
{
  int i;
  double ca,sa,x,y,z;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (PERSP) {
	glFrustum(-0.2,0.2,-0.2,0.2,0.3,300);
  } else {
	glOrtho(-1.0,1.0,-1.0,1.0,-10.0,10.0);
  }

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  if (PERSP) {
	glTranslatef( posx, posy, posz);
  } else {
	glScalef(posz, posz, posz);
	glTranslatef( posx, posy, 0.0);
  }

  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  if (DoLights) {
     glPushMatrix();
     glTranslatef(light0_position[0],light0_position[1],light0_position[2]);
     gluSphere(sphere, 1, 10, 10);
     glPopMatrix();
  }

  glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
  if (DoLights) {
     glPushMatrix();
     glTranslatef(light1_position[0],light1_position[1],light1_position[2]);
     gluSphere(sphere, 1, 10, 10);
     glPopMatrix();
  }

  glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
  if (DoLights) {
     glPushMatrix();
     glTranslatef(light2_position[0],light2_position[1],light2_position[2]);
     gluSphere(sphere, 1, 10, 10);
     glPopMatrix();
  }


  if (RotType) {

	ca =  cos(angincr2*TORAD);
	sa =  sin(angincr2*TORAD);
	for (i=0; i < 3; i++) {
	   y = RR[i][1];
	   z = RR[i][2];
	   RR[i][1] = ca*y - sa*z;
	   RR[i][2] = ca*z + sa*y;
	}

	ca =  cos(angincr1*TORAD);
	sa =  sin(angincr1*TORAD);
	for (i=0; i < 3; i++) {
	   x = RR[i][0];
	   z = RR[i][2];
	   RR[i][0] = ca*x + sa*z;
	   RR[i][2] = ca*z - sa*x;
	}

	glMultMatrixd((const GLdouble *) RR);

  } else {

	glRotatef( ang2, 1.0, 0.0, 0.0 );
	glRotatef( ang1, 0.0, 0.0, 1.0 );
  }
  
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glPushMatrix();

  if (molon) {
      glCallList(theMol[CMols]);
  }

  if (DoLines) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glPopMatrix();

  for (i=0; i<NSurf; i++) {

      if (TRANS) {
         glEnable(GL_BLEND);
         glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	 diffuseColor[i][3] = tr_val;
      } else {
	 diffuseColor[i][3] = 1.0;
         glDisable(GL_BLEND);
      }

      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, &diffuseColor[i][0]);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambientColor);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, &specularColor[i][0]);
      glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100);
      
      glCallList(theSurf[i]);

  }

  glFlush();
  glutSwapBuffers();
}

void
outputEPS(int size, int doSort, char *filename)
{
  GLfloat *feedbackBuffer;
  GLint returned;
  FILE *file;

  feedbackBuffer = (GLfloat*) calloc(size, sizeof(GLfloat));
  glFeedbackBuffer(size, GL_3D_COLOR, feedbackBuffer);
  (void) glRenderMode(GL_FEEDBACK);
  glPushMatrix();
  dispsf();
  glPopMatrix();
  returned = glRenderMode(GL_RENDER);
  if (filename) {
    file = fopen(filename, "w");
    if (file) {
      spewWireFrameEPS(file, doSort, returned, feedbackBuffer, "rendereps");
    } else {
      printf("Could not open %s\n", filename);
    }
  } else {
    /* Helps debugging to be able to see the decode feedback
       buffer as text. */
    printBuffer(returned, feedbackBuffer);
  }
  free(feedbackBuffer);
}


void Reshape(int width, int height)
{

  Vsize = width;
  if (height > width) Vsize = height;

  glViewport(0, 0, Vsize, Vsize);
}


static void Key( unsigned char key, int x, int y )
{
   (void) x;
   (void) y;

   if (moving) return;

   switch (key) {
   case 'x':
      light0_position[0] += 0.5;
      break;
   case 'X':
      light0_position[0] -= 0.5;
      break;
   case 'y':
      light0_position[1] += 0.5;
      break;
   case 'Y':
      light0_position[1] -= 0.5;
      break;
   case 'z':
      light0_position[2] += 0.5;
      break;
   case 'Z':
      light0_position[2] -= 0.5;
      break;
   case 'i':
      light1_position[0] += 0.5;
      break;
   case 'I':
      light1_position[0] -= 0.5;
      break;
   case 'j':
      light1_position[1] += 0.5;
      break;
   case 'J':
      light1_position[1] -= 0.5;
      break;
   case 'k':
      light1_position[2] += 0.5;
      break;
   case 'K':
      light1_position[2] -= 0.5;
      break;
   case '<':
   case ',':
      light2_position[2] += 0.5;
      break;
   case '>':
   case '.':
      light2_position[2] -= 0.5;
      break;
   case 'F':
   case 'f':
      glutFullScreen();
      break;
   case 'T':
   case 't':
      if (TRANS) TRANS = 0;
      else TRANS = 1;
      break;
   case 'P':
   case 'p':
      if (PERSP) {
	PERSP = 0;
	if (poszz > 0.0) posz = (1.0 / poszz);
      } else {
	PERSP = 1;
	posz = poszz*-1.85;
      }
      break;
   case '1':
      if (l1on) {
	l1on = 0;
	glDisable(GL_LIGHT0);
      } else {
	l1on = 1;
	glEnable(GL_LIGHT0);
      }
      break;
   case '2':
      if (l2on) {
	l2on = 0;
	glDisable(GL_LIGHT1);
      } else {
	l2on = 1;
	glEnable(GL_LIGHT1);
      }
      break;
   case '3':
      if (l3on) {
	l3on = 0;
	glDisable(GL_LIGHT2);
      } else {
	l3on = 1;
	glEnable(GL_LIGHT2);
      }
      break;
   case '0':
      if (DoLights) {
	DoLights = 0;
      } else {
	DoLights = 1;
      }
      break;
   case 27:
      exit(0);
      break;
   case 'M':
   case 'm':
      if (molon) {
	molon = 0;
      } else {
	molon = 1;
      }
      break;
   case 'N':
   case 'n':
      CMols++;
      if (CMols >= NMols) CMols = 0;
      break;
   case 'g':
      if (DoFog) {
          DoFog = 0;
          glDisable(GL_FOG);	  
      } else {
         DoFog = 1;
         glEnable(GL_FOG);
         {
             GLfloat fogColor[4] = {0.5, 0.5, 0.5, 1.0};
   
	     if (ColBG != -1) {
                fogColor[0] = clrColor[ColBG][0];
                fogColor[1] = clrColor[ColBG][1];
                fogColor[2] = clrColor[ColBG][2];
             } else {
                fogColor[0] = clrColor[bgcol][0];
                fogColor[1] = clrColor[bgcol][1];
                fogColor[2] = clrColor[bgcol][2];
	     }
             glFogi (GL_FOG_MODE, GL_EXP2);
             glFogfv (GL_FOG_COLOR, fogColor);
             glFogf (GL_FOG_DENSITY, 0.03);
             glHint (GL_FOG_HINT, GL_DONT_CARE);
             glFogf (GL_FOG_START, posz);
             glFogf (GL_FOG_END, posz+5);
          }
       }
      break;
   }
   angincr1 = 0.0;
   angincr2 = 0.0;
   glutPostRedisplay();
}

static void SpecialKey( int key, int x, int y )
{
   (void) x;
   (void) y;

   if (moving) return;

   switch (key) {
   case GLUT_KEY_LEFT:
      yrot -= ROTINCR;
      break;
   case GLUT_KEY_RIGHT:
      yrot += ROTINCR;
      break;
   case GLUT_KEY_UP:
      xrot += ROTINCR;
      break;
   case GLUT_KEY_DOWN:
      xrot -= ROTINCR;
      break;
   case GLUT_KEY_F9:
      save_bmp();
      break;
   case GLUT_KEY_F10:
      save_rgb(rgbfile);
      break;
   case GLUT_KEY_F11:
      save_ps(1);
      break;
   case GLUT_KEY_F12:
      save_pixmap(ppmfile);
      break;
   case GLUT_KEY_PAGE_UP:
      fdens += 0.001;
      glFogf (GL_FOG_DENSITY, fdens);
      break;
   case GLUT_KEY_PAGE_DOWN:
      fdens -= 0.001;
      glFogf (GL_FOG_DENSITY, fdens);
      break;
   default:
      return;
   }
   angincr1 = 0.0;
   angincr2 = 0.0;
   dispsf();
   glutPostRedisplay();
}

static void
motion(int x, int y)
{
 int i;
  if (moving) {
    moving = 0;
    ang1 = ang1 + (x - startx);
    ang2 = ang2 + (y - starty);
/*
    angincr1 = (double) (x - startx) ;
    angincr2 = (double) (y - starty) ;
*/
    for (i=0; i < HISTMAX-1; i++) AngIncrHist1[i] = AngIncrHist1[i+1];
    for (i=0; i < HISTMAX-1; i++) AngIncrHist2[i] = AngIncrHist2[i+1];
    AngIncrHist1[HISTMAX-1] = (double) (x - startx) ;
    AngIncrHist2[HISTMAX-1] = (double) (y - starty) ;
    angincr1 = 0.0;
    for (i=0; i < HISTMAX-1; i++) angincr1 = angincr1 + AngIncrHist1[i];
    angincr1 = MULTH * angincr1 / (double) HISTMAX;
    angincr2 = 0.0;
    for (i=0; i < HISTMAX-1; i++) angincr2 = angincr2 + AngIncrHist2[i];
    angincr2 = MULTH * angincr2 / (double) HISTMAX;

    startx = x;
    starty = y;
    glutPostRedisplay();
    moving = 1;
  }
  if (mmoving) {
    posx = posx + (x - mstartx) / 600.0;
    posy = posy - (y - mstarty) / 600.0;
    mstartx = x;
    mstarty = y;
    angincr1 = 0.0;
    angincr2 = 0.0;
    glutPostRedisplay();
  }
  if (smoving) {
    posz = posz - (y - sstarty) / 60.0;
    if (posz < 0.005 && !PERSP) posz = 0.005;
    sstarty = y;
    angincr1 = 0.0;
    angincr2 = 0.0;
    glutPostRedisplay();
  }
}

static void
mouse(int button, int state, int x, int y)
{
  int modf;

  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_DOWN) {
      modf = glutGetModifiers();
      if (modf & GLUT_ACTIVE_SHIFT) {
         mmoving = 1;
         mstartx = x;
         mstarty = y;
      } else if (modf & GLUT_ACTIVE_CTRL) {
         smoving = 1;
         sstarty = y;
      } else {
         moving = 1;
         startx = x;
         starty = y;
      }
    }
    if (state == GLUT_UP) {
      moving = 0;
      mmoving = 0;
      smoving = 0;
    }
  }
}

void menu_select(int iopt)
{
  switch (iopt) {
  case 1:
    molon = 0;
    break;
  case 2:
    molon = 1;
    break;
  case 3:
    if (TRANS) TRANS = 0;
    else TRANS = 1;
    break;
  case 4:
    if (PERSP) {
	PERSP = 0;
	if (poszz > 0.0) posz = (1.0 / poszz);
    } else {
	PERSP = 1;
	posz = poszz*-1.85;
    }
    break;
  case 5:
    if (l1on) {
	l1on = 0;
	glDisable(GL_LIGHT0);
    } else {
	l1on = 1;
	glEnable(GL_LIGHT0);
    }
    break;
  case 6:
    if (l2on) {
	l2on = 0;
	glDisable(GL_LIGHT1);
    } else {
	l2on = 1;
	glEnable(GL_LIGHT1);
    }
    break;
  case 7:
    if (l3on) {
	l3on = 0;
	glDisable(GL_LIGHT2);
    } else {
	l3on = 1;
	glEnable(GL_LIGHT2);
    }
    break;
  case 8:
    glutFullScreen();
    break;
  case 9:
    if (DoLines) DoLines = 0;
    else DoLines = 1;
    break;
  case 10:
    if (AutoRot) AutoRot = 0;
    else {
       AutoRot = 1;
       angincr2 = 0.0;
       glutPostRedisplay();
       return;
    }
    break;
  case 11:
    exit(0);
    break;
  }

  angincr1 = 0.0;
  angincr2 = 0.0;
  glutPostRedisplay();
}

void bg_select(int iopt)
{

  ColBG = iopt;
  glClearColor(clrColor[iopt][0],clrColor[iopt][1],clrColor[iopt][2],
		clrColor[iopt][3]);
  glutPostRedisplay();
}

void anim_select(int iopt)
{

  switch (iopt) {
  case 1:
	if (Anim) Anim = 0;
	else Anim = 1;
	break;
  case 2:
	CMols++;
	if (CMols >= NMols) CMols = 0;
	break;
  case 3:
	AnimSlow = 0;
	AnimCount = 0;
	break;
  case 4:
	AnimSlow = 2;
	AnimCount = 2;
	break;
  case 5:
	AnimSlow = 5;
	AnimCount = 5;
	break;
  case 6:
	AnimSlow = 10;
	AnimCount = 10;
	break;
  }
  glutPostRedisplay();
}

void capture_select(int iopt)
{

  switch (iopt) {
  case 1:
	save_bmp();
	break;
  case 2:
	save_rgb(rgbfile);
	break;
  case 3:
	save_ps(1);
	break;
  case 4:
	glutSetCursor(GLUT_CURSOR_WAIT);
#ifdef WIN32
        outputEPS(objectComplexity, 1, "mogl.ps");
#else
	outputEPS(objectComplexity, 1, "moldenogl.ps");
#endif
	glutSetCursor(GLUT_CURSOR_INHERIT);
	break;
  case 5:
	save_pixmap(ppmfile);
	break;
  case 6:
	if (DoCap) DoCap = 0;
	else DoCap = 1;
	break;
  }
  glutPostRedisplay();
}

#if defined(VMS) || defined(UNDERSC)
void initog(int *ofil)
#else
#ifdef CRAY
void INITOG(int *ofil)
#else
void initog_(int *ofil)
#endif
#endif

{

    int i,bgmenu,animmenu,capmenu;

    static float lmodel_ambient[] = {0.1, 0.1, 0.1, 0.1};
    static float lmodel_twoside[] = {GL_TRUE};
    static float lmodel_local[] = {GL_FALSE};
    static float light0_ambient[] = {0.1, 0.1, 0.1, 1.0};
    static float light0_diffuse[] = {1.0, 1.0, 1.0, 0.0};
    static float light0_specular[] = {1.0, 1.0, 1.0, 0.0};
    static float light1_ambient[] = {0.1, 0.1, 0.1, 1.0};
    static float light1_diffuse[] = {1.0, 1.0, 1.0, 1.0};
    static float light1_specular[] = {1.0, 1.0, 1.0, 0.0};
    static float light2_ambient[] = {0.3, 0.3, 0.3, 1.0};
    static float light2_diffuse[] = {1.0, 1.0, 1.0, 1.0};
    static float light2_specular[] = {1.0, 1.0, 1.0, 0.0};
    int argc = 2;
#ifdef WIN32
    char *argv[] = {"mogl","yes"};
#else
    char *argv[] = {"moldenogl","yes"};
#endif

    if (*ofil) {
       out = fopen("molden.ogl","w");

       if (out != NULL) {
          dowrt = 1;
       } else {
          dowrt = 0;
       }
    } else {
      dowrt = 0;
    }

    if (dowrt) fprintf(out,"[MOLDENOGL]\n");

#ifdef __CYGWIN__
    if (getenv("DISPLAY") == NULL) putenv("DISPLAY=localhost:0.0");
#endif
    glutInit(&argc, argv);

    glutInitWindowSize( width, height);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );

    if (glutCreateWindow("Molden") == GL_FALSE) {
        exit(1);
    }

    glutReshapeFunc(Reshape);
    glutSpecialFunc(SpecialKey);
    glutKeyboardFunc(Key);
    glutDisplayFunc(dispsf);
    glutMotionFunc(motion);
    glutMouseFunc(mouse);
    glutIdleFunc(idle);

    bgmenu = glutCreateMenu(bg_select);
    for (i=0;i<Nclr;i++) glutAddMenuEntry(clrNames[i],i);

    animmenu = glutCreateMenu(anim_select);

    glutAddMenuEntry("Animation On/Off", 1);
    glutAddMenuEntry("Next Point [n]", 2);
    glutAddMenuEntry("Dont Slow Down", 3);
    glutAddMenuEntry("Slow Down Factor 2", 4);
    glutAddMenuEntry("Slow Down Factor 5", 5);
    glutAddMenuEntry("Slow Down Factor 10", 6);

    capmenu = glutCreateMenu(capture_select);
    glutAddMenuEntry("BMP Pixmap [F9]", 1);
    glutAddMenuEntry("IRIS RGB format [F10]", 2);
    glutAddMenuEntry("Postscript [F11]", 3);
    glutAddMenuEntry("Vector Postscript", 4);
    glutAddMenuEntry("PPM Pixmap [F12]", 5);
    glutAddMenuEntry("RGB file each update", 6);

    glutCreateMenu(menu_select);
    glutAddSubMenu("BackGround Color",bgmenu);
    glutAddSubMenu("Animation",animmenu);
    glutAddSubMenu("Screen Capture",capmenu);
    glutAddMenuEntry("Molecule Off", 1);
    glutAddMenuEntry("Molecule On", 2);
    glutAddMenuEntry("Transparency Toggle", 3);
    glutAddMenuEntry("Perspective Toggle", 4);
    glutAddMenuEntry("Light1 Toggle", 5);
    glutAddMenuEntry("Light2 Toggle", 6);
    glutAddMenuEntry("Light3 Toggle", 7);
    glutAddMenuEntry("Fullscreen", 8);
    glutAddMenuEntry("SurfLines Toggle", 9);
    glutAddMenuEntry("Auto Rotate", 10);
    glutAddMenuEntry("Quit", 11);

    glutAttachMenu(GLUT_RIGHT_BUTTON);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if (RotType) {
	glPushMatrix();
	glLoadMatrixd((const GLdouble *) RR);
	glRotatef( ang2, 1.0, 0.0, 0.0 );
	glRotatef( ang1, 0.0, 0.0, 1.0 );
	glGetDoublev(GL_MODELVIEW_MATRIX,(GLdouble *) RR);
	glPopMatrix();
    }

    glEnable(GL_DEPTH_TEST);

    if (l1on) glEnable(GL_LIGHT0);
    else glDisable(GL_LIGHT0);
    if (l2on) glEnable(GL_LIGHT1);
    else glDisable(GL_LIGHT1);
    if (l3on) glEnable(GL_LIGHT2);
    else glDisable(GL_LIGHT2);

    glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_position);

    glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

    glLightfv(GL_LIGHT2, GL_AMBIENT, light2_ambient);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, light2_diffuse);
    glLightfv(GL_LIGHT2, GL_SPECULAR, light2_specular);
    glLightfv(GL_LIGHT2, GL_POSITION, light2_position);

    glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, lmodel_local);
    glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, lmodel_twoside);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glEnable(GL_LIGHTING);

    glEnable(GL_NORMALIZE);
    glShadeModel(GL_SMOOTH);
    glLineWidth(1.5);
    glClearColor(clrColor[bgcol][0],clrColor[bgcol][1],clrColor[bgcol][2],
		clrColor[bgcol][3]);
    PERSP = perspon;
    glClearIndex(0);
    glClearDepth(1);

    glPushMatrix();
    glLoadIdentity();

    cyl = gluNewQuadric();
    sphere = gluNewQuadric();

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
    int mopt;

    mopt = 0;
    ogmol(r,adjus,natoms,nat,NULL,icol,xsym,ysym,zsym,vdwr,&mopt,NULL,NULL,NULL);
    BuildList(r,cnst,nnpts1,nnpts2,dens);

    NSurf++;
    NMols++;
    glutMainLoop();

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


    AmbAndDiff = 0;
    ogmol(r,adjus,natoms,nat,iatclr,icol,xsym,ysym,zsym,vdwr,mopt,conn,nconn,iconn);

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

    if (dowrt) fclose(out);
    
    NSurf++;
    NMols++;

    if (perspon) {
	posz = poszz*-1.85;
	light2_position[2] = poszz*2.5;
    } else {
	if (poszz > 0.0) {
	   posz = (1.0 / poszz);
	   light2_position[2] = posz;
	}
    }
    glutMainLoop();
}


