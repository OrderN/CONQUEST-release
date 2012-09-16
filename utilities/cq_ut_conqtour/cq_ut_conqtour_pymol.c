/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_pymol.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_pymol.c *
 *
 * NAME
 *  cq_ut_conqtour_pymol.c
 * PURPOSE
 *
 * USES
 *
 * AUTHOR
 *  torralba
 * CREATION DATE
 *  Oct 7, 2010
 * MODIFICATION HISTORY
 *
 * *****/

#include "cq_ut_conqtour.h"

int printPymolScript(index_t *idx, density_t *d, char *file)
{
  //  FILE *fp;
  //  int i, j, k;
  //  double x, y, z;
  //  //   double dx,dy,dz;
  //  //   double absmax;
  //  //   GLuint cut;
  //  double pt[4][3];
  //  double aux1[3], aux2[3];

  if ((idx->fp = fopen(file, "w")) == NULL)
    {
      fprintf(stderr, "I can't create pymol script '%s'\n", file);
      return (1); /* Error */
    }
  if (drawSlice(idx, d)) /* This will write the file */
    {
      fclose(idx->fp);
      idx->fp = NULL; /* Be explicit */
      return (1); /* Error */
    }
  if (ctrl.verbose)
    printf("INFO: Wrote slice to pymol script '%s'\n", file);
  fclose(idx->fp);
  idx->fp = NULL; /* Be explicit */
  return (0);
}

/* All these functions parallel the ones in the slicevis file */
int printPymolScriptHead(index_t *idx, density_t *d)
{
  if (fprintf(idx->fp, "from pymol.cgo import *\n") < 0)
    return (1);
  if (fprintf(idx->fp, "from pymol import cmd\n\n") < 0)
    return (1);
  if (fprintf(idx->fp, "slice = []\n\n") < 0)
    return (1);
  if (fprintf(idx->fp, "slice.extend( [BEGIN, TRIANGLES] +\n") < 0)
    return (1);
  if (fprintf(idx->fp, "[NORMAL, %f, %f, %f] )\n", idx->surfnormal[0], idx->surfnormal[1], idx->surfnormal[2]) < 0)
    return (1);

  return (0); /* No problem */
}

int printPymolScriptTail(index_t *idx, density_t *d)
{
  if (fprintf(idx->fp, "slice.extend( [END] )\n") < 0)
    return (1);
  if (fprintf(idx->fp, "cmd.load_cgo(slice,'slice')\n") < 0)
    return (1);

  return (0); /* No problem */
}

//fprintf(fp, "slice.extend( [VERTEX, %f, %f, %f] )\n", pt[3][0], pt[3][1], pt[3][2]);
//fprintf(fp, "slice.extend( [COLOR, %f, %f, %f] )\n", c, 0.0, 0.0);

int writePymolVertex(index_t *idx, double x, double y, double z)
{
  int i;
  double shift[3] = { 0.0, 0.0, 0.0 };

  if(ctrl.shftcoords)
    for(i = 0; i < 3 ; ++i)
      shift[i] = limits.realshft[i];
  if (fprintf(idx->fp, "slice.extend( [VERTEX, %f, %f, %f] )\n", x + shift[0], y + shift[1], z + shift[2]) < 0)
    return (1);
  return(0);
}

int drawPymolQuad(double **pt, int **c, index_t *idx, density_t *d)
{
  int state;

  state=0;
  colorifyPymolGradient(idx, d, c[3][X], c[3][Y]);
  state += writePymolVertex(idx, pt[3][0], pt[3][1], pt[3][2]);
  //  glVertex3f(pt[3][0], pt[3][1], pt[3][2]);
  colorifyPymolGradient(idx, d, c[2][X], c[2][Y]);
  state += writePymolVertex(idx, pt[2][0], pt[2][1], pt[2][2]);
  //  glVertex3f(pt[2][0], pt[2][1], pt[2][2]);
  colorifyPymolGradient(idx, d, c[0][X], c[0][Y]);
  state += writePymolVertex(idx, pt[0][0], pt[0][1], pt[0][2]);
  //  glVertex3f(pt[0][0], pt[0][1], pt[0][2]);

  colorifyPymolGradient(idx, d, c[3][X], c[3][Y]);
  state += writePymolVertex(idx, pt[3][0], pt[3][1], pt[3][2]);
  //  glVertex3f(pt[3][0], pt[3][1], pt[3][2]);
  colorifyPymolGradient(idx, d, c[0][X], c[0][Y]);
  state += writePymolVertex(idx, pt[0][0], pt[0][1], pt[0][2]);
  //  glVertex3f(pt[0][0], pt[0][1], pt[0][2]);
  colorifyPymolGradient(idx, d, c[1][X], c[1][Y]);
  state += writePymolVertex(idx, pt[1][0], pt[1][1], pt[1][2]);
  //  glVertex3f(pt[1][0], pt[1][1], pt[1][2]);

  if(state) state = 1; /* Normalise the return value */

  return (state);
}

int drawPymolTriangle(double **pt, int **c, index_t *idx, density_t *d)
{
  int state;

  state=0;
  colorifyPymolGradient(idx, d, c[0][X], c[0][Y]);
  state += writePymolVertex(idx, pt[0][0], pt[0][1], pt[0][2]);
  //  glVertex3f(pt[0][0], pt[0][1], pt[0][2]);
  colorifyPymolGradient(idx, d, c[1][X], c[1][Y]);
  state += writePymolVertex(idx, pt[1][0], pt[1][1], pt[1][2]);
  //  glVertex3f(pt[1][0], pt[1][1], pt[1][2]);
  colorifyPymolGradient(idx, d, c[2][X], c[2][Y]);
  state += writePymolVertex(idx, pt[2][0], pt[2][1], pt[2][2]);
  //  glVertex3f(pt[2][0], pt[2][1], pt[2][2]);

  if(state) state = 1; /* Normalise the return value */

  return (state);
}

int drawPymolCorner(double **pt, int x, int y, int side, index_t *idx, density_t *d)
{
  densitycut_t *spt;
  int state = 0;

  spt = idx->dindex[idx->curr];
  colorifyPymolGradient(idx, d, x, y);
  state += writePymolVertex(idx, pt[0][0], pt[0][1], pt[0][2]);
  //  glVertex3f(pt[0][0], pt[0][1], pt[0][2]);
  colorifyPymolGradient2(idx, d, spt->cturn[side]);
  state += writePymolVertex(idx, pt[1][0], pt[1][1], pt[1][2]);
  //  glVertex3f(pt[1][0], pt[1][1], pt[1][2]);
  colorifyPymolGradient(idx, d, x, y + 1);
  state += writePymolVertex(idx, pt[2][0], pt[2][1], pt[2][2]);
  //  glVertex3f(pt[2][0], pt[2][1], pt[2][2]);

  if(state) state = 1; /* Normalise the return value */

  return (state);
}

int drawPymolPeak(double **pt, int y, index_t *idx, density_t *d)
{
  densitycut_t *spt;
  int state = 0;

  spt = idx->dindex[idx->curr];
  colorifyPymolGradient(idx, d, 0, y); /* Note this means the left border color */
  state += writePymolVertex(idx, pt[0][0], pt[0][1], pt[0][2]);
  //  glVertex3f(pt[0][0], pt[0][1], pt[0][2]);
  colorifyPymolGradient2(idx, d, spt->cturn[2]);
  state += writePymolVertex(idx, pt[1][0], pt[1][1], pt[1][2]);
  //  glVertex3f(pt[1][0], pt[1][1], pt[1][2]);
  colorifyPymolGradient(idx, d, spt->pt[X] - 1, y); /* Note this means the right border color */
  state += writePymolVertex(idx, pt[2][0], pt[2][1], pt[2][2]);
  //  glVertex3f(pt[2][0], pt[2][1], pt[2][2]);

  if(state) state = 1; /* Normalise the return value */

  return (state);
}

int colorifyPymolGradient(index_t *idx, density_t *d, int x, int y)
{
  rgbcolor_t rgb;

  getColorFromScale(d, idx->dindex[idx->curr]->data[x][y], &rgb);
  if (fprintf(idx->fp, "slice.extend( [COLOR, %f, %f, %f] )\n", rgb.r / 255., rgb.g / 255., rgb.b / 255.) < 0)
    return (1);
  return (0);
}

int colorifyPymolGradient2(index_t *idx, density_t *d, double v)
{
  double c;
  double tmp;
  double absmax;
  double scale;

  absmax = MAX(fabs(d->nmax),fabs(d->pmax));
  c = pow(fabs(v / absmax), 1. / 4.);
  if (v > 0.0)
    {
      if (fprintf(idx->fp, "slice.extend( [COLOR, %f, %f, %f] )\n", c, 0.0, 0.0) < 0)
        return (1);
    }
  else
    {
      if (fprintf(idx->fp, "slice.extend( [COLOR, %f, %f, %f] )\n", 0.0, 0.0, c) < 0)
        return (1);
    }
  return (0);
}
