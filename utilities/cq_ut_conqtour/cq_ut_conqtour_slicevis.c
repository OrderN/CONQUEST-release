/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_slicevis.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_slicevis.c *
 *
 * NAME
 *  cq_ut_conqtour_slicevis.c
 * PURPOSE
 *  Visual representation of a density slice
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

int drawSlice(index_t *idx, density_t *d)
{
  int i, j, k;
  double **pt;
  double aux[3];
  densitycut_t *spt;
  int inipt, endpt; /* For each line, there is a central loop that must be limited */
  int peak; /* The number of the corner that closes a triangle, when needed */
  int errors = 0;

  pt = (double **) malloc(4 * sizeof(double *));
  for (i = 0; i < 4; ++i)
    pt[i] = (double *) malloc(3 * sizeof(double));
  spt = idx->dindex[idx->curr];
  //   glFlush();
  if (idx->fp == NULL) /* This means we want a visual representation */
    {
      if (glIsList(spt->slice) == GL_TRUE)
        glDeleteLists(spt->slice, 1);
      spt->slice = glGenLists(1);
      glNewList(spt->slice, GL_COMPILE);
      glBegin(GL_TRIANGLES);
      glNormal3f(idx->surfnormal[0], idx->surfnormal[1], idx->surfnormal[2]);
    }
  else /* Otherwise, we want a file representation (as a pymol script) */
    {
      if (printPymolScriptHead(idx, d))
        return (1); /* Error */
    }
  /* Last line : special */
  /* The slice could "overflow"; this doesn't happen because the steps are chosen carefully */
  for (i = 0; i < spt->pt[Y] - 1; ++i)
    {
      /* The following case is special: it usually means "corner of a triangle" */
      /* but sometimes it happens at the last row of some other polygons too */
      if (spt->lim[LEFT][i + 1] > spt->lim[RIGHT][i + 1])
        {
#ifdef DEBUG
          printf("Special case found: %d\n", idx->curr);
#endif
          if (spt->nocorners == 3) /* A triangle; let's try to fix this */
            {
              for (j = 0; j < 3; ++j)
                {
                  if (spt->refcorner[0] != j && spt->refcorner[1] != j)
                    peak = j;
                }
              for (j = 0; j < 3; ++j)
                {
                  aux[j] = spt->rectangle[0][j] + i * spt->yaxis[j];
                  pt[0][j] = aux[j] + spt->xbord[LEFT][i] * spt->xaxis[j];
                  pt[1][j] = spt->corner[peak][j];
                  pt[2][j] = aux[j] + spt->xbord[RIGHT][i] * spt->xaxis[j];
                }
              errors += drawPeak(pt, idx, d, i);
              continue;
            }
          else
            continue;
        }
      /* Each line has to be drawn in five steps */
      /* Step 1. A (usually trapezoidal) quadrangle on the left side */
      for (j = 0; j < 3; ++j)
        {
          /* Move to the origin of the line and define 4 corners */
          aux[j] = spt->rectangle[0][j] + i * spt->yaxis[j];
          pt[0][j] = aux[j] + spt->xbord[LEFT][i] * spt->xaxis[j];
          pt[1][j] = aux[j] + spt->lim[LEFT][i] * spt->xaxis[j];
          pt[2][j] = aux[j] + spt->xbord[LEFT][i + 1] * spt->xaxis[j] + spt->yaxis[j];
          pt[3][j] = aux[j] + spt->lim[LEFT][i + 1] * spt->xaxis[j] + spt->yaxis[j];
          /* We always use 0 for x, because it is were we store the limits */
          /* It is also a signal to drawQuad that this is the left border */
        }
      errors += drawQuad(pt, idx, d, 0, i);
      /* Step 2. A series of triangles, until both this line and the next are fully within the slice*/
      /* There are two possibilities */
      /*   A. The limit of the next row is to the left of the limit in this row */
      /*   B. The limit of the next row is to the right of the limit in this row */
      //printf("CUNNA %d %d\n", spt->lim[LEFT][i], spt->lim[LEFT][i+1]);
      if (spt->lim[LEFT][i] >= spt->lim[LEFT][i + 1]) /* Case A */
        {
          inipt = spt->lim[LEFT][i + 1];
          for (j = 0; j < 3; ++j)
            {
              pt[0][j] = aux[j] + spt->lim[LEFT][i] * spt->xaxis[j];
              pt[1][j] = aux[j] + (spt->lim[LEFT][i + 1] - 1) * spt->xaxis[j] + spt->yaxis[j];
            }
          while (inipt < spt->lim[LEFT][i])
            {
              for (j = 0; j < 3; ++j)
                {
                  pt[1][j] += spt->xaxis[j];
                  pt[2][j] = pt[1][j] + spt->xaxis[j];
                }
              errors += drawTriangle(pt, idx, d, spt->lim[LEFT][i], i, inipt, 1);
              ++inipt;
            }
        }
      else /* Case B */
        {
          inipt = spt->lim[LEFT][i];
          for (j = 0; j < 3; ++j)
            {
              pt[0][j] = aux[j] + (spt->lim[LEFT][i] - 1) * spt->xaxis[j];
              pt[1][j] = aux[j] + spt->lim[LEFT][i + 1] * spt->xaxis[j] + spt->yaxis[j];
            }
          while (inipt < spt->lim[LEFT][i + 1])
            {
              for (j = 0; j < 3; ++j)
                {
                  pt[0][j] += spt->xaxis[j];
                  pt[2][j] = pt[0][j] + spt->xaxis[j];
                }
              errors += drawTriangle(pt, idx, d, spt->lim[LEFT][i + 1], i + 1, inipt, -1);
              ++inipt;
            }
        }
      /* Step 3. A series of triangles, similar to Step 2, but symmetric, */
      /*         on the other side of the polygon */
      /*         NOTE: We do this here, and not after step 4, because */
      /*               this way we can calculate endpt, before drawing */
      /*               the main body of the slice (step 4) */
      if (spt->lim[RIGHT][i] < spt->lim[RIGHT][i + 1]) /* Case B (see step 2) */
        {
          endpt = spt->lim[RIGHT][i];
          for (j = 0; j < 3; ++j)
            {
              pt[0][j] = aux[j] + spt->lim[RIGHT][i] * spt->xaxis[j];
              pt[1][j] = aux[j] + (spt->lim[RIGHT][i] - 1) * spt->xaxis[j] + spt->yaxis[j];
            }
          while (endpt < spt->lim[RIGHT][i + 1])
            {
              for (j = 0; j < 3; ++j)
                {
                  pt[1][j] += spt->xaxis[j];
                  pt[2][j] = pt[1][j] + spt->xaxis[j];
                }
              errors += drawTriangle(pt, idx, d, spt->lim[RIGHT][i], i, endpt, 1);
              ++endpt;
            }
          /* But, unlike the left side, on this side the endpt must return to the limit */
          endpt = spt->lim[RIGHT][i];
        }
      else /* Case A */
        {
          endpt = spt->lim[RIGHT][i + 1];
          for (j = 0; j < 3; ++j)
            {
              pt[0][j] = aux[j] + (spt->lim[RIGHT][i + 1] - 1) * spt->xaxis[j];
              pt[1][j] = aux[j] + spt->lim[RIGHT][i + 1] * spt->xaxis[j] + spt->yaxis[j];
            }
          while (endpt < spt->lim[RIGHT][i])
            {
              for (j = 0; j < 3; ++j)
                {
                  pt[0][j] += spt->xaxis[j];
                  pt[2][j] = pt[0][j] + spt->xaxis[j];
                }
              errors += drawTriangle(pt, idx, d, spt->lim[RIGHT][i + 1], i + 1, endpt, -1);
              ++endpt;
            }
          /* But, unlike the left side, on this side the endpt must return to the limit */
          endpt = spt->lim[RIGHT][i + 1];
        }

      /* Step 4. A series of inner quadrangles, which are actually rectangles */
      for (j = 0; j < 3; ++j)
        pt[0][j] = aux[j] + inipt * spt->xaxis[j];
      //                                          + (spt->lim[LEFT][i]) * spt->xaxis[j];
      //printf("INIPT %d %d\n", inipt, endpt);
      //          for(j=spt->lim[LEFT][i]; j < spt->lim[RIGHT][i]; ++j)
      for (j = inipt; j < endpt; ++j)
        {
          for (k = 0; k < 3; ++k)
            pt[1][k] = pt[0][k] + spt->xaxis[k];
          for (k = 0; k < 3; ++k)
            pt[2][k] = pt[0][k] + spt->yaxis[k];
          for (k = 0; k < 3; ++k)
            pt[3][k] = pt[0][k] + spt->xaxis[k] + spt->yaxis[k];
          /* Draw a quadrangle */
          errors += drawQuad(pt, idx, d, j, i);
          for (k = 0; k < 3; ++k)
            pt[0][k] += spt->xaxis[k];
        }
      /* Step 5. A final quatrangle (usually trapezoidal), to cap the line */
      for (j = 0; j < 3; ++j)
        {
          /* Move to the origin of the line and define 4 corners */
          aux[j] = spt->rectangle[0][j] + i * spt->yaxis[j];
          pt[0][j] = aux[j] + spt->lim[RIGHT][i] * spt->xaxis[j];
          pt[1][j] = aux[j] + spt->xbord[RIGHT][i] * spt->xaxis[j];
          pt[2][j] = aux[j] + spt->lim[RIGHT][i + 1] * spt->xaxis[j] + spt->yaxis[j];
          pt[3][j] = aux[j] + spt->xbord[RIGHT][i + 1] * spt->xaxis[j] + spt->yaxis[j];
          /* We always use 0 for x, because it is were we store the limits */
          /* It is also a signal to drawQuad that this is the left border */
        }
      errors += drawQuad(pt, idx, d, spt->pt[X] - 1, i);
    }
  /* Draw the corners at the two potential turns */
  for (i = LEFT; i <= RIGHT; ++i)
    if (spt->turn[i][LINE] != -1)
      {
        for (j = 0; j < 3; ++j)
          {
            aux[j] = spt->rectangle[0][j] + spt->turn[i][LINE] * spt->yaxis[j];
            pt[0][j] = aux[j] + spt->xbord[i][spt->turn[i][LINE]] * spt->xaxis[j];
            pt[1][j] = spt->corner[spt->turn[i][CORNER]][j];
            pt[2][j] = aux[j] + spt->xbord[i][spt->turn[i][LINE] + 1] * spt->xaxis[j] + spt->yaxis[j];
          }
        /* Remember that the data for borders is in data[0] & data[pt[X]-1] (see steps 1 & 5) */
        errors += drawCorner(pt, idx, d, (i) ? spt->pt[X] - 1 : 0, spt->turn[i][LINE], i);
      }
  if (idx->fp == NULL) /* This means we want a visual representation */
    {
      glEnd();
      glEndList();
    }
  else /* Otherwise, we want a file representation (as a pymol script) */
    {
      if (errors)
        return (1);
      if (printPymolScriptTail(idx, d))
        return (1);
    }
  for (i = 0; i < 4; ++i)
    free(pt[i]);
  free(pt);
  requestRedisplay();
  return 0; /* No problems found */
}

/*
 *  pt -> Pointer to four points at the corners of the QUAD
 */
int drawQuad(double **pt, index_t *idx, density_t *d, int x, int y)
{
  int i;
  int **c; /* This is used for translating between requested x,y and the actual data array indexes */
  densitycut_t *spt;
  int result; /* To take note whether there was an error */

  c = (int **) malloc(4 * sizeof(int *));
  for (i = 0; i < 4; ++i)
    c[i] = (int *) malloc(2 * sizeof(int));
  spt = idx->dindex[idx->curr];
  /* Usually the four corners are adjacent in the data array */
  /* but the borders are special cases */
  if (x == 0) /* This is at the left border */
    {
      c[0][X] = x;
      c[0][Y] = y;
      c[1][X] = spt->lim[LEFT][y];
      c[1][Y] = y;
      c[2][X] = x;
      c[2][Y] = y + 1;
      c[3][X] = spt->lim[LEFT][y + 1];
      c[3][Y] = y + 1;
    }
  else if (x == spt->pt[X] - 1) /* Right border */
    {
      c[0][X] = spt->lim[RIGHT][y];
      c[0][Y] = y;
      c[1][X] = x;
      c[1][Y] = y;
      c[2][X] = spt->lim[RIGHT][y + 1];
      c[2][Y] = y + 1;
      c[3][X] = x;
      c[3][Y] = y + 1;
    }
  else /* Usual rectangle */
    {
      c[0][X] = x;
      c[0][Y] = y;
      c[1][X] = x + 1;
      c[1][Y] = y;
      c[2][X] = x;
      c[2][Y] = y + 1;
      c[3][X] = x + 1;
      c[3][Y] = y + 1;
    }
  if (idx->fp == NULL) /* Visual */
    {
      colorifyGradient(idx, d, c[3][X], c[3][Y]);
      glVertex3f(pt[3][0], pt[3][1], pt[3][2]);
      colorifyGradient(idx, d, c[2][X], c[2][Y]);
      glVertex3f(pt[2][0], pt[2][1], pt[2][2]);
      colorifyGradient(idx, d, c[0][X], c[0][Y]);
      glVertex3f(pt[0][0], pt[0][1], pt[0][2]);

      colorifyGradient(idx, d, c[3][X], c[3][Y]);
      glVertex3f(pt[3][0], pt[3][1], pt[3][2]);
      colorifyGradient(idx, d, c[0][X], c[0][Y]);
      glVertex3f(pt[0][0], pt[0][1], pt[0][2]);
      colorifyGradient(idx, d, c[1][X], c[1][Y]);
      glVertex3f(pt[1][0], pt[1][1], pt[1][2]);

      result = 0; /* Everything ok */
    }
  else /* Pymol script */
    {
      result = drawPymolQuad(pt, c, idx, d);
    }
  for (i = 0; i < 4; ++i)
    free(c[i]);
  free(c);

  return (result);
}

/*
 *  pt -> Pointer to four points at the corners of the QUAD
 *  sign -> Indicates whether the "pivoting point" of the triangles is
 *          in the next line (+1) or in the previous one (-1)
 */
int drawTriangle(double **pt, index_t *idx, density_t *d, int x, int y, int ini, int sign)
{
  int i;
  int **c; /* This is used for translating between requested x,y and the actual data array indexes */
  int result; /* To take note whether there was an error */

  c = (int **) malloc(3 * sizeof(int *));
  for (i = 0; i < 3; ++i)
    c[i] = (int *) malloc(2 * sizeof(int));

  if (sign > 0)
    {
      c[0][X] = x;
      c[0][Y] = y;
      c[1][X] = ini;
      c[1][Y] = y + 1;
      c[2][X] = ini + 1;
      c[2][Y] = y + 1;
    }
  else
    {
      c[0][X] = ini;
      c[0][Y] = y - 1;
      c[1][X] = x;
      c[1][Y] = y;
      c[2][X] = ini + 1;
      c[2][Y] = y - 1;
    }
  if (idx->fp == NULL) /* Visual */
    {
      colorifyGradient(idx, d, c[0][X], c[0][Y]);
      glVertex3f(pt[0][0], pt[0][1], pt[0][2]);
      colorifyGradient(idx, d, c[1][X], c[1][Y]);
      glVertex3f(pt[1][0], pt[1][1], pt[1][2]);
      colorifyGradient(idx, d, c[2][X], c[2][Y]);
      glVertex3f(pt[2][0], pt[2][1], pt[2][2]);

      result = 0;
    }
  else
    {
      result = drawPymolTriangle(pt, c, idx, d);
    }
  for (i = 0; i < 3; ++i)
    free(c[i]);
  free(c);

  return (result);
}

int drawCorner(double **pt, index_t *idx, density_t *d, int x, int y, int side)
{
  densitycut_t *spt;
  int result;

  if (idx->fp == NULL) /* Visual */
    {
      spt = idx->dindex[idx->curr];
      colorifyGradient(idx, d, x, y);
      glVertex3f(pt[0][0], pt[0][1], pt[0][2]);
      colorifyGradient2(idx, d, spt->cturn[side]);
      glVertex3f(pt[1][0], pt[1][1], pt[1][2]);
      colorifyGradient(idx, d, x, y + 1);
      glVertex3f(pt[2][0], pt[2][1], pt[2][2]);

      result = 0;
    }
  else
    {
      result = drawPymolCorner(pt, x, y, side, idx, d);
    }

  return (result);
}

int drawPeak(double **pt, index_t *idx, density_t *d, int y)
{
  densitycut_t *spt;
  int result;

  if (idx->fp == NULL) /* Visual */
    {
      spt = idx->dindex[idx->curr];
      colorifyGradient(idx, d, 0, y); /* Note this means the left border color */
      glVertex3f(pt[0][0], pt[0][1], pt[0][2]);
      colorifyGradient2(idx, d, spt->cturn[2]);
      glVertex3f(pt[1][0], pt[1][1], pt[1][2]);
      colorifyGradient(idx, d, spt->pt[X] - 1, y); /* Note this means the right border color */
      glVertex3f(pt[2][0], pt[2][1], pt[2][2]);

      result = 0;
    }
  else
    {
      result = drawPymolPeak(pt, y, idx, d);
    }

  return (result);
}

void colorifyGradient(index_t *idx, density_t *d, int x, int y)
{
  rgbcolor_t rgb;

  getColorFromScale(d, idx->dindex[idx->curr]->data[x][y], &rgb);
  glColor3f(rgb.r / 255., rgb.g / 255., rgb.b / 255.);
}

void colorifyGradient2(index_t *idx, density_t *d, double v)
{
  double c;
  double tmp;
  double absmax;
  double scale;

  absmax = MAX(fabs(d->nmax),fabs(d->pmax));
  c = pow(fabs(v / absmax), 1. / 4.);
  //   c=pow(fabs(v/absmax),1./1.);
  if (v > 0.0)
    glColor3f(c, 0, 0);
  else
    glColor3f(0, 0, c);
}

