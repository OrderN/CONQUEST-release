/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_math.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_math.c *
 *
 * NAME
 *  cq_ut_conqtour_math.c
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

double trilinearInterpolation(double *pt, density_t *d)
{
  int i;
  double tmp[3];
  double d00, d10, d01, d11;
  double d0, d1;
  int ref[3];

  /* Trilinear interpolation */
  for (i = 0; i < 3; ++i)
    {
      /* Adding a very small number (ZEROTOL2) favors that points that are supposed to be */
      /* at exactly a certain grid point get this as the reference */
      ref[i] = floor(fabs((pt[i] + ZEROTOL2) / d->del[i])); /* Origin of interpolating cell (min=0; max=d->pt-1) */
      if (ref[i] == d->pt[i] - 1)
        --ref[i]; /* Make point on or beyond upper limit one below limit */
      tmp[i] = (pt[i] / d->del[i]) - ref[i]; /* This, we allow to be out of bounds (i.e. extrapolation) */
    }
  d00 = (d->data[ref[0] + 1][ref[1]][ref[2]] - d->data[ref[0]][ref[1]][ref[2]]) * tmp[0] + d->data[ref[0]][ref[1]][ref[2]];
  d10 = (d->data[ref[0] + 1][ref[1] + 1][ref[2]] - d->data[ref[0]][ref[1] + 1][ref[2]]) * tmp[0]
      + d->data[ref[0]][ref[1] + 1][ref[2]];
  d01 = (d->data[ref[0] + 1][ref[1]][ref[2] + 1] - d->data[ref[0]][ref[1]][ref[2] + 1]) * tmp[0] + d->data[ref[0]][ref[1]][ref[2]
      + 1];
  d11 = (d->data[ref[0] + 1][ref[1] + 1][ref[2] + 1] - d->data[ref[0]][ref[1] + 1][ref[2] + 1]) * tmp[0] + d->data[ref[0]][ref[1]
      + 1][ref[2] + 1];
  d0 = (d10 - d00) * tmp[1] + d00;
  d1 = (d11 - d01) * tmp[1] + d01;
  return ((d1 - d0) * tmp[2] + d0);
}

double dotProduct(double *v1, double *v2)
{
  int i;
  double tmp;

  tmp = 0.0;
  for (i = 0; i < 3; ++i)
    tmp += v1[i] * v2[i];
  return (tmp);
}

void crossProduct(double *a, double *b, double *result)
{
  result[0] = a[1] * b[2] - a[2] * b[1];
  result[1] = a[2] * b[0] - a[0] * b[2];
  result[2] = a[0] * b[1] - a[1] * b[0];
}

/* Input angle in degrees */
void rotateVector(GLdouble angle, GLdouble *axis, GLdouble *v)
{
  int i;
  GLdouble unit_axis[3];
  GLdouble copy[3];
  GLdouble mod;

  angle *= DEG_TO_RAD;

  mod = 0;
  for (i = 0; i < 3; ++i) /* Normalise axis, just in case */
    {
      unit_axis[i] = axis[i];
      mod += unit_axis[i] * unit_axis[i];
      copy[i] = v[i];
    }
  mod = sqrt(mod);
  for (i = 0; i < 3; ++i)
    unit_axis[i] /= mod;

  v[0] = (1 + (1 - cos(angle)) * (pow(unit_axis[0], 2) - 1)) * copy[0] + ((1 - cos(angle)) * unit_axis[0] * unit_axis[1]
      - unit_axis[2] * sin(angle)) * copy[1] + (unit_axis[1] * sin(angle) + (1 - cos(angle)) * unit_axis[0] * unit_axis[2])
      * copy[2];
  v[1] = ((1 - cos(angle)) * unit_axis[0] * unit_axis[1] + unit_axis[2] * sin(angle)) * copy[0] + (1 + (1 - cos(angle)) * (pow(
      unit_axis[1], 2) - 1)) * copy[1] + ((1 - cos(angle)) * unit_axis[1] * unit_axis[2] - unit_axis[0] * sin(angle)) * copy[2];
  v[2] = ((1 - cos(angle)) * unit_axis[0] * unit_axis[2] - unit_axis[1] * sin(angle)) * copy[0] + ((1 - cos(angle))
      * unit_axis[1] * unit_axis[2] + unit_axis[0] * sin(angle)) * copy[1] + (1 + (1 - cos(angle)) * (pow(unit_axis[2], 2) - 1))
      * copy[2];
}

/* This function also returns the modulus of the vector before normalization */
double normalize(double *v)
{
  int i;
  double modulus;

  modulus = 0.0;
  for (i = 0; i < 3; ++i)
    modulus += pow(v[i], 2.0);
  modulus = sqrt(modulus);
  for (i = 0; i < 3; ++i)
    v[i] /= modulus;

  return (modulus);
}
