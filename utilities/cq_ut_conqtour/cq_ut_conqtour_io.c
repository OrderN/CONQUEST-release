/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_io.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_io.c *
 *
 * NAME
 *  cq_ut_conqtour_io.c
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

#include <string.h>
#include "cq_ut_conqtour.h"

#define CQ_USING_CQINFO  0x01
#define CQ_USING_LIMITS  0x02
#define CQ_USING_DENSITY 0x04
#define CQ_USING_MAXIMA  0x08

#define MAXUNKNOWN 100 /* Maximum number of unknown atom types that can be found in a file */

typedef struct
{
  int z; /* Atomic number */
  char *symbol; /* Atomic symbol */
  double weight; /* Atomic weight */
} atominfo_t;

static atominfo_t periodictable[] =
  {
    { 1, "H", 1.01 },
    { 2, "He", 4.00 },
    { 3, "Li", 6.94 },
    { 4, "Be", 9.01 },
    { 5, "B", 10.81 },
    { 6, "C", 12.01 },
    { 7, "N", 14.01 },
    { 8, "O", 15.99 },
    { 9, "F", 18.99 },
    { 10, "Ne", 20.18 },
    { 11, "Na", 22.99 },
    { 12, "Mg", 24.31 },
    { 13, "Al", 26.98 },
    { 14, "Si", 28.09 },
    { 15, "P", 30.97 },
    { 16, "S", 32.06 },
    { 17, "Cl", 35.45 },
    { 18, "Ar", 39.95 },
    { 19, "K", 39.10 },
    { 20, "Ca", 40.08 },
    { 21, "Sc", 44.96 },
    { 22, "Ti", 47.87 },
    { 23, "V", 50.94 },
    { 24, "Cr", 51.99 },
    { 25, "Mn", 54.94 },
    { 26, "Fe", 55.85 },
    { 27, "Co", 58.93 },
    { 28, "Ni", 58.69 },
    { 29, "Cu", 63.55 },
    { 30, "Zn", 65.39 },
    { 31, "Ga", 69.72 },
    { 32, "Ge", 72.64 },
    { 33, "As", 74.92 },
    { 34, "Se", 78.96 },
    { 35, "Br", 79.90 },
    { 36, "Kr", 83.80 },
    { 37, "Rb", 85.47 },
    { 38, "Sr", 87.62 },
    { 39, "Y", 88.91 },
    { 40, "Zr", 91.22 },
    { 41, "Nb", 92.91 },
    { 42, "Mo", 95.94 },
    { 43, "Tc", 98.00 },
    { 44, "Ru", 101.07 },
    { 45, "Rh", 102.91 },
    { 46, "Pd", 106.42 },
    { 47, "Ag", 107.87 },
    { 48, "Cd", 112.41 },
    { 49, "In", 114.82 },
    { 50, "Sn", 118.71 },
    { 51, "Sb", 121.76 },
    { 52, "Te", 127.60 },
    { 53, "I", 126.90 },
    { 54, "Xe", 131.29 },
    { 55, "Cs", 132.91 },
    { 56, "Ba", 137.33 },
    { 57, "La", 138.91 },
    { 58, "Ce", 140.12 },
    { 59, "Pr", 140.91 },
    { 60, "Nd", 144.24 },
    { 61, "Pm", 145.00 },
    { 62, "Sm", 150.36 },
    { 63, "Eu", 151.96 },
    { 64, "Gd", 157.25 },
    { 65, "Tb", 158.93 },
    { 66, "Dy", 162.50 },
    { 67, "Ho", 164.93 },
    { 68, "Er", 167.26 },
    { 69, "Tm", 168.93 },
    { 70, "Yb", 173.04 },
    { 71, "Lu", 174.97 },
    { 72, "Hf", 178.49 },
    { 73, "Ta", 180.95 },
    { 74, "W", 183.84 },
    { 75, "Re", 186.21 },
    { 76, "Os", 190.23 },
    { 77, "Ir", 192.22 },
    { 78, "Pt", 195.08 },
    { 79, "Au", 196.97 },
    { 80, "Hg", 200.59 },
    { 81, "Tl", 204.38 },
    { 82, "Pb", 207.20 },
    { 83, "Bi", 208.98 },
    { 84, "Po", 209.00 },
    { 85, "At", 210.00 },
    { 86, "Rn", 222.00 },
    { 87, "Fr", 223.00 },
    { 88, "Ra", 226.00 },
    { 89, "Ac", 227.00 },
    { 90, "Th", 232.04 },
    { 91, "Pa", 231.04 },
    { 92, "U", 238.03 },
    { 93, "Np", 237.00 },
    { 94, "Pu", 244.00 },
    { 95, "Am", 243.00 },
    { 96, "Cm", 247.00 },
    { 97, "Bk", 247.00 },
    { 98, "Cf", 251.00 },
    { 99, "Es", 252.00 },
    { 100, "Fm", 257.00 },
    { 101, "Md", 258.00 },
    { 102, "No", 259.00 },
    { 103, "Lr", 262.00 },
    { 104, "Rf", 261.00 },
    { 105, "Db", 262.00 },
    { 106, "Sg", 266.00 },
    { 107, "Bh", 264.00 },
    { 108, "Hs", 277.00 },
    { 109, "Mt", 268.00 },
    { 0, NULL, 0.0 } };
//static atominfo_t atomweights[] = "|H|1.01|He|4.00|"
//  "Li|6.94|Be|9.01|B|10.81|C|12.01|N|14.01|O|15.99|F|18.99|Ne|20.18|"
//  "Na|22.99|Mg|24.31|Al|26.98|Si|28.09|P|30.97|S|32.06|Cl|35.45|Ar|39.95|"
//  "K|39.10|Ca|40.08|Sc|44.96|Ti|47.87|V|50.94|Cr|51.99|Mn|54.94|Fe|55.85|"
//  "Co|58.93|Ni|58.69|Cu|63.55|Zn|65.39|Ga|69.72|Ge|72.64|As|74.92|Se|78.96|Br|79.90|Kr|83.80|"
//  "Rb|85.47|Sr|87.62|Y|88.91|Zr|91.22|Nb|92.91|Mo|95.94|Tc|98.00|Ru|101.07|"
//  "Rh|102.91|Pd|106.42|Ag|107.87|Cd|112.41|In|114.82|Sn|118.71|Sb|121.76|Te|127.60|I|126.90|Xe|131.29|"
//  "Cs|132.91|Ba|137.33|La|138.91|Ce|140.12|Pr|140.91|Nd|144.24|Pm|145.00|Sm|150.36|Eu|151.96|Gd|157.25|"
//  "Tb|158.93|Dy|162.50|Ho|164.93|Er|167.26|Tm|168.93|Yb|173.04|Lu|174.97|Hf|178.49|Ta|180.95|W|183.84|"
//  "Re|186.21|Os|190.23|Ir|192.22|Pt|195.08|Au|196.97|Hg|200.59|Tl|204.38|Pb|207.20|Bi|208.98|Po|209.00|At|210.00|Rn|222.00|"
//  "Fr|223.00|Ra|226.00|Ac|227.00|Th|232.04|Pa|231.04|U|238.03|Np|237.00|Pu|244.00|Am|243.00|Cm|247.00|Bk|247.00|Cf|251.00|"
//  "Es|252.00|Fm|257.00|Md|258.00|No|259.00|Lr|262.00|Rf|261.00|Db|262.00|Sg|266.00|Bh|264.00|Hs|277.00|Mt|268.00|";

// TODO Pass control and limit structures

// TODO Check errors carefully
int readXplorDensity(char *name, density_t *d)
{
  char line[MAXLIN];
  int di; /* Dummy integer */
  double dd; /* Dummy double */
  int i, j, k;
  int i2, j2, k2; /* Counters for the actual indexes in the loaded matrix */
  int old[3];
  FILE *fp;
  //  int shifted[3]; /* A point after applying for the shift */
  int end[3]; /* The last point within limits */
  int inlimits[3];
  int ptinfile[3];

  for (i = 0; i < 3; ++i)
    old[i] = d->pt[i];
  if ((fp = fopen(name, "r")) == NULL)
    {
      if (ctrl.verbose)
        printf("I can't open charge file '%s'\n", name);
      return (1); /* File not found */
    }

  fgets(line, MAXLIN, fp); /* Ignore 4 lines */
  fgets(line, MAXLIN, fp);
  fgets(line, MAXLIN, fp);
  fgets(line, MAXLIN, fp);
  fgets(line, MAXLIN, fp); /* Read grid points */
  sscanf(line, "%d %d %d %d %d %d %d %d %d", &d->pt[0], &di, &di, &d->pt[1], &di, &di, &d->pt[2], &di, &di);//
  /* Decide whether we are within limits */
  // TODO IMP Check that the griddims are not negative and that there are at least 2 grid points
  // TODO Do this kind of thing in a function, because it is reused in other functions
  for (i = 0; i < 3; ++i)
    {
      /* Strides must be non-negative; if they are not, set to 1 */
      if (limits.strides[i] <= 0)
        limits.strides[i] = 1;
      if (limits.griddims[i] > 0) /* A grid window is given */
        {
          /* But the number of points in the file is not enough: Adjust it */
          if ((double) ((d->pt[i] - limits.gridshft[i]) / limits.strides[i]) < limits.griddims[i])
            {
              if (ctrl.verbose)
                fprintf(stderr, "WARNING: The requested grid limits are too large for the file (coord. %s). "
                  "Adjusting to maximum\n", REPORT_COORDINATE(i));
              limits.griddims[i] = (int) ((double) (d->pt[i] - limits.gridshft[i]) / limits.strides[i]);
              if (limits.griddims[i] < 2)
                {
                  limits.griddims[i] = 2;
                  if (ctrl.verbose)
                    fprintf(stderr, "WARNING: Your limits imply no points will be loaded for coord. %s. "
                      "Adjusting to 2 points. Check your limits before loading\n", REPORT_COORDINATE(i));
                }
            }
        }
      else
        {
          /* There was no grid window: Load as much as possible, taking into account the stride and grid shift*/
          limits.griddims[i] = (int) ((double) (d->pt[i] - limits.gridshft[i]) / limits.strides[i]);
          if (limits.griddims[i] < 2)
            {
              limits.griddims[i] = 2;
              if (ctrl.verbose)
                fprintf(stderr, "WARNING: Your limits imply no points will be loaded for coord. %s. "
                  "Adjusting to 2 points. Check your limits before loading\n", REPORT_COORDINATE(i));
            }

        }
    }
#ifdef DEBUG
  printf("Read grid: %d %d %d\n", limits.griddims[0], limits.griddims[1], limits.griddims[2]);
  printf("Grid points = %d %d %d\n", d->pt[0], d->pt[1], d->pt[2]);
#endif
  fgets(line, MAXLIN, fp); /* Read lattice data */
  sscanf(line, "%lf %lf %lf %lf %lf %lf", &d->lvm[0], &d->lvm[1], &d->lvm[2], &d->angle[0], &d->angle[1], &d->angle[2]);
#ifdef DEBUG
  printf("Lattice = %f %f %f\n", d->lvm[0], d->lvm[1], d->lvm[2]);
#endif
  if (d->angle[0] < 89.9 || d->angle[0] > 90.1 || d->angle[1] < 89.9 || d->angle[1] > 90.1 || d->angle[2] < 89.9 || d->angle[2]
      > 90.1)
    {
      if (ctrl.verbose)
        printf("Sorry, but we only handle orthorhombic cells for the moment\n");
      return (2); /* Wrong file format */
    }
  /* Fix the density and limits structures */
  for (i = 0; i < 3; ++i)
    {
      ptinfile[i] = d->pt[i];
      /* Calculate the spacing between grid points */
      d->del[i] = d->lvm[i] / (d->pt[i] - 1);
      /* If the real shift is 0 at input, set it to match the input grid shift */
      /* Otherwise, leave it as it was (the user "must" have a reason) */
      if (limits.realshft[i] < ZEROTOL)
        {
          limits.realshft[i] = d->del[i] * limits.gridshft[i];
        }
      /* Correct spacing by grid strides */
      d->del[i] *= limits.strides[i];
      /* Recalculate the lattice vectors (this will be the correct value for the loaded data) */
      d->lvm[i] = d->del[i] * (limits.griddims[i] - 1);
      /* Set the real window so that it matches the grid window */
      limits.realdims[i] = d->lvm[i];
      /* And reset the number of density points (to match what will be loaded) */
      d->pt[i] = limits.griddims[i];
    }

  fgets(line, MAXLIN, fp);
  if (strstr(line, "ZYX") == NULL)
    {
      if (ctrl.verbose)
        printf("Only ZXY-ordered files are handled\n");
      return (2); /* Wrong format */
    }
  // TODO Store an informative name for the density
  resetDensity(d, old); /* From here, the previous data are lost (hence, no display) */
#ifdef DEBUG
  printf("Allocated %d %d %d\n", d->pt[0], d->pt[1], d->pt[2]);
#endif
  /* Finally, read the data */
  d->pmin = VERYLARGE;
  d->pmax = ZERO;
  d->nmin = VERYLARGE;
  d->nmax = ZERO;
  i2 = -1; /* Start "window" counter */
  for (i = 0; i < 3; ++i)
    end[i] = limits.gridshft[i] + limits.griddims[i] * limits.strides[i];
  for (i = 0; i < ptinfile[2]; ++i)
    {
      //      shifted[2] = i - limits.gridshft[2];
      if (i >= limits.gridshft[2] && i < end[2] && !((i - limits.gridshft[2]) % limits.strides[2]))
        {
          inlimits[2] = TRUE;
          ++i2;
        }
      else
        inlimits[2] = FALSE;
      fscanf(fp, "%d", &di);
      if (di != i)
        {
          if (ctrl.verbose)
            printf("Block error in xplor file\n");
          return (3); /* Incomplete file */
        }
      j2 = -1;
      for (j = 0; j < ptinfile[1]; ++j)
        {
          //          shifted[1] = j - limits.gridshft[1];
          if (j >= limits.gridshft[1] && j < end[1] && !((j - limits.gridshft[1]) % limits.strides[1]))
            {
              inlimits[1] = TRUE;
              ++j2;
            }
          else
            inlimits[1] = FALSE;
          k2 = -1;
          for (k = 0; k < ptinfile[0]; ++k)
            {
              //              shifted[0] = k - limits.gridshft[0];
              if (k >= limits.gridshft[0] && k < end[0] && !((k - limits.gridshft[0]) % limits.strides[0]))
                {
                  inlimits[0] = TRUE;
                  ++k2;
                }
              else
                inlimits[0] = FALSE;
              //              fscanf(fp, "%lf", &d->data[k][j][i]);
              fscanf(fp, "%lf", &dd);
              //if(d->data[k][j][i] < 0.0) d->data[k][j][i] = 0.0;
              //              if(i >= limits.gridshft[2] && j >= limits.gridshft[1] && k >= limits.gridshft[0] &&
              //                  i <= end[2] && j <= end[1] && k <= end[0] &&
              //                  !((i - limits.gridshft[2]) % limits.strides[2]) &&
              //                  !((j - limits.gridshft[1]) % limits.strides[1]) &&
              //                  !((k - limits.gridshft[0]) % limits.strides[0]))
              if (inlimits[0] && inlimits[1] && inlimits[2])
                {
                  //                  printf("%d %d %d = %e\n", k2,j2,i2, dd);
                  d->data[k2][j2][i2] = dd;
                  //                  if (d->data[k2][j2][i2] > 0.0)
                  //                    {
                  //                      if (d->pmin > d->data[k2][j2][i2])
                  //                        d->pmin = d->data[k2][j2][i2];
                  //                      if (d->pmax < d->data[k2][j2][i2])
                  //                        d->pmax = d->data[k2][j2][i2];
                  //                    }
                  //                  else if (d->data[k2][j2][i2] < 0.0) /* We explicitly ignore 0.0 */
                  //                    {
                  //                      if (d->nmin > fabs(d->data[k2][j2][i2]))
                  //                        d->nmin = fabs(d->data[k2][j2][i2]);
                  //                      if (d->nmax < fabs(d->data[k2][j2][i2]))
                  //                        d->nmax = fabs(d->data[k2][j2][i2]);
                  //                    }
                }
              // I do this for the moment to keep the same color scale for all limits */
              if (dd > 0.0)
                {
                  if (d->pmin > dd)
                    d->pmin = dd;
                  if (d->pmax < dd)
                    d->pmax = dd;
                }
              else if (dd < 0.0) /* We explicitly ignore 0.0 */
                {
                  if (d->nmin > fabs(dd))
                    d->nmin = fabs(dd);
                  if (d->nmax < fabs(dd))
                    d->nmax = fabs(dd);
                }
            }
        }
    }
#ifdef DEBUG
  printf("Pos : Min = %e   Max = %e\n", d->pmin, d->pmax);
  printf("Neg : Min = %e   Max = %e\n", d->nmin, d->nmax);
#endif

  /* Reset some of the grid limits to make them more intuitive */
  for (i = 0; i < 3; ++i)
    {
      limits.gridshft[i] = 0;
      limits.strides[i] = 1;
    }
  return 0; /* Everything was fine */
}

int writeXplorDensity(char *name, density_t *d, limits_t *l)
{
  int i, j;
  int x, y, z;
  FILE *fp;
  int sectsize; /* The size of each "block" of the xplor file */

  if ((fp = fopen(name, "w")) == NULL)
    {
      fprintf(stderr, "I can't write to file '%s'\n", name);
      return 1; /* Error */
    }
  fprintf(fp, "\n%8d\n", 2);
  fprintf(fp, "REMARKS FILENAME=\"%s\"\n", name);
  fprintf(fp, "REMARKS X-PLOR Density map\n");
  /* Decide whether we are within limits */
  // TODO IMP Check that the griddims are not negative and that there are at least 2 grid points
  // TODO Do this kind of thing in a function, because it is reused in other functions
  for (i = 0; i < 3; ++i)
    {
      /* Strides must be non-negative; if they are not, set to 1 */
      if (l->strides[i] <= 0)
        l->strides[i] = 1;
      if (l->griddims[i] > 0) /* A grid window is given */
        {
          /* But the number of points is not enough: Adjust it */
          if ((double) ((d->pt[i] - l->gridshft[i]) / l->strides[i]) < l->griddims[i])
            {
              if (ctrl.verbose)
                fprintf(stderr, "WARNING: The requested grid limits are too large for the file (coord. %s). "
                  "Adjusting to maximum\n", REPORT_COORDINATE(i));
              l->griddims[i] = (int) ((double) (d->pt[i] - l->gridshft[i]) / l->strides[i]);
              if (l->griddims[i] < 2)
                {
                  l->griddims[i] = 2;
                  if (ctrl.verbose)
                    fprintf(stderr, "WARNING: Your limits imply no points will be loaded for coord. %s. "
                      "Adjusting to 2 points. Check your limits before loading\n", REPORT_COORDINATE(i));
                }
            }
        }
      else
        {
          /* There was no grid window: Use as much as possible, taking into account the stride and grid shift*/
          l->griddims[i] = (int) ((double) (d->pt[i] - l->gridshft[i]) / l->strides[i]);
          if (l->griddims[i] < 2)
            {
              l->griddims[i] = 2;
              if (ctrl.verbose)
                fprintf(stderr, "WARNING: Your limits imply no points will be loaded for coord. %s. "
                  "Adjusting to 2 points. Check your limits before loading\n", REPORT_COORDINATE(i));
            }
        }
      /* Make the real window match the grid window */
      l->realdims[i] = d->lvm[i] * (l->griddims[i] - 1) / (d->pt[i] - 1);
    }
  /* Take note of the number of grid points that will be actually written */
  fprintf(fp, "% 8d% 8d% 8d% 8d% 8d% 8d% 8d% 8d% 8d\n", l->griddims[0], 0, l->griddims[0] - 1, l->griddims[1], 0, l->griddims[1]
      - 1, l->griddims[2], 0, l->griddims[2] - 1);
  /* Cell info; the angles are hardcoded bc. of Conquest's limitation to orthorhombic */
  fprintf(fp, "% 12.5E% 12.5E% 12.5E% 12.5E% 12.5E% 12.5E\n", l->realdims[0], l->realdims[1], l->realdims[2], 90.0, 90.0, 90.0);
  fprintf(fp, "ZYX\n");

  sectsize = l->griddims[0] * l->griddims[1];
  z = l->gridshft[2];
  for (i = 0; i < l->griddims[2]; ++i)
    {
      y = l->gridshft[1] - l->strides[1]; /* Prepare for the first y line */
      fprintf(fp, "%8d", i);
      for (j = 0; j < sectsize; ++j)
        {
          if (j % l->griddims[0] == 0) /* This is the beginning of an x line */
            {
              x = l->gridshft[0];
              y += l->strides[1];
            }
          if (j % 6 == 0)
            fprintf(fp, "\n");
          fprintf(fp, "% 12.5E", d->data[x][y][z]);
          x += l->strides[0];
        }
      fprintf(fp, "\n");
      z += l->strides[2];
    }
  fclose(fp);

  /* If the real shift is unset, make it match the grid shift, so that if now we want to */
  /* write a coordinates file, the output will be shifted to match the piece of density we just wrote */
  for (i = 0; i < 3; ++i)
    {
      /* If the real shift is 0 at input, set it to match the input grid shift */
      /* Otherwise, leave it as it was (the user "must" have a reason) */
      if (limits.realshft[i] < ZEROTOL)
        {
          limits.realshft[i] = d->del[i] * limits.gridshft[i];
        }
    }
  /* Reset some of the grid limits to make them more intuitive */
  /* This is a little bit like a principle of "limits can only be used once" */
  for (i = 0; i < 3; ++i)
    {
      limits.gridshft[i] = 0;
      limits.strides[i] = 1;
    }
}

int readCqDensity(density_t *d)
{
  int i, j, k;
  int old[3];
  FILE *bf, *df;
  int blno[3]; /* Number of blocks per dimension (according to blocks file) */
  int procno; /* Number of processes (according to blocks file) */
  char tmpfile[MAXLIN];
  int noblocks; /* Number of blocks in a given density file */
  int blstart[3]; /* Starting point of the block */
  int blidx; /* Block index */
  int x, y, z; /* Indexes to iterate over a block */
  int idx[3]; /* Index of a particular density point (in the files) */
  int end[3]; /* End-point that defines the window to be loaded in memory (along with the shift) */
  int pos[3]; /* Position of a point relative to the shift */
  double val;

  //TODO Check availability here. For the moment, we assume below that everything is ok

  /* Use of limits:
   * -Grid window = Try to load these many points
   * -Grid shift = Skip the first shift points (ignoreing strides)
   * -Strides = Read every stride points
   * -Real window = Ignored
   * -Real shift = Ignored
   */
  // TODO Check that the grid window is not too large
  for (i = 0; i < 3; ++i)
    {
      old[i] = d->pt[i]; /* Store old window (for deallocation) */
      d->pt[i] = limits.griddims[i]; /* Read a specified window */
    }
  resetDensity(d, old);
  if (FLAG_UP(cqinfo.read, CQ_READ_FROM_BLOCKS)) /* Use a blocks file to rearrange data */
    {
      /* Try to open the blocks file */
      bf = fopen(cqinfo.blocksfl, "r");
      if (bf == NULL)
        {
          if (ctrl.verbose)
            fprintf(stderr, "Error opening blocks file '%s'\n", cqinfo.blocksfl);
          SET_FLAG(ctrl.has, CQ_CTRL_HAS_DENSITY, FALSE); /* Failed to load the density */
          return 1;
        }
      /* Read the number of blocks per dimension and check for compatibility */
      for (i = 0; i < 3; ++i)
        {
          fscanf(bf, "%d", &blno[i]);
          if (cqinfo.pt[i] != cqinfo.blockpt[i] * blno[i])
            {
              if (ctrl.verbose)
                {
                  fprintf(stderr, "Error: Blocks for %s in '%s' (%d) do not match the expected value (%d)\n",
                      REPORT_COORDINATE(i), cqinfo.blocksfl, blno[i], cqinfo.pt[i] / cqinfo.blockpt[i]);
                  fprintf(stderr, "This probably means that the files are nor from the same calculation\n");
                }
              SET_FLAG(ctrl.has, CQ_CTRL_HAS_DENSITY, FALSE); /* Failed to load the density */
              return 1;
            }
        }
      /* Read the number of processes */
      fscanf(bf, "%d", &procno);
      if (procno != cqinfo.cores)
        {
          if (ctrl.verbose)
            {
              fprintf(stderr, "Error: Processes in '%s' (%d) do not match the expected value (%d)\n", cqinfo.blocksfl, procno,
                  cqinfo.cores);
              fprintf(stderr, "This probably means that the files are nor from the same calculation\n");
            }
          SET_FLAG(ctrl.has, CQ_CTRL_HAS_DENSITY, FALSE); /* Failed to load the density */
          return 1;
        }
      /* Read the density files */
      //     for(i=0; i < 3; ++i)      end[i] = limits.gridshft[i] + limits.griddims[i];
      d->pmin = VERYLARGE;
      d->pmax = ZERO;
      d->nmin = VERYLARGE;
      d->nmax = ZERO;
      for (i = 0; i < 3; ++i)
        end[i] = limits.griddims[i] * limits.strides[i];
      for (i = 1; i <= procno; ++i)
        {
          makeFilename(tmpfile, &cqinfo, i, EXT_NONE);
          if (ctrl.verbose)
            if (i % 1000 == 0)
              {
                printf("M");
                fflush(stdout);
              }
            else if (i % 500 == 0)
              {
                printf("D");
                fflush(stdout);
              }
            else if (i % 100 == 0)
              {
                printf("C");
                fflush(stdout);
              }
            else if (i % 50 == 0)
              {
                printf("L");
                fflush(stdout);
              }
            else if (i % 10 == 0)
              {
                printf("|");
                fflush(stdout);
              }
            else
              {
                printf("-");
                fflush(stdout);
              }
          df = fopen(tmpfile, "r");
          if (df == NULL)
            {
              if (ctrl.verbose)
                fprintf(stderr, "Error opening density file '%s'\n", tmpfile);
              SET_FLAG(ctrl.has, CQ_CTRL_HAS_DENSITY, FALSE); /* Failed to load the density */
              return 1;
            }
          //TODO Alternative way of reading from bz2 files, so they don't need to be uncompressed
          /* Read number of blocks in this process */
          fscanf(bf, "%*d %d %*d", &noblocks);
          for (j = 0; j < noblocks; ++j) /* Iterate over block indexes */
            {
              if (fscanf(bf, "%*d %d\n", &blidx) != 1)
                {
                  if (ctrl.verbose)
                    {
                      printf("\n");
                      fprintf(stderr, "Error reading block indexes\n");
                    }
                  SET_FLAG(ctrl.has, CQ_CTRL_HAS_DENSITY, FALSE); /* Failed to load the density */
                  return 1;
                }
              /* Starting point */
              blstart[0] = (blidx - 1) / (blno[1] * blno[2]);
              blstart[1] = (blidx - 1 - blstart[0] * blno[1] * blno[2]) / blno[2];
              blstart[2] = (blidx - 1) - blstart[0] * blno[1] * blno[2] - blstart[1] * blno[2];
              //printf("-----%d | %d %d %d\n", blidx, blstart[0], blstart[1], blstart[2]);
              for (z = 0; z < cqinfo.blockpt[2]; ++z)
                {
                  for (y = 0; y < cqinfo.blockpt[1]; ++y)
                    {
                      for (x = 0; x < cqinfo.blockpt[0]; ++x)
                        {
                          idx[0] = blstart[0] * cqinfo.blockpt[0] + x;
                          idx[1] = blstart[1] * cqinfo.blockpt[1] + y;
                          idx[2] = blstart[2] * cqinfo.blockpt[2] + z;
                          //printf("%d %d %d | %d %d %d = ", x,y,z,idx[0], idx[1], idx[2]);
                          if (fscanf(df, "%lf", &val) != 1)
                            {
                              if (ctrl.verbose)
                                {
                                  printf("\n");
                                  fprintf(stderr, "Error reading density data, in file '%s'\n", tmpfile);
                                }
                              SET_FLAG(ctrl.has, CQ_CTRL_HAS_DENSITY, FALSE); /* Failed to load the density */
                              return 1;
                            }
                          val /= pow(BOHR_TO_ANGSTROM, 3.0); /* In memory, units are always electrons/(angstr.)**3 */
                          for (k = 0; k < 3; ++k)
                            pos[k] = idx[k] - limits.gridshft[k];
                          /* Is this point within the limits? */
                          /*
                           if(idx[0] >= limits.gridshft[0] && idx[0] < end[0]
                           && idx[1] >= limits.gridshft[1] && idx[1] < end[1]
                           && idx[2] >= limits.gridshft[2] && idx[2] < end[2])
                           */
                          if (pos[0] >= 0 && pos[0] < end[0] && pos[0] % limits.strides[0] == 0 && pos[1] >= 0 && pos[1] < end[1]
                              && pos[1] % limits.strides[1] == 0 && pos[2] >= 0 && pos[2] < end[2] && pos[2] % limits.strides[2]
                              == 0)
                            {
                              d->data[pos[0] / limits.strides[0]][pos[1] / limits.strides[1]][pos[2] / limits.strides[2]] = val;
                              //printf("Si = %e\n", d->data[idx[0]-limits.gridshft[0]]
                              //                            [idx[1]-limits.gridshft[1]]
                              //                            [idx[2]-limits.gridshft[2]]);
                              if (val > 0.0)
                                {
                                  if (d->pmin > val)
                                    d->pmin = val;
                                  if (d->pmax < val)
                                    d->pmax = val;
                                }
                              else if (val < 0.0) /* We explicitly ignore 0.0 */
                                {
                                  if (d->nmin > fabs(val))
                                    d->nmin = fabs(val);
                                  if (d->nmax < fabs(val))
                                    d->nmax = fabs(val);
                                }
                            }
                        }
                    }
                }
            }
          fclose(df);
        }

      fclose(bf);
#ifdef DEBUG
      printf("Pos : Min = %e   Max = %e\n", d->pmin, d->pmax);
      printf("Neg : Min = %e   Max = %e\n", d->nmin, d->nmax);
#endif
    }
  else
    {
      fprintf(stderr, "Reading density: SORRY: Requested old method, not implemented\n");
      SET_FLAG(ctrl.has, CQ_CTRL_HAS_DENSITY, FALSE); /* Failed to load the density */
      return 1;
    }
  for (i = 0; i < 3; ++i)
    {
      d->lvm[i] = limits.realdims[i];
      d->angle[i] = 90.0;
      d->del[i] = d->lvm[i] / (d->pt[i] - 1);
    }
  /* Resetting of limits:
   * -Grid window = The same as before, unless there was not enough data to read; then, data read
   * -Grid shift = Zero (because this would be used for output and in principle we want to write everything from the memory)
   * -Strides = One (for the same reason: All data would be written)
   * -Real window = Adjusted to actual data in the memory (as always)
   * -Real shift = Untouched (because it could be used to adjust written coordinates)
   */
  // TODO When we checked whether the data was enough to fill the window, adjust the window limits accordingly
  for (i = 0; i < 3; ++i)
    {
      limits.gridshft[i] = 0;
      limits.strides[1] = 1;
    }

  //  SET_FLAG(ctrl.has, CQ_CTRL_HAS_DENSITY, TRUE); /* Signal that we have a density */
  if (ctrl.verbose)
    printf("\n");
  return 0; /* It seems it worked */
}

/* This function uses the real shift (real and grid windows are equivalents) */
/* There is no correction of the limits; just use them as they are */
/* An empty window actually means to load everything */
int readXyzCoordinates(char *name, xyzcoords_t *c)
{
  char line[MAXLIN];
  int i, j;
  int oldnumatoms;
  FILE *fp;
  char dummy[MAXLIN];
  double coord[3];
  int atcount, atnum;

  if ((fp = fopen(name, "r")) == NULL)
    {
      fprintf(stderr, "I can't open xyz coordinates file '%s'\n", name);
      return 1; /* Error : I/O */
    }

  oldnumatoms = c->numatoms;
  fgets(line, MAXLIN, fp);
  if (sscanf(line, "%d", &atnum) != 1) /* Check that it looks like an xyz */
    {
      if (ctrl.verbose)
        {
          fprintf(stderr, "The coordinates file doesn't look like an XYZ file\n");
        }
      return 2; /* Error : File type */
    }
  /* Deallocate coordinates if they were already allocated */
  if (c->xyz != NULL)
    {
      if (ctrl.verbose)
        {
          printf("Destroying previous coordinates\n");
        }
      for (i = 0; i < oldnumatoms; ++i)
        {
          free(c->atname[i]);
          free(c->xyz[i]);
        }
      free(c->atname);
      free(c->xyz);
      free(c->move);
    }
  /* Use limits, which basically means to load only atoms within the window, taking the shift into account */
  /* For this, first count the number of atoms to be loaded */
  atcount = 0;
  fgets(line, MAXLIN, fp); /* Skip title */
  for (i = 0; i < atnum; ++i)
    {
      if (fgets(line, MAXLIN, fp) == NULL)
        {
          fprintf(stderr, "Error reading XYZ file: '%s'\n", name);
          return 3; /* Error : incomplete */
        }
      if (sscanf(line, "%s %lf %lf %lf", dummy, &coord[0], &coord[1], &coord[2]) != 4)
        {
          fprintf(stderr, "Error reading XYZ file: '%s'\n", name);
          return 3; /* Error : incomplete */
        }
      /* Count atoms within limits */
      /* If a real window component is zero, use all atoms */
      if ((limits.realdims[0] < ZEROTOL
          || (coord[0] >= limits.realshft[0] && coord[0] <= limits.realshft[0] + limits.realdims[0])))
        if ((limits.realdims[1] < ZEROTOL || (coord[1] >= limits.realshft[1] && coord[1] <= limits.realshft[1]
            + limits.realdims[1])))
          if ((limits.realdims[2] < ZEROTOL || (coord[2] >= limits.realshft[2] && coord[2] <= limits.realshft[2]
              + limits.realdims[2])))
            ++atcount;
    }
  /* Back to the beginning */
  fseek(fp, 0, SEEK_SET);
  c->numatoms = atcount;
  /* Read the file*/
  fgets(line, MAXLIN, fp); /* Ignore atoms in file (we already know) */
  fgets(line, MAXLIN, fp); /* Ignore title */
  c->atname = (char **) malloc((size_t) c->numatoms * sizeof(char *));
  c->xyz = (double **) malloc((size_t) c->numatoms * sizeof(double *));
  c->move = (char *) malloc((size_t) c->numatoms * sizeof(char));
  for (i = 0; i < c->numatoms; ++i)
    {
      c->atname[i] = (char *) malloc((size_t) MAXLIN * sizeof(char));
      c->xyz[i] = (double *) malloc((size_t) 3 * sizeof(double));
      c->move[i] = 0;
    }
#ifdef DEBUG
  printf("Going to read: %d atoms\n", atcount);
#endif
  /* Now we truly read */
  atcount = 0;
  for (i = 0; i < atnum; ++i)
    {
      if (fgets(line, MAXLIN, fp) == NULL)
        {
          fprintf(stderr, "Error reading XYZ file: '%s'\n", name);
          return 3; /* Error : incomplete */
        }
      if (sscanf(line, "%s %lf %lf %lf", dummy, &coord[0], &coord[1], &coord[2]) != 4)
        {
          fprintf(stderr, "Error reading XYZ file: '%s'\n", name);
          return 3; /* Error : incomplete */
        }
      if ((limits.realdims[0] < ZEROTOL
          || (coord[0] >= limits.realshft[0] && coord[0] <= limits.realshft[0] + limits.realdims[0])))
        if ((limits.realdims[1] < ZEROTOL || (coord[1] >= limits.realshft[1] && coord[1] <= limits.realshft[1]
            + limits.realdims[1])))
          if ((limits.realdims[2] < ZEROTOL || (coord[2] >= limits.realshft[2] && coord[2] <= limits.realshft[2]
              + limits.realdims[2])))
            {
              strcpy(c->atname[atcount], dummy);
              for (j = 0; j < 3; ++j)
                {
                  c->xyz[atcount][j] = coord[j];
                  if (ctrl.shftcoords) /* If this toggle is True, we need to move coordinates "back" */
                    c->xyz[atcount][j] -= limits.realshft[j];
                }
#ifdef DEBUG
              printf("%d %s %8.3f %8.3f %8.3f\n", atcount, c->atname[atcount], c->xyz[atcount][0], c->xyz[atcount][1], c->xyz[atcount][2]);
#endif
              ++atcount;
            }
    }
  if (ctrl.verbose)
    printf("Info: Read %d atoms from '%s'\n", atnum, name);
  return 0; /* No errors */
}

/* NOTE: Conquest-style coordinates are also stored in an xyz structure */
int readCqCoordinates(cqruninfo_t *cq, xyzcoords_t *c)
{
  char line[MAXLIN];
  int i, j;
  int oldnumatoms;
  FILE *fp;
  char dummy[MAXLIN];
  double coord[3];
  int atcount, atnum;
  double lattice[3][3]; /* To temporarily store lattice vectors (independent of the final values, which depends on limits) */
  int countread; /* Used to store a number of read fields */
  int type; /* Atom type, as listed in the cqruninfo_t structure */
  char move[3][MAXLIN]; /* The movement flags (T/F) */

  /* Check that we have the name of the coordinates file */
  /* Instead of using FLAG_UP(cq->have, CQ_HAVE_COORDINATES), which refers to the readability of the file */
  /* use the length of the file name and try to open it, if defined */
  if (!strlen(cq->coordinatesfl))
    {
      if (ctrl.verbose)
        printf("The Conquest coordinates file is not defined\n");
      return 4; /* Error: Undefined file */
    }
  if ((fp = fopen(cq->coordinatesfl, "r")) == NULL)
    {
      if (ctrl.verbose)
        printf("I can't open Conquest coordinates file '%s'\n", cq->coordinatesfl);
      return 1; /* Error : I/O */
    }
  /* Now we know that the file exists and can be opened for reading */
  /* This doesn't mean it's working, but this flag refers to whether we know the file name and it exists */
  SET_FLAG(cqinfo.have, CQ_HAVE_COORDINATES,TRUE);
  /* Check the we known the expected scaling of coords. (absolute vs. fractional) */
  /* If unset, set to default */
  if (!FLAG_UP(cq->have, CQ_HAVE_FRACTIONAL))
    {
      SET_FLAG(cq->have, CQ_HAVE_FRACTIONAL, TRUE); /* The flag will be set in any case after calling this function */
      SET_FLAG(cq->read, CQ_READ_FRACTIONAL, isTrue(CQ_DFLT_FRACTIONAL));
      if (ctrl.verbose)
        printf("WARNING: Setting flag to default: %s = %s\n", CQ_FLAG_FRACTIONAL, REPORT_FLAG(cq->read, CQ_READ_FRACTIONAL));
    }
  oldnumatoms = c->numatoms;
  countread = 0;
  fgets(line, MAXLIN, fp);
  countread += sscanf(line, "%lf %lf %lf", &lattice[0][0], &lattice[0][1], &lattice[0][2]);
  fgets(line, MAXLIN, fp);
  countread += sscanf(line, "%lf %lf %lf", &lattice[1][0], &lattice[1][1], &lattice[1][2]);
  fgets(line, MAXLIN, fp);
  countread += sscanf(line, "%lf %lf %lf", &lattice[2][0], &lattice[2][1], &lattice[2][2]);
  fgets(line, MAXLIN, fp);
  countread += sscanf(line, "%d", &atnum);
  if (countread != 10) /* Check that it looks like a Conquest coordinates file */
    {
      if (ctrl.verbose) /* Perhaps this message should be left to the python function */
        {
          fprintf(stderr, "The coordinates file doesn't look like a Conquest file\n");
        }
      return 2; /* Error : File type */
    }
  /* Check that the simulation cell looks orthorhombic */
  for (i = 0; i < 3; ++i)
    {
      for (j = 0; j < 3; ++j)
        {
          if (i != j)
            if (lattice[i][j] > ZEROTOL)
              {
                if (ctrl.verbose)
                  {
                    printf("WARNING: At present Conquest only handles orthorhombic cells\n");
                    printf("         Yours doesn't seem to be. Non diagonal components of the lattice will be ignored\n");
                  }
              }
        }
    }
  /* Use limits, which basically means to load only atoms within the window, taking the shift into account */
  /* For this, first count the number of atoms to be loaded */
#ifdef DEBUG
  printf("%f %f %f\n", lattice[0][0], lattice[0][1], lattice[0][2]);
  printf("%f %f %f\n", lattice[1][0], lattice[1][1], lattice[1][2]);
  printf("%f %f %f\n", lattice[2][0], lattice[2][1], lattice[2][2]);
  printf("%d\n", atnum);
#endif
  atcount = 0;
  for (i = 0; i < atnum; ++i)
    {
      //      printf("%d\n", i);
      if (fgets(line, MAXLIN, fp) == NULL)
        {
          fprintf(stderr, "Error reading Conquest coordinate file: '%s'\n", cq->coordinatesfl);
          return 3; /* Error : incomplete */
        }
      if (sscanf(line, "%lf %lf %lf %d %s %s %s", &coord[0], &coord[1], &coord[2], &type, move[2], move[1], move[0]) != 7)
        {
          fprintf(stderr, "Error reading Conquest coordinate file: '%s'\n", cq->coordinatesfl);
          return 3; /* Error : incomplete */
        }
      for (j = 0; j < 3; ++j)
        {
          /* Rescale coordinates, if the file is fractional */
          if (FLAG_UP(cq->read, CQ_READ_FRACTIONAL))
            coord[j] *= lattice[j][j] * BOHR_TO_ANGSTROM;
          /* Wrap around, the inefficient way */
          while (coord[j] < 0.0)
            coord[j] += lattice[j][j];
          while (coord[j] > lattice[j][j])
            coord[j] -= lattice[j][j];
          /* Convert coordinates to Angstrom (always in Bohr, even if the user's units are Angs.) */
          coord[j] *= BOHR_TO_ANGSTROM;
        }
      /* Count atoms within limits */
      /* If a real window component is zero, use all atoms */
      if ((limits.realdims[0] < ZEROTOL
          || (coord[0] >= limits.realshft[0] && coord[0] <= limits.realshft[0] + limits.realdims[0])))
        if ((limits.realdims[1] < ZEROTOL || (coord[1] >= limits.realshft[1] && coord[1] <= limits.realshft[1]
            + limits.realdims[1])))
          if ((limits.realdims[2] < ZEROTOL || (coord[2] >= limits.realshft[2] && coord[2] <= limits.realshft[2]
              + limits.realdims[2])))
            {
#ifdef DEBUG
              printf("Read  : %f %f %f | %f %f %f\n", coord[0], coord[1], coord[2], coord[0] / BOHR_TO_ANGSTROM, coord[1]
                  / BOHR_TO_ANGSTROM, coord[2] / BOHR_TO_ANGSTROM);
#endif
              ++atcount;
            }
          else
            {
#ifdef DEBUG
              printf("UnRead: %f %f %f | %f %f %f\n", coord[0], coord[1], coord[2], coord[0] / BOHR_TO_ANGSTROM, coord[1]
                  / BOHR_TO_ANGSTROM, coord[2] / BOHR_TO_ANGSTROM);
#endif
            }
    }
#ifdef DEBUG
  printf("Going to read: %d atoms\n", atcount);
#endif
  if(atcount == 0)
    {
      if(ctrl.verbose)
        printf("Error: No atoms can be read in file '%s' within current window. Aborting\n", cq->coordinatesfl);
      return(3); /* Incomplete file */
    }
  /* Deallocate coordinates if they were already allocated */
  if (c->xyz != NULL)
    {
      if (ctrl.verbose)
        {
          printf("Destroying previous coordinates\n");
        }
      for (i = 0; i < oldnumatoms; ++i)
        {
          free(c->atname[i]);
          free(c->xyz[i]);
        }
      free(c->atname);
      free(c->xyz);
      free(c->move);
    }
  /* Back to the beginning */
  fseek(fp, 0, SEEK_SET);
  c->numatoms = atcount;
  /* Read the file*/
  for (i = 0; i < 4; ++i) /* Ignore the first 4 lines */
    fgets(line, MAXLIN, fp);
  c->atname = (char **) malloc((size_t) c->numatoms * sizeof(char *));
  c->xyz = (double **) malloc((size_t) c->numatoms * sizeof(double *));
  c->move = (char *) malloc((size_t) c->numatoms * sizeof(char));
  for (i = 0; i < c->numatoms; ++i)
    {
      c->atname[i] = (char *) malloc((size_t) MAXLIN * sizeof(char));
      c->xyz[i] = (double *) malloc((size_t) 3 * sizeof(double));
    }
  /* Now we truly read */
  atcount = 0;
  for (i = 0; i < atnum; ++i)
    {
#ifdef DEBUG
      printf("%d\n", i);
#endif
      if (fgets(line, MAXLIN, fp) == NULL)
        {
          fprintf(stderr, "Error reading Conquest coordinate file: '%s'\n", cq->coordinatesfl);
          return 3; /* Error : incomplete */
        }
      if (sscanf(line, "%lf %lf %lf %d %s %s %s", &coord[0], &coord[1], &coord[2], &type, move[2], move[1], move[0]) != 7)
        {
          fprintf(stderr, "Error reading Conquest coordinate file: '%s'\n", cq->coordinatesfl);
          return 3; /* Error : incomplete */
        }
      for (j = 0; j < 3; ++j)
        {
          /* Rescale coordinates, if the file is fractional */
          if (FLAG_UP(cq->read, CQ_READ_FRACTIONAL))
            coord[j] *= lattice[j][j];
          /* Wrap around, the inefficient way */
          while (coord[j] < 0.0)
            coord[j] += lattice[j][j];
          while (coord[j] > lattice[j][j])
            coord[j] -= lattice[j][j];
          /* Convert coordinates to Angstrom (always in Bohr, even if the user's units are Angs.) */
          coord[j] *= BOHR_TO_ANGSTROM;
        }
      /* A few extra checks this time */
      /* 1. Is atom type is within expected bounds? */
      if (type < 1 || type > cq->numspecies)
        {
          fprintf(stderr, "Error reading Conquest coordinate file: '%s'. Wrong species number\n", cq->coordinatesfl);
          return 2; /* Error: File type */
        }
      /* 2. Are movement flags correct? */
      for (j = 0; j < 3; ++j)
        {
          if (strcmp(move[j], "T") && strcmp(move[j], "F")) /* Neither true nor false */
            {
              fprintf(stderr, "Error reading Conquest coordinate file: '%s'. Wrong movement flags\n", cq->coordinatesfl);
              return 2; /* Error: File type */
            }
        }

      /* Count atoms within limits */
      /* If a real window component is zero, use all atoms */
      if ((limits.realdims[0] < ZEROTOL
          || (coord[0] >= limits.realshft[0] && coord[0] <= limits.realshft[0] + limits.realdims[0])))
        if ((limits.realdims[1] < ZEROTOL || (coord[1] >= limits.realshft[1] && coord[1] <= limits.realshft[1]
            + limits.realdims[1])))
          if ((limits.realdims[2] < ZEROTOL || (coord[2] >= limits.realshft[2] && coord[2] <= limits.realshft[2]
              + limits.realdims[2])))
            {
              /* Do the actual reading here */
              c->move[atcount] = 0; /* Reset movements */
              for (j = 0; j < 3; ++j)
                {
                  c->xyz[atcount][j] = coord[j];
                  if (ctrl.shftcoords) /* If this toggle is True, we need to move coordinates "back" */
                    c->xyz[atcount][j] -= limits.realshft[j];
                  if (!strcmp(move[j], "T")) /* Allow movements this atom and component */
                    c->move[atcount] |= (0x01 << j);
                }
              /* Map atom index to atom name */
              strcpy(c->atname[atcount], cq->typename[type - 1]);
#ifdef DEBUG
              printf("%d %s %8.3f %8.3f %8.3f\n", atcount, c->atname[atcount], c->xyz[atcount][0], c->xyz[atcount][1], c->xyz[atcount][2]);
#endif
              ++atcount;
            }
    }
  if (ctrl.verbose)
    printf("Info: Read %d atoms from '%s'\n", c->numatoms, cq->coordinatesfl);
  return 0; /* No errors */
}

/* If a density is not available, pass NULL for the second parameter */
int writeCqCoordinates(cqruninfo_t *cq, density_t *d, xyzcoords_t *c, char *name)
{
  FILE *fp;
  int i, j;
  int atcount;
  double dims[3]; /* Dimensions to be written to the output file */
  double coord[3]; /* Temporary storage for coordinates */
  char moveflags[MAXLIN]; /* String for "movement" flags */
  int cqatomindex; /* A Conquest atom index (the one listed in the input file for each atom type) */
  double max[3]; /* Largest coordinates */
  int knowntypes; /* Total known atom types */
  int knownfound; /* Count of atoms that are listed in our periodic table */
  int unknownfound; /* Count of atoms that are not listed in our periodic table */
  atominfo_t **atomtable; /* Pointers to atomic data */
  char found;
  char using; /* Flag to take note of the cell dimensions to be used*/
  /* IMP NOTE : If the real window is used as the dimensions of the output file,
   *            the real shift will be used to shift the coordinates.
   *            This is to keep atoms within the cell (as far as possible).
   */

#ifdef DEBUG
  if(d != NULL)
    printf("Got a density\n");
  else
    printf("Did not get a density\n");
#endif

  /* Prepare the table of atom types */
  /* 1. Count the types with have in our periodic table */
  knownfound = 0;
  unknownfound = 0;
  knowntypes = 0;
  while (periodictable[knowntypes].symbol != NULL)
    ++knowntypes;
#ifdef DEBUG
  printf("Known types = %d\n", knowntypes);
#endif
  /* 2. Allocate enough pointers */
  atomtable = (atominfo_t **) malloc((knowntypes + MAXUNKNOWN) * sizeof(atominfo_t *));
  for (i = 0; i < knowntypes + MAXUNKNOWN; ++i)
    atomtable[i] = NULL; /* For the moment, the table is empty */

  /* Check the we known the expected scaling of coords. (absolute vs. fractional) */
  /* If unset, set to default */
  if (!FLAG_UP(cq->have, CQ_HAVE_FRACTIONAL))
    {
      SET_FLAG(cq->have, CQ_HAVE_FRACTIONAL, TRUE); /* The flag will be set in any case after calling this function */
      SET_FLAG(cq->read, CQ_READ_FRACTIONAL, isTrue(CQ_DFLT_FRACTIONAL));
      if (ctrl.verbose)
        printf("WARNING: Setting flag to default: %s = %s\n", CQ_FLAG_FRACTIONAL, REPORT_FLAG(cq->read, CQ_READ_FRACTIONAL));
    }
  /* Check that we have appropriate cell dimensions */
  /* First from Conquest info; the user could set these dimensions to any desired value */
  /* The real window will be used anyway to pick atoms, but written dimensions will be the user's */
  if (cq->lvm[0] > ZEROTOL && cq->lvm[1] > ZEROTOL && cq->lvm[0] > ZEROTOL)
    {
      for (i = 0; i < 3; ++i)
        dims[i] = cq->lvm[i];
      using = CQ_USING_CQINFO;
      if (ctrl.verbose)
        printf("Info: Using cell dimensions from the Conquest info structures\n"
          "      If you don't want this, set the CQ cell dimensions to zero, with cq_info_set_cell_dimensions\n");
    }
  else if (limits.realdims[0] > ZEROTOL && limits.realdims[1] > ZEROTOL && limits.realdims[2] > ZEROTOL) /* Using limits */
    {
      for (i = 0; i < 3; ++i)
        dims[i] = limits.realdims[i];
      using = CQ_USING_LIMITS;
      if (ctrl.verbose)
        printf("Info: Using cell dimensions from limits (real window).\n"
          "      If you don't want this, set the CQ cell dimensions, with cq_info_set_cell_dimensions\n");
    }
  else if (d != NULL && d->lvm[0] > ZEROTOL && d->lvm[1] > ZEROTOL && d->lvm[2] > ZEROTOL)
    {
      /* If the previous two checks fail, it means that the user didn't choose explicit cell dimensions */
      /* (using the cq info structure) and that (s)he wants all the atoms, bc. the real window is zero */
      /* In this case, if a density is available, the cell dimensions will be used */
      for (i = 0; i < 3; ++i)
        dims[i] = d->lvm[i];
      using = CQ_USING_DENSITY;
      if (ctrl.verbose)
        printf("Info: Using cell dimensions from the available density.\n"
          "      If you don't want this, set the CQ cell dimensions, with cq_info_set_cell_dimensions,\n"
          "      or choose a window (and maybe shift) with cq_limits_set_real_window (and cq_limits_set_real_shift)");
    }
  else /* If everything else fails, calculate dimensions using the maximum coordinates of the loaded molecule */
    {
      for (i = 0; i < 3; ++i)
        max[i] = 0.0;
      for (i = 0; i < c->numatoms; ++i)
        {
          for (j = 0; j < 3; ++j)
            if (c->xyz[i][j] > max[j])
              max[j] = c->xyz[i][j];
        }
      for (i = 0; i < 3; ++i)
        dims[i] = max[i];
      using = CQ_USING_MAXIMA;
      if (ctrl.verbose)
        printf("Info: Using largest coordinates to decide the cell dimensions\n"
          "      If you don't want this, set the CQ cell dimensions, with cq_info_set_cell_dimensions\n");
    }
  if (ctrl.verbose)
    printf("Info: Writing a file with the folowing dimensions: %f %f %f A (%f %f %f a0)\n", dims[0], dims[1], dims[2], dims[0]
        / BOHR_TO_ANGSTROM, dims[1] / BOHR_TO_ANGSTROM, dims[2] / BOHR_TO_ANGSTROM);
  /* Count the number of atoms to be written */
  /* Also, prepare a table of atom types, to go with the coordinates file */
  atcount = 0;
  for (i = 0; i < c->numatoms; ++i)
    {
      /* If a real window component is zero, use all atoms */
      if ((limits.realdims[0] < ZEROTOL || (c->xyz[i][0] >= limits.realshft[0] && c->xyz[i][0] <= limits.realshft[0]
          + limits.realdims[0])))
        if ((limits.realdims[1] < ZEROTOL || (c->xyz[i][1] >= limits.realshft[1] && c->xyz[i][1] <= limits.realshft[1]
            + limits.realdims[1])))
          if ((limits.realdims[2] < ZEROTOL || (c->xyz[i][2] >= limits.realshft[2] && c->xyz[i][2] <= limits.realshft[2]
              + limits.realdims[2])))
            {
              //              printf("%d = %s\n", i, c->atname[i]);
              /* Check whether we already have this atom listed... */
              /* ...first in the part of the table for known atom types */
              found = NO;
              for (j = 0; j < knownfound; ++j)
                {
                  //                  printf("  %d\n", j);
                  if (atomtable[j] != NULL && !strcmp(c->atname[i], atomtable[j]->symbol)) /* Found */
                    {
                      found = YES;
                      break;
                    }
                }
              if (j == knownfound) /* Not found in the list we have so far */
                {
                  for (j = 0; j < unknownfound; ++j)
                    {
                      //                      printf("U %d\n", j);
                      if (atomtable[knowntypes + j] != NULL && !strcmp(c->atname[i], atomtable[knowntypes + j]->symbol)) /* Found */
                        {
                          found = YES;
                          break;
                        }
                    }
                }
              if (found == NO) /* Not found in our list so far */
                for (j = 0; j < knowntypes; ++j) /* Very inefficient search of the atom name in the periodic table */
                  {
                    //                    printf("    %d\n", j);
                    if (!strcmp(c->atname[i], periodictable[j].symbol)) /* Found */
                      {
                        /* Assign a pointer to the atom information */
                        atomtable[knownfound] = &periodictable[j];
                        ++knownfound;
                        found = YES;
                        break;
                      }
                  }
              if (found == NO) /* Not know */
                {
                  /* Prepare a new structure to indicate that this is unknown */
                  atomtable[knowntypes + unknownfound] = (atominfo_t *) malloc(sizeof(atominfo_t));
                  atomtable[knowntypes + unknownfound]->symbol = (char *) malloc(MAXLIN * sizeof(char));
                  atomtable[knowntypes + unknownfound]->z = 0;
                  strcpy(atomtable[knowntypes + unknownfound]->symbol, c->atname[i]);
                  atomtable[knowntypes + unknownfound]->weight = 0.0;
                  if (ctrl.verbose)
                    {
                      printf("WARNING: Atom type '%s' is not in our periodic table\n", c->atname[i]);
                      if (strlen(c->atname[i]) > 2)
                        printf("WARNING: Atom names should not be longer than 2 characters\n", c->atname[i]);
                    }
                  ++unknownfound;
                }
              ++atcount;
            }
    }
  if (atcount > 0)
    {
      if ((fp = fopen(name, "w")) == NULL)
        {
          fprintf(stderr, "I can't create the coordinates file '%s'\n", name);
          return 1; /* Error : I/O */
        }
      /* Write lattice vectors */
      fprintf(fp, "% 10.6f % 10.6f % 10.6f\n", dims[0] / BOHR_TO_ANGSTROM, 0.0, 0.0);
      fprintf(fp, "% 10.6f % 10.6f % 10.6f\n", 0.0, dims[1] / BOHR_TO_ANGSTROM, 0.0);
      fprintf(fp, "% 10.6f % 10.6f % 10.6f\n", 0.0, 0.0, dims[2] / BOHR_TO_ANGSTROM);
      fprintf(fp, "%6d\n", atcount);
      for (i = 0; i < c->numatoms; ++i)
        {
          /* If a real window component is zero, use all atoms */
          if ((limits.realdims[0] < ZEROTOL || (c->xyz[i][0] >= limits.realshft[0] && c->xyz[i][0] <= limits.realshft[0]
              + limits.realdims[0])))
            if ((limits.realdims[1] < ZEROTOL || (c->xyz[i][1] >= limits.realshft[1] && c->xyz[i][1] <= limits.realshft[1]
                + limits.realdims[1])))
              if ((limits.realdims[2] < ZEROTOL || (c->xyz[i][2] >= limits.realshft[2] && c->xyz[i][2] <= limits.realshft[2]
                  + limits.realdims[2])))
                {
                  //printf("%s -- %d -- ", c->atname[i], cqatomindex);
                  for (j = 0; j < knowntypes + unknownfound; ++j)
                    {
                      if (j == knownfound) /* Jump to the end of the table (i.e. to the "unknown types") */
                        j = knowntypes;
                      if (!strcmp(c->atname[i], atomtable[j]->symbol)) /* Found */
                        {
                          //printf("Found %d -- ", j);
                          if (j < knowntypes)
                            cqatomindex = j + 1; /* CQ indexes start at 1 */
                          else
                            /* This type is not in the periodic table: Set correct index */
                            cqatomindex = j + 1 - knowntypes + knownfound;
                          break;
                        }
                    }
                  //printf("%d\n", cqatomindex);
                  moveflags[0] = '\0';
                  for (j = 0; j < 3; ++j)
                    {
                      if ((c->move[i] >> j) & 0x01)
                        strcat(moveflags, " T");
                      else
                        strcat(moveflags, " F");
                      coord[j] = c->xyz[i][j];
                      coord[j] /= BOHR_TO_ANGSTROM; /* Use bohr */
                      if(using == CQ_USING_LIMITS || ctrl.shftcoords) /* Shift if necessary/requested */
                        {
                          coord[j] -= limits.realshft[j] / BOHR_TO_ANGSTROM;
                        }
                      /* If fractional, rescale the coordinates */
                      if (FLAG_UP(cq->have, CQ_HAVE_FRACTIONAL) && FLAG_UP(cq->read, CQ_READ_FRACTIONAL))
                        {
                          coord[j] /= (dims[j] / BOHR_TO_ANGSTROM);
                        }
                    }
                  fprintf(fp, "%14.6f %14.6f %14.6f %3d %s\n", coord[0], coord[1], coord[2], cqatomindex, moveflags);
                }
        }
    }
  else
    {
      if (ctrl.verbose)
        printf("WARNING: No atoms to write within current limits: File not written\n");
      return 2; /* No atoms to write */
    }
  fclose(fp);
  /* Print the atom table we found */
  printf("----- FYI: You could copy the next few lines into your Conquest_input file -----\n");
  printf("%s %s\n", CQ_FLAG_FRACTIONAL, REPORT_FLAG_SHORT(cq->read, CQ_READ_FRACTIONAL));
  printf("%s %d\n", CQ_FLAG_SPECIES, knownfound + unknownfound);
  printf("%%block %s\n", CQ_FLAG_CHEMBLOCK);
  for (i = 0; i < knownfound; ++i)
    {
      printf("%2d %6.2f  %2s\n", i + 1, atomtable[i]->weight, atomtable[i]->symbol);
    }
  for (i = 0; i < unknownfound; ++i)
    {
      printf("%2d %6.2f  %2s\n", knownfound + i + 1, atomtable[knowntypes + i]->weight, atomtable[knowntypes + i]->symbol);
    }
  printf("%%endblock\n");
  printf("----- FYI END ------------------------------------------------------------------\n");
  /* Free atom table */
  /* Only the "unknown" (not in the periodic table) entries were dynamically allocated */
  for (i = knowntypes; i < knowntypes + MAXUNKNOWN; ++i)
    {
      if (atomtable[i] != NULL)
        {
          free(atomtable[i]->symbol);
          free(atomtable[i]);
        }
    }
  free(atomtable);

  return 0; /* No errors */
}

/* If a density is not available, pass NULL for the second parameter */
int writeXyzCoordinates(cqruninfo_t *cq, density_t *d, xyzcoords_t *c, char *name)
{
  FILE *fp;
  int i, j;
  int atcount;
  double dims[3]; /* Dimensions to be written to the output file */
  double coord[3]; /* Temporary storage for coordinates */
  double max[3]; /* Largest coordinates */
  char using; /* Flag to take note of the cell dimensions to be used*/
  /* IMP NOTE : If the real window is used as the dimensions of the output file,
   *            the real shift will be used to shift the coordinates.
   *            This is to keep atoms within the cell (as far as possible).
   */

#ifdef DEBUG
  if(d != NULL)
    printf("Got a density\n");
  else
    printf("Did not get a density\n");
#endif

  /* "Fractionality" of coordinates is not relevant here */

  /* Check that we have appropriate cell dimensions */
  /* First from Conquest info; the user could set these dimensions to any desired value */
  /* The real window will be used anyway to pick atoms, but written dimensions will be the user's */
  if (cq->lvm[0] > ZEROTOL && cq->lvm[1] > ZEROTOL && cq->lvm[0] > ZEROTOL)
    {
      for (i = 0; i < 3; ++i)
        dims[i] = cq->lvm[i];
      using = CQ_USING_CQINFO;
      if (ctrl.verbose)
        printf("Info: Using cell dimensions from the Conquest info structures\n"
          "      If you don't want this, set the CQ cell dimensions to zero, with cq_info_set_cell_dimensions\n");
    }
  else if (limits.realdims[0] > ZEROTOL && limits.realdims[1] > ZEROTOL && limits.realdims[2] > ZEROTOL) /* Using limits */
    {
      for (i = 0; i < 3; ++i)
        dims[i] = limits.realdims[i];
      using = CQ_USING_LIMITS;
      if (ctrl.verbose)
        printf("Info: Using cell dimensions from limits (real window).\n"
          "      If you don't want this, set the CQ cell dimensions, with cq_info_set_cell_dimensions\n");
    }
  else if (d != NULL && d->lvm[0] > ZEROTOL && d->lvm[1] > ZEROTOL && d->lvm[2] > ZEROTOL)
    {
      /* If the previous two checks fail, it means that the user didn't choose explicit cell dimensions */
      /* (using the cq info structure) and that (s)he wants all the atoms, bc. the real window is zero */
      /* In this case, if a density is available, the cell dimensions will be used */
      for (i = 0; i < 3; ++i)
        dims[i] = d->lvm[i];
      using = CQ_USING_DENSITY;
      if (ctrl.verbose)
        printf("Info: Using cell dimensions from the available density.\n"
          "      If you don't want this, set the CQ cell dimensions, with cq_info_set_cell_dimensions,\n"
          "      or choose a window (and maybe shift) with cq_limits_set_real_window (and cq_limits_set_real_shift)");
    }
  else /* If everything else fails, calculate dimensions using the maximum coordinates of the loaded molecule */
    {
      for (i = 0; i < 3; ++i)
        max[i] = 0.0;
      for (i = 0; i < c->numatoms; ++i)
        {
          for (j = 0; j < 3; ++j)
            if (c->xyz[i][j] > max[j])
              max[j] = c->xyz[i][j];
        }
      for (i = 0; i < 3; ++i)
        dims[i] = max[i];
      using = CQ_USING_MAXIMA;
      if (ctrl.verbose)
        printf("Info: Using largest coordinates to decide the cell dimensions\n"
          "      If you don't want this, set the CQ cell dimensions, with cq_info_set_cell_dimensions\n");
    }
  if (ctrl.verbose)
    printf("Info: Writing a file with the folowing dimensions: %f %f %f A (%f %f %f a0)\n", dims[0], dims[1], dims[2], dims[0]
        / BOHR_TO_ANGSTROM, dims[1] / BOHR_TO_ANGSTROM, dims[2] / BOHR_TO_ANGSTROM);
  /* Count the number of atoms to be written */
  atcount = 0;
  for (i = 0; i < c->numatoms; ++i)
    {
      /* If a real window component is zero, use all atoms */
      if ((limits.realdims[0] < ZEROTOL || (c->xyz[i][0] >= limits.realshft[0] && c->xyz[i][0] <= limits.realshft[0]
          + limits.realdims[0])))
        if ((limits.realdims[1] < ZEROTOL || (c->xyz[i][1] >= limits.realshft[1] && c->xyz[i][1] <= limits.realshft[1]
            + limits.realdims[1])))
          if ((limits.realdims[2] < ZEROTOL || (c->xyz[i][2] >= limits.realshft[2] && c->xyz[i][2] <= limits.realshft[2]
              + limits.realdims[2])))
            {
              ++atcount;
            }
    }
  if (atcount > 0)
    {
      if ((fp = fopen(name, "w")) == NULL)
        {
          fprintf(stderr, "I can't create the coordinates file '%s'\n", name);
          return 1; /* Error : I/O */
        }
      fprintf(fp, "%12d\n", atcount);
      fprintf(fp, "Your comment here\n");
      for (i = 0; i < c->numatoms; ++i)
        {
          /* If a real window component is zero, use all atoms */
          if ((limits.realdims[0] < ZEROTOL || (c->xyz[i][0] >= limits.realshft[0] && c->xyz[i][0] <= limits.realshft[0]
              + limits.realdims[0])))
            if ((limits.realdims[1] < ZEROTOL || (c->xyz[i][1] >= limits.realshft[1] && c->xyz[i][1] <= limits.realshft[1]
                + limits.realdims[1])))
              if ((limits.realdims[2] < ZEROTOL || (c->xyz[i][2] >= limits.realshft[2] && c->xyz[i][2] <= limits.realshft[2]
                  + limits.realdims[2])))
                {
                  for (j = 0; j < 3; ++j)
                    {
                      coord[j] = c->xyz[i][j];
                      if(using == CQ_USING_LIMITS || ctrl.shftcoords) /* Shift if necessary/requested */
                        {
                          coord[j] -= limits.realshft[j];
                        }
                    }
                  fprintf(fp, "%2s %8.3f %8.3f %8.3f\n", c->atname[i], coord[0], coord[1], coord[2]);
                }
        }
    }
  else
    {
      if (ctrl.verbose)
        printf("WARNING: No atoms to write within current limits: File not written\n");
      return 2; /* No atoms to write */
    }
  fclose(fp);

  return 0; /* No errors */
}

int makeFilename(char *name, cqruninfo_t *info, int coreno, char *ext)
{
  int i;
  char format[MAXLIN]; /* Used to create a convenient format for density files, etc. */
  int figures; /* Number of decimal figures necessary to store the number of cores */

  figures = (int) ceil(log10(info->cores));
  if (figures < 3)
    figures = 3; /* Keep backward compatibility */
  if (strlen(info->densityfl) + strlen(ext) + figures > MAXLIN - 1)
    return 1;
  sprintf(format, "%s.%%0%dd%s", info->densityfl, figures, ext);
  sprintf(name, format, coreno);
  return 0;
}
