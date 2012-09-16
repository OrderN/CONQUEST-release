/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_slicelog.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_slicelog.c *
 *
 * NAME
 *  cq_ut_conqtour_slicelog.c
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

void calculateDrawSlice(index_t *idx, box_t *b, density_t *d)
{
  int i, j, k, l;
  double scalefactor;
  double tmp1, tmp2, tmp3;
  int tmpi, tmpi2;
  double t; /* The parameter of the line equation */
  double aux1[3], aux2[3], aux3[3];
  double rectdir[4][3]; /* Directions of rectangle edges */
  double pang[6]; /* Angles "on plane", used to order vertices */
  //   double modulus;
  // UNCOMMENT and remove globals when debugged
  //   int pivot[4];
  //   int minindx, minindy, maxindx, maxindy;
  //   double minx, miny, maxx, maxy;
  double x, y;
  double ang, minang, totang;
  double tmprect[4][3];
  double area;
  double maxpt, tmppt; /* Used to calculate grid point spacing after rotating */
  double uaxis[3][3] =
    {
      { 1.0, 0.0, 0.0 }, /* Used to calculate projections onto orthonormal axes */
      { 0.0, 1.0, 0.0 },
      { 0.0, 0.0, 1.0 } };
  int found; /* Used to count corners on a rectangle edge */
  double xcoord[2]; /* X coordinate of the two reference corner */
  //   double xcoord;
  double ycoord[6]; /* Y coordinate of polygon corners in rectangle coordinates */
  //   double tolerance;        /* Used to find intersecting vertices... */
  //   int round;               /* ...being more tolerant if numerical error makes us fail */
  densitycut_t *spt;

  if (idx->dindex[idx->curr] != NULL)
    return; /* This slice is ready */
  /* If the slice is not ready, allocate and prepare it */
  idx->dindex[idx->curr] = (densitycut_t *) malloc((size_t) sizeof(densitycut_t));
  spt = idx->dindex[idx->curr]; /* For more readable code and, potentially, faster access */
  /* Find vertices of the polygon defined by the intersection */
  /* of the plane with the cell box */
  scalefactor = 1.0; /* Used to correct the in-plane point in case it went out of the box */
  spt->nocorners = 0;
  while(spt->nocorners < 3)
    {
      /* First, define an "inplane" point, based on the origin */
      for (i = 0; i < 3; ++i)
        spt->inplane[i] = idx->inplane[i] + scalefactor * ((int) idx->curr - (int) idx->orig) * idx->surfnormal[i] * idx->dz;
#ifdef DEBUG
      printf("Inplane = %e | %d %d | %f %f %f\n", scalefactor, idx->curr, idx->orig,
          spt->inplane[0],
          spt->inplane[1],
          spt->inplane[2]);
#endif
      spt->nocorners = 0;
      for (i = 0; i < 12; ++i)
        {
          tmp1 = dotProduct(idx->surfnormal, b->edge[i]);
          /* Check that the edge and plane are not coplanar */
          /* If they are, we just ignore */
          if (tmp1 < -ZEROTOL || tmp1 > ZEROTOL)
            {
              for (j = 0; j < 3; ++j)
                aux1[j] = spt->inplane[j] - b->refpoint[i][j];
              tmp2 = dotProduct(idx->surfnormal, aux1);
              t = tmp2 / tmp1;
              if (t >= 0.0 && t <= 1.0)
                {
                  for (j = 0; j < 3; ++j)
                    spt->corner[spt->nocorners][j] = b->edge[i][j] * t + b->refpoint[i][j];
#ifdef DEBUG
                  printf("Corner = %d %f | %f %f %f | %f %f %f\n", i, t, spt->corner[spt->nocorners][0], spt->corner[spt->nocorners][1], spt->corner[spt->nocorners][2], b->edge[i][0], b->edge[i][1], b->edge[i][2]);
#endif
                  ++spt->nocorners;
                }
            }
        }
      scalefactor -= ZEROTOL;
    }
#ifdef DEBUG
              printf("No. Corners = %d\n", spt->nocorners);
#endif
  /* Before we continue, it is convenient to have the vertices ordered, */
  /* so they form a convex polynomial. To achieve this, change to 2D */
  /* coordinates on the intersecting plane, with origin at spt->inplane   */
  /* 1. Prepare the basis, using one arbitrary vertex to define x and the normal */
  for (i = 0; i < 3; ++i)
    aux1[i] = spt->corner[0][i] - spt->inplane[i];
  crossProduct(aux1, idx->surfnormal, aux2);
  normalize(aux1);
  normalize(aux2);
  //printf("%f %f %f | %f %f %f\n", aux1[0], aux1[1], aux1[2], aux2[0], aux2[1], aux2[2]);
  /* 2. Now calculate angles to the first vertex */
  /*    We do it like this, instead of a simple dot product, */
  /*    because acos is limited to one quadrant */
  pang[0] = 0.0; /* For first vertex, angle with itself */
  for (i = 1; i < spt->nocorners; ++i)
    {
      for (j = 0; j < 3; ++j)
        aux3[j] = spt->corner[i][j] - spt->inplane[j];
      tmp1 = dotProduct(aux1, aux3); /* On-plane x-coord */
      tmp2 = dotProduct(aux2, aux3); /* On-plane y-coord */
      pang[i] = RAD_TO_DEG * atan2(tmp2, tmp1);
      //printf("%f %f %f    %f %f\n", aux3[0], aux3[1], aux3[2], tmp1, tmp2);
      //printf("A %f\n", pang[i]);
    }
  /* Finally, reorder vertices (simple bubble) */
  for (i = 0; i < spt->nocorners - 1; ++i)
    {
      for (j = 0; j < spt->nocorners - 1 - i; ++j)
        {
          if (pang[j] > pang[j + 1]) /* Swap... */
            {
              tmp1 = pang[j]; /* ...angles... */
              pang[j] = pang[j + 1];
              pang[j + 1] = tmp1;
              for (k = 0; k < 3; ++k) /* ...and coordinates */
                {
                  tmp1 = spt->corner[j][k];
                  spt->corner[j][k] = spt->corner[j + 1][k];
                  spt->corner[j + 1][k] = tmp1;
                }
            }
        }
    }
  //for(i=0; i < spt->nocorners; ++i)  printf("O %f\n", pang[i]);

  /* Find the minimum area enclosing rectangle, using rotating calipers */
  /* Step 1: Take any pair of consecutive points; use the line they define */
  /*         as one axis, and the perpendicular as another; calculate four */
  /*         extreme vertices of the polygon */
  for (i = 0; i < 3; ++i)
    aux1[i] = spt->corner[1][i] - spt->corner[0][i]; /* x-axis, on plane */
  //   crossProduct(aux1, spt->inplane, aux2);
  crossProduct(aux1, idx->surfnormal, aux2); /* y-axis, on plane */
  normalize(aux1);
  normalize(aux2);
  /* Using, arbitrarily, the first vertex as origin, get on plane coords of all others */
  /* and choose extremes */
  minindx = minindy = maxindx = maxindy = 0;
  minx = miny = maxx = maxy = 0.0;
  for (i = 1; i < spt->nocorners; ++i)
    {
      /* Vector position of the vertex */
      for (j = 0; j < 3; ++j)
        aux3[j] = spt->corner[i][j] - spt->corner[0][j];
      /* Coordinates in this system */
      x = dotProduct(aux1, aux3);
      y = dotProduct(aux2, aux3);
#ifdef DEBUG
      printf("%d %f %f\n", i, x, y);
#endif
      /* Update of extremes */
      if (x < minx)
        {
          minindx = i;
          minx = x;
        }
      if (x > maxx)
        {
          maxindx = i;
          maxx = x;
        }
      if (y < miny)
        {
          minindy = i;
          miny = y;
        }
      if (y > maxy)
        {
          maxindy = i;
          maxy = y;
        }
    }
#ifdef DEBUG
  printf("Minmax : %d %d %d %d | %f %f %f %f\n", minindx, maxindx, minindy, maxindy, minx, maxx, miny, maxy);
#endif
  /* Step 2: Calculate the enclosing rectangle and its area */
  for (i = 0; i < 3; ++i)
    {
      /* Remember the axes */
      spt->xaxis[i] = aux1[i];
      spt->yaxis[i] = aux2[i];
      /* Take note of the enclosing rectangle (in world coordinates) */
      spt->rectangle[0][i] = spt->corner[0][i] + minx * spt->xaxis[i] + miny * spt->yaxis[i];
      spt->rectangle[1][i] = spt->corner[0][i] + maxx * spt->xaxis[i] + miny * spt->yaxis[i];
      spt->rectangle[2][i] = spt->corner[0][i] + maxx * spt->xaxis[i] + maxy * spt->yaxis[i];
      spt->rectangle[3][i] = spt->corner[0][i] + minx * spt->xaxis[i] + maxy * spt->yaxis[i];
      /* Prepare for area */
      aux1[i] = spt->rectangle[1][i] - spt->rectangle[0][i];
      aux2[i] = spt->rectangle[3][i] - spt->rectangle[0][i];
    }

  crossProduct(aux1, aux2, aux3);
#ifdef DEBUG
  printf("%f %f %f | %f %f %f | %f %f %f\n",
      aux1[0], aux1[1], aux1[2], aux2[0], aux2[1], aux2[2], aux1[2], aux3[0], aux3[1], aux3[2]);
#endif
  spt->area = 0.0;
  for (i = 0; i < 3; ++i)
    spt->area += pow(aux3[i], 2.0);
  spt->area = sqrt(spt->area);

#ifdef DEBUG
  printf("Area = %f | %f %f %f | %f %f %f | %f %f %f\n", spt->area,
      spt->rectangle[0][0], spt->rectangle[0][1], spt->rectangle[0][2],
      spt->rectangle[1][0], spt->rectangle[1][1], spt->rectangle[1][2],
      spt->rectangle[3][0], spt->rectangle[3][1], spt->rectangle[3][2]);
#endif

  /* To simplify handling, get directions of the rectangle edges */
  for (i = 0; i < 4; ++i)
    {
      for (k = 0; k < 3; ++k)
        {
          rectdir[i][k] = spt->rectangle[(i + 1) % 4][k] - spt->rectangle[i][k];
        }
    }
  /* Step 3: Find angles between the sides of the rectangle and the next edges of the polygon */
  /*         Then rotate all sides of the rectangle (the calipers) by the smallest angle */
  /*         Recalculate rectangle and area. Repeat until all edges of the polygon have been visited */
  pivot[0] = minindy; /* Initial pivotal points (by index to corners) */
  pivot[1] = maxindx;
  pivot[2] = maxindy;
  pivot[3] = minindx;
  totang = 0.0;
  while (totang < 90.0)
    {
      minang = 365.0;
      for (i = 0; i < 4; ++i) /* Loop over edges of enclosing rectangle */
        {
          for (k = 0; k < 3; ++k)
            {
              //             aux1[k] = spt->rectangle[(i+1)%4][k] - spt->rectangle[i][k];
              aux1[k] = rectdir[i][k];
              aux2[k] = spt->corner[(pivot[i] + 1) % spt->nocorners][k] - spt->corner[pivot[i]][k];
            }
          normalize(aux1);
          normalize(aux2);
          tmp1 = dotProduct(aux1, aux2);
          if (tmp1 > 1.0)
            tmp1 = 1.0; /* acos will give NaN if > 1.0 bc. of numerical error */
          ang = acos(tmp1);
          if (ang < ZEROTOL2) /* Colinear edge: Move to next edge */
            {
#ifdef DEBUG
              printf("%d era cero\n", i);
#endif
              ++pivot[i];
              pivot[i] = pivot[i] % spt->nocorners;
              for (k = 0; k < 3; ++k)
                {
                  //                aux1[k] = spt->rectangle[(i+1)%4][k] - spt->rectangle[i][k];
                  aux2[k] = spt->corner[(pivot[i] + 1) % spt->nocorners][k] - spt->corner[pivot[i] % spt->nocorners][k];
                }
              //             normalize(aux1);
              normalize(aux2);
              tmp1 = dotProduct(aux1, aux2);
              if (tmp1 > 1.0)
                tmp1 = 1.0; /* acos will give NaN if > 1.0 bc. of numerical error */
              ang = acos(tmp1);
            }
          ang *= RAD_TO_DEG;
#ifdef DEBUG
          printf("Ang %d : %f\n", i, ang);
#endif
          if (ang < minang)
            minang = ang;
        }
      /* Rotate the calipers by the angle; at least one of them */
      /* will be colinear to one edge of the polygon */
      for (i = 0; i < 4; ++i)
        rotateVector(-minang, idx->surfnormal, rectdir[i]);
      totang += minang;
#ifdef DEBUG
      printf("Tot %f %f\n", totang, minang);
#endif

      /* Recalculate rectangle vertices */
      for (i = 0; i < 4; ++i)
        {
          tmp1 = (rectdir[i][0] * (spt->corner[pivot[(i + 1) % 4]][1] - spt->corner[pivot[i]][1]) - rectdir[i][1]
              * (spt->corner[pivot[(i + 1) % 4]][0] - spt->corner[pivot[i]][0])) / (rectdir[i][1] * rectdir[(i + 1) % 4][0]
              - rectdir[i][0] * rectdir[(i + 1) % 4][1]);
          for (j = 0; j < 3; ++j)
            {
              tmprect[i][j] = tmp1 * rectdir[(i + 1) % 4][j] + spt->corner[pivot[(i + 1) % 4]][j];
            }
        }
      /* Recalculate area */
      for (i = 0; i < 3; ++i)
        {
          aux1[i] = tmprect[1][i] - tmprect[0][i];
          aux2[i] = tmprect[3][i] - tmprect[0][i];
        }
      crossProduct(aux1, aux2, aux3);
      area = 0.0;
      for (i = 0; i < 3; ++i)
        area += pow(aux3[i], 2.0);
      area = sqrt(area);
#ifdef DEBUG
      printf("A %e | %f %f %f | %f %f %f | %f %f %f\n", area,
          tmprect[0][0], tmprect[0][1], tmprect[0][2], tmprect[1][0], tmprect[1][1], tmprect[1][2],
          tmprect[3][0], tmprect[3][1], tmprect[3][2]);
#endif
      /* Take note of the new rectangle if the area is smaller, */
      /* but not if there was numerical error */
      if (area < spt->area && area > ZEROTOL2)
        {
          spt->area = area;
          for (i = 0; i < 4; ++i)
            {
              for (j = 0; j < 3; ++j)
                {
                  spt->rectangle[i][j] = tmprect[i][j];
                }
            }
        }
    }

  /* --------------Allocate data and populate matrix by interpolation---------------- */
  // FIXME THIS COMMENT IS WRONG
  /* 1. Find the axes of the rectangle, in such a way that */
  /*    two vertices are on the x axis. This is achieved by calculating */
  /*    vertex coordinates in the rectangle basis */
  /*    and rotating if only one vertex is found */
  /*    (Ideally, one would one to have the longest side as x axis, */
  /*    but this is not always possible with the stated constraint */
  for (i = 0; i < 3; ++i)
    {
      aux1[i] = spt->rectangle[1][i] - spt->rectangle[0][i];
      aux2[i] = spt->rectangle[3][i] - spt->rectangle[0][i];
    }
  tmp1 = normalize(aux1); /* Lengths before normalization, in tmp vars */
  tmp2 = normalize(aux2);
  //   if(tmp1 > tmp2)
  //    {
  for (i = 0; i < 3; ++i)
    {
      spt->xaxis[i] = aux1[i];
      spt->yaxis[i] = aux2[i];
    }
  //    }
  //   else
  //    {
  //      for(i=0; i < 3; ++i)
  //       {
  //         spt->xaxis[i]=aux2[i];
  //         spt->yaxis[i]=aux1[i];
  //       }
  //      tmp3=tmp1;       /* Swap lengths of rectangle sides */
  //      tmp1=tmp2;
  //      tmp2=tmp3;
  //    }
  /* Find two corners in first rectangle edge */
  found = 0;
  //   round = 0;
  //   tolerance = ZEROTOL2;
  while (found < 2)
    {
      for (i = 0; i < spt->nocorners; ++i)
        {
          for (j = 0; j < 3; ++j)
            aux1[j] = spt->corner[i][j] - spt->rectangle[0][j];
          //         if((ycoord[i]=dotProduct(aux1,spt->yaxis)) < tolerance)    /* This is one */
          if ((ycoord[i] = dotProduct(aux1, spt->yaxis)) < ZEROTOL2) /* This is one */
            {
              xcoord[found] = dotProduct(aux1, spt->xaxis);
              spt->refcorner[found] = i;
              ++found;
              if (found == 2)
                {
                  /* The corners must be ordered so that refcorner[0] is "to the left". This must be checked */
                  /* Essentially, they always are, except if one is 0 and the other is the no of corners */
                  if (spt->refcorner[0] == 0 && spt->refcorner[1] == spt->nocorners - 1) /* Swap */
                    {
                      tmp3 = spt->refcorner[0];
                      spt->refcorner[0] = spt->refcorner[1];
                      spt->refcorner[1] = tmp3;
                      tmp3 = xcoord[0];
                      xcoord[0] = xcoord[1];
                      xcoord[1] = tmp3;
                    }
                  if (xcoord[0] > xcoord[1]) /* We need to invert to orientation of the cornes (clockwise to anti, etc) */
                    {
#ifdef DEBUG
                      printf("Invert %e %e\n", xcoord[0], xcoord[1]);
#endif
                      for (k = 0; k < spt->nocorners / 2; ++k)
                        {
                          for (l = 0; l < 3; ++l)
                            {
                              tmpi = spt->refcorner[0];
                              tmp3 = spt->corner[(tmpi - k) % spt->nocorners][l];
                              tmpi2 = spt->refcorner[1];
                              spt->corner[(tmpi - k) % spt->nocorners][l] = spt->corner[tmpi2 + k][l];
                              spt->corner[tmpi2 + k][l] = tmp3;
                            }
                        }
                    }
                  // The following line, commented out to be able to diagnose error
                  //               break;   /* We have two corners; no need for more */
                }
            }
#ifdef DEBUG
          printf("C %e %e\n",
              dotProduct(aux1,spt->xaxis),
              dotProduct(aux1,spt->yaxis));
#endif
        }
      // THIS SHOULD NEVER, REPEAT, NEVER HAPPEN !!!!!
      if (found > 2)
        {
          fprintf(stderr, "Error: too many corners found\n");
          exit(-1);
        }
      if (found < 2) /* Bad luck; only one corner on that edge */
        {
#ifdef DEBUG
          printf("Failed!\n"); // Not really; just this one time
#endif
          //          ++round;
          /* Rotate rectangle (=> also change sign of axis vectors) and try again */
          /* This should work in a maximum of 3 rotations */
          for (i = 0; i < 3; ++i)
            {
              /* Rotate vertices (equivalent to 90 degrees) */
              tmp3 = spt->rectangle[0][i];
              for (j = 1; j < 4; ++j)
                spt->rectangle[j - 1][i] = spt->rectangle[j][i];
              spt->rectangle[3][i] = tmp3;
              /* The operation below implements the necessary rotation (and reflection) of axes */
              tmp3 = spt->xaxis[i];
              spt->xaxis[i] = spt->yaxis[i];
              spt->yaxis[i] = -tmp3;
              tmp3 = tmp1; /* Swap lengths of rectangle sides */
              tmp1 = tmp2;
              tmp2 = tmp3;
            }
          found = 0; /* Second time should work... */
          //          if(!(round%2))  tolerance *= 1.5;     /* ...but if doesn't, be more tolerant */
        }
    }
#ifdef DEBUG
  printf("Winners = %d %d\n", spt->refcorner[0], spt->refcorner[1]);
#endif

  /* 2. Decide the steps in each direction (relative to the rectangle edges) */
  /*    The algorithm is as follows: The major and minor edges of the rectangle are projected */
  /*    onto the x, y and z axes and the number of points that would correspond to that projection, */
  /*    using the spacing of the original grid is calculated. Each rectangle edge gets the spacing */
  /*    of the cell axis that produces most points. This guarantees that all intervals are */
  /*    sampled at least once and avoids visual discontinuities */
  /* Choose x axis spacing (on rectangle) */
  maxpt = 0.0;
  for (i = 0; i < 3; ++i)
    {
      tmppt = tmp1 * fabs(dotProduct(spt->xaxis, uaxis[i])) / d->del[i];
#ifdef DEBUG
      printf("TM %f\n", tmppt);
#endif
      if (tmppt > maxpt)
        {
          maxpt = tmppt;
          spt->dx = d->del[i];
        }
    }
  /* Choose y axis spacing (on rectangle) */
  maxpt = 0.0;
  for (i = 0; i < 3; ++i)
    {
      tmppt = tmp2 * fabs(dotProduct(spt->yaxis, uaxis[i])) / d->del[i];
#ifdef DEBUG
      printf("TM %f\n", tmppt);
#endif
      if (tmppt > maxpt)
        {
          maxpt = tmppt;
          spt->dy = d->del[i];
        }
    }
#ifdef DEBUG
  printf("INTER %f %f\n", spt->dx, spt->dy);
#endif
  /* 3. Calculate number of interpolation points per dimension */
  /*    This is the number of intervals, plus two points to treat */
  /*    the borders as a special case (possibly redundant) */
  /*    Then, readjust the intervals so that the number of grid points is exact */
  spt->pt[X] = 2 + (int) floor(tmp1 / spt->dx);
  spt->pt[Y] = 2 + (int) floor(tmp2 / spt->dy);
#ifdef DEBUG
  printf("Puntos = %d %d | %e %e | %e %e\n", spt->pt[X], spt->pt[Y],
      tmp1, tmp2, spt->dx, spt->dy);
#endif
  spt->dx = tmp1 / (spt->pt[X] - 1);
  spt->dy = tmp2 / (spt->pt[Y] - 1);
#ifdef DEBUG
  printf("Puntos = %d %d | %e %e | %e %e\n", spt->pt[X], spt->pt[Y],
      tmp1, tmp2, spt->dx, spt->dy);
#endif
  /* 4. Rescale internal axes accordingly (so that the vectors represent exactly one step) */
  for (i = 0; i < 3; ++i)
    {
      spt->xaxis[i] *= spt->dx;
      spt->yaxis[i] *= spt->dy;
    }
  /* 5. Allocate space for data and border info. */
  spt->xbord[LEFT] = (double *) malloc((size_t) spt->pt[Y] * sizeof(double));
  spt->xbord[RIGHT] = (double *) malloc((size_t) spt->pt[Y] * sizeof(double));
  spt->lim[LEFT] = (int *) malloc((size_t) spt->pt[Y] * sizeof(int));
  spt->lim[RIGHT] = (int *) malloc((size_t) spt->pt[Y] * sizeof(int));
  spt->data = (double **) malloc((size_t) spt->pt[X] * sizeof(double *));
  for (i = 0; i < spt->pt[X]; ++i)
    {
      spt->data[i] = (double *) malloc((size_t) spt->pt[Y] * sizeof(double));
    }
//  /* 6. Prepare openGL list for display */
//  spt->slice = glGenLists(1); /* Reserve list for precompiled openGL object */
//#ifdef DEBUG
//  printf("Slice = %d\n", spt->slice);
//#endif
//
//  /* 7. Interpolate data on the intersecting surface and make list for display */
  /* 6. Interpolate data on the intersecting surface and make list for display */
  populateData(idx, b, d);

  /* The point was to draw the slice, so do it*/
  drawSlice(idx, d);
}

void populateData(index_t *idx, box_t *b, density_t *d)
{
  int i, j, k;
  int tmp1;
  double pt[3];
  double tmp[3];
  double d00, d10, d01, d11, d0, d1;
  int ref[3];
  int jumpline;
  int lcorner, rcorner; /* Indexes of corners were relevant polygons edges start */
  double ldir[3], rdir[3]; /* Directions of those edges */
  double lproj[2], rproj[2]; /* Projetions on the rectangle (i.e. rectangle coordinates) */
  double llen, rlen;
  densitycut_t *spt; /* Pointer to the slice (this mainly makes the code more readable) */
  int outside;

  spt = idx->dindex[idx->curr];
  /* Directions of relevant (intersenting) edges */
  lcorner = spt->refcorner[0];
  rcorner = spt->refcorner[1];
  for (i = 0; i < 3; ++i)
    {
#ifdef DEBUG
      printf("Aqui %d %d %d %d\n", lcorner, rcorner,
          (lcorner-1+spt->nocorners) % spt->nocorners,
          (rcorner+1) % spt->nocorners);
#endif
      ldir[i] = spt->corner[(lcorner - 1 + spt->nocorners) % spt->nocorners][i] - spt->corner[lcorner][i];
      rdir[i] = spt->corner[(rcorner + 1) % spt->nocorners][i] - spt->corner[rcorner][i];
    }

  /* Reset info on "turns" */
  spt->turn[LEFT][LINE] = -1; /* It means unset */
  spt->turn[LEFT][CORNER] = -1;
  spt->turn[RIGHT][LINE] = -1;
  spt->turn[RIGHT][CORNER] = -1;

//  resetIndex(NULL,NULL,idx,d,b);
//  exit(0);
  ///* The last line has to be treated differently */
  //   for(i=0; i < spt->pt[Y]-1; ++i)
  for (i = 0; i < spt->pt[Y]; ++i)
    {
      jumpline = NO;
      spt->lim[LEFT][i] = 0;
      spt->lim[RIGHT][i] = spt->pt[X] - 1;
      /* Move to next line, over the surface of the rectangle */
      for (j = 0; j < 3; ++j)
        pt[j] = spt->rectangle[0][j] + spt->yaxis[j] * i;

#ifdef DEBUG
      printf("%d ------------------------\n", i);
#endif
      spt->xbord[LEFT][i] = findBorder(pt, spt, ldir, lcorner, &outside);
      if (outside == YES) /* This means we past a corner, so we need to change the reference */
        {
          lcorner = (lcorner - 1 + spt->nocorners) % spt->nocorners;
          for (j = 0; j < 3; ++j)
            {
              ldir[j] = spt->corner[(lcorner - 1 + spt->nocorners) % spt->nocorners][j] - spt->corner[lcorner][j];
            }
          //printf("LCOR %d %d | %f %f %f\n", lcorner, (lcorner-1+spt->nocorners) % spt->nocorners, ldir[0], ldir[1], ldir[2]);
          spt->xbord[LEFT][i] = findBorder(pt, spt, ldir, lcorner, &outside);
          /* We allow to be outside only if it is the last row; in that case, just ignore */
          if (outside && i != spt->pt[Y] - 1)
            {
              fprintf(stderr, "Confusing edges: This is probably a bug\n");
              exit(-1);
            }
          else /* We have identified a "turn", i.e. a change of edge */
            {
              spt->turn[LEFT][LINE] = i - 1; /* We use the previous line */
              spt->turn[LEFT][CORNER] = lcorner;
            }
        }
      spt->xbord[RIGHT][i] = findBorder(pt, spt, rdir, rcorner, &outside);
      if (outside == YES)
        {
          rcorner = (rcorner + 1) % spt->nocorners;
          for (j = 0; j < 3; ++j)
            {
              rdir[j] = spt->corner[(rcorner + 1) % spt->nocorners][j] - spt->corner[rcorner][j];
            }
          spt->xbord[RIGHT][i] = findBorder(pt, spt, rdir, rcorner, &outside);
          /* We allow to be outside only if it is the last row; in that case, just ignore */
          if (outside && i != spt->pt[Y] - 1)
            {
              fprintf(stderr, "Confusing edges\n");
              exit(-1);
            }
          else
            {
              spt->turn[RIGHT][LINE] = i - 1;
              spt->turn[RIGHT][CORNER] = rcorner;
            }
        }
      //      spt->lim[LEFT][i]=(int)ceil(spt->xbord[LEFT][i]);
      spt->lim[LEFT][i] = (int) ceil(spt->xbord[LEFT][i]);
      spt->lim[RIGHT][i] = (int) floor(spt->xbord[RIGHT][i]);
#ifdef DEBUG
      printf("Bint %d %e | %d %e\n", spt->lim[LEFT][i], spt->xbord[LEFT][i], spt->lim[RIGHT][i], spt->xbord[RIGHT][i]);
#endif
      for (j = 0; j < spt->pt[X]; ++j)
        spt->data[j][i] = 0.0;
      for (j = 0; j < 3; ++j)
        tmp[j] = pt[j];
      //      for(j=0; j < spt->pt[X]; ++j)
      /* First, interpolate the borders and store the result at the borders of the rectangle */
      /* This way, if the limit and the borders are at the same point, we only store once */
      /* It also simplifies border plotting, because border information is always at the same place */
      for (j = 0; j < 3; ++j)
        pt[j] = tmp[j] + spt->xaxis[j] * spt->xbord[LEFT][i];
      spt->data[0][i] = trilinearInterpolation(pt, d);
      for (j = 0; j < 3; ++j)
        pt[j] = tmp[j] + spt->xaxis[j] * spt->xbord[RIGHT][i];
      spt->data[spt->pt[X] - 1][i] = trilinearInterpolation(pt, d);
      /* Then, interpolate points on the regular grid */
      //      for(j=spt->lim[LEFT][i]+1; j < spt->lim[RIGHT][i]; ++j)
      for (j = 0; j < 3; ++j)
        pt[j] = tmp[j] + spt->xaxis[j] * spt->lim[LEFT][i];
      for (j = spt->lim[LEFT][i]; j <= spt->lim[RIGHT][i]; ++j)
        {
          /* Trilinear interpolation */
          spt->data[j][i] = trilinearInterpolation(pt, d);
          /* Move to next point on x axis, over the surface of the rectangle */
          for (k = 0; k < 3; ++k)
            pt[k] += spt->xaxis[k];
        }
      /* Finally, interpolate the corners that participate in turns */
      for (j = LEFT; j <= RIGHT; ++j)
        if (spt->turn[j][CORNER] != -1)
          spt->cturn[j] = trilinearInterpolation(spt->corner[spt->turn[j][CORNER]], d);
      /* If the slice is a triangle, we are likely to need to interpolate the last corner explicitly */
      if (spt->nocorners == 3)
        {
          for (j = 0; j < 3; ++j)
            {
#ifdef DEBUG
              printf("Pasa %d %d %d\n", j, spt->refcorner[0], spt->refcorner[1]);
#endif
              if (spt->refcorner[0] != j && spt->refcorner[1] != j)
                {
                  tmp1 = j;
                  spt->cturn[2] = trilinearInterpolation(spt->corner[tmp1], d);
                }
            }
        }
    }
#ifdef DEBUG
  printf("Turns : %d %d %d %d\n",
      spt->turn[LEFT][LINE], spt->turn[LEFT][CORNER],
      spt->turn[RIGHT][LINE], spt->turn[RIGHT][CORNER]);
#endif
  /*
   for(i=0; i < spt->pt[Y]; ++i)
   printf("Lim %d %d %d\n", i, spt->lim[LEFT][i], spt->lim[RIGHT][i]);
   */
}

/* Parameters:
 pt  -> The starting point of the row of which we want to calculate the borders
 spt -> Pointer to the slice that contains the data
 dir -> Direction of the edge that intersects the current row;  the intersection is the border of the polygon
 Note: This vector has to have the length of the edge, so that when expressed in parametric units
 the parameter is btw 0 and 1, > 1 meaning that the next edge should be used
 corner -> Index of the corner from which the edge starts
 flag   -> This will be true if the edge is intersected outside of the data cell
 Returns the border's coordinate (in rectangle space)
 */
double findBorder(double *pt, densitycut_t *spt, double *dir, int corner, int *flag)
{
  /* Rectangle coordinates */
  double c[2]; /* Projection of the corner onto the rectangle */
  double d[2]; /* Projection of dir onto the rectangle */
  double p[2]; /* Projection of row origin (pt) onto the rectangle */

  double l1, l2; /* Solution of the intersection in parametric form */

  /* Project direction of the edge onto the rectangle */
  c[0] = dotProduct(spt->corner[corner], spt->xaxis) / spt->dx;
  c[1] = dotProduct(spt->corner[corner], spt->yaxis) / spt->dy;
  d[0] = dotProduct(dir, spt->xaxis) / spt->dx;
  d[1] = dotProduct(dir, spt->yaxis) / spt->dy;
  p[0] = dotProduct(pt, spt->xaxis) / spt->dx;
  p[1] = dotProduct(pt, spt->yaxis) / spt->dy;
  //printf("%e %e | %e %e | %e %e | %e %e\n", c[0], c[1], d[0], d[1], p[0], p[1], spt->dx, spt->dy);

  /* Calculate intersection */
  l1 = (d[0] * (p[1] - c[1]) - d[1] * (p[0] - c[0])) / (d[1] * spt->dx);
  if (fabs(d[0]) > ZEROTOL2)
    l2 = (spt->dx * l1 + p[0] - c[0]) / d[0];
  else
    l2 = (p[1] - c[1]) / d[1];
  if (l2 > 1.0 + ZEROTOL2)
    *flag = YES;
  else
    *flag = NO;

#ifdef DEBUG
  printf("Bor %e %e\n", l1, l2);
#endif
  return l1;
}

