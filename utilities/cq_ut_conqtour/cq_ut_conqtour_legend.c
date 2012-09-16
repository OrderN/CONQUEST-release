/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_legend.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_legend.c *
 *
 * NAME
 *  cq_ut_conqtour_legend.c
 * PURPOSE
 *  Functions to create a legend for the density
 * USES
 *  
 * AUTHOR
 *  torralba
 * CREATION DATE
 *  Feb 15, 2011
 * MODIFICATION HISTORY
 *
 * *****/

#include <string.h>
#include "cq_ut_conqtour.h"

void makeLegend(GLuint *legend, ctrl_t *c, density_t *d)
{
  int i, j, n;
  int borderstep; /* Use to create a gradation for borders thicker than one pixel */
  int tmp;
  legend_t *l;
  int pos; /* Position within the image array */
  int size; /* Size of the image */
  char *image; /* The legend bar as an RGBA image (i.e. RGB with transparency) */
  rgbcolor_t rgb; /* A color to be displayed */
  double val;
  GLdouble windowpos[4]; /* Inintial window position */
  double totalw, totalh; /* Total width and height, including borders and ticks */
  //  scalemarks_t marks;
  //  int bsx, bsy; /* Border shifts (i.e. thickness) */

  l = &(d->legend);
  //  l->bt = 2; /* Border thickness */
  //  l->tl = 4; /* Tick length */
  totalw = l->w + 2 * l->bt + l->tl;
  totalh = l->h + 2 * l->bt;
  size = totalw * totalh * 4;
  image = (char *) malloc((size_t) size * sizeof(char));
  memset(image, 0, size);
  allocateMarks(&(l->marks), l->h);

  if (glIsList(*legend) == GL_TRUE)
    glDeleteLists(*legend, 1);
  *legend = glGenLists(1);
  glNewList(*legend, GL_COMPILE_AND_EXECUTE);
  glAlphaFunc(GL_GEQUAL, 0.5);
  glEnable(GL_ALPHA_TEST);
  decideMarks(&(l->marks), c, d);
  /* Make border first */
  borderstep = (int) (255. / (l->bt + 1));
  for (i = 0; i < l->bt; ++i)
    {
      drawRectangle(image, totalw, totalh, i, i, l->w + 2 * l->bt - i - 1, l->h + 2 * l->bt - i - 1, i * borderstep);
    }
  //  printf("Row %e %e %d %d\n", totalw, totalh, l->h, image);
  /* Now, fill the scale with colors, plus add ticks and labels */
  for (i = 0; i < l->h; ++i)
    {
      //printf("%d\n", i);
      getColorFromScale(d, l->marks.val[i], &rgb);
      for (j = 0; j < l->w; ++j)
        {
          pos = ((i + l->bt) * totalw + j + l->bt) * 4;
          if (rgb.good == YES)
            {
              /* Color */
              image[pos] = rgb.r;
              image[pos + 1] = rgb.g;
              image[pos + 2] = rgb.b;
              image[pos + 3] = 255; /* Opaque */
            }
          else /* This color is out of range : Represent as a striped pattern */
            {
              if ((i + j) % 4 == 0)
                {
                  image[pos] = 0;
                  image[pos + 1] = 0;
                  image[pos + 2] = 0;
                  image[pos + 3] = 255; /* Opaque */
                }
              else
                {
                  image[pos] = 128;
                  image[pos + 1] = 128;
                  image[pos + 2] = 128;
                  image[pos + 3] = 255; /* Opaque */
                }
            }
        }
      if (l->marks.tick[i] == YES)
        {
          for (j = l->w + 2 * l->bt; j < totalw; ++j)
            {
              pos = ((i + l->bt) * totalw + j) * 4;
              /* Tick */
              image[pos] = 0;
              image[pos + 1] = 0;
              image[pos + 2] = 0;
              image[pos + 3] = 255; /* Opaque */
              pos = ((i + l->bt - 1) * totalw + j) * 4;
              /* Shadow */
              image[pos] = 128;
              image[pos + 1] = 128;
              image[pos + 2] = 128;
              image[pos + 3] = 255; /* Opaque */
            }
        }
    }
  //  printf("Row %f %f %p\n", totalw, totalh, image);
  //  for(i = 0 ; i < totalw * totalh; ++i)
  //    image[i] = 0;
  glDrawPixels(totalw, totalh, GL_RGBA, GL_UNSIGNED_BYTE, image);
  glEndList();
  free(image);
}

void drawRectangle(char *image, int sx, int sy, int x0, int y0, int x1, int y1, int color)
{
  int i, j;
  int pos;

  for (i = y0; i <= y1; ++i)
    {
      pos = (i * sx + x0) * 4;
      image[pos] = color;
      image[pos + 1] = color;
      image[pos + 2] = color;
      image[pos + 3] = 255; /* Opaque */
    }
  for (i = x0; i <= x1; ++i)
    {
      pos = (y1 * sx + i) * 4;
      image[pos] = color;
      image[pos + 1] = color;
      image[pos + 2] = color;
      image[pos + 3] = 255; /* Opaque */
    }
  for (i = y1; i >= y0; --i)
    {
      pos = (i * sx + x1) * 4;
      image[pos] = color;
      image[pos + 1] = color;
      image[pos + 2] = color;
      image[pos + 3] = 255; /* Opaque */
    }
  for (i = x1; i >= x0; --i)
    {
      pos = (y0 * sx + i) * 4;
      image[pos] = color;
      image[pos + 1] = color;
      image[pos + 2] = color;
      image[pos + 3] = 255; /* Opaque */
    }
}

void allocateMarks(scalemarks_t *s, int p)
{
  s->pixels = p;
  s->val = (double *) malloc(p * sizeof(double));
  s->tick = (char *) malloc(p * sizeof(char));
  s->label = (char *) malloc(p * sizeof(char));
  if (s->val == NULL || s->tick == NULL || s->label == NULL)
    {
      fprintf(stderr, "Error allocating scale marks\n");
      exit(-1);
    }
}

void deallocateMarks(scalemarks_t *s)
{
  if (s->val != NULL)
    free(s->val);
  if (s->tick != NULL)
    free(s->tick);
  if (s->label != NULL)
    free(s->label);
}

void deallocateIntervals(density_t *d)
{
  int i;

  if (d->legend.nointer)
    {
      if (d->legend.val != NULL)
        free(d->legend.val);
      if (d->legend.col != NULL)
        free(d->legend.col);
    }
}

/**** Marks (ticks and labels) ****/
void decideMarks(scalemarks_t *marks, ctrl_t *c, density_t *d)
{
  int i;

  /* Get the real value for this row */
  switch (FLAG_UP(d->legend.method, CQ_LEGEND_TICKS_MASK))
    {
  case CQ_LEGEND_TICKS_LINEAR:
    makeLinearMarks(marks, c, d);
    break;
  case CQ_LEGEND_TICKS_SIGNABSLOG:
    //    makeLinearMarks(marks, c, d);
    makeSignedAbsoluteLogMarks(marks, c, d);
    break;
  default:
    if (c->verbose)
      printf("WARNING: Using default method for legend ticks\n");
    SET_MULTIBIT_FLAG(c->legend.method, CQ_LEGEND_DFLT_METHOD_TICKS, CQ_LEGEND_TICKS_MASK);
    decideMarks(marks, c, d);
    }
}

void makeLinearMarks(scalemarks_t *m, ctrl_t *c, density_t *d)
{
  int i;
  double step; /* Value step from row to row of the scale */
  int row; /* A given pixel row in the scale */
  double tmpval;

  //  d->legend.tickinterval = 0.1;
  step = (d->legend.max - d->legend.min) / d->legend.h;
  //printf("St %e %e\n", step, d->legend.tickinterval);
  for (i = 0; i < m->pixels; ++i)
    {
      /* Just a linear interpolation */
      m->val[i] = d->legend.min + i * step;
      m->tick[i] = NO; /*For the moment */
      m->label[i] = NO;
    }
  /* Use origin to mark first tick */
  row = (int) ((d->legend.tickorigin - d->legend.min) / step);
  if (row < 0 || row >= m->pixels)
    {
      fprintf(stderr, "Error creating scale: Origin out of bounds = %d\n", row);
      exit(-1);
    }
  m->val[row] = d->legend.tickorigin; /* Use a clear number as a more readable approximation */
  m->label[row] = YES;
  m->tick[row] = YES;
  if (d->legend.tickinterval < ZEROTOL) /* The user did not provide an interval : Calculate one using no. ticks suggestion */
    {
      d->legend.tickinterval = (d->legend.max - d->legend.min) / (d->legend.noticks - 1);
      //      printf("In %e\n", d->legend.tickinterval);
    }
  d->legend.noticks = 1; /* Count the actual number of ticks */
  /* Mark "positives" (i.e. above the origin) */
  tmpval = d->legend.tickorigin + d->legend.tickinterval;
  while (tmpval <= d->legend.max)
    {
      row = (int) ((tmpval - d->legend.min) / step);
      if (row >= 0 && row < m->pixels)
        {
          m->val[row] = tmpval;
          m->label[row] = YES;
          m->tick[row] = YES;
          ++d->legend.noticks;
        }
      else
        {
          if (c->verbose)
            {
              printf("WARNING: Tick out of bounds (pixel = %d) ignored\n", row);
            }
        }
      tmpval += d->legend.tickinterval;
    }
  /* Mark "negatives" (i.e. below the origin) */
  tmpval = d->legend.tickorigin - d->legend.tickinterval;
  while (tmpval >= d->legend.min)
    {
      row = (int) ((tmpval - d->legend.min) / step);
      if (row >= 0 && row < m->pixels)
        {
          m->val[row] = tmpval;
          m->label[row] = YES;
          m->tick[row] = YES;
          ++d->legend.noticks;
        }
      else
        {
          if (c->verbose)
            {
              printf("WARNING: Tick out of bounds (pixel = %d) ignored\n", row);
            }
        }
      tmpval -= d->legend.tickinterval;
    }
  //  m->tick[0] = YES;
}

void makeSignedAbsoluteLogMarks(scalemarks_t *m, ctrl_t *c, density_t *d)
{
  int i, j;
  //  int nointer;
  double step; /* Value step from row to row of the scale */
  int row; /* A given pixel row in the scale */
  double tmpval;
  //  int ordpmin, ordpmax, ordnmin, ordnmax; /* Minimum and maximum positive and negative magnitude orders */
  double ordpmin, ordpmax, ordnmin, ordnmax; /* Minimum and maximum positive and negative magnitude orders */
  int negpt, pospt; /* Scale rows used to represent negative and positive values, respectively */
  int negrow, posrow; /* Rows where the negative and positive values start, i.e. lowest absolute values */
  double exponent; /* Exponent to be represented */
  int sign;
  double tmporig; /* Used to keep tick origins */

  /* Calculate the number of intervals and step per scale row */
  if (d->pmin >= VERYLARGE) /* This indicates that no negative values are present */
    {
      ordpmin = 0.0;
      ordpmax = 0.0;
    }
  else
    {
      ordpmin = log10(d->pmin);
      ordpmax = log10(d->pmax);
    }
  if (d->nmin >= VERYLARGE)
    {
      ordnmin = 0.0;
      ordnmax = 0.0;
    }
  else
    {
      ordnmin = log10(d->nmin); /* nmin and nmax are already in absolute value */
      ordnmax = log10(d->nmax);
    }
  //  nointer = 1;
  //  if ((ordnmax - ordnmin) > 0)
  //    nointer += 1 + (ordnmax - ordnmin);
  //  if ((ordpmax - ordpmin) > 0)
  //    nointer += 1 + (ordpmax - ordpmin);
  step = (ordnmax - ordnmin + ordpmax - ordpmin) / (d->legend.h - CQ_LEGEND_POINTS_BETWEEN_SIGNS - 2);
  sign = -1;
  //printf("Step %e\n", step);
  //  exit(-1);
  negpt = (int) ((double) (ordnmax - ordnmin)) / step;
  pospt = (int) ((double) (ordpmax - ordpmin)) / step;
  exponent = ordnmax;
  step = -step; /* Go toward more negative exponents */
  //printf("%e %e %e %e %e\n", exponent, ordnmax, ordnmin, ordpmin, ordpmax);
  for (i = 0; i < m->pixels; ++i)
    {
      //printf("Ex %d %e\n", i, exponent);
      m->val[i] = sign * pow(10.0, exponent);
      exponent += step;
      if (exponent < ordnmin - ZEROTOL2 && step < 0.0) /* We went into the dark zone!! No data here: show stripes to indicate */
        {
          negrow = i;
          //          printf("Ex-------- %e\n", exponent);
          if (i + CQ_LEGEND_POINTS_BETWEEN_SIGNS >= m->pixels)
            {
              fprintf(stderr, "Error creating scale: Tick out of bounds\n");
              exit(-1);
            }
          for (j = i; j <= i + CQ_LEGEND_POINTS_BETWEEN_SIGNS; ++j)
            {
              m->val[j] = 0.0;
              if (j - i == ceil((double) CQ_LEGEND_POINTS_BETWEEN_SIGNS / 2.0))
                {
                  m->tick[j] = YES;
                  m->label[j] = NO; /* Don't label zero */
                }
              else
                {
                  m->tick[j] = NO;
                  m->label[j] = NO;
                }
            }
          exponent = ordpmin;
          step = -step; /* Go toward more positive exponents */
          sign = 1; /* And produce positive numbers */
          i += CQ_LEGEND_POINTS_BETWEEN_SIGNS;
          posrow = i;
        }
      m->tick[i] = NO; /*For the moment */
      m->label[i] = NO;
    }

  /* Decide the tick interval */
  /* IMPORTANT: In this case the tick interval is interpreted as the EXPONENT interval, no the linear value */
  if (d->legend.tickinterval < ZEROTOL) /* The user did not provide an interval : Calculate one */
    {
      d->legend.tickinterval = ((ordpmax - ordpmin) + (ordnmax - ordnmin)) / (d->legend.noticks - 1);
      //printf("In2 %e | %e %e %e %e | %e\n", d->legend.tickinterval, d->pmax, d->pmin, d->nmax, d->nmin, step);
    }
  //d->legend.tickinterval = 1;
  /**** Ticks for positive data ****/
  if (d->pmin < VERYLARGE) /* "If we have data": pmin = VERYLARGE if no data was read from file */
    {
      /* Decide positive origin */
      if (d->legend.tickorigin < d->pmin + ZEROTOL || d->legend.tickorigin > d->pmax) /* Out of range */
        tmporig = ceil(log10(d->pmin)); /* First whole order of magnitude above the positive minimum */
      else
        tmporig = log10(d->legend.tickorigin);
      //printf("Or+ %e\n", tmporig);
      //exit(-1);
      d->legend.noticks = 1; /* Count the actual number of ticks */
      /* Marks above the positive origin) */
      //tmporig = -6;
      tmpval = tmporig;
      while (tmpval <= ordpmax)
        {
          row = 1 + posrow + (int) ((tmpval - ordpmin) / step);
          if (row < 0 || row >= m->pixels)
            {
              fprintf(stderr, "Error creating scale: Tick out of bounds\n");
              exit(-1);
            }
          m->val[row] = pow(10.0, tmpval);
          m->label[row] = YES;
          m->tick[row] = YES;
          ++d->legend.noticks;
          tmpval += d->legend.tickinterval;
        }
      /* Marks below the positive origin */
      tmpval = tmporig - d->legend.tickinterval;
      while (tmpval >= ordpmin)
        {
          row = 1 + posrow + (int) ((tmpval - ordpmin) / step);
          if (row < 0 || row >= m->pixels)
            {
              fprintf(stderr, "Error creating scale: Tick out of bounds\n");
              exit(-1);
            }
          m->val[row] = pow(10.0, tmpval);
          m->label[row] = YES;
          m->tick[row] = YES;
          ++d->legend.noticks;
          tmpval -= d->legend.tickinterval;
        }
    }
  /**** Ticks for negative data ****/
  if (d->nmin < VERYLARGE) /* "If we have data": nmin = VERYLARGE if no data was read from file */
    {
      //printf("OOO\n");
      /* Decide negative origin */
      if (d->legend.tickorigin < d->nmin + ZEROTOL || d->legend.tickorigin > d->nmax) /* Out of range */
        tmporig = ceil(log10(d->nmin)); /* First whole order of magnitude above (in abs. val.) the negative minimum */
      else
        tmporig = log10(d->legend.tickorigin);
      //printf("Or++ %e\n", tmporig);
      //exit(-1);
      /* Marks below the negative origin */
      //tmporig = -6;
      tmpval = tmporig;
      while (tmpval <= ordnmax)
        {
          row = negrow - 1 - (int) ((tmpval - ordnmin) / step);
          if (row < 0 || row >= m->pixels)
            {
              fprintf(stderr, "Error creating scale: Tick out of bounds\n");
              exit(-1);
            }
          m->val[row] = -pow(10.0, tmpval);
          m->label[row] = YES;
          m->tick[row] = YES;
          ++d->legend.noticks;
          tmpval += d->legend.tickinterval;
        }
      /* Marks above the negative origin */
      tmpval = tmporig - d->legend.tickinterval;
      while (tmpval >= ordnmin)
        {
          row = negrow - 1 - (int) ((tmpval - ordnmin) / step);
          if (row < 0 || row >= m->pixels)
            {
              fprintf(stderr, "Error creating scale: Tick out of bounds\n");
              exit(-1);
            }
          m->val[row] = -pow(10.0, tmpval);
          m->label[row] = YES;
          m->tick[row] = YES;
          ++d->legend.noticks;
          tmpval -= d->legend.tickinterval;
        }
      //  printf("Finis\n");
    }
}

/**** Color interval methods ****/

void makeLinearColorIntervals(ctrl_t *c, density_t *d)
{
  int i;
  double colorinterval;

  deallocateIntervals(d); /* Avoid memory leaks */
  /* To create color intervals, use hint for number of ticks */
  /* NOTE: The suggested number of ticks is not necessarily the same as the actual number */
  /*       This will depend on the ticks origin and the ticks interval */
  d->legend.nointer = c->legend.nointer; /* Use the control-defined value */
  colorinterval = (d->legend.max - d->legend.min) / d->legend.nointer;
  d->legend.val = (double *) malloc((size_t) (1 + d->legend.nointer) * sizeof(double));
  d->legend.col = (rgbcolor_t *) malloc((size_t) (1 + d->legend.nointer) * sizeof(rgbcolor_t));
  if (d->legend.val == NULL || d->legend.col == NULL)
    {
      fprintf(stderr, "Error allocating legend data\n");
      exit(-1);
    }
  for (i = 0; i <= d->legend.nointer; ++i)
    {
      d->legend.val[i] = d->legend.min + (colorinterval) * i;
      //printf("CI %e\n", d->legend.val[i]);
    }
}

void makeSignedAbsoluteLogColorIntervals(ctrl_t *c, density_t *d)
{
  int i, j;
  int ordpmin, ordpmax, ordnmin, ordnmax; /* Minimum and maximum positive and negative magnitude orders */

  deallocateIntervals(d); /* Avoid memory leaks */
  /* To create color intervals, use the orders of magnitude spanned by the data */
  /* This is almost a logarithmic scale, except that negative numbers are also included */
  if (d->pmin >= VERYLARGE) /* This indicates that no negative values are present */
    {
      ordpmin = 0;
      ordpmax = 0;
    }
  else
    {
      ordpmin = (int) floor(log10(d->pmin));
      ordpmax = (int) ceil(log10(d->pmax));
    }
  if (d->nmin >= VERYLARGE)
    {
      ordnmin = 0;
      ordnmax = 0;
    }
  else
    {
      ordnmin = (int) floor(log10(d->nmin)); /* nmin and nmax are already in absolute value */
      ordnmax = (int) ceil(log10(d->nmax));
    }
  d->legend.nointer = 1;
  if ((ordnmax - ordnmin) > 0)
    d->legend.nointer += 1 + (ordnmax - ordnmin);
  if ((ordpmax - ordpmin) > 0)
    d->legend.nointer += 1 + (ordpmax - ordpmin);
  //  printf("%d | %e %e %e\n", d->legend.nointer, d->pmin, d->nmin, VERYLARGE - 1);
  d->legend.val = (double *) malloc((size_t) (1 + d->legend.nointer) * sizeof(double));
  d->legend.col = (rgbcolor_t *) malloc((size_t) (1 + d->legend.nointer) * sizeof(rgbcolor_t));
  if (d->legend.val == NULL || d->legend.col == NULL)
    {
      fprintf(stderr, "Error allocating legend data\n");
      exit(-1);
    }
  j = 0;
  if ((ordnmax - ordnmin) > 0)
    for (i = ordnmax; i >= ordnmin; --i)
      {
        d->legend.val[j] = -pow(10.0, i);
        //        printf("CI %e\n", d->legend.val[j]);
        ++j;
      }
  d->legend.val[j] = 0.0;
  //  printf("CI %e\n", d->legend.val[j]);
  ++j;
  if ((ordpmax - ordpmin) > 0)
    for (i = ordpmin; i <= ordpmax; ++i)
      {
        d->legend.val[j] = pow(10.0, i);
        //        printf("CI %e\n", d->legend.val[j]);
        ++j;
      }
}

/**** Color methods ****/
void makeSaturatedRedBlueColors(ctrl_t *c, density_t *d)
{
  int i;
  int border; /* Point at which sign of the interval changes */
  //  double low, high; /* Lowest and highest values */
  //  double middle; /* The value of the "border point" (from blue to red); usually 0, but it could be low or high */

  if (d->legend.nointer == 0)
    {
      if (c->verbose)
        printf("WARNING: makeSaturatedRedBlueColors call ignored: Undefined color intervals\n");
      return;
    }
  /* Ordered intervals assumed, of course */
  /* Find the point where sign changes */
  border = 0;
  while (d->legend.val[border] < ZEROTOL && border < d->legend.nointer)
    ++border;
  for (i = 0; i <= d->legend.nointer; ++i)
    {
      if (i < border)
        {
          d->legend.col[i].r = 0;
          //          d->legend.col[i].g = (int) (255 * ((double)(border - i) / (double)border));
          d->legend.col[i].g = (int) (255 * ((double) (i) / (double) border));
          d->legend.col[i].b = 255;
        }
      else
        {
          d->legend.col[i].r = 255;
          //          d->legend.col[i].g = (int) (255 * ((double)(i - border) / (double)(d->legend.nointer - border)));
          d->legend.col[i].g = (int) (255 * ((double) (d->legend.nointer - i) / (double) (d->legend.nointer - border)));
          d->legend.col[i].b = 0;
        }

    }
  //  low = d->legend.val[0];
  //  high = d->legend.val[d->legend.nointer];
  //  if (low > 0.0)
  //    {
  //      middle = low;
  //    }
  //  else if (high < 0.0)
  //    {
  //      middle = high;
  //    }
  //  else
  //    middle = 0.0;
  //  for (i = 0; i <= d->legend.nointer; ++i)
  //    {
  //      if (d->legend.val[i] < 0.0)
  //        {
  //          d->legend.col[i].r = 0;
  //          d->legend.col[i].g = (int) 255 * (d->legend.val[i] - middle) / (low - middle);
  //          d->legend.col[i].b = 255;
  //        }
  //      else
  //        {
  //          d->legend.col[i].r = 255;
  //          d->legend.col[i].g = (int) 255 * (d->legend.val[i] - middle) / (high - middle);
  //          d->legend.col[i].b = 0;
  //        }
  //    }
}

void makeRedBlueColors(ctrl_t *c, density_t *d)
{
  int i;
  int border; /* Point at which sign of the interval changes */
  //  double low, high; /* Lowest and highest values */
  //  double middle; /* The value of the "border point" (from blue to red); usually 0, but it could be low or high */

  if (d->legend.nointer == 0)
    {
      if (c->verbose)
        printf("WARNING: makeSaturatedRedBlueColors call ignored: Undefined color intervals\n");
      return;
    }

  if (d->legend.nointer == 0)
    {
      if (c->verbose)
        printf("WARNING: makeSaturatedRedBlueColors call ignored: Undefined color intervals\n");
      return;
    }
  /* Ordered intervals assumed, of course */
  /* Find the point where sign changes */
  border = 0;
  while (d->legend.val[border] < ZEROTOL && border < d->legend.nointer)
    ++border;
  for (i = 0; i <= d->legend.nointer; ++i)
    {
      if (i < border)
        {
          d->legend.col[i].r = 0;
          d->legend.col[i].g = 0;
          d->legend.col[i].b = (int) (255 * ((double) (border - i) / (double) border));
        }
      else
        {
          d->legend.col[i].r = (int) (255 * ((double) (i - border) / (double) (d->legend.nointer - border)));
          d->legend.col[i].g = 0;
          d->legend.col[i].b = 0;
        }

    }
}

void makeRainbowColors(ctrl_t *c, density_t *d)
{
  double h; /* The parameters of the HSV color model */
  double s = 1.0;
  double v = 1.0;
  double chroma;
  double hp;
  double step;
  double half, fraction;
  int integer;
  int i;

  if (d->legend.nointer == 0)
    {
      if (c->verbose)
        printf("WARNING: makeSaturatedRedBlueColors call ignored: Undefined color intervals\n");
      return;
    }
  chroma = v * s;
  step = 360.0 / d->legend.nointer;
  for (i = 0; i <= d->legend.nointer; ++i)
    {
      hp = (d->legend.nointer - i) * step / 60.0; /* d->legend.nointer - i assigns hot colors to most positive values */
      half = hp / 2.0;
      integer = (int) half; /* Integer part of the number */
      /* In the following, half - entero is the fractional part of half */
      /* The fractional part of hp / 2.0 times 2.0 gives "hp (real) modulus 2" */
      /* The rest produces a saw with amplitude = chroma and period 2 (equivalent to 120 degrees) */
      fraction = chroma * (1.0 - fabs(2.0 * (half - integer) - 1.0));

      /* Finally, transform to RGB */
      fraction += v - chroma;
      if (hp >= 0.0 && hp < 1.0)
        {
          d->legend.col[i].r = (int) 255 * v;
          d->legend.col[i].g = (int) 255 * fraction;
          d->legend.col[i].b = (int) 255 * (v - chroma);
        }
      else if (hp >= 1.0 && hp < 2.0)
        {
          d->legend.col[i].r = (int) 255 * fraction;
          d->legend.col[i].g = (int) 255 * v;
          d->legend.col[i].b = (int) 255 * (v - chroma);
        }
      else if (hp >= 2.0 && hp < 3.0)
        {
          d->legend.col[i].r = (int) 255 * (v - chroma);
          d->legend.col[i].g = (int) 255 * v;
          d->legend.col[i].b = (int) 255 * fraction;
        }
      else if (hp >= 3.0 && hp < 4.0)
        {
          d->legend.col[i].r = (int) 255 * (v - chroma);
          d->legend.col[i].g = (int) 255 * fraction;
          d->legend.col[i].b = (int) 255 * v;
        }
      else if (hp >= 4.0 && hp < 5.0)
        {
          d->legend.col[i].r = (int) 255 * fraction;
          d->legend.col[i].g = (int) 255 * (v - chroma);
          d->legend.col[i].b = (int) 255 * v;
        }
      else if (hp >= 5.0 && hp <= 6.0)
        {
          d->legend.col[i].r = (int) 255 * v;
          d->legend.col[i].g = (int) 255 * (v - chroma);
          d->legend.col[i].b = (int) 255 * fraction;
        }
      //      /* Use hot colors for most positive values */
      //      d->legend.col[i].r = 255 - d->legend.col[i].r;
      //      d->legend.col[i].g = 255 - d->legend.col[i].g;
      //      d->legend.col[i].b = 255 - d->legend.col[i].b;
      //      printf("%d %f %f %f | %d %d %d\n", i, hp, fraction, fraction / 255, d->legend.col[i].r, d->legend.col[i].g, d->legend.col[i].b);
    }
  //exit(-1);
}

void getColorFromScale(density_t *d, double val, rgbcolor_t *c)
{
  int min, mid, max; /* Extremes to play "guess a number" with the intervals */
  int smin, sval, smax; /* Signs of interval extremes and user value */
  double stepdiff; /* Interval value */
  double mydiff; /* Distance between our point and the lower end of the interval */

  /* Is this color within range? */
  //  printf("%e | %e %e | %e %e\n", val, d->legend.min, d->legend.max, -d->nmin, d->pmin);
  if (val < d->legend.min || val > d->legend.max || (val > -d->nmin && val < d->pmin))
    {
      c->good = NO;
      c->r = 0;
      c->g = 0;
      c->b = 0;
      return;
    }
  /* First, find the relevant interval */
  min = 0;
  max = d->legend.nointer;
  while (max - min > 1)
    {
      mid = (max + min) / 2;
      //      printf("jui %d %d %d | %e %e %e\n", min, max, mid, d->legend.val[min], d->legend.val[max], val);
      if (val > d->legend.val[mid])
        min = mid;
      else
        max = mid;
    }
  //  printf("---------------jui %d %d %d | %e %e %e\n", min, max, mid, d->legend.val[min], d->legend.val[max], val);
  switch (FLAG_UP(d->legend.method, CQ_LEGEND_INTERPOL_MASK))
    {
  case CQ_LEGEND_INTERPOL_CONSTANT:
    memcpy(c, &d->legend.col[max], sizeof(rgbcolor_t));
    break;
  case CQ_LEGEND_INTERPOL_LINEAR:
    mydiff = (val - d->legend.val[min]);
    stepdiff = (d->legend.val[max] - d->legend.val[min]);
    c->r = (int) (d->legend.col[min].r + mydiff * ((double) (d->legend.col[max].r - d->legend.col[min].r)) / stepdiff);
    c->g = (int) (d->legend.col[min].g + mydiff * ((double) (d->legend.col[max].g - d->legend.col[min].g)) / stepdiff);
    c->b = (int) (d->legend.col[min].b + mydiff * ((double) (d->legend.col[max].b - d->legend.col[min].b)) / stepdiff);
    //    printf("%d %d %d\n", c->r, c->g, c->b);
    break;
  case CQ_LEGEND_INTERPOL_SIGNABSLOG:
    smin = SIGN(d->legend.val[min]);
    smax = SIGN(d->legend.val[max]);
    sval = SIGN(val);
    if (smin != smax || smin != sval)
      {
        if (ctrl.verbose && d->legend.warn)
          {
            d->legend.warn = FALSE;
            fprintf(stderr, "WARNING: You asked for log color interpolation within a color interval that crosses zero.\n"
                "         Expect trouble: The color scale is likely to be wrong\n");
          }
      }
    mydiff = (log10(fabs(val)) - log10(fabs(d->legend.val[min])));
    stepdiff = (log10(fabs(d->legend.val[max])) - log10(fabs(d->legend.val[min])));
    /* This is in fact a linear interpolation of the exponents */
    c->r = (int) (d->legend.col[min].r + mydiff * ((double) (d->legend.col[max].r - d->legend.col[min].r)) / stepdiff);
    c->g = (int) (d->legend.col[min].g + mydiff * ((double) (d->legend.col[max].g - d->legend.col[min].g)) / stepdiff);
    c->b = (int) (d->legend.col[min].b + mydiff * ((double) (d->legend.col[max].b - d->legend.col[min].b)) / stepdiff);
    //    printf("%d %d %d\n", c->r, c->g, c->b);
    break;
  default:
    if (ctrl.verbose)
      printf("WARNING: Using default method for interpolating color\n");
    SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_DFLT_METHOD_INTERPOL, CQ_LEGEND_INTERPOL_MASK);
    SET_MULTIBIT_FLAG(d->legend.method, CQ_LEGEND_DFLT_METHOD_INTERPOL, CQ_LEGEND_INTERPOL_MASK);
    getColorFromScale(d, val, c); /* Try again */
    return;
    }
  c->good = YES;
}
