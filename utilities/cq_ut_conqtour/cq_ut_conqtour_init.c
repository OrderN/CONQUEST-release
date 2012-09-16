/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_init.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_init.c *
 *
 * NAME
 *  cq_ut_conqtour_init.c
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
#include <string.h>
#include <ctype.h>
#include <stdarg.h>

static double atominfo[MAXATOMS][4] =
  {
    { 0.75, 0.4, 0.4, 0.4 },
    { 1.00, 0.0, 0.3, 0.0 },
    { 1.00, 0.0, 0.0, 0.3 },
    { 1.00, 0.3, 0.0, 0.0 },
    { 1.00, 0.3, 0.0, 0.3 } };

void init(int argc, char **argv, index_t *idx, density_t *d)
{
  int i, j, k;
  double x, y, z;

  density.data = NULL;
  setupControl(argc, argv, &ctrl);
  resetLimits(&limits);
  processCommandLine(argc, argv);

  /* Glut initialization */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(600, 450);
  // TODO Remove initial position
#ifdef DEBUG
  glutInitWindowPosition(800, 0);
#endif
  glutCreateWindow("Conqtour");
  checkWindowVisibility(ctrl.has); /* Hide the window if no density is available and will check periodically until one is */

  /* Set up callbacks */
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(specialKeys);
  glutMouseFunc(mouse);
  glutMotionFunc(mouseMotion);

  /* Prepare view*/
  iluminate();
  resetView(&ctrl);
  glClearColor(1.0, 1.0, 1.0, 0.0); /* White background */
  glEnable(GL_DEPTH_TEST);

  /* Prepare visual objects */
  idx->dindex = NULL;
  if (ctrl.has & CQ_CTRL_HAS_DENSITY)
    {
      setBox(&density, &box);
      prepareBoundingBox(&boundingbox);
      prepareLegend(&ctrl, &density);
      resetIndex(NULL, NULL, idx, d, &box);
    }
  prepareAxes(&axes);
  if (ctrl.has & CQ_CTRL_HAS_COORDS)
    {
      molecule = glGenLists(1);
      prepareMol(molecule, 1.0);
    }

  // TODO These should be in ctrl_t
  modifiers = 0;
  rendermode = GL_RENDER;
  glSelectBuffer(PICKBUFSIZE, pickbuffer);

  /* Find available information about a Conquest run */
  resetConquestInfo(&cqinfo);
  //strcpy(cqinfo.inputfl,"uuuuuuuuuuuuuu76235472");
  //strcpy(cqinfo.outputfl,"76235472");
  //strcpy(cqinfo.blocksfl,"76235472");
  //strcpy(cqinfo.densityfl,"thmkiycbSFs");
  //cqinfo.numspecies=7349;
  //limits.gridshft[0]=100;

  //limits.strides[0] = 3;
  //limits.strides[1] = 2;
  //limits.strides[2] = 1;
  //limits.gridshft[0]=216;
  //limits.gridshft[1]=250;
  //limits.griddims[0]=216;
  //limits.griddims[1]=250;
  /*
   limits.griddims[2]=250;
   limits.griddims[0]=108;
   limits.griddims[1]=125;
   limits.griddims[2]=125;
   */
  //limits.gridshft[1]=100;
  //limits.gridshft[2]=100;
  findConquestAvailability(ctrl.verbose);
  //printf("--------------------------------------\n");
  //resetConquestInfo(&cqinfo);
  //findConquestAvailability(ctrl.verbose);
}

void processCommandLine(int argc, char **argv)
{
  int i;
  int processOptions = YES;
  char *option;

  initializeQueue(&scripts);
  for (i = 1; i < argc; ++i)
    {
      checkArgumentLength(argv[i]);
      if (argv[i][0] == '-' && processOptions == YES) /* Option */
        {
          if (argv[i][1] == '-')
            {
              if (argv[i][2] == '\0')
                {
                  processOptions = NO;
                  continue;
                }
            }
          option = &argv[i][1];
          if (!strcmp(option, "v") || !strcmp(option, "-verbose"))
            {
              ctrl.verbose = VERBOSE;
            }
          else if (!strcmp(option, "d") || !strcmp(option, "-density"))
            {
              ++i;
              if (i == argc)
                showUsageAndExit(); /* Insufficient arguments */
              else
                {
                  checkArgumentLength(argv[i]);
                  if (!readXplorDensity(argv[i], &density))
                    ctrl.has |= CQ_CTRL_HAS_DENSITY;
                  else
                    ctrl.has &= ~CQ_CTRL_HAS_DENSITY;
                }
            }
          else if (!strcmp(option, "c") || !strcmp(option, "-coordinates"))
            {
              ++i;
              if (i == argc)
                showUsageAndExit(); /* Insufficient arguments */
              else
                {
                  // TODO Check whether it is xyz or some other format
                  checkArgumentLength(argv[i]);
                  if (!readXyzCoordinates(argv[i], &coords))
                    {
                      ctrl.has |= CQ_CTRL_HAS_COORDS;
                      ctrl.has |= CQ_CTRL_HAS_COORDS_XYZ;
                      ctrl.display |= CQ_CTRL_SHOW_MOLECULE;
                    }
                  else
                    {
                      ctrl.has &= ~CQ_CTRL_HAS_COORDS;
                      ctrl.has &= ~CQ_CTRL_HAS_COORDS_XYZ;
                      ctrl.display &= ~CQ_CTRL_SHOW_MOLECULE;
                    }
                }
            }
          else
            showUsageAndExit();
        }
      else if (argv[i][0] != '-' || processOptions == NO) /* Python script */
        {
          //printf("Python file: %s\n", argv[i]);
          checkArgumentLength(argv[i]);
          pushNameToQueue(&scripts, argv[i]);
        }
    }
}

void checkArgumentLength(char *arg)
{
  if (strlen(arg) > MAXLIN - 1)
    {
      fprintf(stderr, "Argument is too long\n");
      exit(-1);
    }
}

void showUsageAndExit(void)
{
  fprintf(stderr, "Usage: %s [-v|--verbose] "
    "[-d|--density <xplor file>] "
    "[-c|--coordinates <xyz file>] "
    "[python scripts...]\n", PROGRAM_NAME);
  exit(-1);
}

void initializeQueue(namequeue_t *q)
{
  q->head = NULL;
  q->tail = NULL;
}

void pushNameToQueue(namequeue_t *q, char *name)
{
  namequeuenode_t *tmp;

  tmp = q->tail;
  q->tail = malloc(sizeof(namequeuenode_t));
  if (q->tail == NULL)
    {
      fprintf(stderr, "Error allocating node for name queue\n");
      exit(-1);
    }
  q->tail->next = NULL;
  if (q->head == NULL)
    q->head = q->tail; /* First node */
  if (tmp != NULL)
    tmp->next = q->tail; /* Previous tail points to new tail */
  strcpy(q->tail->name, name); /* Copy name (the input could destroyed any time) */
}

char *popNameFromQueue(namequeue_t *q)
{
  namequeuenode_t *tmp;
  char *name;

  name = malloc(MAXLIN * sizeof(char));
  if (name == NULL)
    {
      fprintf(stderr, "Error allocating name\n");
      exit(-1);
    }
  if (q->head == NULL)
    return NULL; /* Empty queue */
  strcpy(name, q->head->name);
  tmp = q->head;
  q->head = q->head->next; /* Point to next node */
  if (q->head == NULL) /* Empty queue */
    q->tail = NULL;
  free(tmp); /* Free former head node */

  return name; /* The caller must free !!! */
}

void prepareBoundingBox(GLuint *bb)
{
  int i;

  /* Create a representation of the bounding box */
  *bb = glGenLists(1);
  glNewList(*bb, GL_COMPILE);
  glBegin(GL_LINES);
  glColor3f(0.0, 0.0, 0.0);
  for (i = 0; i < 12; ++i)
    {
      glVertex3f(box.refpoint[i][0], box.refpoint[i][1], box.refpoint[i][2]);
      glVertex3f(box.refpoint[i][0] + box.edge[i][0], box.refpoint[i][1] + box.edge[i][1], box.refpoint[i][2] + box.edge[i][2]);
    }
  glEnd();
  glEndList();
}

void prepareAxes(GLuint *a)
{
  /* Create a simple representation of the axes */
  *a = glGenLists(1);
  glNewList(*a, GL_COMPILE);
  glBegin(GL_LINES);
  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(ctrl.center[0], ctrl.center[1], ctrl.center[2]);
  glVertex3f(density.lvm[0], ctrl.center[1], ctrl.center[2]);
  glColor3f(0.0, 1.0, 0.0);
  glVertex3f(ctrl.center[0], ctrl.center[1], ctrl.center[2]);
  glVertex3f(ctrl.center[0], density.lvm[1], ctrl.center[2]);
  glColor3f(0.0, 0.0, 1.0);
  glVertex3f(ctrl.center[0], ctrl.center[1], ctrl.center[2]);
  glVertex3f(ctrl.center[0], ctrl.center[1], density.lvm[2]);
  glEnd();
  glEndList();
}

void setBox(density_t *d, box_t *b)
{
  /* Prepare box in a way convenient to intersect it with a plane */
  b->refpoint[0][0] = 0.0;
  b->refpoint[0][1] = 0.0;
  b->refpoint[0][2] = 0.0;
  b->edge[0][0] = d->lvm[0];
  b->edge[0][1] = 0.0;
  b->edge[0][2] = 0.0;
  b->refpoint[1][0] = d->lvm[0];
  b->refpoint[1][1] = 0.0;
  b->refpoint[1][2] = 0.0;
  b->edge[1][0] = 0.0;
  b->edge[1][1] = d->lvm[1];
  b->edge[1][2] = 0.0;
  b->refpoint[2][0] = d->lvm[0];
  b->refpoint[2][1] = d->lvm[1];
  b->refpoint[2][2] = 0.0;
  b->edge[2][0] = 0.0;
  b->edge[2][1] = 0.0;
  b->edge[2][2] = d->lvm[2];
  b->refpoint[3][0] = d->lvm[0];
  b->refpoint[3][1] = 0.0;
  b->refpoint[3][2] = 0.0;
  b->edge[3][0] = 0.0;
  b->edge[3][1] = 0.0;
  b->edge[3][2] = d->lvm[2];
  b->refpoint[4][0] = 0.0;
  b->refpoint[4][1] = 0.0;
  b->refpoint[4][2] = 0.0;
  b->edge[4][0] = 0.0;
  b->edge[4][1] = 0.0;
  b->edge[4][2] = d->lvm[2];
  b->refpoint[5][0] = 0.0;
  b->refpoint[5][1] = 0.0;
  b->refpoint[5][2] = d->lvm[2];
  b->edge[5][0] = d->lvm[0];
  b->edge[5][1] = 0.0;
  b->edge[5][2] = 0.0;
  b->refpoint[6][0] = d->lvm[0];
  b->refpoint[6][1] = 0.0;
  b->refpoint[6][2] = d->lvm[2];
  b->edge[6][0] = 0.0;
  b->edge[6][1] = d->lvm[1];
  b->edge[6][2] = 0.0;
  b->refpoint[7][0] = 0.0;
  b->refpoint[7][1] = 0.0;
  b->refpoint[7][2] = d->lvm[2];
  b->edge[7][0] = 0.0;
  b->edge[7][1] = d->lvm[1];
  b->edge[7][2] = 0.0;
  b->refpoint[8][0] = 0.0;
  b->refpoint[8][1] = 0.0;
  b->refpoint[8][2] = 0.0;
  b->edge[8][0] = 0.0;
  b->edge[8][1] = d->lvm[1];
  b->edge[8][2] = 0.0;
  b->refpoint[9][0] = 0.0;
  b->refpoint[9][1] = d->lvm[1];
  b->refpoint[9][2] = 0.0;
  b->edge[9][0] = 0.0;
  b->edge[9][1] = 0.0;
  b->edge[9][2] = d->lvm[2];
  b->refpoint[10][0] = 0.0;
  b->refpoint[10][1] = d->lvm[1];
  b->refpoint[10][2] = d->lvm[2];
  b->edge[10][0] = d->lvm[0];
  b->edge[10][1] = 0.0;
  b->edge[10][2] = 0.0;
  b->refpoint[11][0] = 0.0;
  b->refpoint[11][1] = d->lvm[1];
  b->refpoint[11][2] = 0.0;
  b->edge[11][0] = d->lvm[0];
  b->edge[11][1] = 0.0;
  b->edge[11][2] = 0.0;

  b->corner[0][0] = 0.0;
  b->corner[0][1] = 0.0;
  b->corner[0][2] = 0.0;
  b->corner[1][0] = d->lvm[0];
  b->corner[1][1] = 0.0;
  b->corner[1][2] = 0.0;
  b->corner[2][0] = 0.0;
  b->corner[2][1] = d->lvm[1];
  b->corner[2][2] = 0.0;
  b->corner[3][0] = 0.0;
  b->corner[3][1] = 0.0;
  b->corner[3][2] = d->lvm[2];
  b->corner[4][0] = d->lvm[0];
  b->corner[4][1] = d->lvm[1];
  b->corner[4][2] = 0.0;
  b->corner[5][0] = d->lvm[0];
  b->corner[5][1] = 0.0;
  b->corner[5][2] = d->lvm[2];
  b->corner[6][0] = 0.0;
  b->corner[6][1] = d->lvm[1];
  b->corner[6][2] = d->lvm[2];
  b->corner[7][0] = d->lvm[0];
  b->corner[7][1] = d->lvm[1];
  b->corner[7][2] = d->lvm[2];

}

void prepareMol(GLuint molecule, double scale)
{
  int i, j;
  int index;
  double red, green, blue;

  if (scale > 1.0)
    scale = 1.0;
  if (scale < 0.3)
    scale = 0.3;

  //printf("%d\n", coords.numatoms);
  glNewList(molecule, GL_COMPILE);
  //        glPushMatrix();
  for (i = 0; i < coords.numatoms; ++i)
    {
      //printf("Atom\n");
      if (!strcmp(coords.atname[i], "H"))
        index = 0;
      else if (!strcmp(coords.atname[i], "C"))
        index = 1;
      else if (!strcmp(coords.atname[i], "N"))
        index = 2;
      else if (!strcmp(coords.atname[i], "O"))
        index = 3;
      else
        index = 4;

      //if(glutGetModifiers() & GLUT_ACTIVE_SHIFT)  index=4;
      //if(modifiers & GLUT_ACTIVE_SHIFT)  index=4;

      //if(i > 3149) scale = 0.3;

      red = atominfo[index][1];
      green = atominfo[index][2];
      blue = atominfo[index][3];
      if (ctrl.selectedatoms)
        {
          for (j = 0; j < ctrl.selectedatoms; ++j)
            if (ctrl.selection[j] == i)
              {
                red = 0.7;
                green = 0.7;
                blue = 0.0;
              }
        }
      glPushMatrix();
      glColor3f(red, green, blue);
      glTranslatef(coords.xyz[i][0], coords.xyz[i][1], coords.xyz[i][2]);
      glLoadName(molecule + i); /* Name atoms for picking*/
      glutSolidSphere(atominfo[index][0] * scale, 20, 20);
      glPopMatrix();
      //           glTranslatef(-coords.xyz[i][0], -coords.xyz[i][1], -coords.xyz[i][2]);
    }
  //        glPopMatrix();
  glEndList();
}

void prepareLegend(ctrl_t *c, density_t *d)
{
  d->legend.x = c->legend.x;
  d->legend.y = c->legend.y;
  d->legend.w = c->legend.w;
  d->legend.h = c->legend.h;
  d->legend.bt = c->legend.bt;
  d->legend.tl = c->legend.tl;
  d->legend.tickorigin = c->legend.tickorigin;
  d->legend.tickinterval = c->legend.tickinterval;
  d->legend.noticks = c->legend.noticks;
  d->legend.method = c->legend.method;
  d->legend.max = (d->pmin <= d->pmax) ? d->pmax : -d->nmin; /* If the condition is false, all points are negative */
  d->legend.min = (d->nmin <= d->nmax) ? -d->nmax : d->pmin; /* If the condition is false, all points are positive */
  strcpy(d->legend.format, c->legend.format);
  d->legend.warn = c->legend.warn;
  /* Correct origin if out of range */
  if (d->legend.tickorigin < d->legend.min || d->legend.tickorigin > d->legend.max)
    {
      if (d->legend.min > CQ_LEGEND_DFLT_TICKORIGIN) /* If the default is also out of range (under), use minimum datum */
        {
          d->legend.tickorigin = d->legend.min;
        }
      else if (d->legend.max < CQ_LEGEND_DFLT_TICKORIGIN) /* If it is out of range (over), use maximum */
        {
          d->legend.tickorigin = d->legend.max;
        }
      else /* Otherwise, use the default */
        {
          d->legend.tickorigin = CQ_LEGEND_DFLT_TICKORIGIN;
        }
    }
  /* Make color intervals */
  switch (FLAG_UP(d->legend.method, CQ_LEGEND_INTERVAL_MASK))
    {
  case CQ_LEGEND_INTERVAL_LINEAR:
    makeLinearColorIntervals(c, d);
    break;
  case CQ_LEGEND_INTERVAL_SIGNABSLOG:
    makeSignedAbsoluteLogColorIntervals(c, d);
    break;
  default:
    if (c->verbose)
      printf("WARNING: Using default method to create color intervals\n");
    SET_MULTIBIT_FLAG(c->legend.method, CQ_LEGEND_DFLT_METHOD_INTERVAL, CQ_LEGEND_INTERVAL_MASK);
    prepareLegend(c, d); /* Try again */
    return; /* What comes next should have been done already */
    }
  /* Make colors */
  switch (FLAG_UP(d->legend.method, CQ_LEGEND_COLOR_MASK))
    {
  case CQ_LEGEND_COLOR_SATREDBLUE:
    makeSaturatedRedBlueColors(c, d);
    break;
  case CQ_LEGEND_COLOR_REDBLUE:
    makeRedBlueColors(c, d);
    break;
  case CQ_LEGEND_COLOR_RAINBOW:
    makeRainbowColors(c, d);
    break;
  default:
    if (c->verbose)
      printf("WARNING: Using default method to create color scale\n");
    SET_MULTIBIT_FLAG(c->legend.method, CQ_LEGEND_DFLT_METHOD_COLOR, CQ_LEGEND_COLOR_MASK);
    prepareLegend(c, d); /* Try again */
    return; /* What comes next should have been done already */
    }
  makeLegend(&legend, c, d);
  //  printf("Done\n");
}

void resetDensity(density_t *d, int *old)
{
  int i, j;

  /* Deallocate the density if it was already allocated (e.g. loading new density) */
  if (d->data != NULL)
    {
      for (i = 0; i < old[0]; ++i)
        {
          for (j = 0; j < old[1]; ++j)
            {
              //printf("%d %d\n", i, j);
              if (d->data[i][j] != NULL)
                free(d->data[i][j]);
            }
          if (d->data[i] != NULL)
            free(d->data[i]);
        }
      free(d->data);
    }
  //printf("==%d | %d %d %d\n", d->data, d->pt[0], d->pt[1], d->pt[2]);

  /* Allocate space for the density */
  /* NOTE: It would be more efficient to have indices in the order ZYX */
  /*       but we can afford having a more intuitive order */
  d->data = (double ***) malloc((size_t) d->pt[0] * sizeof(double **));
  for (i = 0; i < d->pt[0]; ++i)
    {
      d->data[i] = (double **) malloc((size_t) d->pt[1] * sizeof(double *));
      for (j = 0; j < d->pt[1]; ++j)
        {
          d->data[i][j] = (double *) malloc((size_t) d->pt[2] * sizeof(double));
        }
    }

}

/* NULL normal means z-axis; NULL inplane means center of cell */
void resetIndex(GLdouble *normal, GLdouble *inplane, index_t *idx, density_t *d, box_t *b)
{
  GLdouble aux1[3], aux2[3], aux3[3];
  int i, j;
  double zmin, zmax, zcoord;

  //   pthread_mutex_lock(&ctrl.mutex);

  /* Deallocate previous index, if it existed */
  if (idx->dindex != NULL)
    {
      for (i = 0; i < idx->cuts; ++i)
        {
          if (idx->dindex[i] != NULL)
            {
#ifdef DEBUG
              printf("Deleting %d %d\n", idx->dindex[i]->slice, idx->cuts);
#endif
              if (glIsList(idx->dindex[i]->slice) == GL_TRUE)
                glDeleteLists(idx->dindex[i]->slice, 1);
              if (idx->dindex[i]->data != NULL)
                {
                  for (j = 0; j < idx->dindex[i]->pt[0]; ++j)
                    {
                      //printf("oo %d %d\n", i, j);
                      free(idx->dindex[i]->data[j]);
                    }
                  free(idx->dindex[i]->data);
                  free(idx->dindex[i]->lim[LEFT]);
                  free(idx->dindex[i]->lim[RIGHT]);
                  free(idx->dindex[i]->xbord[LEFT]);
                  free(idx->dindex[i]->xbord[RIGHT]);
                }
              free(idx->dindex[i]);
            }
        }
      free(idx->dindex);
    }

  /* Prepare normal vector and inplane point */
  if (normal == NULL)
    {
      aux1[0] = 0.0;
      aux1[1] = 0.0;
      aux1[2] = 1.0;
    }
  else
    {
      for (i = 0; i < 3; ++i)
        aux1[i] = normal[i];
      normalize(aux1); /* Just in case it is not */
    }
  if (inplane == NULL)
    {
      for (i = 0; i < 3; ++i)
        aux2[i] = d->lvm[i] / 2.0;
    }
  else
    {
      if (inplane[0] > 0.0 && inplane[0] < d->lvm[0] && inplane[1] > 0.0 && inplane[1] < d->lvm[1] && inplane[2] > 0.0
          && inplane[2] < d->lvm[2])
        {
          if (ctrl.verbose)
            {
              fprintf(stderr, "WARNING: Requested in-plane point for slice is out of bounds\n");
              fprintf(stderr, "         Using center of the box\n");
            }
          for (i = 0; i < 3; ++i)
            aux2[i] = inplane[i];
        }
      else
        {
          for (i = 0; i < 3; ++i)
            aux2[i] = d->lvm[i] / 2.0;
        }
    }

  /* Calculate number of cuts and initial current cut */
  zmin = zmax = 0.0;
  for (i = 0; i < 8; ++i)
    {
      for (j = 0; j < 3; ++j)
        aux3[j] = b->corner[i][j] - aux2[j];
      zcoord = dotProduct(aux1, aux3);
#ifdef DEBUG
      printf("z %d %f | %f %f %f\n", i, zcoord, b->corner[i][0], b->corner[i][1], b->corner[i][2]);
#endif
      if (zcoord > zmax)
        zmax = zcoord;
      else if (zcoord < zmin)
        zmin = zcoord;
    }

  // IMPORTANT NOTE : The question of the increments for arbitrary normal must be solved
  //                  For the moment, use the values deducted from input file
  idx->dz = d->del[2];

  idx->cuts = (int) floor(fabs(zmin / d->del[2]));
  idx->curr = idx->cuts; /* This is the value of the initial slice */
  idx->cuts += 1; /* Count the initial (central) slice */
  idx->cuts += (int) floor(fabs(zmax / d->del[2]));
#ifdef DEBUG
  printf("Extremes = %f %f %d %d | %d %d\n", zmin, zmax, idx->cuts, idx->curr,
      (int)floor(fabs(zmin/d->del[2])), (int)floor(fabs(zmax/d->del[2])));
#endif

  /* Allocate enough space for the index */
  idx->dindex = (densitycut_t **) malloc((size_t) idx->cuts * sizeof(densitycut_t *));
  for (i = 0; i < idx->cuts; ++i)
    idx->dindex[i] = NULL;
  //   idx->dindex[idx->curr] = (densitycut_t *)malloc((size_t)sizeof(densitycut_t));

  /* Establish origin for the whole index */
  for (i = 0; i < 3; ++i)
    {
      idx->surfnormal[i] = aux1[i];
      idx->inplane[i] = aux2[i];
      idx->orig = idx->curr;
    }
  idx->fp = NULL; /* This signals that slices will drawn visually, rather than to a file (e.g. pymol script) */

  // TODO Move out of this function, and perhaps split in Calculate and Draw
  /* Now, create the initial slice */
  calculateDrawSlice(idx, b, d);

  //   pthread_mutex_unlock(&ctrl.mutex);
  /*
   dindex[currentcut]->surfnormal[0]=0.0;
   dindex[currentcut]->surfnormal[1]=0.0;
   dindex[currentcut]->surfnormal[2]=1.0;
   dindex[currentcut]->inplane[0]=d->lvm[0]/2.0;
   dindex[currentcut]->inplane[1]=d->lvm[1]/2.0;
   dindex[currentcut]->inplane[2]=d->lvm[2]/2.0;
   */
}

void resetConquestInfo(cqruninfo_t *info)
{
  int i;

  if (info->typename != NULL) /* Free type names, if they exist */
    {
      for (i = 0; i < info->numspecies; ++i)
        free(info->typename[i]);
      free(info);
    }
  memset(info, 0, sizeof(cqruninfo_t));
}

void resetLimits(limits_t *l)
{
  memset(l, 0, sizeof(limits_t));
  l->strides[0] = l->strides[1] = l->strides[2] = 1;
}

/* Find out as much information as possible about a Conquest run */
// TODO Modify as necessary to accommodate the pyCQ user's interface, most likely by breaking into smaller functions
// TODO Describe side effects in details, e.g. which flags will be always reset from available info
int findConquestAvailability(int verbose)
{
  int i, j;
  FILE *fci; /* Conquest input file */
  FILE *fco; /* Conquest output file */
  char value[MAXLIN];
  char tmpfile[MAXLIN]; /* Used to build CQ filenames */
  char units[MAXLIN]; /* Units read from the output as a string */
  int problemcnt; /* Count of files that cannot be accessed */
  double conv; /* Conversion factor for distance units */

  if (verbose)
    {
      printf("-------------------------------------------------------------------\n");
      printf("-------------------------------------------------------------------\n");
      printf("----                 KNOWN CONQUEST PARAMETERS                 ----\n");
      printf("----   Note: this is NOT what you choose in the command line   ----\n");
      printf("----         and it is only used if you load densities and/or  ----\n");
      printf("----         coordinates from the pyCQ prompt                  ----\n");
      printf("-------------------------------------------------------------------\n");
      printf("-------------------------------------------------------------------\n");
    }
  /* Use defaults only if the values are not available, i.e. given by the user */

  /* Decide name of Conquest input file */
  if (!strlen(cqinfo.inputfl))
    {
      strcpy(cqinfo.inputfl, CQ_DFLT_INPUT); /* Start with the default input file */
      if (verbose)
        printf("Info: Setting input file to default: %s\n", CQ_DFLT_INPUT);
    }
  else
    {
      if (verbose)
        printf("Info: Keeping predefined input file value: %s\n", cqinfo.inputfl);
    }
  /* Have Conquest input file name !! Try to open it */
  fci = fopen(cqinfo.inputfl, "r");
  SET_FLAG(cqinfo.have, CQ_HAVE_INPUT, (fci != NULL));
  if (fci == NULL)
    {
      if (verbose)
        printf("WARNING: CQ input file not found for reading: %s\n", cqinfo.inputfl);
      if (!strlen(cqinfo.outputfl))
        {
          strcpy(cqinfo.outputfl, CQ_DFLT_OUTPUT);
          if (verbose)
            printf("Info: Setting output file to default: %s\n", CQ_DFLT_OUTPUT);
        }
      else
        {
          if (verbose)
            printf("Info: Keeping predefined output file value: %s\n", cqinfo.outputfl);
        }
    }
  if (FLAG_UP(cqinfo.have, CQ_HAVE_INPUT))
    {
      if (!strlen(cqinfo.outputfl)) /* If not user-specified, search in input */
        {
          if (!getFlag(fci, CQ_FLAG_OUTPUTFILE, cqinfo.outputfl)) /* Input specifies output */
            {
              strcpy(cqinfo.outputfl, CQ_DFLT_OUTPUT);
              if (verbose)
                printf("Info: Setting output file to default: %s\n", CQ_DFLT_OUTPUT);
            }
          else if (verbose)
            printf("Info: Using %s to set output file name: %s\n", cqinfo.inputfl, cqinfo.outputfl);
        }
      else
        {
          if (verbose)
            printf("Info: Keeping predefined output file value: %s\n", cqinfo.outputfl);
        }
    }
  /* Have Conquest output file name !! Try to open it */
  fco = fopen(cqinfo.outputfl, "r");
  SET_FLAG(cqinfo.have, CQ_HAVE_OUTPUT, (fco != NULL));
  if (verbose && fco == NULL)
    printf("WARNING: CQ output file not found for reading: %s\n", cqinfo.outputfl);
  if (strlen(cqinfo.coordinatesfl)) /* If the user-specified an input file, keep it */
    {
      if (verbose)
        printf("Info: Keeping predefined coordinates file value: %s\n", cqinfo.coordinatesfl);
    }
  else if (FLAG_UP(cqinfo.have, CQ_HAVE_INPUT))
    {
      if (!getFlag(fci, CQ_FLAG_COORDFILE, cqinfo.coordinatesfl)) /* Input specifies output? */
        {
          strcpy(cqinfo.coordinatesfl, CQ_DFLT_COORDINATES);
          if (verbose)
            printf("Info: Setting coordinates file to default: %s\n", CQ_DFLT_COORDINATES);
        }
      else if (verbose)
        printf("Info: Using %s to set coordinates file name: %s\n", cqinfo.inputfl, cqinfo.coordinatesfl);
    }
  else
    {
      strcpy(cqinfo.coordinatesfl, CQ_DFLT_COORDINATES);
      if (verbose)
        printf("Info: Setting coordinates file to default: %s\n", CQ_DFLT_COORDINATES);
    }
  /* Have coordinates file name !! Check that it can be read */
  /* Note: access returns 0 if the file can be opened with the requested permissions */
  SET_FLAG(cqinfo.have, CQ_HAVE_COORDINATES,!access(cqinfo.coordinatesfl, R_OK));
  if (verbose && !FLAG_UP(cqinfo.have, CQ_HAVE_COORDINATES))
    printf("WARNING: CQ coordinates file not found for reading: %s\n", cqinfo.coordinatesfl);

  /* CQ_HAVE_FRACTIONAL is used to give the user a mechanism to preserve the value of CQ_FLAG_FRACTIONAL */
  if (FLAG_UP(cqinfo.have, CQ_HAVE_FRACTIONAL)) /* The value of cqinfo.read is meaningful and was previously set */
    {
      if (verbose)
        printf("Info: Keeping value of flag: %s = %s\n", CQ_FLAG_FRACTIONAL, REPORT_FLAG(cqinfo.read, CQ_READ_FRACTIONAL));
    }
  else /* Still don't know whether CQ coordinates should be considered fractional (in case we have them) */
    {
      SET_FLAG(cqinfo.have, CQ_HAVE_FRACTIONAL, TRUE); /* The flag will be set in any case after calling this function */
      /* Get information from input, if available */
      if (FLAG_UP(cqinfo.have, CQ_HAVE_INPUT) && getFlag(fci, CQ_FLAG_FRACTIONAL, value))
        {
          SET_FLAG(cqinfo.read, CQ_READ_FRACTIONAL, isTrue(value));
          if (verbose)
            printf("Info: Setting flag: %s = %s\n", CQ_FLAG_FRACTIONAL, REPORT_FLAG(cqinfo.read, CQ_READ_FRACTIONAL));
        }
      else
        {
          SET_FLAG(cqinfo.read, CQ_READ_FRACTIONAL, isTrue(CQ_DFLT_FRACTIONAL));
          if (verbose)
            printf("Info: Setting flag to default: %s = %s\n", CQ_FLAG_FRACTIONAL, REPORT_FLAG(cqinfo.read, CQ_READ_FRACTIONAL));
        }
    }
  /* Try to get the number of species in the run */
  if (cqinfo.numspecies != 0) /* We already know; keep the value */
    {
      if (verbose)
        printf("Info: Keeping value of flag: %s = %d\n", CQ_FLAG_SPECIES, cqinfo.numspecies);
    }
  else
    {
      if (FLAG_UP(cqinfo.have, CQ_HAVE_INPUT) && getFlag(fci, CQ_FLAG_SPECIES, value)) /* Input specifies the number of species */
        {
          /* We don't need to check for errno, etc. in the conversion because getting 0 IS an error */
          cqinfo.numspecies = (int) strtol(value, NULL, 10);
          if (cqinfo.numspecies == 0)
            {
              if (verbose)
                printf("WARNING: Invalid number of species in %s\n", cqinfo.inputfl);
            }
          else if (verbose)
            printf("Info: Number of species = %d\n", cqinfo.numspecies);
        }
      else /* Set the number of species */
        {
          cqinfo.numspecies = 0;
          if (verbose)
            printf("WARNING: The number of species is unknown\n");
        }
    }
  /* Try to read the names of each species */
  if (cqinfo.numspecies > 0 && cqinfo.typename != NULL)
    {
      if (verbose)
        {
          printf("Info: Keeping atom types, as follows:\n");
          for (i = 0; i < cqinfo.numspecies; ++i)
            {
              printf("Info:  Species %3d %s\n", i + 1, cqinfo.typename[i]);
            }
        }
    }
  else if (FLAG_UP(cqinfo.have, CQ_HAVE_INPUT) && cqinfo.numspecies > 0)
    {
      if (findBlockLabel(fci, "ChemicalSpeciesLabel"))
        {
          cqinfo.typename = (char **) malloc((size_t) cqinfo.numspecies * sizeof(char *));
          for (i = 0; i < cqinfo.numspecies; ++i)
            {
              cqinfo.typename[i] = (char *) malloc((size_t) MAXLIN * sizeof(char));
              if (fgets(value, MAXLIN - 1, fci) && sscanf(value, "%*d %*f %s", cqinfo.typename[i]) == 1)
                {
                  if (verbose)
                    printf("Info:  Species %3d = %s\n", i + 1, cqinfo.typename[i]);
                }
              else
                {
                  for (j = 0; j <= i; ++j)
                    free(cqinfo.typename[j]);
                  free(cqinfo.typename);
                  cqinfo.typename = NULL;
                  if (verbose)
                    printf("WARNING: Error reading species. Deleting the whole block\n");
                  break;
                }
            }
        }
      else
        {
          if (verbose)
            printf("WARNING: Atom type names could not be found in %s\n", cqinfo.inputfl);
        }
    }
  else
    {
      if (verbose)
        printf("WARNING: It is not possible to find atom type names from %s\n", cqinfo.inputfl);
    }
  /* If the distance unit is set, we keep it. (This is to preserve the user's previous choice) */
  if (FLAG_UP(cqinfo.units, CQ_UNITS_MASKDIST)) /* The value of cqinfo.units is meaningful and was previously set */
    {
      if (verbose)
        printf("Info: Keeping value of flag: %s = %s\n", CQ_FLAG_DISTUNITS, REPORT_DISTUNITS(cqinfo.units));
    }
  else /* Still don't know what the distance units are */
    {
      /* Get information from input, if available */
      if (FLAG_UP(cqinfo.have, CQ_HAVE_INPUT) && getFlag(fci, CQ_FLAG_DISTUNITS, value))
        {
          if (isEqualStr(value, CQ_TAG_BOHR1) || isEqualStr(value, CQ_TAG_BOHR2))
            {
              SET_FLAG(cqinfo.units, CQ_UNITS_MASKDIST, FALSE); /* Clear the bits for distance units */
              SET_FLAG(cqinfo.units, CQ_UNITS_BOHR, TRUE); /* Set units to bohr */
              if (verbose)
                printf("Info: Setting flag: %s = %s\n", CQ_FLAG_DISTUNITS, REPORT_DISTUNITS(cqinfo.units));
            }
          else if (isEqualStr(value, CQ_TAG_ANGSTROM))
            {
              SET_FLAG(cqinfo.units, CQ_UNITS_MASKDIST, FALSE); /* Clear the bits for distance units */
              SET_FLAG(cqinfo.units, CQ_UNITS_ANGSTROM, TRUE); /* Set units to angstrom */
              if (verbose)
                printf("Info: Setting flag: %s = %s\n", CQ_FLAG_DISTUNITS, REPORT_DISTUNITS(cqinfo.units));
            }
          else
            {
              SET_FLAG(cqinfo.units, CQ_UNITS_MASKDIST, FALSE); /* Clear the bits for distance units */
              SET_FLAG(cqinfo.units, CQ_UNITS_DFLTDIST, TRUE); /* Set units to default */
              if (verbose)
                printf("Info: Setting flag to default: %s = %s\n", CQ_FLAG_DISTUNITS, REPORT_DISTUNITS(cqinfo.units));
            }
        }
      else
        {
          SET_FLAG(cqinfo.units, CQ_UNITS_MASKDIST, FALSE); /* Clear the bits for distance units */
          SET_FLAG(cqinfo.units, CQ_UNITS_DFLTDIST, TRUE); /* Set units to default */
          if (verbose)
            printf("Info: Setting flag to default: %s = %s\n", CQ_FLAG_DISTUNITS, REPORT_DISTUNITS(cqinfo.units));
        }
    }
  /* CQ_HAVE_FROM_BLOCKS is set if the block method (from file vs old) is already known */
  if (FLAG_UP(cqinfo.have, CQ_HAVE_FROM_BLOCKS))
    {
      if (verbose)
        printf("Info: Keeping value of flag: %s = %s\n", CQ_FLAG_FROM_BLOCKS, REPORT_FLAG(cqinfo.read, CQ_READ_FROM_BLOCKS));
    }
  else /* Still don't know whether blocks should be read from file or calculated by the old method */
    {
      SET_FLAG(cqinfo.have, CQ_HAVE_FROM_BLOCKS, TRUE); /* The flag will be set in any case after calling this function */
      /* This information cannot be found in the input (as far as I know) */
      /* Note that the CQ flag Grid.ReadBlocks is not exactly equivalent to this, since that one usually defaults to false */
      /* Here, by default, a file is read; only the user can specify otherwise to force the old method */
      // TODO Write a pyCQ function to allow the user to do this
      SET_FLAG(cqinfo.read, CQ_READ_FROM_BLOCKS, isTrue(CQ_DFLT_FROM_BLOCKS));
      if (verbose)
        printf("Info: Setting flag to default: %s = %s\n", CQ_FLAG_FROM_BLOCKS, REPORT_FLAG(cqinfo.read, CQ_READ_FROM_BLOCKS));
    }
  /* Decide about the partitioning method */
  if (FLAG_UP(cqinfo.have, CQ_HAVE_HILBERT)) /* Partitioning method has already been chosen */
    {
      if (verbose)
        printf("Info: Keeping value of flag: %s = %s\n", CQ_FLAG_PARTITIONER, REPORT_PARTITIONER(cqinfo.read, CQ_READ_HILBERT));
    }
  else
    {
      SET_FLAG(cqinfo.have, CQ_HAVE_HILBERT, TRUE); /* After this, the partitioning method is chosen */
      if (FLAG_UP(cqinfo.have, CQ_HAVE_INPUT))
        {
          if (!getFlag(fci, CQ_FLAG_PARTITIONER, value)) /* Input specifies partitioner method */
            {
              SET_FLAG(cqinfo.read, CQ_READ_HILBERT, isEqualStr(CQ_HILBERT_NAME,CQ_DFLT_PARTITIONER));
              if (verbose)
                printf("Info: Setting flag to default: %s = %s\n", CQ_FLAG_PARTITIONER, CQ_DFLT_PARTITIONER);
              //           strcpy(cqinfo.blocksfl, CQ_DFLT_BLOCKS);
              //           printf("Info: Setting blocks file to default: %s\n", CQ_DFLT_BLOCKS);
            }
          else /* We have partitioning method flag */
            {
              SET_FLAG(cqinfo.read, CQ_READ_HILBERT, isEqualStr(value,CQ_HILBERT_NAME));
              if (verbose)
                printf("Info: Setting flag: %s = %s\n", CQ_FLAG_PARTITIONER, value);
            }
        }
      else
        {
          SET_FLAG(cqinfo.read, CQ_READ_HILBERT, isEqualStr(CQ_HILBERT_NAME,CQ_DFLT_PARTITIONER));
          if (verbose)
            printf("Info: Setting flag to default: %s = %s\n", CQ_FLAG_PARTITIONER, CQ_DFLT_PARTITIONER);
          //        strcpy(cqinfo.blocksfl, CQ_DFLT_BLOCKS);
          //        printf("Info: Setting blocks file to default: %s\n", CQ_DFLT_BLOCKS);
        }
    }
  /* Decide the name of the blocks file */
  if (!strlen(cqinfo.blocksfl)) /* If not user-specified, search in input */
    {
      if (FLAG_UP(cqinfo.have, CQ_HAVE_INPUT))
        {
          if (!getFlag(fci, CQ_FLAG_BLOCKSFILE, cqinfo.blocksfl)) /* CQ Input-specified? */
            {
              if (FLAG_UP(cqinfo.read, CQ_READ_HILBERT)) /* No flag in CQ input: Use defaults */
                {
                  strcpy(cqinfo.blocksfl, CQ_DFLT_BLOCKS_HILBERT);
                  if (verbose)
                    printf("Info: Setting blocks file to default for Hilbert partitioning: %s\n", CQ_DFLT_BLOCKS_HILBERT);
                }
              else
                {
                  strcpy(cqinfo.blocksfl, CQ_DFLT_BLOCKS_RASTER);
                  if (verbose)
                    printf("Info: Setting blocks file to default for non-Hilbert partitioning: %s\n", CQ_DFLT_BLOCKS_RASTER);
                }
            }
          else if (verbose)
            printf("Info: Using %s to set blocks file name: %s\n", cqinfo.inputfl, cqinfo.blocksfl);
        }
      else /* No CQ input: Use defaults */
        {
          if (FLAG_UP(cqinfo.read, CQ_READ_HILBERT))
            {
              strcpy(cqinfo.blocksfl, CQ_DFLT_BLOCKS_HILBERT);
              if (verbose)
                printf("Info: Setting blocks file to default for Hilbert partitioning: %s\n", CQ_DFLT_BLOCKS_HILBERT);
            }
          else /* We still keep this, in case the default changes */
            {
              strcpy(cqinfo.blocksfl, CQ_DFLT_BLOCKS_RASTER);
              if (verbose)
                printf("Info: Setting blocks file to default for non-Hilbert partitioning: %s\n", CQ_DFLT_BLOCKS_RASTER);
            }
        }
    }
  else if (verbose) /* User-defined, already present before the call to this function */
    printf("Info: Keeping predefined blocks file value: %s\n", cqinfo.blocksfl);
  /* Have blocks file name !! Check that it can be read */
  /* Note: access returns 0 if the file can be opened with the requested permissions */
  SET_FLAG(cqinfo.have, CQ_HAVE_BLOCKS,!access(cqinfo.blocksfl, R_OK));
  if (verbose && !FLAG_UP(cqinfo.have, CQ_HAVE_BLOCKS))
    printf("WARNING: CQ blocks file not found for reading: %s\n", cqinfo.blocksfl);
  /* Decide the name of the density files */
  if (!strlen(cqinfo.densityfl)) /* If not user-specified, search in input */
    {
      strcpy(cqinfo.densityfl, CQ_DFLT_DENSITY);
      if (verbose)
        printf("Info: Setting density file prefix to default: %s\n", CQ_DFLT_DENSITY);
    }
  else /* User-defined, already present before the call to this function */
    {
      if (verbose)
        printf("Info: Keeping predefined density file prefix value: %s\n", cqinfo.densityfl);
    }
  /* From here, try to collect information from the output file */
  if (cqinfo.cores > 0) /* Value, still unknown */
    {
      if (verbose)
        printf("Info: Keeping value of number of processes = %d\n", cqinfo.cores);
    }
  else if (!FLAG_UP(cqinfo.have, CQ_HAVE_OUTPUT)) /* Therefore, if we don't have the output, we cannot check (without guessing) */
    {
      if (verbose)
        printf("WARNING: CQ output file missing: Unknown number of processes\n");
    }
  else
    {
      /* Have density file name stub. Check whether we can open the files for reading */
      /* First, we need to figure out the number of processes of the job */
      if (scanPattern(fco, CQ_PATTERN_CORES, &cqinfo.cores) != 1)
        {
          if (verbose)
            printf("WARNING: Number of processes not found in the CQ output file\n");
          cqinfo.cores = 0; /* This indicates that the no of cores is not known */
        }
      else
        {
          if (verbose)
            printf("Info: Number of processes = %d\n", cqinfo.cores);
          /* Check whether density files are accessible for reading */
          problemcnt = 0;
          for (i = 1; i <= cqinfo.cores; ++i)
            {
              if (makeFilename(tmpfile, &cqinfo, i, EXT_NONE))
                {
                  if (verbose)
                    printf("WARNING: Density file names are too long and cannot be used\n");
                  problemcnt = -1; /* We have problems, but is not that files cannot be accessed */
                  break;
                }
              /* Remember: access returns 0 if access is OK */
              if (access(tmpfile, R_OK))
                ++problemcnt;
            }
          // TODO Check availability also if keeping stub name
          if (problemcnt == 0)
            {
              SET_FLAG(cqinfo.have, CQ_HAVE_DENSITY, TRUE);
              if (verbose)
                printf("Info: All density files are available for reading\n");
            }
          else
            {
              SET_FLAG(cqinfo.have, CQ_HAVE_DENSITY, FALSE);
              if (problemcnt > 0)
                if (verbose)
                  printf("WARNING: %d density files not available for reading\n", problemcnt);
            }
        }
    }
  if (cqinfo.pt[0] != 0 && cqinfo.pt[1] != 0 && cqinfo.pt[2] != 0)
    {
      if (verbose)
        printf("Info: Keeping value of number of grid points = %d %d %d\n", cqinfo.pt[0], cqinfo.pt[1], cqinfo.pt[2]);
    }
  else if (!FLAG_UP(cqinfo.have, CQ_HAVE_OUTPUT))
    {
      if (verbose)
        printf("WARNING: CQ output file missing: Unknown number of grid points\n");
    }
  else
    {
      /* Try to read the number of grid points */
      if (scanPattern(fco, CQ_PATTERN_GRIDPTS, &cqinfo.pt[0], &cqinfo.pt[1], &cqinfo.pt[2]) != 3)
        {
          if (verbose)
            printf("WARNING: Number of grid points not found in the CQ output file\n");
          cqinfo.pt[0] = 0; /* This indicates that the no of grid points is not known */
          cqinfo.pt[1] = 0;
          cqinfo.pt[2] = 0;
        }
      else
        {
          if (verbose)
            printf("Info: Number of grid points = %d %d %d\n", cqinfo.pt[0], cqinfo.pt[1], cqinfo.pt[2]);
        }
    }
  if (cqinfo.blockpt[0] != 0 && cqinfo.blockpt[1] != 0 && cqinfo.blockpt[2] != 0) /* Value, still unknown */
    {
      if (verbose)
        printf("Info: Keeping value of number of grid points per block = %d %d %d\n", cqinfo.blockpt[0], cqinfo.blockpt[1],
            cqinfo.blockpt[2]);
    }
  else if (!FLAG_UP(cqinfo.have, CQ_HAVE_OUTPUT))
    {
      if (verbose)
        printf("WARNING: CQ output file missing: Unknown number of grid points per block\n");
    }
  else
    {
      /* Try to read the number of grid points per block */
      if (scanPattern(fco, CQ_PATTERN_BLOCKPTS, &cqinfo.blockpt[0], &cqinfo.blockpt[1], &cqinfo.blockpt[2]) != 3)
        {
          if (verbose)
            printf("WARNING: Number of grid points per block not found in the CQ output file\n");
          cqinfo.blockpt[0] = 0; /* This indicates that the no of grid points is not known */
          cqinfo.blockpt[1] = 0;
          cqinfo.blockpt[2] = 0;
        }
      else
        {
          if (verbose)
            printf("Info: Number of grid points per block = %d %d %d\n", cqinfo.blockpt[0], cqinfo.blockpt[1], cqinfo.blockpt[2]);
        }
    }
  if (FLAG_UP(cqinfo.units, CQ_UNITS_BOHR))
    conv = BOHR_TO_ANGSTROM;
  else
    conv = 1.0;
  if (cqinfo.lvm[0] > ZEROTOL && cqinfo.lvm[1] > ZEROTOL && cqinfo.lvm[2] > ZEROTOL) /* Value, already known */
    {
      if (verbose)
        printf("Info: Keeping value of dimensions = %f %f %f %s\n", cqinfo.lvm[0] / conv, cqinfo.lvm[1] / conv, cqinfo.lvm[2]
            / conv, REPORT_DISTUNITS(cqinfo.units));
    }
  else if (!FLAG_UP(cqinfo.have, CQ_HAVE_OUTPUT))
    {
      if (verbose)
        printf("WARNING: CQ output file missing: Unknown dimensions\n");
    }
  else
    {
      /* Try to read the dimensions of the cell from the output (this is not necessarily what will be read in) */
      if (scanPattern(fco, CQ_PATTERN_BOX, &cqinfo.lvm[0], &cqinfo.lvm[1], &cqinfo.lvm[2], units) != 4)
        {
          if (verbose)
            printf("WARNING: Dimensions of the simulation box not found in the CQ output file\n");
          cqinfo.lvm[0] = 0.0; /* This indicates that the dimensions are not known */
          cqinfo.lvm[1] = 0.0;
          cqinfo.lvm[2] = 0.0;
        }
      else
        {
          /* The dimensions are always stored in Angstrom */
          for (i = 0; i < 3; ++i)
            {
              cqinfo.lvm[i] *= conv;
            }
          if (verbose)
            {
              /* The dimensions are reported in user's units */
              printf("Info: Dimensions (according to %s): a = %f b = %f c = %f %s\n", cqinfo.outputfl, cqinfo.lvm[0] / conv,
                  cqinfo.lvm[1] / conv, cqinfo.lvm[2] / conv, units);
              if ((FLAG_UP(cqinfo.units, CQ_UNITS_BOHR) && !(isEqualStr(units, CQ_TAG_BOHR1) || isEqualStr(units, CQ_TAG_BOHR2)))
                  || (FLAG_UP(cqinfo.units, CQ_UNITS_ANGSTROM) && !isEqualStr(units, CQ_TAG_ANGSTROM)))
                printf("WARNING: Input- and output-file dimensions do not agree\n");
            }
        }
    }
  for (i = 0; i < 3; ++i)
    if (!limits.strides[i]) /* Not initialized */
      {
        limits.strides[i] = 1; /* By default, read all grid points*/
      }
  /* In order to establish the limits, we need to know both the no. of grid points and the dimensions */
  /* Otherwise, we keep whatever values we have (from the user) */
  if (cqinfo.lvm[0] > ZEROTOL && cqinfo.lvm[1] > ZEROTOL && cqinfo.lvm[2] > ZEROTOL && cqinfo.pt[0] > 0 && cqinfo.pt[1] > 0
      && cqinfo.pt[2] > 0)
    {
      /* First, check that the user-provided grid shifts and dimensions are reasonable */
      /* If they are, keep them; if not, reset them and warn the user */
      for (i = 0; i < 3; ++i)
        {
          /* Grid limits */
          if (limits.gridshft[i] >= cqinfo.pt[i] - limits.strides[i])
            {
              if (verbose)
                printf("WARNING: User-provided grid shift is not valid for coord. %s; resetting\n", REPORT_COORDINATE(i));
              limits.gridshft[i] = 0;
            }
          if (limits.griddims[i] == 0) /* The value is not set yet; do it, from shift to cell boundary */
            {
              limits.griddims[i] = (int) ceil((cqinfo.pt[i] - limits.gridshft[i]) / (double) limits.strides[i]);
            }
          if (limits.gridshft[i] + limits.griddims[i] * limits.strides[i] > cqinfo.pt[i])
            {
              if (verbose)
                printf("WARNING: User-provided grid window is not valid for coord. %s; setting to maximum\n",
                    REPORT_COORDINATE(i));
              limits.griddims[i] = (int) ceil((cqinfo.pt[i] - limits.gridshft[i]) / (double) limits.strides[i]);
            }
          /* Real limits */
          /* The shifts need not be exactly on grid points; in fact, they can be anything, */
          /* since real shifts are used to shift the whole data e.g. on output */
          /* The dimensions (loaded window), on the other hand, must match the spacing, including stride */
          if (limits.realshft[i] < ZEROTOL)
            limits.realshft[i] = 0.0; /* Make sure the number is zero */
          //          limits.realdims[i] = conv * (limits.griddims[i] - 1) * limits.strides[i] * (cqinfo.lvm[i] / (cqinfo.pt[i] - 1));
          /* Expect stored values in Angstrom */
          limits.realdims[i] = (limits.griddims[i] - 1) * limits.strides[i] * (cqinfo.lvm[i] / (cqinfo.pt[i] - 1));
        }
    }
  if (verbose)
    {
      /* Even though units are stored in Angstrom, reports are in the units of the CQ files */
      printf("Info: Grid limits: Shifts = %d %d %d\n", limits.gridshft[0], limits.gridshft[1], limits.gridshft[2]);
      printf("Info: Grid limits: Window = %d %d %d\n", limits.griddims[0], limits.griddims[1], limits.griddims[2]);
      printf("Info: Real limits: Shifts = %8.3f %8.3f %8.3f %s\n", limits.realshft[0] / conv, limits.realshft[1] / conv,
          limits.realshft[2] / conv, REPORT_DISTUNITS(cqinfo.units));
      printf("Info: Real limits: Window = %8.3f %8.3f %8.3f %s\n", limits.realdims[0] / conv, limits.realdims[1] / conv,
          limits.realdims[2] / conv, REPORT_DISTUNITS(cqinfo.units));
      printf("Info: Strides = %d %d %d\n", limits.strides[0], limits.strides[1], limits.strides[2]);
    }
}

int getFlag(FILE *fp, char *flag, char *value)
{
  char line[MAXLIN], tmp[MAXLIN], tmp2[MAXLIN];
  char *occurrence;
  int i, found = 0;

  /* Change flag to lowercase */
  for (i = 0; flag[i] != '\0' && i < MAXLIN - 1; ++i)
    tmp[i] = tolower(flag[i]);
  tmp[i] = '\0';
  /* Rewind file to the beginning, in case it is not */
  fseek(fp, 0L, SEEK_SET);
  while (!found && fgets(line, MAXLIN, fp) != NULL)
    {
      /* Line to lowercase (ignoring comments) */
      for (i = 0; line[i] != '\0' && line[i] != '#' && i < MAXLIN - 1; ++i)
        tmp2[i] = tolower(line[i]);
      tmp2[i] = '\0';
      if ((occurrence = strstr(tmp2, tmp)) != NULL)
        {
          //        sscanf(occurrence, "%*s %s", value);
          // To get the value, we still use the original line (i.e. the case is not ignored)
          // The index is calculated as difference of pointers
          sscanf(&line[occurrence - tmp2], "%*s %s", value);
          return 1; /* Found */
        }
    }
  value[0] = '\0';
  return 0; /* Not found */
}

int findBlockLabel(FILE *fp, char *label)
{
  char line[MAXLIN], tmp[MAXLIN], tmp2[MAXLIN];
  int end;
  int i;

  /* Start the label by the block prefix */
  strcpy(tmp, CQ_LABEL_BLOCK);
  end = strlen(tmp);
  /* Change flag to lowercase */
  for (i = 0; label[i] != '\0' && i < MAXLIN - 1; ++i)
    tmp[end + i] = tolower(label[i]);
  tmp[end + i] = '\0';
  /* Rewind file to the beginning, in case it is not */
  fseek(fp, 0L, SEEK_SET);
  while (fgets(line, MAXLIN, fp) != NULL)
    {
      /* Line to lowercase (ignoring comments) */
      for (i = 0; line[i] != '\0' && line[i] != '#' && i < MAXLIN - 1; ++i)
        tmp2[i] = tolower(line[i]);
      tmp2[i] = '\0';
      if (strstr(tmp2, tmp) != NULL)
        return 1; /* The next line is the block data */
    }
  return 0; /* Not found */
}

int isTrue(char *value)
{
  int i;
  char tmp[MAXLIN];

  for (i = 0; value[i] != '\0' && i < MAXLIN - 1; ++i)
    tmp[i] = tolower(value[i]);
  tmp[i] = '\0';
  if (!strcmp(tmp, "true") || !strcmp(tmp, "t") || !strcmp(tmp, ".true."))
    return TRUE;
  else
    return FALSE; /* Anything that is not true is false (i.e. no complaint) */
}

int isEqualStr(char *value1, char *value2)
{
  int i;
  char tmp1[MAXLIN], tmp2[MAXLIN];

  for (i = 0; value1[i] != '\0' && i < MAXLIN - 1; ++i)
    tmp1[i] = tolower(value1[i]);
  tmp1[i] = '\0';
  for (i = 0; value2[i] != '\0' && i < MAXLIN - 1; ++i)
    tmp2[i] = tolower(value2[i]);
  tmp2[i] = '\0';

  if (!strcmp(tmp1, tmp2))
    return TRUE;
  else
    return FALSE;
}

/* This function is similar to fscanf, except that it searches within lines
 * (not necessarily at the beginning). The pattern can span several lines, separated by '\n'.
 * NOTE: The pattern finishes with the last conversion specifier
 *       (indicated by %, like in printf). In other words, constant strings
 *       after the last % specifier are ignored.
 */
int scanPattern(FILE *fp, char *pattern, ...)
{
  char lpat[MAXLIN]; /* Local copy of the pattern */
  int len; /* Length of the pattern */
  int expected; /* Number of expected items */
  int read; /* Number of items correctly read */
  va_list ap;
  char *block; /* A block of lines that may contain the pattern */
  char *start; /* Pointer to the actual start of the pattern within the block */
  char *ptpat; /* Pointer within the pattern */
  int lines; /* Number of lines to be considered in the block */
  char prefix[MAXLIN]; /* The first part of the pattern, used for searching */
  char *tmp1, *tmp2; /* Auxiliary pointers to the end of the prefix */
  char *prxend; /* The actual pointer to the end of the prefix */
  int i;

  /* Make a local copy of the pattern, truncating if necessary */
  memcpy(lpat, pattern, ((len = strlen(pattern)) >= MAXLIN) ? (len = MAXLIN - 1) : len);
  lpat[len] = '\0'; /* Terminate string */
  if (len == 0)
    return 0; /* NULL pattern -> No matches */
  //printf("%d %d\n\n%s\n\n%s\n", len, strlen(pattern), pattern, lpat);
  read = 0;
  lines = 1;
  va_start(ap, pattern);
  block = (char *) malloc((size_t) MAXLIN); /* Start by assuming it is only one line... */
  ptpat = lpat; /* Start at the beginning */
  while ((ptpat = strstr(ptpat + 1, "\n")) != NULL) /* ...but resized if necessary */
    {
      ++lines;
      block = realloc(block, lines * MAXLIN);
    }
  //printf("%d\n", lines);
  /* Create a prefix to search. This is either the fragment before the first % character,
   * or the first line if there is no % in it
   */
  tmp1 = strstr(pattern, "%");
  tmp2 = strstr(pattern, "\n");
  if (tmp1 != NULL)
    {
      if (tmp2 != NULL)
        {
          if (tmp1 < tmp2)
            prxend = tmp1;
          else
            prxend = tmp2;
        }
      else
        prxend = tmp1;
    }
  else
    {
      if (tmp2 != NULL)
        prxend = tmp2;
      else
        prxend = pattern + strlen(pattern); /* No '%' and no '\n' */
    }
  /* Prepare the prefix */
  for (i = 0; i < prxend - pattern; ++i)
    prefix[i] = pattern[i];
  prefix[prxend - pattern] = '\0'; /* Terminate the string */
  //printf("P = --%s--\n", prefix);

  /* Move to the start of the file, in case we are not there */
  fseek(fp, 0L, SEEK_SET);
  while (fgets(block, MAXLIN, fp) != NULL)
    {
      //printf("%s", block);
      if ((start = strstr(block, prefix)) != NULL) /* Prefix found in this line */
        {
          /* Read the remaining lines of the pattern */
          for (i = 0; i < lines - 1; ++i)
            {
              //printf("Block building:\n%s\n", block);
              if (fgets(&block[strlen(block)], MAXLIN, fp) == NULL) /* Error reading pattern : Stop */
                {
                  free(block);
                  va_end(ap);
                  return -1;
                }
            }
          //printf("Block found:\n%s\n", start);
          /* This is the time to (try and) get the data from the pattern and finish */
          read = vsscanf(start, pattern, ap);
          return (read);
        }
    }
  free(block);
  va_end(ap);
}
