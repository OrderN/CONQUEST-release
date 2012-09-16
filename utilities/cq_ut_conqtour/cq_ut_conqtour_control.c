/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_control.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_control.c *
 *
 * NAME
 *  cq_ut_conqtour_control.c
 * PURPOSE
 *
 * USES
 *
 * AUTHOR
 *  torralba
 * CREATION DATE
 *  Oct 21, 2010
 * MODIFICATION HISTORY
 *
 * *****/

#include <string.h>
#include "cq_ut_conqtour.h"

// TODO Even though ctrl is currently global, it should be a parameter to functions. For the moment. try to do it as much as possible

void setupControl(int argc, char **argv, ctrl_t *c)
{
  //   memset(c, 0, sizeof(ctrl_t));           /* Initialise the structure */
  pthread_mutex_init(&c->mutex, NULL); /* Set up a mutex for safe thread control */
  pthread_cond_init(&c->cond, NULL); /* Set up a condition to test display state */
  c->display |= CQ_CTRL_SHOW_BOX;
  c->display |= CQ_CTRL_SHOW_LEGEND;
  c->display &= ~CQ_CTRL_SHOW_AXES;
  c->display &= ~CQ_CTRL_SHOW_RECTANGLE;
  c->has &= ~CQ_CTRL_HAS_DENSITY;
  c->shftcoords = FALSE;
  c->spherescale = 1.0;
  c->selectedatoms = 0;
  c->pid = getpid();
  c->redisplays = 0;
  memset(&(c->legend), 0, sizeof(legend_t));
  c->legend.x = CQ_LEGEND_DFLT_X;
  c->legend.y = CQ_LEGEND_DFLT_Y;
  c->legend.w = CQ_LEGEND_DFLT_W;
  c->legend.h = CQ_LEGEND_DFLT_H;
  c->legend.bt = CQ_LEGEND_DFLT_BORDER_THICKNESS;
  c->legend.tl = CQ_LEGEND_DFLT_TICK_LENGTH;
  c->legend.tickorigin = CQ_LEGEND_DFLT_TICKORIGIN;
  c->legend.noticks = CQ_LEGEND_DFLT_APPROX_TICKS; /* Used to suggest a number of ticks */
  c->legend.nointer = CQ_LEGEND_DFLT_INTERVALS;
  c->legend.warn = TRUE;
  strcpy(c->legend.format, CQ_LEGEND_DFLT_FORMAT);
  SET_MULTIBIT_FLAG(c->legend.method, CQ_LEGEND_DFLT_METHOD_TICKS, CQ_LEGEND_TICKS_MASK);
  SET_MULTIBIT_FLAG(c->legend.method, CQ_LEGEND_DFLT_METHOD_INTERPOL, CQ_LEGEND_INTERPOL_MASK);
  SET_MULTIBIT_FLAG(c->legend.method, CQ_LEGEND_DFLT_METHOD_INTERVAL, CQ_LEGEND_INTERVAL_MASK);
  SET_MULTIBIT_FLAG(c->legend.method, CQ_LEGEND_DFLT_METHOD_COLOR, CQ_LEGEND_COLOR_MASK);
}

void resetView(ctrl_t *c)
{
  int i;

  /* Calculate center of the cell and radius of inscribing sphere*/
  c->size = 0.0;
  for (i = 0; i < 3; ++i)
    {
      c->center[i] = density.lvm[i] / 2.0;
      c->size += pow(c->center[i], 2.0);
    }
  c->size = sqrt(c->size);
  /* Define initial orientation and distance */
  c->up[0] = 0.0;
  c->up[1] = 1.0;
  c->up[2] = 0.0;
  c->normal[0] = 0.0;
  c->normal[1] = 0.0;
  c->normal[2] = 1.0;
  /* Go a little further (to have a reasonably low field) */
  c->distance = 5.0 * ctrl.size;
  c->field = 2.0 * RAD_TO_DEG * atan2(ctrl.size, ctrl.distance);
#ifdef DEBUG
  printf("Dist %f %f | %f %f %f\n", c->distance, c->field, c->center[0], c->center[1], c->center[2]);
#endif
  /* Make sure that the clipping planes don't hide anything */
  c->near = 1.5; //0.7;
  c->far = 1000.0;
}

/* The difference between resetView and setView is that the former recalculates view,
 * whereas the latter uses the values present in the control structure
 * (except, at present, the field, which is calculated, and near and far,
 * which are hardcoded - I know, I know :-))
 */
void setView(ctrl_t *c)
{
  int i;

//  /* Calculate the radius of inscribing sphere and the field that spans it */
//  c->size = 0.0;
//  for (i = 0; i < 3; ++i)
//    {
//      c->size += pow(density.lvm[i] / 2.0, 2.0);
//    }
//  c->size = sqrt(c->size);
//  c->field = 2.0 * RAD_TO_DEG * atan2(ctrl.size, ctrl.distance);
  /* Make sure that the clipping planes don't hide anything */
  c->near = 1.5; //0.7;
  c->far = 1000.0;
}

void reorientRecenter(ctrl_t *c)
{
  int i;

  for (i = 0; i < 3; ++i)
    {
      ctrl.normal[i] = idx.surfnormal[i];
      ctrl.up[i] = idx.dindex[idx.curr]->rectangle[3][i] - idx.dindex[idx.curr]->rectangle[0][i];
    }
  normalize(ctrl.up);
  recenter(c);
}

void recenter(ctrl_t *c) /* And redisplay view */
{
  int i;

  for (i = 0; i < 3; ++i)
    ctrl.center[i] = density.lvm[i] / 2.0;
  requestRedisplay();
}

