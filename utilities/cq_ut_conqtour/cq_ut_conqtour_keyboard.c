/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_keyboard.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_keyboard.c *
 *
 * NAME
 *  cq_ut_conqtour_keyboard.c
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

#include <pthread.h>
extern pthread_t pythread;

void keyboard(unsigned char key, int x, int y)
{
  int i;

  //modifiers = glutGetModifiers();
  //printf("Teclas %d\n", modifiers & GLUT_ACTIVE_SHIFT);
  pthread_mutex_lock(&ctrl.mutex);
  /* Ignore keys (except ESC) if the density is not loaded */
  if (!FLAG_UP(ctrl.has, CQ_CTRL_HAS_DENSITY) && key != 27)
    {
      pthread_mutex_unlock(&ctrl.mutex);
      return;
    }
  switch (key)
    {
  case 'i':
    writePpmImage("screenshot.ppm");
  case 'a':
    if (idx.curr > 0)
      --idx.curr;
    calculateDrawSlice(&idx, &box, &density);
    requestRedisplay();
    break;
  case 's':
    if (idx.curr < idx.cuts - 1)
      ++idx.curr;
    calculateDrawSlice(&idx, &box, &density);
    requestRedisplay();
    break;
  case 'r':
    ctrl.display ^= CQ_CTRL_SHOW_RECTANGLE;
    requestRedisplay();
    break;
  case 'R':
    resetView(&ctrl);
    requestRedisplay();
    break;
  case 'b':
    ctrl.display ^= CQ_CTRL_SHOW_BOX;
    requestRedisplay();
    break;
  case 'x':
    ctrl.display ^= CQ_CTRL_SHOW_AXES;
    requestRedisplay();
    break;
  case 'm':
    ctrl.display ^= CQ_CTRL_SHOW_MOLECULE;
    requestRedisplay();
    break;
  case 'l':
    ctrl.display ^= CQ_CTRL_SHOW_LEGEND;
    requestRedisplay();
    break;
  case 'p':
    printPymolScript(&idx, &density, "slice.py");
    break;
  case 'C': /* Reorient to see plane from "above", and recenter */
    reorientRecenter(&ctrl);
    break;
  case 'c': /* Recenter (but don't rescale) */
    recenter(&ctrl);
    break;
  case 27: /* ESC */
    printf("\n");
    exit(0);
    }
  pthread_mutex_unlock(&ctrl.mutex);
}

void specialKeys(int key, int x, int y)
{
  pthread_mutex_lock(&ctrl.mutex);
  modifiers = glutGetModifiers();
  switch (key)
    {
  case GLUT_KEY_UP:
    if (idx.curr > 0)
      --idx.curr;
    calculateDrawSlice(&idx, &box, &density);
    requestRedisplay();
    break;
  case GLUT_KEY_DOWN:
    if (idx.curr < idx.cuts - 1)
      ++idx.curr;
    calculateDrawSlice(&idx, &box, &density);
    requestRedisplay();
    break;
  case GLUT_KEY_LEFT:
    if (ctrl.spherescale > 0.3)
      {
        ctrl.spherescale -= 0.05;
        prepareMol(molecule, ctrl.spherescale);
        requestRedisplay();
      }
    break;
  case GLUT_KEY_RIGHT:
    if (ctrl.spherescale < 1.0)
      {
        ctrl.spherescale += 0.05;
        prepareMol(molecule, ctrl.spherescale);
        requestRedisplay();
      }
    requestRedisplay();
    break;
    }
  pthread_mutex_unlock(&ctrl.mutex);
}

