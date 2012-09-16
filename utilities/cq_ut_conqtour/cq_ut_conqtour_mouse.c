/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_mouse.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_mouse.c *
 *
 * NAME
 *  cq_ut_conqtour_mouse.c
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

void mouse(int button, int state, int x, int y)
{
  int hits;

  modifiers = glutGetModifiers();
  //printf("Teclas %d\n", modifiers & GLUT_ACTIVE_SHIFT);
  /* Reset the selection (if necessary) */
  //   if(state == GLUT_DOWN && ctrl.selectedatoms == 3)      /* This way, reset after click */
  pthread_mutex_lock(&ctrl.mutex);
  if (ctrl.selectedatoms == 3) /* This way, reset immediately */
    {
      ctrl.selectedatoms = 0;
      prepareMol(molecule, ctrl.spherescale);
      requestRedisplay();
    }
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
      glutTimerFunc(300, resetClickCount, 0);
      ++clicks;
      ctrl.x = x;
      ctrl.y = y;
      if (clicks == 1) /* We can do everything with just one button */
        if (modifiers & GLUT_ACTIVE_SHIFT)
          ctrl.action = TRANSLATE;
        else if (modifiers & GLUT_ACTIVE_CTRL)
          ctrl.action = SCALE;
        else
          ctrl.action = ROTATE;
      else if (clicks > 1)
        {
#ifdef DEBUG
          printf("Double click\n");
#endif
          rendermode = GL_SELECT;
          glRenderMode(GL_SELECT);
          /* We need to call display directly (and it must contain locks for indirect calls) */
          /* Therefore, we unlock and relock */
          pthread_mutex_unlock(&ctrl.mutex);
          display();
          pthread_mutex_lock(&ctrl.mutex);
          //         requestRedisplay();
          rendermode = GL_RENDER;
          hits = glRenderMode(GL_RENDER);
          ctrl.display |= CQ_CTRL_SHOW_SHAPE;
          requestRedisplay();
          //         reshape(ctrl.w, ctrl.h);
#ifdef DEBUG
          printf("H = %d\n", hits);
#endif
          processHits(hits);
        }
    }
  else if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
    {
      ctrl.y = y;
      ctrl.action = SCALE;
    }
  else if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN)
    {
      ctrl.x = x;
      ctrl.y = y;
      ctrl.action = TRANSLATE;
    }
  else
    ctrl.action = NONE;
  pthread_mutex_unlock(&ctrl.mutex);
}

void mouseMotion(int x, int y)
{
  //   int angle,
  int i;
  int scrctr[2]; /* Screen center */
  double ang; /* Angle between 2 vectors */
  double axis[3]; /* Axis perpendicular to 2 vectors */
  double aux1[3];
  double aux2[3];
  double dtmp, tmp[3]; /* For temporary calculations */

  double radius;
  double radtrack; /* Radius of virtual trackball (in pixels) */
  double z; /* Auxiliary coordinate on virtual trackball */
  double lenfield; /* Width of view field */

  pthread_mutex_lock(&ctrl.mutex);
  //printf("Modif %d\n", modifiers & GLUT_ACTIVE_SHIFT);
  if (ctrl.action == NONE) /* Unlock and go */
    {
      pthread_mutex_unlock(&ctrl.mutex);
      return;
    }
  else if (ctrl.action == TRANSLATE)
    {
      lenfield = ctrl.distance * sin(ctrl.field * DEG_TO_RAD);
      /* We have the screen y, but we need x; find it */
      crossProduct(ctrl.normal, ctrl.up, axis);
      for (i = 0; i < 3; ++i)
        {
          ctrl.center[i] += lenfield * (x - ctrl.x) * axis[i] / ctrl.w + lenfield * (y - ctrl.y) * ctrl.up[i] / ctrl.w;
        }
      ctrl.x = x;
      ctrl.y = y;
    }
  else if (ctrl.action == ROTATE)
    {
      //printf("%d %d | %d %d || %f %f %f | %f %f %f\n",
      //  x,y, ctrl.x, ctrl.y, ctrl.up[0], ctrl.up[1], ctrl.up[2], ctrl.normal[0], ctrl.normal[1], ctrl.normal[2]);
      scrctr[0] = ctrl.w / 2;
      scrctr[1] = ctrl.h / 2;
      /* Calculate radius from center (in pixels, on screen plane) */
      radius = pow((x - scrctr[0]), 2.0) + pow((y - scrctr[1]), 2.0);
      radius = sqrt(radius);
      radtrack = 0.95 * (ctrl.h / 2.0);

      /* Control works in two different ways */
      /*  1. Inside a certain circunference around the center, */
      /*     movements are similar to rotating a world map with one's hand */
      /*  2. Outside of that circunference, movements are rotations */
      /*     around the normal to the screen */
      /* The radius is given by 95% half of the window heigth */
      if (radius > 0.95 * radtrack)
        {
          ang = atan2(y - scrctr[1], x - scrctr[0]) - atan2(ctrl.y - scrctr[1], ctrl.x - scrctr[0]);
          rotateVector(RAD_TO_DEG * ang, ctrl.normal, ctrl.up);
        }
      else
        {
          /* Here, we asume we are touching a 3D sphere */
          /* Therefore, from the 2D screen coordinates, we must calculate the third */
          /* using the radius of the virtual circunference */
          aux2[0] = x - scrctr[0];
          aux2[1] = -(y - scrctr[1]); /* The negative changes from screen coords. to (up, normal, upnormal) coords. */
          z = sqrt(pow(radtrack, 2.0) - pow(aux2[0], 2.0) - pow(aux2[1], 2.0));
          aux2[2] = z;

          /* It is possible, that the control point is outside of the virtual sphere, */
          /* while the current point lies inside. This is checked here */
          aux1[0] = ctrl.x - scrctr[0];
          aux1[1] = -(ctrl.y - scrctr[1]);
          ctrl.z = pow(radtrack, 2.0) - pow(aux1[0], 2.0) - pow(aux1[1], 2.0);
          if (ctrl.z > 0.0)
            {
              /* Now that we know that we are within the virtual sphere, */
              /* get third component of auxiliary vector */
              ctrl.z = sqrt(ctrl.z);
              aux1[2] = ctrl.z;

              /* Angle between ctrl and current vectors*/
              // TODO Use functions to calculate angle
              ang = 0.0;
              for (i = 0; i < 3; ++i)
                ang += aux1[i] * aux2[i];
              ang /= sqrt(pow(aux1[0], 2.0) + pow(aux1[1], 2.0) + pow(aux1[2], 2.0));
              ang /= sqrt(pow(aux2[0], 2.0) + pow(aux2[1], 2.0) + pow(aux2[2], 2.0));
              ang = RAD_TO_DEG * acos(ang);

              /* Get the third basis component */
              crossProduct(ctrl.up, ctrl.normal, ctrl.upnormal);

              /* Axis perpendicular to ctrl and current vectors */
              /* (going from ctrl to current) */
              crossProduct(aux2, aux1, axis);

              /* Before we use the axis, express it in world coordinates */
              tmp[0] = axis[0] * ctrl.upnormal[0] + axis[1] * ctrl.up[0] + axis[2] * ctrl.normal[0];
              tmp[1] = axis[0] * ctrl.upnormal[1] + axis[1] * ctrl.up[1] + axis[2] * ctrl.normal[1];
              tmp[2] = axis[0] * ctrl.upnormal[2] + axis[1] * ctrl.up[2] + axis[2] * ctrl.normal[2];
              for (i = 0; i < 3; ++i)
                axis[i] = tmp[i];

              /* Finally, rotate up and eye (normal) vectors around the axis */
              /* The goal is to make the ctrl vector go to the current vector */
              /* This, of course, rotates the whole scene, represented by up and normal */
              rotateVector(ang, axis, ctrl.up);
              rotateVector(ang, axis, ctrl.normal);
            }
          //          else printf("Negative\n");
        }
      ctrl.x = x;
      ctrl.y = y;
    }
  else if (ctrl.action == SCALE)
    {
      dtmp = ctrl.distance;
      ctrl.distance += 0.5 * (ctrl.y - y);
      if (ctrl.distance < 0.0)
        ctrl.distance = dtmp;
      ctrl.y = y;
    }
  requestRedisplay();
  pthread_mutex_unlock(&ctrl.mutex);
}

void resetClickCount(int v)
{
  /* This is to control double clicks */
  clicks = 0;
}

// TODO Check that the selected atoms are not out of bounds
void processHits(int h)
{
  int i, j, index;
  GLuint names; /* Number of names for a hit */
  GLuint item; /* A particular hit */
  GLuint zmin, zmax; /* Distances to objects (relative to observer) */
  GLuint nearestz, closest; /* Closest z and object selected */
  double aux1[3], aux2[3], aux3[3]; /* Auxiliary vectors */
  //    double modulus;

  if (h == 0) /* No selection; reset ongoing selection if any */
    {
      ctrl.selectedatoms = 0;
      prepareMol(molecule, ctrl.spherescale);
      requestRedisplay();
      return;
    }

  // Find the nearest object. This will be selected.
  index = 0;
  nearestz = 0xffffffff; /* Furthest possible */
  for (i = 0; i < h; ++i)
    {
      names = pickbuffer[index++];
      zmin = pickbuffer[index++];
      zmax = pickbuffer[index++];
#ifdef DEBUG
      printf(" Hit %d : %u %u %u\n", i, names, zmin, zmax);
#endif
      /* By design, names = 1, but we keep general */
      for (j = 0; j < names; ++j)
        {
          item = pickbuffer[index++];
#ifdef DEBUG
          printf("   Name %d : %u\n", j, item);
#endif
          if (item == -1)
            return; /* This is NOT a selection, but an error */
          if (zmin < nearestz)
            {
              nearestz = zmin;
              closest = item;
            }
        }
    }
#ifdef DEBUG
  printf("     Closest object = %u %d\n", closest, molecule);
#endif
  if (closest >= molecule && closest < molecule + coords.numatoms) /* A true atom selection */
    {
      /* Add to selection list */
      ctrl.selection[ctrl.selectedatoms] = closest - molecule;
      printf("\nPicked atom = %d (%10.6f, %10.6f, %10.6f) \n", ctrl.selection[ctrl.selectedatoms],
          coords.xyz[ctrl.selection[ctrl.selectedatoms]][0], coords.xyz[ctrl.selection[ctrl.selectedatoms]][1],
          coords.xyz[ctrl.selection[ctrl.selectedatoms]][2]);
      ++ctrl.selectedatoms;
      prepareMol(molecule, ctrl.spherescale);
      requestRedisplay();
      /* We select three atoms at a time (to define a plane) */
      if (ctrl.selectedatoms == 3)
        {
          // TODO Make this message appear in a subwindow of the graphical window
          printf("Info: Three atoms picked = %d %d %d\npyCQ> ", ctrl.selection[0], ctrl.selection[1], ctrl.selection[2]);
          /* Recalculate normal (slice and control) */
          for (i = 0; i < 3; ++i)
            {
              aux1[i] = coords.xyz[ctrl.selection[1]][i] - coords.xyz[ctrl.selection[0]][i];
              aux2[i] = coords.xyz[ctrl.selection[2]][i] - coords.xyz[ctrl.selection[1]][i];
            }
          //         crossProduct(aux1, aux2, densitycut.surfnormal);
          crossProduct(aux1, aux2, aux3);
          /*
           modulus=0.0;
           for(i=0; i < 3; ++i)   modulus += pow(densitycut.surfnormal[i], 2.0);
           modulus = sqrt(modulus);
           */
          //         if(modulus < ZEROTOL)  /* Probably parallel vectors : Reset */
          //         if(normalize(densitycut.surfnormal) < ZEROTOL)  /* Probably parallel vectors : Reset */
          if (normalize(aux3) < ZEROTOL) /* Probably parallel vectors : Reset */
            {
              ctrl.selectedatoms = 0;
              prepareMol(molecule, ctrl.spherescale);
              requestRedisplay();
              return;
            }
          for (i = 0; i < 3; ++i)
            {
              //            densitycut.surfnormal[i] /= modulus;
              //            ctrl.normal[i] = densitycut.surfnormal[i];
              /* Prepare view for display */
              ctrl.normal[i] = aux3[i];
            }
          /* Recalculate up vector (control), pointing opposite to the bisection */
          // NOTE: THIS IS NOT DONE LIKE THIS ANYMORE
          //       Now we use one of the rectangle sides (see below)
          //         modulus=0.0;
          /*
           for(i=0; i < 3; ++i)
           {
           ctrl.up[i]=aux1[i]-aux2[i];
           //            modulus += pow(ctrl.up[i], 2.0);
           }
           //         for(i=0; i < 3; ++i)  ctrl.up[i] /= modulus;
           normalize(ctrl.up);
           */

          /* Take note of a point on the new slice (the second atom of the selection) */
          //         for(i=0; i < 3; ++i)    dindex[currentcut]->inplane[i] = coords.xyz[ctrl.selection[1]][i];

          /* Reset the index, using second atom of the selection as point "inplane" */
          resetIndex(aux3, coords.xyz[ctrl.selection[1]], &idx, &density, &box);

          //         calculateDrawSlice(idx.dindex[idx.curr], &box);
          calculateDrawSlice(&idx, &box, &density);

          /* Set new up vector (control) */
          for (i = 0; i < 3; ++i)
            ctrl.up[i] = idx.dindex[idx.curr]->rectangle[3][i] - idx.dindex[idx.curr]->rectangle[0][i];
          normalize(ctrl.up);

          prepareMol(molecule, ctrl.spherescale);
          requestRedisplay();
        }
    }
  else
    {
      ctrl.selectedatoms = 0;
      prepareMol(molecule, ctrl.spherescale);
      requestRedisplay();
    }
}

