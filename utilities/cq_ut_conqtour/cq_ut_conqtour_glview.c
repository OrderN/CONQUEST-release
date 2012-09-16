/* -----------------------------------------------------------------------------
 * $Id: $
 * ----------------------------------------------------------------------------
 * File cq_ut_conqtour_glview.c
 * ----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_glview.c *
 *
 * NAME
 *  cq_ut_conqtour_glview.c
 * PURPOSE
 *  General display of the openGL scene
 * USES
 *
 * AUTHOR
 *  torralba
 * CREATION DATE
 *  Oct 7, 2010
 * MODIFICATION HISTORY
 *
 * *****/

#include <time.h>
#include <string.h>
#include "cq_ut_conqtour.h"

extern pid_t conqtourPid;

void display(void)
{
  int viewport[4];
  int i, j, k;
  int w, h;
  int jumpline;
  double pt[3];
  static char message[] = "Charge density not loaded - Disabled display";
  char *m;
  char label[MAXLIN];
  int posx, posy;

  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  //printf("Display\n");
  /* Run pending commands before redisplaying */
  runCommandQueue();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  /* If there is no density, don't display anything and return */
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      /* In principle, the window will be hiden at this point, but just in case, show a message and try to hide */
      glColor3f(0.0, 0.0, 0.7);
      glWindowPos2i((ctrl.w - strlen(message) * 10) / 2, ctrl.h / 2);
      for (m = message; *m != '\0'; ++m)
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, *m);
      glutSwapBuffers();
      if (ctrl.redisplays) /* This message counts as displayed */
        {
          ctrl.redisplays = 0;
          pthread_cond_signal(&ctrl.cond); /* Signal the python thread to continue */
        }
      glutHideWindow(); /* Immediately hide the window... */
      checkWindowVisibility(ctrl.has); /* ...and set a timer that periodically checks whether it should be shown */
//      glutShowWindow();
      pthread_mutex_unlock(&ctrl.mutex);
      return;
    }
  /* This would normally be done in shape, but it is not called when loading
   * a density file from the python interpreter (because the window is not touched).
   * We use for picking as well
   */
  if (ctrl.display & CQ_CTRL_SHOW_SHAPE)
    {
      w = glutGet(GLUT_WINDOW_WIDTH);
      h = glutGet(GLUT_WINDOW_HEIGHT);
      glViewport(0, 0, (GLsizei) w, (GLsizei) h);
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluPerspective(ctrl.field, (GLfloat) w / (GLfloat) h, ctrl.near, ctrl.far);
      /* Take note of window size */
      ctrl.w = w;
      ctrl.h = h;
      ctrl.display &= ~CQ_CTRL_SHOW_SHAPE; /* Not needed anymore */
    }

  if (rendermode == GL_SELECT)
    {
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      /* Get current portview */
      glGetIntegerv(GL_VIEWPORT, viewport);
      gluPickMatrix(ctrl.x, ctrl.h - ctrl.y, PICKTOL, PICKTOL, viewport);
      gluPerspective(ctrl.field, (GLfloat) ctrl.w / (GLfloat) ctrl.h, ctrl.near, ctrl.far);
    }

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(ctrl.center[0] + ctrl.normal[0] * ctrl.distance, ctrl.center[1] + ctrl.normal[1] * ctrl.distance, ctrl.center[2]
      + ctrl.normal[2] * ctrl.distance, ctrl.center[0], ctrl.center[1], ctrl.center[2], ctrl.up[0], ctrl.up[1], ctrl.up[2]);

  /* Not strictly necessary, since it would be ignored in GL_RENDER */
  if (rendermode == GL_SELECT)
    {
      glInitNames();
      glPushName(-1); /* This is the first, unusable name in the name stack */
    }
  if (rendermode == GL_RENDER)
    {
      if (ctrl.display & CQ_CTRL_SHOW_BOX)
        glCallList(boundingbox);
      if (ctrl.display & CQ_CTRL_SHOW_AXES)
        glCallList(axes);
      //      glCallList(idx.dindex[idx.curr]->slice);
      //printf("%d\n", idx.dindex[idx.curr]->slice);
    }
  /* Draw and name the density plane, so we can discard selections behind it */
  /* Next line, ignored in GL_RENDER mode */
  //   glLoadName(idx.dindex[idx.curr]->slice);
  glLoadName(molecule + coords.numatoms + 1);
  glCallList(idx.dindex[idx.curr]->slice);
  if (ctrl.display & CQ_CTRL_SHOW_MOLECULE)
    {
      glCallList(molecule);
    }

  /*
   for(i=0; i < idx.dindex[idx.curr]->nocorners ; ++i)
   {
   glPushMatrix();
   glColor3f(0.0, 0.0, 0.7);
   glColor3f(0.0,0.0,i*0.2);
   //printf("%e %e %e\n", densitycut.corner[i][0], densitycut.corner[i][1], densitycut.corner[i][2]);
   glTranslatef(idx.dindex[idx.curr]->corner[i][0], idx.dindex[idx.curr]->corner[i][1], idx.dindex[idx.curr]->corner[i][2]);
   glutSolidSphere(0.3, 20, 20);
   glPopMatrix();
   }
   */
  /*
   glPushMatrix();
   glColor3f(1.0,0.0,0.0);
   glTranslatef(densitycut.corner[pivot[0]][0], densitycut.corner[pivot[0]][1], densitycut.corner[pivot[0]][2]);
   glutWireCube(1.5);
   glPopMatrix();
   glPushMatrix();
   glColor3f(1.0,0.0,0.0);
   glTranslatef(densitycut.corner[pivot[1]][0], densitycut.corner[pivot[1]][1], densitycut.corner[pivot[1]][2]);
   glutWireCube(1.7);
   glPopMatrix();
   glPushMatrix();
   glColor3f(1.0,0.0,0.0);
   glTranslatef(densitycut.corner[pivot[2]][0], densitycut.corner[pivot[2]][1], densitycut.corner[pivot[2]][2]);
   glutWireCube(1.9);
   glPopMatrix();
   glPushMatrix();
   glColor3f(1.0,0.0,0.0);
   glTranslatef(densitycut.corner[pivot[3]][0], densitycut.corner[pivot[3]][1], densitycut.corner[pivot[3]][2]);
   glutWireCube(2.1);
   glPopMatrix();
   */

  if (rendermode == GL_RENDER)
    {
      if (ctrl.display & CQ_CTRL_SHOW_RECTANGLE)
        {
          /**/
          glPushMatrix();
          glColor3f(0.2, 0.0, 0.0);
          glTranslatef(idx.dindex[idx.curr]->rectangle[0][0], idx.dindex[idx.curr]->rectangle[0][1],
              idx.dindex[idx.curr]->rectangle[0][2]);
          //glutSolidSphere(0.3, 20, 20);
          glutWireCube(1.5);
          glPopMatrix();
          glPushMatrix();
          glColor3f(0.4, 0.0, 0.0);
          glTranslatef(idx.dindex[idx.curr]->rectangle[1][0], idx.dindex[idx.curr]->rectangle[1][1],
              idx.dindex[idx.curr]->rectangle[1][2]);
          //glutSolidSphere(0.3, 20, 20);
          glutWireCube(1.7);
          glPopMatrix();
          glPushMatrix();
          glColor3f(0.6, 0.0, 0.0);
          glTranslatef(idx.dindex[idx.curr]->rectangle[2][0], idx.dindex[idx.curr]->rectangle[2][1],
              idx.dindex[idx.curr]->rectangle[2][2]);
          //glutSolidSphere(0.3, 20, 20);
          glutWireCube(1.9);
          glPopMatrix();
          glPushMatrix();
          glColor3f(0.8, 0.0, 0.0);
          glTranslatef(idx.dindex[idx.curr]->rectangle[3][0], idx.dindex[idx.curr]->rectangle[3][1],
              idx.dindex[idx.curr]->rectangle[3][2]);
          //glutSolidSphere(0.3, 20, 20);
          glutWireCube(2.1);
          glPopMatrix();
          /**/
          glBegin(GL_LINE_LOOP);
          glColor3f(1.0, 0.0, 1.0);
          glVertex3f(idx.dindex[idx.curr]->rectangle[0][0], idx.dindex[idx.curr]->rectangle[0][1],
              idx.dindex[idx.curr]->rectangle[0][2]);
          glVertex3f(idx.dindex[idx.curr]->rectangle[1][0], idx.dindex[idx.curr]->rectangle[1][1],
              idx.dindex[idx.curr]->rectangle[1][2]);
          glVertex3f(idx.dindex[idx.curr]->rectangle[2][0], idx.dindex[idx.curr]->rectangle[2][1],
              idx.dindex[idx.curr]->rectangle[2][2]);
          glVertex3f(idx.dindex[idx.curr]->rectangle[3][0], idx.dindex[idx.curr]->rectangle[3][1],
              idx.dindex[idx.curr]->rectangle[3][2]);
          glEnd();
        }

      glBegin(GL_LINES);
      glColor3f(0.0, 0.0, 1.0);
      for (i = 0; i < idx.dindex[idx.curr]->nocorners; ++i)
        {
          glVertex3f(idx.dindex[idx.curr]->corner[i][0], idx.dindex[idx.curr]->corner[i][1], idx.dindex[idx.curr]->corner[i][2]);
          //   if(i+1 < idx.dindex[idx.curr]->nocorners)
          glVertex3f(idx.dindex[idx.curr]->corner[(i + 1) % idx.dindex[idx.curr]->nocorners][0], idx.dindex[idx.curr]->corner[(i
              + 1) % idx.dindex[idx.curr]->nocorners][1],
              idx.dindex[idx.curr]->corner[(i + 1) % idx.dindex[idx.curr]->nocorners][2]);
          /**
           else
           glVertex3f(idx.dindex[idx.curr]->corner[0][0],
           idx.dindex[idx.curr]->corner[0][1],
           idx.dindex[idx.curr]->corner[0][2]);
           **/

        }
      glEnd();

      /*
       glBegin(GL_POLYGON);
       glColor3f(0.0,0.0,0.2);
       for(i=0; i < idx.dindex[idx.curr]->nocorners ; ++i)
       {
       glVertex3f(idx.dindex[idx.curr]->corner[i][0], idx.dindex[idx.curr]->corner[i][1], idx.dindex[idx.curr]->corner[i][2]);
       }
       glEnd();
       */
      /*
       //glRasterPos2i(0, 0);            // Rasterizes in world coordinates (i.e. rotates with objects, etc.)
       glAlphaFunc(GL_GEQUAL, 0.5);      // Rasterizes at fixed window positions (does not rotate; useful to show a scale or logo)
       glEnable(GL_ALPHA_TEST);
       glWindowPos2i(0, 0);
       glDrawPixels(logo_width, logo_height, GL_RGBA, GL_UNSIGNED_BYTE, logo_image);
       */
//      unsigned char image[10*30*3];
//      char color;
//      int x,y,n;
//
//      memset(image, 100, 10*30*3);
//      for(y=0; y < 10; ++y)
//        {
//          for(x=0; x < 30; ++x)
//            {
//              for(n=0; n < 3; ++n)
//                {
//                  if(n==0)
//                    {
//                      color = 0;
//                    }
//                  else if(n == 1)
//                    {
//                      color = 256-x;
//                    }
//                  else
//                    {
//                      color = 200 + y;
//                    }
//                  image[(y*30+x)*3 + n] = color;//((l*10+m)*3 + n)%256;
//                  //printf("%d\n", image[(y*10+x)*3 + n]);
//                }
//            }
//        }
//      glWindowPos2i(100, 100);
//      glDrawPixels(10, 30, GL_RGB, GL_UNSIGNED_BYTE, image);
//      char label[] = "-1e-4";
//      glColor3f(0.0, 0.0, 0.0);
//      //glColor3f(0.5, 0.5, 0.5);
//      glWindowPos2i(140, 100);
//      for(x = 0; x < strlen(label); ++x)
//        {
//          glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, label[x]);
//        }
//      glColor3f(0.5, 0.5, 0.5);
//      //glColor3f(0.0, 0.0, 0.0);
//      glWindowPos2i(141, 99);
//      for(x = 0; x < strlen(label); ++x)
//        {
//          glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, label[x]);
//        }
      if (ctrl.display & CQ_CTRL_SHOW_LEGEND)
        {
//          int windowpos[4];

          glWindowPos2i(density.legend.x, ctrl.h - density.legend.y - density.legend.h - 2 * density.legend.bt);
//          glRasterPos2i(density.legend.x, ctrl.h - density.legend.y - density.legend.h - 2);
//          glRasterPos2i(0,0);
//          glGetIntegerv(GL_CURRENT_RASTER_POSITION, windowpos);
//          printf("%d %d %d %d\n", windowpos[0], windowpos[1], windowpos[2], windowpos[3]);
//          glRasterPos2i(density.legend.x, ctrl.h - density.legend.y - density.legend.h - 2);
          glCallList(legend);
          /* The labels are better done here, not in the command list */
          for(i = 0; i < density.legend.marks.pixels; ++i)
            {
              if(density.legend.marks.label[i] == YES) /* No label without tick :-P */
                {
                  sprintf(label, density.legend.format, density.legend.marks.val[i]);
                  glColor3f(0.0, 0.0, 0.0);
                  posx = density.legend.x + density.legend.w + 2 * density.legend.bt + density.legend.tl + 2;
                  posy = ctrl.h - density.legend.y - density.legend.h - density.legend.bt + i - 3; /* This 3, related to font */
                  glWindowPos2i(posx, posy);
                  for(j = 0; j < strlen(label); ++j)
                    {
                      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, label[j]);
                    }
                  glColor3f(0.5, 0.5, 0.5);
                  glWindowPos2i(posx + 1, posy - 1);
                  for(j = 0; j < strlen(label); ++j)
                    {
                      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, label[j]);
                    }
                }
            }

//
//          int i, j, n;
//          legend_t *l;
//          int pos; /* Position within the image array */
//          int size; /* Size of the image */
//          unsigned char *image; /* The legend bar as an RGBA image (i.e. RGB with transparency) */
////          unsigned char image[30000]; /* The legend bar as an RGBA image (i.e. RGB with transparency) */
//
//          l = &(density.legend);
//          size = (l->w + 3) * (l->h + 3) * 4;
//          printf("Size %d %d %d %d %d\n", size, l->x, l->y, l->w, l->h);
//          image = (unsigned char *)malloc((size_t)size * sizeof(unsigned char));
//          memset(image, 255, size);
//
//          for (i = 0; i < l->h + 3; ++i)
//            {
//              for(j = 0; j < l->w + 3; ++j)
//                {
//                  pos = (i * (l->w + 3) + j) * 4;
//                  if( i < l->h + 2 && j < l->w + 2)
//                    {
//                      if(i == 0 || i == l->h + 1 || j == 0 || j == l->w + 1)
//                        {
//                          /* Black border */
//                          image[pos] = 0;
//                          image[pos + 1] = 0;
//                          image[pos + 2] = 0;
//                          image[pos + 3] = 255; /* Opaque */
//                        }
//                      else
//                        {
//                          /* Color */
//                          image[pos] = 255; /* FOR THE MOMENT, red */
//                          image[pos + 1] = 0;
//                          image[pos + 2] = 0;
//                          image[pos + 3] = 255;
//                        }
//                    }
//                  else /* Last column and row : Draw a grey shadow */
//                    {
//                      if( i != 0 && j != 0)
//                        {
//                          /* Grey */
//                          image[pos] = 128;
//                          image[pos + 1] = 128;
//                          image[pos + 2] = 128;
//                          image[pos + 3] = 255;
//                        }
//                    }
//                }
//            }
//          printf("%d %d\n", l->x - 1, l->y - 1);
//          glWindowPos2i(l->x - 1, l->y - 1);
////          glWindowPos2i(100, 100);
//          glDrawPixels(l->w + 3, l->h + 3, GL_RGBA, GL_UNSIGNED_BYTE, image);
////          glDrawPixels(10, 30, GL_RGB, GL_UNSIGNED_BYTE, image);
        }
    }

  /*
   glPointSize(4);
   glBegin(GL_POINTS);
   for(i=0; i < idx.dindex[idx.curr]->pt[1]; ++i)
   {
   jumpline=NO;
   for(j=0; j < 3; ++j)   pt[j] = idx.dindex[idx.curr]->rectangle[0][j]
   + idx.dindex[idx.curr]->yaxis[j] * i;
   for(j=0; j < idx.dindex[idx.curr]->pt[0]; ++j)
   {
   if( pt[0] < -ZEROTOL2
   || pt[1] < -ZEROTOL2
   || pt[2] < -ZEROTOL2
   || pt[0] > density.lvm[0] + ZEROTOL2
   || pt[1] > density.lvm[1] + ZEROTOL2
   || pt[2] > density.lvm[2] + ZEROTOL2 )
   {
   glColor3f(0.0,0.0,1.0);
   glVertex3f(pt[0], pt[1], pt[2]);
   if(jumpline == NO)   ;
   else              {  break; }
   }
   else
   {
   glColor3f(0.5,0.5,0.5);
   glVertex3f(pt[0], pt[1], pt[2]);
   jumpline = YES;
   }
   for(k=0; k < 3; ++k)   pt[k] += idx.dindex[idx.curr]->xaxis[k];
   }
   }
   glEnd();
   */

  if (rendermode == GL_RENDER)
    glutSwapBuffers();
  //printf("R %d\n", ctrl.redisplays);
  if (ctrl.redisplays)
    {
      ctrl.redisplays = 0;
      pthread_cond_signal(&ctrl.cond); /* Signal the python thread to continue */
    }
  pthread_mutex_unlock(&ctrl.mutex); /* Be safe */
}

void reshape(int w, int h)
{
  int i;

  pthread_mutex_lock(&ctrl.mutex);
  //printf("Shape\n");
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(ctrl.field, (GLfloat) w / (GLfloat) h, ctrl.near, ctrl.far);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  /* The up vector is the minus y axis, so that the origin is at the upper-left corner */
  gluLookAt(ctrl.center[0] + ctrl.normal[0] * ctrl.distance, ctrl.center[1] + ctrl.normal[1] * ctrl.distance, ctrl.center[2]
      + ctrl.normal[2] * ctrl.distance, ctrl.center[0], ctrl.center[1], ctrl.center[2], ctrl.up[0], ctrl.up[1], ctrl.up[2]);
  /* Take note of window size */
  ctrl.w = w;
  ctrl.h = h;
  //printf("Shape\n");
  pthread_mutex_unlock(&ctrl.mutex);
}

void iluminate(void)
{
  float ambient[] =
    { 0.0, 0.0, 0.0, 1.0 };
  float diffuse[] =
    { 1.0, 1.0, 1.0, 1.0 };
  float specular[] =
    { 1.0, 1.0, 1.0, 1.0 };
  float specular2[] =
    { 0.3, 0.3, 0.3, 1.0 };
  float emission[] =
    { 0.0, 0.0, 0.0, 1.0 };
  float position0[] =
    { 1.0, 0.5, 0.0, 0.0 };
  float position1[] =
    { -1.0, 0.5, 0.0, 0.0 };
  float lmodel_ambient[] =
    { 0.85, 0.85, 0.85, 1.0 };

  glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
  glLightfv(GL_LIGHT0, GL_POSITION, position0);
  glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT1, GL_POSITION, position1);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);

  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular2);
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emission);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 40);
}

