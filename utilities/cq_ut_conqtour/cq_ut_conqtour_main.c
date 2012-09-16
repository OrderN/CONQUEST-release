/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_main.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_main.c *
 *
 * NAME
 *  cq_ut_conqtour_main.c
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

pthread_mutex_t ini_mutex;

int main(int argc, char **argv)
{
  pthread_t pythread;
  struct sigaction action;
  sigset_t mask;

  /* Create thread for python interpreter */
  pthread_mutex_init(&ini_mutex, NULL); /* Initialize mutex and condition to signal */

  pthread_mutex_lock(&ini_mutex); /* Keep the python thread waiting until openGL is initialized */

  // TODO Check error after creating the python thread
  pthread_create(&pythread, NULL, pythonConsole, 0); /* Create the thread and wait until initialization is complete */
  atexit(pythreadFinalize);

  /* Set handle for SIGUSR1, for communication from the Python thread */
  sigfillset (&mask); /* Enable all signals */
  action.sa_handler = awakeMainThread;
  action.sa_mask = mask;
  action.sa_flags = 0;
  sigaction(SIGUSR1, &action, NULL);

  /* Initialize Conqtour and OpenGL */
  init(argc, argv, &idx, &density);

  /* Initialization is complete (as far as we can tell): unlock the python thread */
  pthread_mutex_unlock(&ini_mutex);

  /* Enter Glut main loop */
  glutMainLoop();
}
