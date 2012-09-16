/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_pythread.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_pythread.c *
 *
 * NAME
 *  cq_ut_conqtour_pythread.c
 * PURPOSE
 *  Create and handle a thread where a python interpreter will run
 * USES
 *
 * AUTHOR
 *  torralba
 * CREATION DATE
 *  Oct 19, 2010
 * MODIFICATION HISTORY
 *
 * *****/

#include <Python.h>
#include "cq_ut_conqtour.h"

extern pthread_mutex_t ini_mutex;
extern pthread_cond_t ini_cond;

void *pythonConsole(void *id)
{
  long tid;
  char *name;
  FILE *fs;

  tid = (long) id;

  /* We must wait until the main thread has initialized openGL */
  /* Otherwise, executing a script could lead to crash */
  pthread_mutex_lock(&ini_mutex);
  //printf("Locked: worker\n");
  //    pthread_cond_wait(&ini_cond, &ini_mutex);
  //printf("Inied -1\n");
  pthread_mutex_unlock(&ini_mutex);

  //printf("Inied 0\n");
  Py_Initialize(); /* First, initialize python */
  //printf("Inied 1\n");
  init_conqtour_pymodule(); /* Make conqtour functions accessible... */
  //printf("Inied 2\n");
  PyRun_SimpleString("from conqtour import *"); /* ...and load them directly */

  //system("sleep 1");

  /* Execute the scripts requested in the command line */
  while (scripts.head != NULL)
    {
      name = popNameFromQueue(&scripts);
      if ((fs = fopen(name, "r")) == NULL)
        {
          fprintf(stderr, "Python script '%s' could not be opened\n", name);
        }
      printf("Running python script: %s\n", name);
      PyRun_SimpleFile(fs, name);
      printf("Done python script: %s\n", name);
      free(name);
      fclose(fs);
    }

  /* Check whether stdin is a terminal */
  /* If it is not, close the thread and continue without python */
  /* We accept scripts via command line, but not via redirection */
  // TODO Run everything in a single call to PyRun_SimpleString and use try...except, so that if readline is not available it continues smoothly
  if (isatty(fileno(stdin)))
    {
      /* Ready to start work in python interpreter */
      PyRun_SimpleString("import readline"); /* Support editing and history */
      PyRun_SimpleString("import rlcompleter"); /* Support autocompletion */
      PyRun_SimpleString("readline.parse_and_bind('tab: complete')");
      PyRun_SimpleString("import sys"); /* Change the prompt so we know where we are */
      PyRun_SimpleString("sys.ps1='pyCQ> '");
      PyRun_SimpleString("sys.ps2='..... '");
      PyRun_SimpleString("del sys"); /* Do not accumulate unnecessary imported modules */

      //PyRun_SimpleString("execfile('rotate.py')");
      /* Start listening to user's input */
#ifdef CQ_INFINITE_PYTHON
      for(;;) /* Swallow Ctrl+D (i.e. EOF) */
#endif
      PyRun_AnyFile(stdin, "cq_ut_conqtour Python terminal");
      /* The following, necessary in case we remove to infinite loop */
      Py_Finalize();
      exit(0);
    }
  else
    {
      printf("Python interpreter disabled: stdin is not a terminal\n");
      Py_Finalize();
      pthread_exit(NULL);
    }
}

// TODO Thread safety in all calls to common routines from the python interface

void pythreadFinalize(void)
{
  /* Probably python should be finalized properly here */

  /* Make sure that the shell produces echo on return */
  /* This is necessary because we are using readline */
  /* and sometimes it doesn't reset the initial echo state */
  /* But if stdin is not a terminal (e.g. is a file), ignore */
  if (isatty(fileno(stdin)))
    system("stty echo");

  /* More cleanup could be done here */
}

/* Schedule a redisplay of the view */
void requestRedisplay(void)
{
  //printf("Request\n");
  ++ctrl.redisplays;
  /*
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    glutHideWindow();
  else
    glutShowWindow();
  */
  /* Strictly speaking, this should not be called from the python thread, but it is the unavoidable minimum */
  glutPostRedisplay();
  kill(ctrl.pid, SIGUSR1); /* Awake the main thread; display should be called and the command queue emptied */
}

