/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_glthread.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_glthread.c *
 *
 * NAME
 *  cq_ut_conqtour_glthread.c
 * PURPOSE
 *  Functions to be executed only from the main (openGL) thread,
 *  but normally requested from the python thread. These functions
 *  (service functions) are used to avoid calling openGL functions from the python thread,
 *  which results in undefined behaviour
 * USES
 *  
 * AUTHOR
 *  torralba
 * CREATION DATE
 *  Feb 3, 2011
 * MODIFICATION HISTORY
 *
 * *****/

#include "cq_ut_conqtour.h"
#include <stdarg.h>
#include <time.h>
#include <string.h>

#define DELAY 100

int srvPrintNumbers(comm_t *c);

/*** These functions are always called from the GLUT "thread" (main thread) ***/

void awakeMainThread(int v) /*  Executed as a result of a signal SIGUSR1 */
{
  ; /* Do nothing, just awake the thread */
}

/* This function must be called from the GL thread and if the window is hidden, */
/* which happens when there is no density, sets a GLUT timer that checks until */
/* a density is available */
void checkWindowVisibility(int has)
{
  /* Ignore the past value and use the actual control */
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      glutHideWindow();
      glutTimerFunc(DELAY, checkWindowVisibility, ctrl.has);
    }
  else
    glutShowWindow();
}

/* The function runCommandQueue is responsible for executing (mainly) openGL-containing service functions */
/* that are requested from the python thread. It must be called from a thread-safe and non-reentrant point. */
/* For this reason, it is used within display() */
void runCommandQueue(void)
{
//  static working = NO;
  cmdqueuenode_t *node;

  //  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  static unsigned int count = 0;
  ++count;
  //printf("**************Count: %d %d\n", count, pycmds.head);

  //glutPostRedisplay();
  //  if (working)
  //    {
  //      printf("Returning\n");
  //      return;
  //    }
  //  working = YES;
  while ((node = popCmdFromQueue(&pycmds)) != NULL)
    {
//      printf("--N %d\n", node);
      /* Execute the requested function */
      node->cmd(node->params);
      /* Free memory */
      deallocateCommunicator(node->params);
      free(node);
    }
//  pthread_mutex_unlock(&ctrl.mutex);
//  printf("** 1\n");
//  glutPostRedisplay();
//  printf("** 2\n");
//  working = NO;
}

void waitForRedisplay(ctrl_t *c)
{
  struct timespec ts;
  struct timeval tp;

  gettimeofday(&tp, NULL);

  /* Convert from timeval to timespec */
  ts.tv_sec  = tp.tv_sec;
  ts.tv_nsec = tp.tv_usec * 1000;
  ts.tv_sec += WAITING_TIME;


//  clock_gettime(CLOCK_REALTIME, &ts);
//  ts.tv_sec += WAITING_TIME;
  if (ctrl.redisplays > 0)
    {
//      printf("Waiting\n");
      pthread_cond_timedwait(&ctrl.cond, &ctrl.mutex, &ts);
//      printf("Pass\n");
    }
}

/**** Command queue-handling functions ****/
void initializeCommands(cmdqueue_t *q)
{
  q->head = NULL;
  q->tail = NULL;
}

/* Add a new command to the tail of the queue */
void pushCmdToQueue(cmdqueue_t *q, int(*c)(comm_t *), comm_t *p)
{
  cmdqueuenode_t *tmp;

  tmp = q->tail;
  q->tail = malloc(sizeof(cmdqueuenode_t));
  if (q->tail == NULL)
    {
      fprintf(stderr, "Error allocating node for command queue\n");
      exit(-1);
    }
  q->tail->next = NULL;
  if (q->head == NULL)
    {
      q->head = q->tail; /* First node */
      //      printf("First node: %d\n", q->head);
    }
  if (tmp != NULL)
    tmp->next = q->tail; /* Previous tail points to new tail */
  //writeQueue(q);
  /* Fill in the structure */
  q->tail->cmd = c; /* The command */
  q->tail->params = p; /* The parameters */
}

/* Get a command from the head of the queue */
/* This function returns a node but doesn't deallocate it. */
/* That is the responsibility of the caller */
cmdqueuenode_t *popCmdFromQueue(cmdqueue_t *q)
{
  cmdqueuenode_t *node; /* Node (Command + Parameters) to be returned */

  if (q->head == NULL)
    return NULL; /* Empty queue */
  node = q->head;
  q->head = q->head->next; /* Point to next node */
  if (q->head == NULL) /* Empty queue */
    q->tail = NULL;
  //  printf("Next: %d\n", q->head);

  return node; /* The caller must free !!! */
}

void writeQueue(cmdqueue_t *q)
{
  int i;
  cmdqueuenode_t *node;

  printf("H T %d %d\n", (int)q->head, (int)q->tail);
  node = q->head;
  if (node != NULL)
    while (node->next != NULL)
      {
        printf("-C %d\n", (int)node);
        node = q->head->next;
        printf("-N %d\n", (int)node);
      }
  printf("End\n");
}

/**** Communicator-handling functions: Used to create the parameters of a service function ****/
comm_t *allocateCommunicator(int num, ...)
{
  int i;
  comm_t *c; /* The communicator to be created */

  va_list args;
  va_start(args, num);

  /* Allocate space for the communicator */
  c = (comm_t *) malloc((size_t) sizeof(comm_t));
  c->num = num;
  if(num == 0)
    {
      va_end(args);
      return c;
    }
  c->p = (void **) malloc((size_t) num * sizeof(void *));
  c->length = (commlen_t *) malloc((size_t) num * sizeof(commlen_t));
  c->size = (size_t *) malloc((size_t) num * sizeof(size_t));
  c->dealloc = (char *) malloc((size_t) num * sizeof(char));
  if (c == NULL || c->p == NULL || c->length == NULL || c->size == NULL)
    {
      fprintf(stderr, "Error allocating communicator\n");
      exit(-1);
    }
  /* The user provides two parameters per array of data, or NULL for an actual pointer: */
  /* 1. Length of the array, i.e. number of elements (type: commlen_t); if this is NULL, an actual pointer is intended */
  /* 2. Size of each element of the array (type: size_t); not given for actual pointers */
  for (i = 0; i < num; ++i)
    {
      c->length[i] = va_arg(args, commlen_t);
      if ((void *) c->length[i] != NULL) /* An array */
        {
          c->size[i] = va_arg(args, size_t);
          /* Allocate space for the array */
          c->p[i] = (void *) malloc((size_t) c->length[i] * c->size[i]);
          //          printf("Allocated %d\n", c->length[i] * c->size[i]);
          if (c->p[i] == NULL)
            {
              fprintf(stderr, "Error allocating data for communicator\n");
              exit(-1);
            }
          c->dealloc[i] = TRUE; /* The space allocated for c->p[i] must be deallocated in the future */
        }
      else /* An actual pointer */
        {
//          printf("No,no\n");
          c->p[i] = NULL; /* The user will assign this pointer after allocation */
          c->size[i] = (size_t) NULL; /* Not an array; the p field is used as an actual pointer */
          c->dealloc[i] = FALSE;
        }
    }
  va_end(args);
  return c;
}

void deallocateCommunicator(comm_t *c)
{
  int i;

//  printf("Fuera commun\n");
  for (i = 0; i < c->num; ++i)
    {
      if(c->dealloc[i]) /* Deallocate only when necessary */
        {
//          printf("Fuera %d\n", i);
          free(c->p[i]);
        }
    }
  if(c->num > 0)
    {
      free(c->p);
      free(c->length);
      free(c->size);
      free(c->dealloc);
    }
  //  printf("Ahora\n");
  free(c);
  //  printf("Ya\n");
}

/**** Service functions ****/
int srvPrintNumbers(comm_t *c)
{
  int i, j;
  int *ip;
  double *dp;

  //  printf("Entra: %d\n", c->num);
  //  for (i = 0; i < c->num; ++i)
  //    printf("L %d\n", c->length[i]);
  //  return 1;
  for (i = 0; i < c->num; ++i)
    {
      printf("--%d\n", i);
      ip = c->p[i]; /* A cast from (void *) to (int *) is necessary */
      for (j = 0; j < c->length[i]; ++j)
        {
          if (i % 2 == 0)
            {
              ip = c->p[i]; /* A cast from (void *) to (int *) is necessary */
              printf("%d ", ip[j]);
            }
          else
            {
              dp = c->p[i]; /* A cast from (void *) to (int *) is necessary */
              printf("%f ", dp[j]);
            }
        }
      printf("\n");
    }
  return 0;
}

int srvResetIndex(comm_t *c)
{
  if(c->num != 5)
    {
      fprintf(stderr, "Error: Call to srvResetIndex requires 5 parameters; %d given\n", c->num);
      exit(-1);
    }
  resetIndex(c->p[0], c->p[1], c->p[2], c->p[3], c->p[4]);

  return 0;
}

int srvWritePpmImage(comm_t *c)
{
  if(c->num != 1)
    {
      fprintf(stderr, "Error: Call to srvWritePpmImage requires 1 parameter; %d given\n", c->num);
      exit(-1);
    }
  printf("%s\n", (char *)c->p[0]);
  writePpmImage(c->p[0]);

  return 0;
}

int srvWritePymolSlice(comm_t *c)
{
  if(c->num != 3)
    {
      fprintf(stderr, "Error: Call to srvWritePymolSlice requires 3 parameter; %d given\n", c->num);
      exit(-1);
    }
  printPymolScript(c->p[0], c->p[1], c->p[2]);

  return 0;
}

//int srvWriteCqCoordinates(comm_t *c)
//{
//  if(c->num != 4)
//    {
//      fprintf(stderr, "Error: Call to srvWriteCqCoordinates requires 4 parameters; %d given\n", c->num);
//      exit(-1);
//    }
//  writeCqCoordinates(c->p[0], c->p[1], c->p[2], c->p[3]);
//
//  return 0;
//}

//int srvWriteXplorDensity(comm_t *c)
//{
//  if(c->num != 3)
//    {
//      fprintf(stderr, "Error: Call to srvWriteXplorDensity requires 3 parameters; %d given\n", c->num);
//      exit(-1);
//    }
////  printf("%s\n", (char *)c->p[0]);
//  writeXplorDensity(c->p[0], c->p[1], c->p[2]);
//
//  return 0;
//}

int srvRemakeCell(comm_t *c)
{
  density_t *d;

  if(c->num != 3)
    {
      fprintf(stderr, "Error: Call to srvRemakeCell requires 3 parameters; %d given\n", c->num);
      exit(-1);
    }
  /* Delete the gl list for the previous box, if any, and recreate box */
  if (glIsList(boundingbox) == GL_TRUE)
    glDeleteLists(boundingbox, 1);
  prepareBoundingBox(&boundingbox);
  /* Delete the gl list for the previous axes and recreate*/
  if (glIsList(axes) == GL_TRUE)
    glDeleteLists(axes, 1);
  prepareAxes(&axes);
  /* Reset the legend*/
  d = c->p[1];
  memset(&(d->legend), 0, sizeof(legend_t));
  prepareLegend(&ctrl, c->p[1]); /* Passing density */
  resetIndex(NULL, NULL, c->p[0], c->p[1], c->p[2]); /* Passing index, density and box */

  return 0;
}

int srvDeleteMoleculeList(comm_t *c)
{
  if(c->num != 0)
    {
      fprintf(stderr, "Error: Call to srvDeleteMoleculeList doesn't take any parameters\n");
      exit(-1);
    }
  /* Delete the gl list for the previous box, if any, and recreate box */
  if (glIsList(molecule) == GL_TRUE)
    glDeleteLists(molecule, 1);

  return 0;
}

int srvCreateMoleculeList(comm_t *c)
{
  if(c->num != 0)
    {
      fprintf(stderr, "Error: Call to srvCreateMoleculeList doesn't take any parameters\n");
      exit(-1);
    }
  /* Delete the gl list for the previous box, if any, and recreate box */
  molecule = glGenLists(1); /* Prepare visuals */
  prepareMol(molecule, 1.0);

  return 0;
}

int srvCalculateDrawSlice(comm_t *c)
{
  if(c->num != 3)
    {
      fprintf(stderr, "Error: Call to srvCalculateDrawSlice requires 3 parameters; %d given\n", c->num);
      exit(-1);
    }
  calculateDrawSlice(c->p[0], c->p[1], c->p[2]);

  return 0;
}

int srvResetRedraw(comm_t *c)
{
  int i;
  index_t *pidx;
  ctrl_t *pc;

  if(c->num != 6)
    {
      fprintf(stderr, "Error: Call to srvResetRedraw requires 6 parameters; %d given\n", c->num);
      exit(-1);
    }
  resetIndex(c->p[0], c->p[1], c->p[2], c->p[3], c->p[4]);
  calculateDrawSlice(c->p[2], c->p[4], c->p[3]);

  pidx = c->p[2];
  pc = c->p[5];
  /* Set new up vector (control) */
  for (i = 0; i < 3; ++i)
    pc->up[i] = pidx->dindex[pidx->curr]->rectangle[3][i] - pidx->dindex[pidx->curr]->rectangle[0][i];
  normalize(pc->up);
  prepareMol(molecule, pc->spherescale);

  return 0;
}

int srvRecreateLegend(comm_t *c)
{
  int i;
  ctrl_t *pc;
  density_t *pd;
  index_t *pi;
  int oldcurr; /* Used to keep current slice number */

  if(c->num != 3)
    {
      fprintf(stderr, "Error: Call to srvRecreateLegend requires 3 parameters; %d given\n", c->num);
      exit(-1);
    }
  pc = c->p[0];
  pd = c->p[1];
  pi = c->p[2];

  prepareLegend(pc, pd);
  oldcurr = pi->curr;
  /* Redraw all the (previously created) slices */
  for(i = 0; i < pi->cuts; ++i)
    {
      pi->curr = i;
      if(pi->dindex[i] != NULL)
        drawSlice(pi, pd);
    }
  pi->curr = oldcurr;
  return 0;
}
