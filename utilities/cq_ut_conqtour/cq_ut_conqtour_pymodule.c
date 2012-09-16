/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour_pymodule.c
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour_pymodule.c *
 *
 * NAME
 *  cq_ut_conqtour_pymodule.c
 * PURPOSE
 *
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

static PyObject *ConqtourError;
static PyObject *ConqtourWarning;

/**** Lists to make legend labels more accesible from the python interpreter ****/

/* List of legend (collective) methods. */
/* Each entry refers to a whole set of ticks, interval and interpolation (but not color) methods */
static char *conqtourLegendCollectiveMethods[] =
  { "linear", "log", NULL };
static PyObject *conqtourListLegendCollectiveMethods;

/* List of legend tick methods */
static char *conqtourLegendTickMethods[] =
  { "linear", "log", NULL };
static PyObject *conqtourListLegendTickMethods;

/* List of legend interval methods */
static char *conqtourLegendIntervalMethods[] =
  { "linear", "log", NULL };
static PyObject *conqtourListLegendIntervalMethods;

/* List of legend interpolation methods */
static char *conqtourLegendInterpolationMethods[] =
  { "linear", "log", "constant", NULL };
static PyObject *conqtourListLegendInterpolationMethods;

/* List of legend color methods */
static char *conqtourLegendColorMethods[] =
  { "redblue", "rb", "redgreenblue", "rgb", "rainbow", NULL };
static PyObject *conqtourListLegendColorMethods;

/* Names of the dictionaries within the python interpreter */
#define conqtourListLegendCollectiveMethodsName "CQ_LEGEND_METHOD"
#define conqtourListLegendTickMethodsName "CQ_LEGEND_TICKS"
#define conqtourListLegendIntervalMethodsName "CQ_LEGEND_INTERVALS"
#define conqtourListLegendInterpolationMethodsName "CQ_LEGEND_INTERPOLATION"
#define conqtourListLegendColorMethodsName "CQ_LEGEND_COLORS"

/* NOTE: "Redisplays before acting" are a safety mechanism to guarantee the expected view */

#ifdef DEBUG
static PyObject *conqtour_0(PyObject *self, PyObject *args)
  {
    int i;
    comm_t *params; /* Parameters of the service function */
    int p0[3] =
      { 14, 2, 7};
    double p1[5] =
      { 0.1, 0.3, -0.67, 99.4, 0.0};
    int p2[2] =
      { 21, 2};
    double p3[4] =
      { 0.76, 0.2, 0.8, -1.5};
    int *ip;
    double *dp;
    cmdqueuenode_t *node;
    static count=0;

    ++count;
    printf("---------------Empieza %d\n", count);
    pthread_mutex_lock(&ctrl.mutex); /* Be safe */
    requestRedisplay(); /* Force redisplay before acting... */
    waitForRedisplay(&ctrl);

    params = allocateCommunicator(4, COMM_DATA(3,int), COMM_DATA(5,double), COMM_DATA(2,int), COMM_DATA(4,double));
    ip = params->p[0]; /* Cast of the void pointer is always necessary */
    for (i = 0; i < 3; ++i)
    ip[i] = p0[i];
    dp = params->p[1]; /* Cast of the void pointer is always necessary */
    for (i = 0; i < 5; ++i)
    dp[i] = p1[i];
    ip = params->p[2]; /* Cast of the void pointer is always necessary */
    for (i = 0; i < 2; ++i)
    ip[i] = p2[i];
    dp = params->p[3]; /* Cast of the void pointer is always necessary */
    for (i = 0; i < 4; ++i)
    dp[i] = p3[i];
    //  writeQueue(&pycmds);
    pushCmdToQueue(&pycmds, srvPrintNumbers, params);
    //  while ((node = popCmdFromQueue(&pycmds)) != NULL)
    //    {
    //      /* Execute the requested function */
    //      node->cmd(node->params);
    //      /* Free memory */
    //      deallocateCommunicator(node->params);
    //      printf("Node\n");
    //      free(node);
    //      writeQueue(&pycmds);
    //      printf("Node-end\n");
    //    }

    //  int (*fg)(comm_t *);
    //  fg = srvPrintNumbers;
    //  printf("Hey %d %d %d\n", pycmds.head->cmd, fg, pycmds.head->params->num);
    requestRedisplay();
    waitForRedisplay(&ctrl);
    printf("---------------Termina %d\n", count);

    pthread_mutex_unlock(&ctrl.mutex);
    Py_INCREF(Py_None);
    return Py_None;
  }
#endif

static PyObject *conqtour_get_current_slice(PyObject *self, PyObject *args)
{
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  return Py_BuildValue("i", idx.curr);
}

static PyObject *conqtour_get_number_atoms(PyObject *self, PyObject *args)
{
  if (!(ctrl.has & CQ_CTRL_HAS_COORDS))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires coordinates");
      return NULL;
    }
  return Py_BuildValue("i", coords.numatoms);
}

static PyObject *conqtour_get_number_slices(PyObject *self, PyObject *args)
{
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  return Py_BuildValue("i", idx.cuts);
}

static PyObject *conqtour_get_slice_plane(PyObject *self, PyObject *args)
{
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  return Py_BuildValue("[[ddd][ddd]]", idx.surfnormal[0], idx.surfnormal[1], idx.surfnormal[2], idx.dindex[idx.curr]->inplane[0],
      idx.dindex[idx.curr]->inplane[1], idx.dindex[idx.curr]->inplane[2]);
}

static PyObject *conqtour_info_find(PyObject *self, PyObject *args)
{
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  findConquestAvailability(ctrl.verbose);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_info_set_cores(PyObject *self, PyObject *args)
{
  int cores;

  if (!PyArg_ParseTuple(args, "i", &cores))
    return NULL;
  if (cores >= 0) /* Zero is used to reset the number of cores */
    {
      pthread_mutex_lock(&ctrl.mutex); /* Be safe */
      cqinfo.cores = cores;
      pthread_mutex_unlock(&ctrl.mutex);
    }
  else
    {
      PyErr_SetString(ConqtourError, "The number of CQ-run cores must be non-negative\n");
      return NULL;
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_info_set_fractional_coordinates(PyObject *self, PyObject *args)
{
  PyObject *bool;

  if (!PyArg_ParseTuple(args, "O", &bool))
    return NULL;
  if (PyObject_IsTrue(bool))
    {
      pthread_mutex_lock(&ctrl.mutex); /* Be safe */
      SET_FLAG(cqinfo.have, CQ_HAVE_FRACTIONAL, TRUE);
      SET_FLAG(cqinfo.read, CQ_READ_FRACTIONAL, TRUE);
      pthread_mutex_unlock(&ctrl.mutex);
    }
  else
    {
      pthread_mutex_lock(&ctrl.mutex); /* Be safe */
      SET_FLAG(cqinfo.have, CQ_HAVE_FRACTIONAL, TRUE);
      SET_FLAG(cqinfo.read, CQ_READ_FRACTIONAL, FALSE);
      pthread_mutex_unlock(&ctrl.mutex);
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_info_set_grid_points(PyObject *self, PyObject *args)
{
  int pt[3];

  if (!PyArg_ParseTuple(args, "(iii)", &pt[0], &pt[1], &pt[2]))
    return NULL;
  if (pt[0] >= 0 && pt[1] >= 0 && pt[2] >= 0) /* Zero in any component is used to reset */
    {
      pthread_mutex_lock(&ctrl.mutex); /* Be safe */
      cqinfo.pt[0] = pt[0];
      cqinfo.pt[1] = pt[1];
      cqinfo.pt[2] = pt[2];
      pthread_mutex_unlock(&ctrl.mutex);
    }
  else
    {
      PyErr_SetString(ConqtourError, "The number of CQ-run grid points must be non-negative\n");
      return NULL;
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_info_set_block_points(PyObject *self, PyObject *args)
{
  int pt[3];

  if (!PyArg_ParseTuple(args, "(iii)", &pt[0], &pt[1], &pt[2]))
    return NULL;
  if (pt[0] >= 0 && pt[1] >= 0 && pt[2] >= 0) /* Zero in any component is used to reset */
    {
      pthread_mutex_lock(&ctrl.mutex); /* Be safe */
      cqinfo.blockpt[0] = pt[0];
      cqinfo.blockpt[1] = pt[1];
      cqinfo.blockpt[2] = pt[2];
      pthread_mutex_unlock(&ctrl.mutex);
    }
  else
    {
      PyErr_SetString(ConqtourError, "The number of CQ-run points per block must be non-negative\n");
      return NULL;
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_info_set_cell_dimensions(PyObject *self, PyObject *args)
{
  double lvm[3]; /* Lattice vector modules */

  if (!PyArg_ParseTuple(args, "(ddd)", &lvm[0], &lvm[1], &lvm[2]))
    return NULL;
  if (lvm[0] > ZEROTOL && lvm[1] > ZEROTOL && lvm[2] > ZEROTOL)
    {
      pthread_mutex_lock(&ctrl.mutex); /* Be safe */
      cqinfo.lvm[0] = lvm[0];
      cqinfo.lvm[1] = lvm[1];
      cqinfo.lvm[2] = lvm[2];
      pthread_mutex_unlock(&ctrl.mutex);
    }
  else if ((lvm[0] > -ZEROTOL && lvm[0] <= ZEROTOL) || (lvm[1] > -ZEROTOL && lvm[1] <= ZEROTOL) || (lvm[2] > -ZEROTOL && lvm[2]
      <= ZEROTOL))
    {
      pthread_mutex_lock(&ctrl.mutex); /* Be safe */
      cqinfo.lvm[0] = 0.0;
      cqinfo.lvm[1] = 0.0;
      cqinfo.lvm[2] = 0.0;
      pthread_mutex_unlock(&ctrl.mutex);
    }
  else
    {
      PyErr_SetString(ConqtourError, "The cell dimensions in a CQ run must be non-negative\n");
      return NULL;
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_info_set_density_stub(PyObject *self, PyObject *args)
{
  char *file;
  comm_t *params;

  if (!PyArg_ParseTuple(args, "s", &file))
    return NULL;
  if (strlen(file) > MAXLIN - 1)
    {
      PyErr_SetString(ConqtourError, "The file name is too long\n");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  strcpy(cqinfo.densityfl, file);
  pthread_mutex_unlock(&ctrl.mutex);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_info_set_coordinates_file(PyObject *self, PyObject *args)
{
  char *file;
  comm_t *params;

  if (!PyArg_ParseTuple(args, "s", &file))
    return NULL;
  if (strlen(file) > MAXLIN - 1)
    {
      PyErr_SetString(ConqtourError, "The file name is too long\n");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  strcpy(cqinfo.coordinatesfl, file);
  /* Note: access returns 0 if the file can be opened with the requested permissions */
  SET_FLAG(cqinfo.have, CQ_HAVE_COORDINATES,!access(cqinfo.coordinatesfl, R_OK));
  pthread_mutex_unlock(&ctrl.mutex);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_info_set_input_file(PyObject *self, PyObject *args)
{
  char *file;
  comm_t *params;

  if (!PyArg_ParseTuple(args, "s", &file))
    return NULL;
  if (strlen(file) > MAXLIN - 1)
    {
      PyErr_SetString(ConqtourError, "The file name is too long\n");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  strcpy(cqinfo.inputfl, file);
  pthread_mutex_unlock(&ctrl.mutex);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_info_set_blocks_file(PyObject *self, PyObject *args)
{
  char *file;
  comm_t *params;

  if (!PyArg_ParseTuple(args, "s", &file))
    return NULL;
  if (strlen(file) > MAXLIN - 1)
    {
      PyErr_SetString(ConqtourError, "The file name is too long\n");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  strcpy(cqinfo.blocksfl, file);
  pthread_mutex_unlock(&ctrl.mutex);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_info_set_output_file(PyObject *self, PyObject *args)
{
  char *file;
  comm_t *params;

  if (!PyArg_ParseTuple(args, "s", &file))
    return NULL;
  if (strlen(file) > MAXLIN - 1)
    {
      PyErr_SetString(ConqtourError, "The file name is too long\n");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  strcpy(cqinfo.outputfl, file);
  pthread_mutex_unlock(&ctrl.mutex);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_limits_get(PyObject *self, PyObject *args)
{
  return Py_BuildValue("[[iii][iii][iii][ddd][ddd]]", limits.griddims[0], limits.griddims[1], limits.griddims[2],
      limits.gridshft[0], limits.gridshft[1], limits.gridshft[2], limits.strides[0], limits.strides[1], limits.strides[2],
      limits.realdims[0], limits.realdims[1], limits.realdims[2], limits.realshft[0], limits.realshft[1], limits.realshft[2]);
}

static PyObject *conqtour_limits_set(PyObject *self, PyObject *args)
{
  GLdouble d[3];
  int i, j;
  //  double conv;
  double maxdim;
  PyObject *list[5];
  int len[5];
  char error[MAXLIN];
  PyObject *datum;
  int tmpi[3][3];
  double tmpd[2][3];

  if (!PyArg_ParseTuple(args, "(OOOOO)", &list[0], &list[1], &list[2], &list[3], &list[4]))
    return NULL;
  /* Take note of our references */
  Py_INCREF(list[0]);
  Py_INCREF(list[1]);
  Py_INCREF(list[2]);
  Py_INCREF(list[3]);
  Py_INCREF(list[4]);
  /* Check that the lists have the right format */
  for (i = 0; i < 5; ++i)
    {
      len[i] = PyList_Size(list[i]);
      if (len[i] != 3)
        {
          Py_DECREF(list[0])
            ;
          Py_DECREF(list[1])
            ;
          Py_DECREF(list[2])
            ;
          Py_DECREF(list[3])
            ;
          Py_DECREF(list[4])
            ;
          sprintf(error, "Sublist %d is not of the required length (3)\n", i);
          PyErr_SetString(PyExc_TypeError, error);
          return NULL;
        }
      for (j = 0; j < 3; ++j)
        {
          datum = PyList_GetItem(list[i], j);
          Py_INCREF(datum); /*We are using it */
          if (i < 3) /* First three sublists are lists of integer */
            {
              if (!PyInt_Check(datum))
                {
                  Py_DECREF(list[0])
                    ;
                  Py_DECREF(list[1])
                    ;
                  Py_DECREF(list[2])
                    ;
                  Py_DECREF(list[3])
                    ;
                  Py_DECREF(list[4])
                    ;
                  Py_DECREF(datum)
                    ;
                  sprintf(error, "Element %d of sublist %d is not an integer\n", j, i);
                  PyErr_SetString(PyExc_TypeError, error);
                  return NULL;
                }
              tmpi[i][j] = (int) PyInt_AsLong(datum);
              //printf("%d\n", tmpi[i][j]);
              Py_DECREF(datum)
                ;
            }
          else /* The rest are lists of doubles */
            {
              if (!PyFloat_Check(datum) && !PyInt_Check(datum))
                {
                  Py_DECREF(list[0])
                    ;
                  Py_DECREF(list[1])
                    ;
                  Py_DECREF(list[2])
                    ;
                  Py_DECREF(list[3])
                    ;
                  Py_DECREF(list[4])
                    ;
                  Py_DECREF(datum)
                    ;
                  sprintf(error, "Element %d of sublist %d is not a float\n", j, i);
                  PyErr_SetString(PyExc_TypeError, error);
                  return NULL;
                }
              tmpd[i - 3][j] = PyFloat_AsDouble(datum);
              //printf("%f\n", tmpd[i-3][j]);
              Py_DECREF(datum)
                ;
            }
        }
    }
  Py_DECREF(list[0])
    ;
  Py_DECREF(list[1])
    ;
  Py_DECREF(list[2])
    ;
  Py_DECREF(list[3])
    ;
  Py_DECREF(list[4])
    ;

  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  /* Now that we have the data, assign */
  for (i = 0; i < 3; ++i)
    {
      limits.griddims[i] = tmpi[0][i];
      limits.gridshft[i] = tmpi[1][i];
      limits.strides[i] = tmpi[2][i];
      limits.realdims[i] = tmpd[0][i];
      limits.realshft[i] = tmpd[1][i];
    }
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_limits_set_grid_shift(PyObject *self, PyObject *args)
{
  GLint s[3];
  int i;
  //  double conv;

  if (!PyArg_ParseTuple(args, "(iii)", &s[0], &s[1], &s[2]))
    return NULL;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  for (i = 0; i < 3; ++i)
    {
      if (s[i] >= 0)
        limits.gridshft[i] = s[i]; /* Negative values keep the current shift*/
      if (cqinfo.pt[i]) /* We have some information about CQ grid */
        {
          if (limits.gridshft[i] >= cqinfo.pt[i] - limits.strides[i])
            {
              if (ctrl.verbose)
                printf("WARNING: User-provided grid shift is not valid for coord. %s; resetting\n", REPORT_COORDINATE(i));
              limits.gridshft[i] = 0;
            }
          if (limits.gridshft[i] + limits.griddims[i] * limits.strides[i] > cqinfo.pt[i])
            {
              if (ctrl.verbose)
                printf("WARNING: User-provided grid window is too large (coord. %s); setting to maximum\n", REPORT_COORDINATE(i));
              limits.griddims[i] = (int) ceil((cqinfo.pt[i] - limits.gridshft[i]) / (double) limits.strides[i]);
            }
          //          /* The real dimensions must much the grid dimensions; correct if possible */
          //          if (FLAG_UP(cqinfo.units, CQ_UNITS_BOHR))
          //            conv = BOHR_TO_ANGSTROM;
          //          else
          //            conv = 1.0;
          if (cqinfo.lvm[i] > ZEROTOL)
            limits.realdims[i] = (limits.griddims[i] - 1) * limits.strides[i] * cqinfo.lvm[i] / (cqinfo.pt[i] - 1);
          //            limits.realdims[i] = conv * (limits.griddims[i] - 1) * limits.strides[i] * cqinfo.lvm[i] / (cqinfo.pt[i] - 1);
        }
    }
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_limits_set_grid_stride(PyObject *self, PyObject *args)
{
  GLint s[3];
  int i;
  //  double conv;

  if (!PyArg_ParseTuple(args, "(iii)", &s[0], &s[1], &s[2]))
    return NULL;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  for (i = 0; i < 3; ++i)
    {
      if (s[i] >= 1)
        limits.strides[i] = s[i];
      if (cqinfo.pt[i]) /* We have some information about CQ grid */
        {
          if (limits.strides[i] > cqinfo.pt[i] - 1) /* We need at least two points when we load */
            {
              if (ctrl.verbose)
                printf("WARNING: The stride for coord. %s is too large; setting to maximum\n", REPORT_COORDINATE(i));
              limits.strides[i] = cqinfo.pt[i] - 1;
            }
          if (limits.gridshft[i] > cqinfo.pt[i] - limits.strides[i] - 1)
            {
              if (ctrl.verbose)
                printf("WARNING: User-provided grid shift is not valid for coord. %s; setting to maximum\n", REPORT_COORDINATE(i));
              limits.gridshft[i] = cqinfo.pt[i] - limits.strides[i] - 1;
            }
          if (limits.griddims[i] > 1 + (cqinfo.pt[i] - limits.gridshft[i]) / limits.strides[i])
            {
              if (ctrl.verbose)
                printf("WARNING: User-provided grid memory points is too large for coord. %s; setting to maximum\n",
                    REPORT_COORDINATE(i));
              limits.griddims[i] = (int) ceil((cqinfo.pt[i] - limits.gridshft[i]) / (double) limits.strides[i]);
            }
          /* The real dimensions must much the grid dimensions; correct if possible */
          //          if (FLAG_UP(cqinfo.units, CQ_UNITS_BOHR))
          //            conv = BOHR_TO_ANGSTROM;
          //          else
          //            conv = 1.0;
          if (cqinfo.lvm[i] > ZEROTOL)
            limits.realdims[i] = (limits.griddims[i] - 1) * limits.strides[i] * cqinfo.lvm[i] / (cqinfo.pt[i] - 1);
          //            limits.realdims[i] = conv * (limits.griddims[i] - 1) * limits.strides[i] * cqinfo.lvm[i] / (cqinfo.pt[i] - 1);
        }
    }
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_limits_set_grid_window(PyObject *self, PyObject *args)
{
  GLint p[3];
  int i;
  //  double conv;

  if (!PyArg_ParseTuple(args, "(iii)", &p[0], &p[1], &p[2]))
    return NULL;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  for (i = 0; i < 3; ++i)
    {
      /*
       if(p[i]==1)        p[i]=2;
       else if(p[i] < 1)  p[i]=limits.griddims[i];
       limits.griddims[i] = p[i];
       */
      if (p[i] == 0 || p[i] == 1) /* Reset the value; negative values are ignored and hence keep the current value */
        {
          limits.griddims[i] = 0; /* This is a special case: the user is asking for a reset */
          limits.realdims[i] = 0.0;
        }
      if (p[i] > 1) /* Hence, p[i] == 1 keeps the previous value */
        limits.griddims[i] = p[i]; /* Must be at least 2; otherwise, keep old value */
      if (cqinfo.pt[i]) /* We have some information about CQ grid */
        {
          if (limits.gridshft[i] + limits.griddims[i] * limits.strides[i] > cqinfo.pt[i])
            {
              if (ctrl.verbose)
                printf("WARNING: User-provided grid window is not valid for coord. %s; setting to maximum\n",
                    REPORT_COORDINATE(i));
              limits.griddims[i] = (int) ceil((cqinfo.pt[i] - limits.gridshft[i]) / (double) limits.strides[i]);
            }
          //          /* The real dimensions must much the grid dimensions; correct if possible */
          //          if (FLAG_UP(cqinfo.units, CQ_UNITS_BOHR))
          //            conv = BOHR_TO_ANGSTROM;
          //          else
          //            conv = 1.0;
          if (cqinfo.lvm[i] > ZEROTOL)
            {
              if (limits.griddims[i] > 1) /* Ignore if the call is a reset */
                limits.realdims[i] = (limits.griddims[i] - 1) * limits.strides[i] * cqinfo.lvm[i] / (cqinfo.pt[i] - 1);
            }
          //            limits.realdims[i] = conv * (limits.griddims[i] - 1) * limits.strides[i] * cqinfo.lvm[i] / (cqinfo.pt[i] - 1);
        }
    }
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_limits_set_real_shift(PyObject *self, PyObject *args)
{
  GLdouble s[3];
  int i;

  if (!PyArg_ParseTuple(args, "(ddd)", &s[0], &s[1], &s[2]))
    return NULL;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  for (i = 0; i < 3; ++i)
    {
      limits.realshft[i] = s[i]; /* Any shift is valid; this is used to shift the output */
    }
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_limits_set_real_shift_from_grid(PyObject *self, PyObject *args)
{
  int i;

  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  for (i = 0; i < 3; ++i)
    {
      if (cqinfo.pt[i] && cqinfo.lvm[i]) /* We have some grid information */
        {
          limits.realshft[i] = limits.gridshft[i] * cqinfo.lvm[i] / (cqinfo.pt[i] - 1);
        }
      else
        limits.realshft[i] = 0.0; /* No info:  reset */
    }
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_limits_set_grid_shift_from_real(PyObject *self, PyObject *args)
{
  int i;

  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  for (i = 0; i < 3; ++i)
    {
      if (cqinfo.pt[i] && cqinfo.lvm[i]) /* We have some grid information */
        {
          limits.gridshft[i] = limits.realshft[i] / (cqinfo.lvm[i] / (cqinfo.pt[i] - 1));
          /* Make the real shift match exactly the grid shift; this is the only function that enforces this */
          limits.realshft[i] = limits.gridshft[i] * cqinfo.lvm[i] / (cqinfo.pt[i] - 1);
          /* Check that the window is not too large */
          if (limits.gridshft[i] + limits.griddims[i] * limits.strides[i] > cqinfo.pt[i])
            {
              if (ctrl.verbose)
                printf("WARNING: User-provided grid window is not valid for coord. %s; setting to maximum\n",
                    REPORT_COORDINATE(i));
              limits.griddims[i] = (int) ceil((cqinfo.pt[i] - limits.gridshft[i]) / (double) limits.strides[i]);
              if (cqinfo.lvm[i] > ZEROTOL)
                {
                  if (limits.griddims[i] > 1) /* Ignore if the call is a reset */
                    limits.realdims[i] = (limits.griddims[i] - 1) * limits.strides[i] * cqinfo.lvm[i] / (cqinfo.pt[i] - 1);
                }
            }
        }
      else
        limits.gridshft[i] = 0; /* No info:  reset */
    }
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_limits_set_real_window(PyObject *self, PyObject *args)
{
  GLdouble d[3];
  int i;
  //  double conv;
  double maxdim;

  if (!PyArg_ParseTuple(args, "(ddd)", &d[0], &d[1], &d[2]))
    return NULL;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  for (i = 0; i < 3; ++i)
    {
      if (d[i] > -ZEROTOL && d[i] < ZEROTOL) /* This is a reset of a window component (both real and grid) */
        {
          limits.realdims[i] = 0.0; /* For the moment, use the requested value, or keep the old one if neg.  */
          limits.griddims[i] = 0; /* The grid window (component) is also reset, to match the real window */
        }
      if (d[i] > ZEROTOL)
        limits.realdims[i] = d[i]; /* For the moment, use the requested value, or keep the old one if neg.  */
      if (cqinfo.pt[i] && cqinfo.lvm[i]) /* CQ grid spacing is known */
        {
          //          if (FLAG_UP(cqinfo.units, CQ_UNITS_BOHR))
          //            conv = BOHR_TO_ANGSTROM;
          //          else
          //            conv = 1.0;
          /* The maximum allowed window depends on the grid shift, not the real shift */
          //          maxdim = conv * cqinfo.lvm[i] * (1.0 - limits.gridshft[i] / (cqinfo.pt[i] - 1));
          maxdim = cqinfo.lvm[i] * (1.0 - limits.gridshft[i] / (cqinfo.pt[i] - 1));
          if (limits.realdims[i] > maxdim + ZEROTOL2)
            {
              if (ctrl.verbose)
                printf("WARNING: Real window is too large for coordinate %s; setting to maximum\n", REPORT_COORDINATE(i));
              limits.realdims[i] = maxdim;
              printf("M %f\n", maxdim);
            }
          /* The grid and real windows must agree */
          /* First, calculate the grid window, based on the current real one */
          //          limits.griddims[i] = 1 + (int) floor((limits.realdims[i] + ZEROTOL2) * (cqinfo.pt[i] - 1) / (conv * cqinfo.lvm[i]
          //              * limits.strides[i]));
          if (limits.realdims[i] > ZEROTOL)
            {
              limits.griddims[i] = 1 + (int) floor((limits.realdims[i] + ZEROTOL2) * (cqinfo.pt[i] - 1) / (cqinfo.lvm[i]
                  * limits.strides[i]));
              /* Finally, adjust the real window to match the grid */
              //          limits.realdims[i] = conv * (limits.griddims[i] - 1) * limits.strides[i] * cqinfo.lvm[i] / (cqinfo.pt[i] - 1);
              limits.realdims[i] = (limits.griddims[i] - 1) * limits.strides[i] * cqinfo.lvm[i] / (cqinfo.pt[i] - 1);
            }
          else
            limits.griddims[i] = 0;
        }
    }
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_view_recenter(PyObject *self, PyObject *args)
{
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex);
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  recenter(&ctrl);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);

  return Py_BuildValue("[ddd]", ctrl.center[0], ctrl.center[1], ctrl.center[2]);
}

static PyObject *conqtour_view_set_center(PyObject *self, PyObject *args)
{
  GLdouble c[3];
  int i;

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "(ddd)", &c[0], &c[1], &c[2]))
    return NULL;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  for (i = 0; i < 3; ++i)
    ctrl.center[i] = c[i];
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_view_rotate(PyObject *self, PyObject *args)
{
  GLdouble angle;

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "d", &angle))
    return NULL;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  rotateVector(angle, ctrl.up, ctrl.normal);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_view_rotate_reslice(PyObject *self, PyObject *args)
{
  GLdouble angle;
  int i;
  double inplane[3];
  comm_t *params;
  double *dp;
  //  index_t *ip;
  //  densitycut_t *dcp;
  //  box_t *bp;

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "d", &angle))
    return NULL;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  waitForRedisplay(&ctrl);

  rotateVector(angle, ctrl.up, ctrl.normal);
  for (i = 0; i < 3; ++i)
    inplane[i] = idx.dindex[idx.curr]->inplane[i];
  params = allocateCommunicator(5, COMM_DATA(3,double), COMM_DATA(3,double), COMM_PTR, COMM_PTR, COMM_PTR);
  dp = params->p[0]; /* Cast of the void pointer to double */
  for (i = 0; i < 3; ++i) /* First parameter : Normal to plane */
    dp[i] = ctrl.normal[i];
  dp = params->p[1]; /* Cast of the void pointer */
  for (i = 0; i < 3; ++i) /* Second parameter : In-plane point */
    dp[i] = inplane[i];
  params->p[2] = &idx; /* Third parameter : Pointer to index */
  params->p[3] = &density; /* Fourth parameter : Pointer to density cut */
  params->p[4] = &box; /* Fith parameter : Pointer to cell box */
  /* The next line is equivalent to calling resetIndex(ctrl.normal, inplane, &idx, &density, &box); */
  pushCmdToQueue(&pycmds, srvResetIndex, params);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_write_ppm_image(PyObject *self, PyObject *args)
{
  char *file;
  comm_t *params;

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "s", &file))
    return NULL;
  if (strlen(file) > MAXLIN - 1)
    {
      PyErr_SetString(ConqtourError, "The file name is too long\n");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  params = allocateCommunicator(1, COMM_DATA(MAXLIN,char));
  strcpy(params->p[0], file); /* First parameter : File name */
  pushCmdToQueue(&pycmds, srvWritePpmImage, params); /* Equivalent to calling writePpmImage(file); */
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_write_pymol_slice(PyObject *self, PyObject *args)
{
  char *file;
  comm_t *params;

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "s", &file))
    return NULL;
  if (strlen(file) > MAXLIN - 1)
    {
      PyErr_SetString(ConqtourError, "The file name is too long\n");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  params = allocateCommunicator(3, COMM_PTR, COMM_PTR, COMM_DATA(MAXLIN,char));
  params->p[0] = &idx; /* First parameter : Pointer to index */
  params->p[1] = &density; /* Second parameter : Pointer to density */
  strcpy(params->p[2], file); /* Third parameter : File name */
  pushCmdToQueue(&pycmds, srvWritePymolSlice, params);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_write_cq_coordinates(PyObject *self, PyObject *args)
{
  char *file;
  int cqerrno;
  char cqerr[MAXLIN], reason[MAXLIN - 100];

  if (!(ctrl.has & CQ_CTRL_HAS_COORDS))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires coordinates");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "s", &file))
    return NULL;
  if (strlen(file) > MAXLIN - 1)
    {
      PyErr_SetString(ConqtourError, "The file name is too long\n");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    cqerrno = writeCqCoordinates(&cqinfo, NULL, &coords, file);
  else
    cqerrno = writeCqCoordinates(&cqinfo, &density, &coords, file);
  if (cqerrno)
    {
      switch (cqerrno)
        {
      case 1:
        strcpy(reason, "Cannot create file");
        break;
      case 2:
        strcpy(reason, "No atoms to write within window");
        break;
      default:
        strcpy(reason, "Unidentified error");
        }
      sprintf(cqerr, "Problem writing coordinates: %s", reason);
      PyErr_SetString(PyExc_IOError, cqerr);
      pthread_mutex_unlock(&ctrl.mutex); /* Don't forget to unlock !!! */
      return NULL;
    }
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_write_xyz_coordinates(PyObject *self, PyObject *args)
{
  char *file;
  int cqerrno;
  char cqerr[MAXLIN], reason[MAXLIN - 100];

  if (!(ctrl.has & CQ_CTRL_HAS_COORDS))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires coordinates");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "s", &file))
    return NULL;
  if (strlen(file) > MAXLIN - 1)
    {
      PyErr_SetString(ConqtourError, "The file name is too long\n");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    cqerrno = writeXyzCoordinates(&cqinfo, NULL, &coords, file);
  else
    cqerrno = writeXyzCoordinates(&cqinfo, &density, &coords, file);
  if (cqerrno)
    {
      switch (cqerrno)
        {
      case 1:
        strcpy(reason, "Cannot create file");
        break;
      case 2:
        strcpy(reason, "No atoms to write within window");
        break;
      default:
        strcpy(reason, "Unidentified error");
        }
      sprintf(cqerr, "Problem writing coordinates: %s", reason);
      PyErr_SetString(PyExc_IOError, cqerr);
      pthread_mutex_unlock(&ctrl.mutex); /* Don't forget to unlock !!! */
      return NULL;
    }
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_write_xplor_density(PyObject *self, PyObject *args)
{
  char *file;

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "s", &file))
    return NULL;
  if (strlen(file) > MAXLIN - 1)
    {
      PyErr_SetString(ConqtourError, "The file name is too long\n");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  writeXplorDensity(file, &density, &limits);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_load_xplor_density(PyObject *self, PyObject *args)
{
  char *file;
  comm_t *params;
  char reason[MAXLIN];
  int cqerrno;

  if (!PyArg_ParseTuple(args, "s", &file))
    return NULL;
  if (strlen(file) > MAXLIN - 1)
    {
      PyErr_SetString(ConqtourError, "The file name is too long\n");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  if (cqerrno = readXplorDensity(file, &density)) /* If there was an error */
    {
      switch (cqerrno)
        {
      case 1:
        strcpy(reason, "File not found");
        break;
      case 2:
        strcpy(reason, "Not the right format?");
        break;
      case 3:
        strcpy(reason, "Incomplete file");
        ctrl.has &= ~CQ_CTRL_HAS_DENSITY; /* At this point, previous density has been destroyed */
        requestRedisplay();
        waitForRedisplay(&ctrl);
        break;
      default:
        strcpy(reason, "Unidentified error");
        }
      PyErr_SetString(PyExc_IOError, reason);
      pthread_mutex_unlock(&ctrl.mutex); /* Don't forget to unlock !!! */
      return NULL;
    }
  ctrl.display |= CQ_CTRL_SHOW_SHAPE;
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  resetView(&ctrl);
  setBox(&density, &box);
  params = allocateCommunicator(3, COMM_PTR, COMM_PTR, COMM_PTR);
  params->p[0] = &idx; /* Third parameter : Pointer to index */
  params->p[1] = &density; /* Fourth parameter : Pointer to density cut */
  params->p[2] = &box; /* Fith parameter : Pointer to cell box */
  pushCmdToQueue(&pycmds, srvRemakeCell, params);
  ctrl.has |= CQ_CTRL_HAS_DENSITY;
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_load_cq_density(PyObject *self, PyObject *args)
{
  comm_t *params;

  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  if (readCqDensity(&density)) /* If there was an error */
    {
      PyErr_SetString(PyExc_IOError, "The density could not be loaded using CQ info");
      pthread_mutex_unlock(&ctrl.mutex); /* Don't forget to unlock !!! */
      return NULL;
    }
  ctrl.display |= CQ_CTRL_SHOW_SHAPE;
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  resetView(&ctrl);
  setBox(&density, &box);
  ctrl.has |= CQ_CTRL_HAS_DENSITY;
  params = allocateCommunicator(3, COMM_PTR, COMM_PTR, COMM_PTR);
  params->p[0] = &idx; /* Third parameter : Pointer to index */
  params->p[1] = &density; /* Fourth parameter : Pointer to density cut */
  params->p[2] = &box; /* Fith parameter : Pointer to cell box */
  pushCmdToQueue(&pycmds, srvRemakeCell, params);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_load_cq_coordinates(PyObject *self, PyObject *args)
{
  char cqerr[MAXLIN], reason[MAXLIN - 100];
  int cqerrno;
  comm_t *params, *params2;

  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  /* We can read directly (without using the command queue) because we don't need OpenGL (yet) */
  if (cqerrno = readCqCoordinates(&cqinfo, &coords)) /* If there was an error */
    {
      switch (cqerrno)
        {
      case 1:
        strcpy(reason, "File not found");
        break;
      case 2:
        strcpy(reason, "Not a Conquest coordinates file/Incomplete information about it?\n"
          "(E.g. Unknown list of atom types)");
        break;
      case 3:
        strcpy(reason, "Incomplete data");
        ctrl.has &= ~CQ_CTRL_HAS_COORDS; /* At this point, previous coords. have been destroyed */
        ctrl.has &= ~CQ_CTRL_HAS_COORDS_XYZ;
        ctrl.display &= ~CQ_CTRL_SHOW_MOLECULE;
        params = allocateCommunicator(0);
        pushCmdToQueue(&pycmds, srvDeleteMoleculeList, params);
        //        if (glIsList(molecule) == GL_TRUE)
        //          glDeleteLists(molecule, 1);
        requestRedisplay();
        waitForRedisplay(&ctrl);
        break;
      case 4:
        strcpy(reason, "Unspecified Conquest coordinates file");
        break;
      default:
        strcpy(reason, "Unidentified error");
        }
      sprintf(cqerr, "Problem with coordinates file: %s", reason);
      PyErr_SetString(PyExc_IOError, cqerr);
      pthread_mutex_unlock(&ctrl.mutex); /* Don't forget to unlock !!! */
      return NULL;
    }
  ctrl.has |= CQ_CTRL_HAS_COORDS; /* We have coordinates !! */
  ctrl.has |= CQ_CTRL_HAS_COORDS_XYZ; /* This is a Conquest file, but the information is stored in an XYZ structure */
  ctrl.display |= CQ_CTRL_SHOW_MOLECULE;
  /* If there was a GL list, destroy */
  params = allocateCommunicator(0);
  pushCmdToQueue(&pycmds, srvDeleteMoleculeList, params);
  //  if (glIsList(molecule) == GL_TRUE)
  //    glDeleteLists(molecule, 1);
  params2 = allocateCommunicator(0);
  pushCmdToQueue(&pycmds, srvCreateMoleculeList, params2);
  //  molecule = glGenLists(1); /* Prepare visuals */
  //  prepareMol(molecule, 1.0);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_load_xyz_coordinates(PyObject *self, PyObject *args)
{
  char *file;
  char cqerr[MAXLIN], reason[MAXLIN - 100];
  int cqerrno;
  comm_t *params, *params2;

  if (!PyArg_ParseTuple(args, "s", &file))
    return NULL;
  if (strlen(file) > MAXLIN - 1)
    {
      PyErr_SetString(ConqtourError, "The file name is too long\n");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  if (cqerrno = readXyzCoordinates(file, &coords)) /* If there was an error */
    {
      switch (cqerrno)
        {
      case 1:
        strcpy(reason, "File not found");
        break;
      case 2:
        strcpy(reason, "Not an XYZ file?");
        break;
      case 3:
        strcpy(reason, "Incomplete file");
        ctrl.has &= ~CQ_CTRL_HAS_COORDS; /* At this point, previous coords. have been destroyed */
        ctrl.has &= ~CQ_CTRL_HAS_COORDS_XYZ;
        ctrl.display &= ~CQ_CTRL_SHOW_MOLECULE;
        params = allocateCommunicator(0);
        pushCmdToQueue(&pycmds, srvDeleteMoleculeList, params);
        //        if (glIsList(molecule) == GL_TRUE)
        //          glDeleteLists(molecule, 1);
        requestRedisplay();
        waitForRedisplay(&ctrl);
        break;
      default:
        strcpy(reason, "Unidentified error");
        }
      sprintf(cqerr, "Problem with coordinates file: %s", reason);
      PyErr_SetString(PyExc_IOError, cqerr);
      pthread_mutex_unlock(&ctrl.mutex); /* Don't forget to unlock !!! */
      return NULL;
    }
  ctrl.has |= CQ_CTRL_HAS_COORDS; /* We have coordinates !! */
  ctrl.has |= CQ_CTRL_HAS_COORDS_XYZ;
  ctrl.display |= CQ_CTRL_SHOW_MOLECULE;
  /* If there was a GL list, destroy */
  params = allocateCommunicator(0);
  pushCmdToQueue(&pycmds, srvDeleteMoleculeList, params);
  //  if (glIsList(molecule) == GL_TRUE)
  //    glDeleteLists(molecule, 1);
  params2 = allocateCommunicator(0);
  pushCmdToQueue(&pycmds, srvCreateMoleculeList, params2);
  //  molecule = glGenLists(1); /* Prepare visuals */
  //  prepareMol(molecule, 1.0);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_view_get(PyObject *self, PyObject *args)
{
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  return Py_BuildValue("[[ddd][ddd][ddd][dddd]]", ctrl.center[0], ctrl.center[1], ctrl.center[2], ctrl.up[0], ctrl.up[1],
      ctrl.up[2], ctrl.normal[0], ctrl.normal[1], ctrl.normal[2], ctrl.distance, ctrl.field, ctrl.near, ctrl.far);
}

static PyObject *conqtour_view_set(PyObject *self, PyObject *args)
{
  int i, j, k;
  double dummy;
  PyObject *list[4];
  char error[MAXLIN];
  int explen[4] =
    { 3, 3, 3, 4 }; /* Expected lengths of each input sublist */
  int len; /* Actual length of the input sublist */
  PyObject *datum;
  double tmpd[4][4]; /* Four lists with a maximum of 4 elements */

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "(OOOO)", &list[0], &list[1], &list[2], &list[3]))
    return NULL;
  /* Take note of our references */
  for (i = 0; i < 4; ++i)
    Py_INCREF(list[i]);
  /* Check that the lists have the right format */
  for (i = 0; i < 4; ++i)
    {
      len = PyList_Size(list[i]);
      if (len != explen[i])
        {
          for (j = 0; j < 4; ++j)
            Py_DECREF(list[j])
              ;
          sprintf(error, "Sublist %d is not of the required length (%d)\n", i, explen[i]);
          PyErr_SetString(PyExc_TypeError, error);
          return NULL;
        }
      for (j = 0; j < explen[i]; ++j)
        {
          datum = PyList_GetItem(list[i], j);
          Py_INCREF(datum); /*We are using it */
          if (!PyFloat_Check(datum) && !PyInt_Check(datum)) /* Must be float type (integer is ok) */
            {
              for (k = 0; k < 4; ++k)
                Py_DECREF(list[k])
                  ;
              Py_DECREF(datum)
                ;
              sprintf(error, "Element %d of sublist %d is not a float\n", j, i);
              PyErr_SetString(PyExc_TypeError, error);
              return NULL;
            }
          tmpd[i][j] = PyFloat_AsDouble(datum);
          //printf("%d %d %f\n", i, j, tmpd[i][j]);
          Py_DECREF(datum)
            ;
        }
    }
  for (i = 0; i < 4; ++i)
    Py_DECREF(list[i])
      ;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  /* Assign values, now they have been checked */
  for (i = 0; i < 3; ++i)
    {
      ctrl.center[i] = tmpd[0][i];
      ctrl.up[i] = tmpd[1][i];
      ctrl.normal[i] = tmpd[2][i];
    }
  ctrl.distance = tmpd[3][0];
  ctrl.field = tmpd[3][1];
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  setView(&ctrl);
  ctrl.has |= CQ_CTRL_HAS_DENSITY;
  requestRedisplay();
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_view_set_reference_plane(PyObject *self, PyObject *args)
{
  int i, j, k;
  PyObject *list[2];
  char error[MAXLIN];
  int explen[2] =
    { 3, 3 }; /* Expected lengths of each input sublist */
  int len; /* Actual length of the input sublist */
  PyObject *datum;
  double tmpd[2][3]; /* Temporary storage, to be used while we check the input */
  comm_t *params;
  double *dp;

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "(OO)", &list[0], &list[1]))
    return NULL;
  /* Take note of our references */
  for (i = 0; i < 2; ++i)
    Py_INCREF(list[i]);
  /* Check that the lists have the right format */
  for (i = 0; i < 2; ++i)
    {
      len = PyList_Size(list[i]);
      if (len != explen[i])
        {
          for (j = 0; j < 2; ++j)
            Py_DECREF(list[j])
              ;
          sprintf(error, "Sublist %d is not of the required length (%d)\n", i, explen[i]);
          PyErr_SetString(PyExc_TypeError, error);
          return NULL;
        }
      for (j = 0; j < explen[i]; ++j)
        {
          datum = PyList_GetItem(list[i], j);
          Py_INCREF(datum); /*We are using it */
          if (!PyFloat_Check(datum) && !PyInt_Check(datum)) /* Must be float type (integer is ok) */
            {
              for (k = 0; k < 2; ++k)
                Py_DECREF(list[k])
                  ;
              Py_DECREF(datum)
                ;
              sprintf(error, "Element %d of sublist %d is not a float\n", j, i);
              PyErr_SetString(PyExc_TypeError, error);
              return NULL;
            }
          tmpd[i][j] = PyFloat_AsDouble(datum);
          //printf("%d %d %f\n", i, j, tmpd[i][j]);
          Py_DECREF(datum)
            ;
        }
    }
  for (i = 0; i < 2; ++i)
    Py_DECREF(list[i])
      ;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);

  /* Prepare the normal for the new index */
  for (i = 0; i < 3; ++i)
    ctrl.normal[i] = tmpd[0][i];
  params = allocateCommunicator(6, COMM_DATA(3,double), COMM_DATA(3,double), COMM_PTR, COMM_PTR, COMM_PTR, COMM_PTR);
  dp = params->p[0]; /* Cast of the void pointer to double */
  for (i = 0; i < 3; ++i) /* First parameter : Normal to plane */
    dp[i] = ctrl.normal[i];
  dp = params->p[1]; /* Cast of the void pointer */
  for (i = 0; i < 3; ++i) /* Second parameter : In-plane point */
    dp[i] = tmpd[1][i];
  params->p[2] = &idx; /* Third parameter : Pointer to index */
  params->p[3] = &density; /* Fourth parameter : Pointer to density cut */
  params->p[4] = &box; /* Fifth parameter : Pointer to cell box */
  params->p[5] = &ctrl; /* Sixth parameter : Pointer to control */
  pushCmdToQueue(&pycmds, srvResetRedraw, params);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_view_set_slice(PyObject *self, PyObject *args)
{
  GLuint plane;
  comm_t *params;

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "i", &plane))
    return NULL;
  if (plane < 0 || plane > idx.cuts - 1) /* Outside of the allowed range */
    {
      PyErr_SetString(ConqtourError, "Requested slice out of bounds");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  idx.curr = plane;
  params = allocateCommunicator(3, COMM_PTR, COMM_PTR, COMM_PTR);
  params->p[0] = &idx; /* Third parameter : Pointer to index */
  params->p[1] = &box; /* Fith parameter : Pointer to cell box */
  params->p[2] = &density; /* Fourth parameter : Pointer to density cut */
  pushCmdToQueue(&pycmds, srvCalculateDrawSlice, params); /*Equivalent to calling calculateDrawSlice(&idx, &box, &density); */
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_legend_set_label_format(PyObject *self, PyObject *args)
{
  char *format;
  char errmsg[MAXLIN];

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "s", &format))
    return NULL;
  if (strlen(format) > MAXLIN - 1)
    {
      PyErr_SetString(ConqtourError, "The format string is too long\n");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);

  /* Force the scale ticks to be recalculated using default values */
  strcpy(ctrl.legend.format, format);
  strcpy(density.legend.format, format);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_legend_get_label_format(PyObject *self, PyObject *args)
{
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  return Py_BuildValue("s", ctrl.legend.format);
}

/* This function is for collective (i.e. ticks + interval + interpolation) methods */
static PyObject *conqtour_legend_set_method(PyObject *self, PyObject *args)
{
  int i;
  char *method;
  comm_t *params;

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "s", &method))
    return NULL;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);

  /* Force the scale ticks to be recalculated using default values */
  ctrl.legend.tickinterval = 0.0;
  ctrl.legend.tickorigin = 0.0;
  if (!strcmp(method, "linear"))
    {
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_TICKS_LINEAR, CQ_LEGEND_TICKS_MASK);
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_INTERVAL_LINEAR, CQ_LEGEND_INTERVAL_MASK);
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_INTERPOL_LINEAR, CQ_LEGEND_INTERPOL_MASK);
    }
  else if (!strcmp(method, "log"))
    {
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_TICKS_SIGNABSLOG, CQ_LEGEND_TICKS_MASK);
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_INTERVAL_SIGNABSLOG, CQ_LEGEND_INTERVAL_MASK);
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_INTERPOL_SIGNABSLOG, CQ_LEGEND_INTERPOL_MASK);
    }
  else
    {
      if (ctrl.verbose)
        printf("WARNING: Method identifier, not in the dictionary. Using defaults\n");
      /* Use default */
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_DFLT_METHOD_TICKS, CQ_LEGEND_TICKS_MASK);
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_DFLT_METHOD_INTERVAL, CQ_LEGEND_INTERVAL_MASK);
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_DFLT_METHOD_INTERPOL, CQ_LEGEND_INTERPOL_MASK);
    }
  /* Add service command to the queue */
  params = allocateCommunicator(3, COMM_PTR, COMM_PTR, COMM_PTR);
  params->p[0] = &ctrl; /* First parameter : Pointer to control */
  params->p[1] = &density; /* Second parameter : Pointer to density cut */
  params->p[2] = &idx; /* Third parameter : Pointer to index */
  pushCmdToQueue(&pycmds, srvRecreateLegend, params);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

/* This function changes the method to create ticks only. */
static PyObject *conqtour_legend_set_tick_method(PyObject *self, PyObject *args)
{
  int i;
  char *method;
  comm_t *params;

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "s", &method))
    return NULL;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);

  /* Force the scale ticks to be recalculated using default values */
  ctrl.legend.tickinterval = 0.0;
  ctrl.legend.tickorigin = 0.0;
  if (!strcmp(method, "linear"))
    {
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_TICKS_LINEAR, CQ_LEGEND_TICKS_MASK);
    }
  else if (!strcmp(method, "log"))
    {
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_TICKS_SIGNABSLOG, CQ_LEGEND_TICKS_MASK);
    }
  else
    {
      if (ctrl.verbose)
        printf("WARNING: Method identifier, not in the dictionary. Using defaults\n");
      /* Use default */
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_DFLT_METHOD_TICKS, CQ_LEGEND_TICKS_MASK);
    }
  /* Add service command to the queue */
  params = allocateCommunicator(3, COMM_PTR, COMM_PTR, COMM_PTR);
  params->p[0] = &ctrl; /* First parameter : Pointer to control */
  params->p[1] = &density; /* Second parameter : Pointer to density cut */
  params->p[2] = &idx; /* Third parameter : Pointer to index */
  pushCmdToQueue(&pycmds, srvRecreateLegend, params);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

/* This function changes the method to interpolate colors within each interval */
static PyObject *conqtour_legend_set_interpolation_method(PyObject *self, PyObject *args)
{
  int i;
  char *method;
  comm_t *params;

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "s", &method))
    return NULL;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);

  if (!strcmp(method, "linear"))
    {
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_INTERPOL_LINEAR, CQ_LEGEND_INTERPOL_MASK);
    }
  else if (!strcmp(method, "log"))
    {
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_INTERPOL_SIGNABSLOG, CQ_LEGEND_INTERPOL_MASK);
    }
  else if (!strcmp(method, "constant"))
    {
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_INTERPOL_CONSTANT, CQ_LEGEND_INTERPOL_MASK);
    }
  else
    {
      if (ctrl.verbose)
        printf("WARNING: Method identifier, not in the dictionary. Using defaults\n");
      /* Use default */
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_DFLT_METHOD_INTERPOL, CQ_LEGEND_INTERPOL_MASK);
    }
  /* Add service command to the queue */
  params = allocateCommunicator(3, COMM_PTR, COMM_PTR, COMM_PTR);
  params->p[0] = &ctrl; /* First parameter : Pointer to control */
  params->p[1] = &density; /* Second parameter : Pointer to density cut */
  params->p[2] = &idx; /* Third parameter : Pointer to index */
  pushCmdToQueue(&pycmds, srvRecreateLegend, params);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

/* This function changes the method to create color intervals only. */
static PyObject *conqtour_legend_set_interval_method(PyObject *self, PyObject *args)
{
  int i;
  char *method;
  comm_t *params;

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "s", &method))
    return NULL;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);

  if (!strcmp(method, "linear"))
    {
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_INTERVAL_LINEAR, CQ_LEGEND_INTERVAL_MASK);
    }
  else if (!strcmp(method, "log"))
    {
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_INTERVAL_SIGNABSLOG, CQ_LEGEND_INTERVAL_MASK);
    }
  else
    {
      if (ctrl.verbose)
        printf("WARNING: Method identifier, not in the dictionary. Using defaults\n");
      /* Use default */
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_DFLT_METHOD_INTERVAL, CQ_LEGEND_INTERVAL_MASK);
    }
  /* Add service command to the queue */
  params = allocateCommunicator(3, COMM_PTR, COMM_PTR, COMM_PTR);
  params->p[0] = &ctrl; /* First parameter : Pointer to control */
  params->p[1] = &density; /* Second parameter : Pointer to density cut */
  params->p[2] = &idx; /* Third parameter : Pointer to index */
  pushCmdToQueue(&pycmds, srvRecreateLegend, params);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

/* This function changes the method to assign colors to intervals */
static PyObject *conqtour_legend_set_color_method(PyObject *self, PyObject *args)
{
  int i;
  char *method;
  comm_t *params;

  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "s", &method))
    return NULL;
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);

  if (!strcmp(method, "redgreenblue") || !strcmp(method, "rgb"))
    {
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_COLOR_SATREDBLUE, CQ_LEGEND_COLOR_MASK);
    }
  else if (!strcmp(method, "redblue") || !strcmp(method, "rb"))
    {
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_COLOR_REDBLUE, CQ_LEGEND_COLOR_MASK);
    }
  else if (!strcmp(method, "rainbow"))
    {
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_COLOR_RAINBOW, CQ_LEGEND_COLOR_MASK);
    }
  else
    {
      if (ctrl.verbose)
        printf("WARNING: Method identifier, not in the dictionary. Using defaults\n");
      /* Use default */
      SET_MULTIBIT_FLAG(ctrl.legend.method, CQ_LEGEND_DFLT_METHOD_COLOR, CQ_LEGEND_COLOR_MASK);
    }
  /* Add service command to the queue */
  params = allocateCommunicator(3, COMM_PTR, COMM_PTR, COMM_PTR);
  params->p[0] = &ctrl; /* First parameter : Pointer to control */
  params->p[1] = &density; /* Second parameter : Pointer to density cut */
  params->p[2] = &idx; /* Third parameter : Pointer to index */
  pushCmdToQueue(&pycmds, srvRecreateLegend, params);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

// TODO Perhaps return atom in plane??
// TODO Check that the atoms are within the 3D box
static PyObject *conqtour_view_pick_atoms_slice(PyObject *self, PyObject *args)
{
  int i;
  GLuint s[3]; /* Selected atoms */
  char errmsg[MAXLIN];
  double aux1[3], aux2[3], aux3[3];
  comm_t *params;
  double *dp;

  if (!(ctrl.has & CQ_CTRL_HAS_COORDS) || !(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires coordinates and 3D data");
      return NULL;
    }
  if (!PyArg_ParseTuple(args, "iii", &s[0], &s[1], &s[2]))
    return NULL;
  for (i = 0; i < 3; ++i)
    if (s[i] < 0 || s[i] >= coords.numatoms)
      {
        sprintf(errmsg, "Atom selection out of bounds : %d", s[i]);
        PyErr_SetString(ConqtourError, errmsg);
        return NULL;
      }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);

  // TODO Most of what follows is common to processHits. Write a function

  /* Recalculate normal (slice and control) */
  for (i = 0; i < 3; ++i)
    {
      aux1[i] = coords.xyz[s[1]][i] - coords.xyz[s[0]][i];
      aux2[i] = coords.xyz[s[2]][i] - coords.xyz[s[1]][i];
    }
  crossProduct(aux1, aux2, aux3);
  if (normalize(aux3) < ZEROTOL) /* Probably parallel vectors : Reset */
    {
      PyErr_SetString(ConqtourWarning, "Picked atoms are colinear. Ignored");
      return NULL;
    }
  for (i = 0; i < 3; ++i)
    {
      /* Prepare view for display */
      ctrl.normal[i] = aux3[i];
    }
  params = allocateCommunicator(6, COMM_DATA(3,double), COMM_DATA(3,double), COMM_PTR, COMM_PTR, COMM_PTR, COMM_PTR);
  dp = params->p[0]; /* Cast of the void pointer to double */
  for (i = 0; i < 3; ++i) /* First parameter : Normal to plane */
    dp[i] = aux3[i];
  /* Reset the index, using second atom of the selection as point "inplane" */
  dp = params->p[1]; /* Cast of the void pointer */
  for (i = 0; i < 3; ++i) /* Second parameter : In-plane point */
    dp[i] = coords.xyz[s[1]][i];
  params->p[2] = &idx; /* Third parameter : Pointer to index */
  params->p[3] = &density; /* Fourth parameter : Pointer to density cut */
  params->p[4] = &box; /* Fifth parameter : Pointer to cell box */
  params->p[5] = &ctrl; /* Sixth parameter : Pointer to control */
  pushCmdToQueue(&pycmds, srvResetRedraw, params);
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_toggle_axes(PyObject *self, PyObject *args)
{
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  ctrl.display ^= CQ_CTRL_SHOW_AXES;
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  if (FLAG_UP(ctrl.display, CQ_CTRL_SHOW_AXES))
    Py_RETURN_TRUE;
else    Py_RETURN_FALSE;
  }

static PyObject *conqtour_toggle_box(PyObject *self, PyObject *args)
{
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  ctrl.display ^= CQ_CTRL_SHOW_BOX;
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  if (FLAG_UP(ctrl.display, CQ_CTRL_SHOW_BOX))
    Py_RETURN_TRUE;
else    Py_RETURN_FALSE;
  }

static PyObject *conqtour_toggle_coordinates_shift(PyObject *self, PyObject *args)
{
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  if (ctrl.shftcoords)
    ctrl.shftcoords = FALSE;
  else
    ctrl.shftcoords = TRUE;
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  if (ctrl.shftcoords)
    Py_RETURN_TRUE;
else    Py_RETURN_FALSE;
  }

static PyObject *conqtour_toggle_legend(PyObject *self, PyObject *args)
{
  if (!(ctrl.has & CQ_CTRL_HAS_DENSITY))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires 3D data");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  ctrl.display ^= CQ_CTRL_SHOW_LEGEND;
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  if (FLAG_UP(ctrl.display, CQ_CTRL_SHOW_LEGEND))
    Py_RETURN_TRUE;
else    Py_RETURN_FALSE;
  }

static PyObject *conqtour_toggle_molecule(PyObject *self, PyObject *args)
{
  if (!(ctrl.has & CQ_CTRL_HAS_COORDS))
    {
      PyErr_SetString(ConqtourWarning, "Function disabled : Requires coordinates");
      return NULL;
    }
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  ctrl.display ^= CQ_CTRL_SHOW_MOLECULE;
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  if (FLAG_UP(ctrl.display, CQ_CTRL_SHOW_MOLECULE))
    Py_RETURN_TRUE;
else    Py_RETURN_FALSE;
  }

static PyObject *conqtour_toggle_verbosity(PyObject *self, PyObject *args)
{
  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  requestRedisplay(); /* Force redisplay before acting... */
  waitForRedisplay(&ctrl);
  ++ctrl.verbose; /* For the moment, there is only one level of verbosity, but keep it general */
  if (ctrl.verbose > MAXVERBOSITY)
    ctrl.verbose = NO;
  requestRedisplay();
  waitForRedisplay(&ctrl);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_BuildValue("i", ctrl.verbose);
}

static PyObject *conqtour_info_reset(PyObject *self, PyObject *args)
{
  int i;

  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  if (cqinfo.typename != NULL)
    {
      for (i = 0; i < cqinfo.numspecies; ++i)
        {
          if (cqinfo.typename[i] != NULL)
            free(cqinfo.typename[i]);
        }
      free(cqinfo.typename);
    }
  memset(&cqinfo, 0, sizeof(cqruninfo_t));
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *conqtour_limits_reset(PyObject *self, PyObject *args)
{
  int i;

  pthread_mutex_lock(&ctrl.mutex); /* Be safe */
  resetLimits(&limits);
  pthread_mutex_unlock(&ctrl.mutex);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef
    conqtourMethods[] =
          {
#ifdef DEBUG
                { "cq_0", conqtour_0, METH_NOARGS, "cq_0()\n\n"
                  "Test command queue\n\n"},
#endif
                { "cq_get_current_slice", conqtour_get_current_slice, METH_NOARGS, "cq_get_current_slice()\n\n"
                  "Return the internal index for the current slice.\n\n"
                  "Requires: 3D data\n\n"
                  "This can be used to reset the slice with cq_view_set_slice" },
                { "cq_get_number_atoms", conqtour_get_number_atoms, METH_NOARGS, "cq_get_number_atoms()\n\n"
                  "Requires: Coordinates\n\n"
                  "Return the number of atoms in our coordinates" },
                { "cq_get_number_slices", conqtour_get_number_slices, METH_NOARGS, "cq_get_number_slices()\n\n"
                  "Requires: 3D data\n\n"
                  "Return the maximum index for a slice along the present normal" },
                { "cq_get_slice_plane", conqtour_get_slice_plane, METH_NOARGS, "cq_get_slice_plane()\n\n"
                  "Requires: 3D data\n\n"
                  "Returns the normal along which slices are calculates and a point in the plane of the current slice\n"
                  "as a list of the form [[normal],[point]].\n\n"
                  "The function to set a new direction and plane for slices is cq_view_set_reference_plane\n" },
                { "cq_info_find", conqtour_info_find, METH_NOARGS, "cq_info_find()\n\n"
                  "Try to get information about a Conquest run from the files in the current directory\n"
                  "Note that it will keep known information (provided by the user or from a previous call)\n"
                  "If fresh information is needed, cq_info_reset should be call before cq_info_find\n\n"
                  "Returns None" },
                { "cq_info_reset", conqtour_info_reset, METH_NOARGS, "cq_info_reset()\n\n"
                  "Force a reset of the known information about a Conquest run.\n"
                  "A consequence is that cq_info_find will do its best to gather information\n"
                  "from available Conquest files in the current directory\n\n"
                  "Returns None" },
                { "cq_info_set_block_points", conqtour_info_set_block_points, METH_VARARGS,
                    "cq_info_set_block_points(points)\n\n"
                      "Set the number of points per block (as a tuple) associated to a Conquest run.\n"
                      "To unset the number, so that cq_info_find will update it, use 0 for any component\n\n"
                      "Returns None" },
                { "cq_info_set_blocks_file", conqtour_info_set_blocks_file, METH_VARARGS, "cq_info_set_blocks_file(file)\n\n"
                  "Set the blocks file name associated to a Conquest run.\n"
                  "To unset the value, use an empty string\n\n"
                  "Returns None" },
                { "cq_info_set_cell_dimensions", conqtour_info_set_cell_dimensions, METH_VARARGS,
                    "cq_info_set_cell_dimensions(dimensions)\n\n"
                      "Set the cell dimensions, i.e. lattice vector modules (as a tuple) associated to a Conquest run.\n"
                      "This is NOT the actual value that will be used 'in memory', which depends on the limits\n"
                      "To unset the value, so that cq_info_find will update it, use 0.0 for any component\n\n"
                      "Returns None" },
                { "cq_info_set_coordinates_file", conqtour_info_set_coordinates_file, METH_VARARGS,
                    "cq_info_set_coordinates_file(file)\n\n"
                      "Set the coordinates file name associated to a Conquest run.\n"
                      "To unset the value, use an empty string\n\n"
                      "Returns None" },
                { "cq_info_set_cores", conqtour_info_set_cores, METH_VARARGS, "cq_info_set_cores(cores)\n\n"
                  "Set the number of cores (MPI processes) associated to a Conquest run.\n"
                  "To unset the number, so that cq_info_find will update it, use 0\n\n"
                  "Returns None" },
                { "cq_info_set_density_stub", conqtour_info_set_density_stub, METH_VARARGS,
                    "cq_info_set_density_stub(file stub)\n\n"
                      "Set the stub for the density files associated to a Conquest run.\n"
                      "To unset the value, use an empty string\n\n"
                      "Returns None" },
                { "cq_info_set_fracional_coordinates", conqtour_info_set_fractional_coordinates, METH_VARARGS,
                    "cq_info_set_fractional_coordinates(True/False)\n\n"
                      "Set whether Conquest coordinates are to be interpreted as fractional or not.\n"
                      "Any true (e.g. True, 1) or false (e.g. False, 0) python objects can be used as the argument.\n\n"
                      "Returns None" },
                { "cq_info_set_grid_points", conqtour_info_set_grid_points, METH_VARARGS, "cq_info_set_grid_points(points)\n\n"
                  "Set the number of grid points (as a tuple) associated to a Conquest run.\n"
                  "This is NOT the number of grid points of the data in memory or to be read\n"
                  "but the actual number that was used in the calculation.\n"
                  "To unset the number, so that cq_info_find will update it, use 0 for any component\n\n"
                  "Returns None" },
                { "cq_info_set_input_file", conqtour_info_set_input_file, METH_VARARGS, "cq_info_set_input_file(file)\n\n"
                  "Set the input-parameters file name associated to a Conquest run.\n"
                  "To unset the value, use an empty string\n\n"
                  "Returns None" },
                { "cq_info_set_output_file", conqtour_info_set_output_file, METH_VARARGS, "cq_info_set_output_file(file)\n\n"
                  "Set the output file name associated to a Conquest run.\n"
                  "To unset the value, use an empty string\n\n"
                  "Returns None" },
                { "cq_legend_get_label_format", conqtour_legend_get_label_format, METH_NOARGS, "cq_legend_get_label_format()\n\n"
                  "Get the legend label format in C printf style\n"
                  "Requires: 3D Data\n\n"
                  "Returns None" },
                { "cq_legend_set_color_method", conqtour_legend_set_color_method, METH_VARARGS,
                    "cq_legend_set_color_method(method)\n\n"
                      "Change the method used to assign colors to the color intervals of the legend\n"
                      "The method is a string identifier. Options are listed in "conqtourListLegendColorMethodsName"\n\n"
                    "    Example:\n"
                    "       pyCQ> "conqtourListLegendColorMethodsName" # List the valid method identifiers\n"
                    "       ['redblue', 'rb', 'rainbow'] # This may be different in your version of Conqtour\n"
                    "       pyCQ> cq_legend_set_color_method('rb')    # Change the method to create color intervals and redisplay\n"
                    "       pyCQ> cq_legend_set_color_method("conqtourListLegendColorMethodsName"[0])    # Alternatively\n\n"
                    "Requires: 3D data\n\n"
                    "Returns None" },
                {
                    "cq_legend_set_interpolation_method",
                    conqtour_legend_set_interpolation_method,
                    METH_VARARGS,
                    "cq_legend_set_interpolation_method(method)\n\n"
                      "Change the method used to interpolate colors within each color interval of the legend\n"
                      "The method is a string identifier. Options are listed in "conqtourListLegendInterpolationMethodsName"\n\n"
                    "    Example:\n"
                    "       pyCQ> "conqtourListLegendInterpolationMethodsName" # List the valid method identifiers\n"
                    "       ['linear', 'log', 'constant'] # This may be different in your version of Conqtour\n"
                    "       pyCQ> cq_legend_set_interpolation_method('linear')    # Change the legend interpolation method and redisplay\n"
                    "       pyCQ> cq_legend_set_interpolation_method("conqtourListLegendInterpolationMethodsName"[0])    # Alternatively\n\n"
                    "Requires: 3D data\n\n"
                    "Returns None" },
                {
                    "cq_legend_set_interval_method",
                    conqtour_legend_set_interval_method,
                    METH_VARARGS,
                    "cq_legend_set_interval_method(method)\n\n"
                      "Change the method used to create color intervals of the legend\n"
                      "The method is a string identifier. Options are listed in "conqtourListLegendIntervalMethodsName"\n\n"
                    "    Example:\n"
                    "       pyCQ> "conqtourListLegendIntervalMethodsName" # List the valid method identifiers\n"
                    "       ['linear', 'log'] # This may be different in your version of Conqtour\n"
                    "       pyCQ> cq_legend_set_interval_method('linear')    # Change the method to create color intervals and redisplay\n"
                    "       pyCQ> cq_legend_set_interval_method("conqtourListLegendIntervalMethodsName"[0])    # Alternatively\n\n"
                    "Requires: 3D data\n\n"
                    "Returns None" },
                { "cq_legend_set_label_format", conqtour_legend_set_label_format, METH_VARARGS,
                    "cq_legend_set_label_format(format)\n\n"
                      "Set the legend label format in C printf style, according to the format string\n"
                      "Requires: 3D Data\n\n"
                      "Returns None" },
                { "cq_legend_set_method", conqtour_legend_set_method, METH_VARARGS, "cq_legend_set_method(method)\n\n"
                  "Change the (collective) method to display colors and represent a legend\n"
                  "The method is a string identifier. Options are listed in "conqtourListLegendCollectiveMethodsName"\n\n"
                "    Example:\n"
                "       pyCQ> "conqtourListLegendCollectiveMethodsName"  # List the valid method identifiers\n"
                "       ['linear', 'log'] # This may be different in your version of Conqtour\n"
                "       pyCQ> cq_legend_set_method('linear')    # Change the legend method and redisplay\n"
                "       pyCQ> cq_legend_set_method("conqtourListLegendCollectiveMethodsName"[0])    # Alternatively\n\n"
                "NOTE: Ticks will be recalculated using default values\n\n"
                "Requires: 3D data\n\n"
                "Returns None" },
                { "cq_legend_set_tick_method", conqtour_legend_set_tick_method, METH_VARARGS,
                    "cq_legend_set_tick_method(method)\n\n"
                      "Change the tick method to display colors and represent a legend\n"
                      "The method is a string identifier. Options are listed in "conqtourListLegendTickMethodsName"\n\n"
                    "    Example:\n"
                    "       pyCQ> "conqtourListLegendTickMethodsName" # List the valid method identifiers\n"
                    "       ['linear', 'log'] # This may be different in your version of Conqtour\n"
                    "       pyCQ> cq_legend_set_tick_method('linear')    # Change the legend tick method and redisplay\n"
                    "       pyCQ> cq_legend_set_tick_method("conqtourListLegendTickMethodsName"[0])    # Alternatively\n\n"
                    "NOTE: Ticks will be recalculated using default values\n\n"
                    "Requires: 3D data\n\n"
                    "Returns None" },
                { "cq_limits_get", conqtour_limits_get, METH_NOARGS, "cq_limits_get()\n\n"
                  "Return the limits used to import and export data.\n"
                  "These limits do not necessarily agree with in-memory data\n\n"
                  "The form of the returned list is [[grid window],[grid shift],[strides],[real window],[real shift]].\n\n" },
                { "cq_limits_reset", conqtour_limits_reset, METH_NOARGS, "cq_limits_reset()\n\n"
                  "Reset the limits used to import and export data.\n\n"
                  "These limits will be all zero, except for the strides, which will be 1\n\n" },
                { "cq_limits_set", conqtour_limits_set, METH_VARARGS, "cq_limits_set(limits)\n\n"
                  "Set the limits used to import and export data.\n"
                  "These limits do not necessarily agree with in-memory data\n\n"
                  "The form of the input list is [[grid window],[grid shift],[strides],[real window],[real shift]].\n\n" },
                { "cq_limits_set_grid_shift", conqtour_limits_set_grid_shift, METH_VARARGS,
                    "cq_limits_set_grid_shift(shifts)\n\n"
                      "Change the shift of the grid; used when reading a density\n"
                      "Use negative values to keep the current shift of a particular coordinate\n\n"
                      "Requires: Nothing, but checks of bounds are only possible when CQ info is available\n\n"
                      "Returns None" },
                { "cq_limits_set_grid_shift_from_real", conqtour_limits_set_grid_shift_from_real, METH_NOARGS,
                    "cq_limits_set_grid_shift_from_real()\n\n"
                      "Set grid shift from the current real shift and available CQ information\n"
                      "If there is no CQ information, it will reset the grid shift to zero\n"
                      "To choose an arbitrary value, use cq_limits_set_grid_shift\n\n"
                      "Returns None" },
                { "cq_limits_set_grid_stride", conqtour_limits_set_grid_stride, METH_VARARGS,
                    "cq_limits_set_grid_stride(strides)\n\n"
                      "Change the stride of the grid; used when reading or writing a density\n"
                      "Use 0 to keep the current stride for a coordinate\n\n"
                      "Requires: Nothing, but checks of bounds are only possible when CQ info is available\n\n"
                      "Returns None" },
                { "cq_limits_set_grid_window", conqtour_limits_set_grid_window, METH_VARARGS,
                    "cq_limits_set_grid_window(points)\n\n"
                      "Change the number of points to be stored in memory; used when reading a density\n"
                      "It must be at least 2 per coordinate, but use zero or one to reset a component (to zero)\n"
                      "and negative numbers to keep the current value\n\n"
                      "Requires: Nothing, but checks of bounds are only possible when CQ info is available\n\n"
                      "Returns None" },
                { "cq_limits_set_real_shift", conqtour_limits_set_real_shift, METH_VARARGS, "cq_limits_set_real_shift(shift)\n\n"
                  "Set real shift (in Angstrom) to any desired value; used to shift output\n"
                  "The grid shift is used to read in and out actual data, so the\n"
                  "real shift is used only to shift output coordinates\n"
                  "To use grid window information to set the real shift,\n"
                  "use conqtour_limits_set_real_shift_from_grid\n\n"
                  "Returns None" },
                { "cq_limits_set_real_shift_from_grid", conqtour_limits_set_real_shift_from_grid, METH_NOARGS,
                    "cq_limits_set_real_shift_from_grid()\n\n"
                      "Set real shift (in Angstrom) from the grid shift and available CQ information\n"
                      "If there is no CQ information, it will reset the real shift to zero\n"
                      "To choose an arbitrary value, use cq_limits_set_real_shift\n\n"
                      "Returns None" },
                { "cq_limits_set_real_window", conqtour_limits_set_real_window, METH_VARARGS,
                    "cq_limits_set_real_window(window)\n\n"
                      "Set real window (in Angstrom) to a specified value.\n"
                      "Use zero to reset a component and negative numbers to keep the current value\n\n"
                      "Requires: Nothing, but checks of bounds are only possible when CQ info is available\n\n"
                      "Returns None" },
                { "cq_load_cq_density", conqtour_load_cq_density, METH_NOARGS, "cq_load_cq_density()\n\n"
                  "Try to load a density using Conquest information, collected from your directory\n"
                  "Returns None" },
                { "cq_load_xplor_density", conqtour_load_xplor_density, METH_VARARGS, "cq_load_xplor_density(file)\n\n"
                  "Open and load a pre-existing Xplor file\n"
                  "Returns None" },
                { "cq_load_cq_coordinates", conqtour_load_cq_coordinates, METH_VARARGS, "cq_load_cq_coordinates()\n\n"
                  "Open and load a Conquest coordinates file (if enough information is available)\n"
                  "Returns None" },
                { "cq_load_xyz_coordinates", conqtour_load_xyz_coordinates, METH_VARARGS, "cq_load_xyz_coordinates(file)\n\n"
                  "Open and load an xyz structure file\n"
                  "Returns None" },
                { "cq_toggle_axes", conqtour_toggle_axes, METH_NOARGS, "cq_toggle_axes()\n\n"
                  "Show/hide the axes (a very primitive representation at the center of the cell)\n\n"
                  "Requires: 3D Data\n\n"
                  "Returns: current state" },
                { "cq_toggle_box", conqtour_toggle_box, METH_NOARGS, "cq_toggle_box()\n\n"
                  "Show/hide the cell box\n\n"
                  "Requires: 3D Data\n\n"
                  "Returns: current state" },
                { "cq_toggle_coordinates_shift", conqtour_toggle_coordinates_shift, METH_NOARGS,
                    "cq_toggle_coordinates_shift()\n\n"
                      "Change the shifting method for shifting coordinates. If this toggle is false,\n"
                      "coordinates are read literally; if true, a shift is applied to them when reading/writing\n\n"
                      "Requires: Nothing\n\n"
                      "Returns: current state" },
                { "cq_toggle_legend", conqtour_toggle_legend, METH_NOARGS, "cq_toggle_legend()\n\n"
                  "Show/hide the legend\n\n"
                  "Requires: 3D Data\n\n"
                  "Returns: current state" },
                { "cq_toggle_molecule", conqtour_toggle_molecule, METH_NOARGS, "cq_toggle_molecule()\n\n"
                  "Show/hide the molecule\n\n"
                  "Requires: Coordinates\n\n"
                  "Returns: current state" },
                { "cq_toggle_verbosity", conqtour_toggle_verbosity, METH_NOARGS, "cq_toggle_verbosity()\n\n"
                  "Be/stop being verbose\n\n"
                  "Returns the verbosity level" },
                { "cq_view_get", conqtour_view_get, METH_NOARGS, "cq_view_get()\n\n"
                  "Return the view parameters in a way that can be used to call cq_view_set\n\n"
                  "    Example:\n"
                  "       pyCQ> view=cq_view_get()   # Now change the view using the mouse\n"
                  "       pyCQ> cq_view_set(view)    # This restores the previous view\n\n"
                  "The form of the returned list is [[center],[up],[normal],[distance,field,near,far]].\n"
                  "At present, 'near' and 'far' are ignored by cq_view_set, \n"
                  "but this could change in the future\n\n"
                  "Requires: 3D data\n\n" },
                { "cq_view_pick_atoms_slice", conqtour_view_pick_atoms_slice, METH_VARARGS,
                    "cq_view_pick_atoms_slice(atom1,atom2,atom3)\n\n"
                      "Change the direction of the slices using three atoms.\n"
                      "A slice index is recalculated and a slice going through the atoms is shown.\n\n"
                      "Requires: Coordinates and 3D data\n\n"
                      "Returns None" },
                { "cq_view_recenter", conqtour_view_recenter, METH_NOARGS, "cq_view_recenter()\n\n"
                  "Recenter view\n\n"
                  "Requires: 3D data\n\n"
                  "Returns the new center as a list" },
                { "cq_view_rotate", conqtour_view_rotate, METH_VARARGS, "cq_view_rotate(angle)\n\n"
                  "Rotate the view (including the slice) 'angle' degrees around the up vector\n\n"
                  "Requires: 3D data\n\n"
                  "Returns None" },
                { "cq_view_rotate_reslice", conqtour_view_rotate_reslice, METH_VARARGS, "cq_view_rotate_reslice(angle)\n\n"
                  "Rotate view 'angle' degrees, but keep slice parallel to screen (non-tracking mode)\n\n"
                  "Requires: 3D data\n\n"
                  "Returns None" },
                { "cq_view_set", conqtour_view_set, METH_VARARGS,
                    "cq_view_set([center][up][normal][distance,field,near,far])\n\n"
                      "Set the view to the parameters passed in a tuple of the form returned by cq_view_get\n\n"
                      "    Example:\n"
                      "       pyCQ> view=cq_view_get()   # Now change the view using the mouse\n"
                      "       pyCQ> cq_view_set(view)    # This restores the previous view\n\n"
                      "The form of the list is [[center],[up],[normal],[distance,field,near,far]].\n"
                      "At present, 'near' and 'far' are ignored, but this could change in the future\n\n"
                      "Requires: 3D data\n\n" },
                { "cq_view_set_center", conqtour_view_set_center, METH_VARARGS, "cq_view_set_center(coordinates)\n\n"
                  "Change the center to the coordinates given in the parameters (as a list or tuple)\n\n"
                  "Requires: 3D data\n\n"
                  "Returns None" },
                { "cq_view_set_reference_plane", conqtour_view_set_reference_plane, METH_VARARGS,
                    "cq_view_set_reference_plane([[normal],[point in plane]])\n\n"
                      "Redraw the slice so that it is on the plane defined by the input parameters\n\n"
                      "    Example:\n"
                      "       pyCQ> cq_view_set_reference_plane([[0,1,0],[7,7,7]])   # Make a slice on the XZ plane\n\n"
                      "If the point is outside of the bounding box, the center of the box will be used\n\n"
                      "To get the current normal and point, use cq_get_slice_plane\n\n"
                      "Requires: 3D data\n\n"
                      "Returns None" },
                { "cq_view_set_slice", conqtour_view_set_slice, METH_VARARGS, "cq_view_set_slice(slice)\n\n"
                  "Change the slice to the 'slice' cut along the present normal\n\n"
                  "Requires: 3D data\n\n"
                  "Returns None" },
                { "cq_write_cq_coordinates", conqtour_write_cq_coordinates, METH_VARARGS, "cq_write_cq_coordinates(file)\n\n"
                  "Write a coordinates file in Conquest format (see manual for many more details).\n\n"
                  "Requires: Coordinates\n\n"
                  "Returns None" },
                { "cq_write_ppm_image", conqtour_write_ppm_image, METH_VARARGS, "cq_write_ppm_image(file)\n\n"
                  "Write a PPM image to the specified file\n\n"
                  "Requires: 3D data\n\n"
                  "Returns None" },
                { "cq_write_pymol_slice", conqtour_write_pymol_slice, METH_VARARGS, "cq_write_pymol_slice(file)\n\n"
                  "Write the current slice to a pymol script\n\n"
                  "Requires: 3D data\n\n"
                  "Returns None" },
                { "cq_write_xplor_density", conqtour_write_xplor_density, METH_VARARGS, "cq_write_xplor_density(file)\n\n"
                  "Write an xplor file from the 3D data, using available limits\n"
                  "Some of them could change as a result of this call, if they exceed available data\n\n"
                  "Requires: 3D data\n\n"
                  "Returns None" },
                { "cq_write_xyz_coordinates", conqtour_write_xyz_coordinates, METH_VARARGS, "cq_write_xyz_coordinates(file)\n\n"
                  "Write a coordinates file in XYZ format (see manual for many more details).\n\n"
                  "Requires: Coordinates\n\n"
                  "Returns None" },
                { NULL, NULL, 0, NULL } /* Sentinel */
          };

void attachListToModule(PyObject *module, PyObject *list, char **entries, char *name)
{
  int i;

  for (i = 0; entries[i] != NULL; ++i)
    ; /* First, count the number of elements */
  list = PyList_New(i);
  for (i = 0; entries[i] != NULL; ++i)
    {
      PyList_SetItem(list, i, Py_BuildValue("s", entries[i]));
    }
  PyModule_AddObject(module, name, list);
}

PyMODINIT_FUNC init_conqtour_pymodule(void)
{
  PyObject *module;
  int i;

  module = Py_InitModule("conqtour", conqtourMethods);
  if (module == NULL)
    return;

  /**** Lists ****/
  /* Legend collective methods */
  attachListToModule(module, conqtourListLegendCollectiveMethods, conqtourLegendCollectiveMethods,
      conqtourListLegendCollectiveMethodsName);
  /* Legend tick methods */
  attachListToModule(module, conqtourListLegendTickMethods, conqtourLegendTickMethods, conqtourListLegendTickMethodsName);
  /* Legend interval methods */
  attachListToModule(module, conqtourListLegendIntervalMethods, conqtourLegendIntervalMethods,
      conqtourListLegendIntervalMethodsName);
  /* Legend interpolation methods */
  attachListToModule(module, conqtourListLegendInterpolationMethods, conqtourLegendInterpolationMethods,
      conqtourListLegendInterpolationMethodsName);
  /* Legend color methods */
  attachListToModule(module, conqtourListLegendColorMethods, conqtourLegendColorMethods, conqtourListLegendColorMethodsName);

  /**** Python exceptions ****/
  ConqtourError = PyErr_NewException("conqtour.error", NULL, NULL);
  Py_INCREF(ConqtourError);
  PyModule_AddObject(module, "error", ConqtourError);
  ConqtourWarning = PyErr_NewException("conqtour.warning", NULL, NULL);
  Py_INCREF(ConqtourWarning);
  PyModule_AddObject(module, "warning", ConqtourWarning);
}
