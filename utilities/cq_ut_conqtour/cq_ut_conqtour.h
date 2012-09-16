/* -----------------------------------------------------------------------------
 * $Id: $
 * -----------------------------------------------------------------------------
 * File cq_ut_conqtour.h
 * -----------------------------------------------------------------------------
 *
 * ***** Conquest/utilities/cq_ut_conqtour.h *
 *
 * NAME
 *  cq_ut_conqtour.h
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

#ifndef CQ_UT_CONQTOUR_H_
#define CQ_UT_CONQTOUR_H_

#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <sys/types.h>
#include <unistd.h>
#include <pthread.h>

//#define ARRAY_ZYX
//#define TRIANGLES

#define PROGRAM_NAME "cq_ut_conqtour"

#define MAXLIN 2056
#define DEG_TO_RAD        0.017453292
#define RAD_TO_DEG        57.29577951
#define BOHR_TO_ANGSTROM  0.529177249

#define ROTATE 0
#define SCALE 1
#define TRANSLATE 2
#define NONE 3

#define YES 1
#define NO 0
#define TRUE 1
#define FALSE 0
#define VERBOSE 1
#define MAXVERBOSITY 1
#define SILENT 0

#define SIDE 0
#define LINE 0
#define CORNER 1
#define LEFT 0
#define RIGHT 1
#define X 0
#define Y 1

#define MAXATOMS 5

#define MAX(a,b) ((a > b)?(a):(b))
#define SIGN(a) ((a > 0.0)?1:-1)

#define PICKTOL 1
#define PICKBUFSIZE 1024
#define ZEROTOL 1.0e-16
#define ZEROTOL2 0.00001      /* This tolerance is key in many calculations of the slice */
#define VERYLARGE 1.0e+300
#define ZERO 0.0

#ifndef WAITING_TIME
#define WAITING_TIME   1
#endif
// TODO TEMPORARY : Remove after debugging
int minindx, minindy, maxindx, maxindy;
double minx, miny, maxx, maxy;
int pivot[4];

/***** Information required to load from Conquest files *****/

#define CQ_UNITS_MASKDIST      0x07  /* Three bits reserved for distance units */
#define CQ_UNITS_BOHR          0x01  /* Either bohr or angstrom */
#define CQ_UNITS_ANGSTROM      0x02
#define CQ_UNITS_DFLTDIST      CQ_UNITS_BOHR
// TODO Change flags for files to XXXXFL (to improve readability)
#define CQ_HAVE_INPUT          0x01  /* CQ input file is available */
#define CQ_HAVE_OUTPUT         0x02  /* CQ output file is available */
#define CQ_HAVE_DENSITY        0x04  /* CQ density files are available */
#define CQ_HAVE_COORDINATES    0x08  /* CQ coordinates are available */
#define CQ_HAVE_BLOCKS         0x10  /* CQ blocks file is available */
#define CQ_HAVE_FROM_BLOCKS    0x20  /* Unset if the block method is unknown */
#define CQ_HAVE_FRACTIONAL     0x40  /* Unset if the kind of coordinate file (fractional or not) is not known */
#define CQ_HAVE_HILBERT        0x80  /* Unset if the partition method is not known */
#define CQ_READ_FROM_BLOCKS    0x01  /* Use a list of blocks; otherwise, use old format */
#define CQ_READ_FRACTIONAL     0x02  /* Coordinates are fractional (used also for output) */
#define CQ_READ_HILBERT        0x04  /* Hilbert partitioning was used */
#define CQ_READ_BZ2            0x08  /* Density files in bzip2 format */
#define CQ_READ_GZ             0x10  /* Density files in gunzip format */

#define CQ_LABEL_BLOCK    "%block "

/* Default names that Conquest understands; can be redefined at compile time */
#ifndef CQ_DFLT_INPUT
#define CQ_DFLT_INPUT     "Conquest_input"
#endif
#ifndef CQ_DFLT_OUTPUT
#define CQ_DFLT_OUTPUT    "Conquest_out"
#endif
#ifndef CQ_DFLT_MAKEBLK
#define CQ_DFLT_MAKEBLK   "make_blk.dat"
#endif
#ifndef CQ_DFLT_MAKEBLK_HILBERT
#define CQ_DFLT_MAKEBLK_HILBERT   "hilbert_make_blk.dat"
#endif
#ifndef CQ_DFLT_COORDINATES
#define CQ_DFLT_COORDINATES       "coordinates.in"   /* This is not from Conquest, which has no default */
#endif
#ifndef CQ_DFLT_DENSITY
#define CQ_DFLT_DENSITY           "chden"            /* But note this is hardcoded in CQ :-( */
#endif
#ifndef CQ_DFLT_BLOCKS_RASTER
#define CQ_DFLT_BLOCKS_RASTER     "make_blk.dat"
#endif
#ifndef CQ_DFLT_BLOCKS_HILBERT
#define CQ_DFLT_BLOCKS_HILBERT    "hilbert_make_blk.dat"
#endif
#ifndef CQ_DFLT_BLOCKS
#define CQ_DFLT_BLOCKS            CQ_DFLT_BLOCKS_HILBERT
#endif
#ifndef CQ_DFLT_FROM_BLOCKS
#define CQ_DFLT_FROM_BLOCKS       "True"
#endif
#ifndef CQ_DFLT_FRACTIONAL
#define CQ_DFLT_FRACTIONAL        "True"
#endif
#ifndef CQ_DFLT_PARTITIONER
#define CQ_DFLT_PARTITIONER       "Hilbert"
#endif

// The next one cannot be changed by the user: "Hilbert" is "Hilbert" is "Hilbert"!!
#define CQ_HILBERT_NAME           "Hilbert"

/* Define input flags in a way that can be overridden at compile time*/
#ifndef CQ_FLAG_OUTPUTFILE
#define CQ_FLAG_OUTPUTFILE   "IO.OutputFile"
#endif
#ifndef CQ_FLAG_COORDFILE
#define CQ_FLAG_COORDFILE    "IO.Coordinates"
#endif
#ifndef CQ_FLAG_BLOCKSFILE
#define CQ_FLAG_BLOCKSFILE   "IO.Blocks"
#endif
#ifndef CQ_FLAG_FROM_BLOCKS
#define CQ_FLAG_FROM_BLOCKS  "ReadBlocksFromFile"
#endif
#ifndef CQ_FLAG_FRACTIONAL
#define CQ_FLAG_FRACTIONAL   "IO.FractionalAtomicCoords"
#endif
#ifndef CQ_FLAG_DISTUNITS
#define CQ_FLAG_DISTUNITS    "General.DistanceUnits"
#endif
#ifndef CQ_FLAG_PARTITIONER
#define CQ_FLAG_PARTITIONER  "General.PartitionMethod"
#endif
#ifndef CQ_FLAG_SPECIES
#define CQ_FLAG_SPECIES      "General.NumberOfSpecies"
#endif
#ifndef CQ_FLAG_CHEMBLOCK
#define CQ_FLAG_CHEMBLOCK    "ChemicalSpeciesLabel"
#endif

/* These tags are for names of units and the like */
#ifndef CQ_TAG_BOHR
#define CQ_TAG_BOHR          "a0"
#endif
#ifndef CQ_TAG_BOHR1
#define CQ_TAG_BOHR1         CQ_TAG_BOHR
#endif
#ifndef CQ_TAG_BOHR2
#define CQ_TAG_BOHR2         "bohr"
#endif
#ifndef CQ_TAG_ANGSTROM
#define CQ_TAG_ANGSTROM      "A"
#endif

#ifndef CQ_PATTERN_CORES
#define CQ_PATTERN_CORES     "The calculation will be performed on %d processors"
#endif
#ifndef CQ_PATTERN_GRIDPTS
#define CQ_PATTERN_GRIDPTS   "The number of cell grid points in each direction is :\n"\
                             " %d cell grid points along x\n"\
                             " %d cell grid points along y\n"\
                             " %d cell grid points along z"
#endif
#ifndef CQ_PATTERN_BLOCKPTS
#define CQ_PATTERN_BLOCKPTS  "The number of cell grid points in each block is :\n"\
                             " %d cell grid points along x\n"\
                             " %d cell grid points along y\n"\
                             " %d cell grid points along z"
#endif
#ifndef CQ_PATTERN_BOX
#define CQ_PATTERN_BOX       "The simulation box has the following dimensions\n"\
                             " a = %lf b = %lf c = %lf %s"
#endif

#define EXT_NONE ""
#define EXT_BZ2  ".bz2"

/***** Macros to handle flags, etc. *****/
#define SET_FLAG(variable,flag,value)       ((value != 0)?(variable |= flag):(variable &= ~flag))
#define SET_MULTIBIT_FLAG(variable,value,mask)       (variable = (variable & ~mask) | (value & mask))
#define FLAG_UP(variable,flag)              (variable & flag)
#define REPORT_FLAG(variable,flag)          ((variable & flag)?"True":"False")
#define REPORT_FLAG_SHORT(variable,flag)    ((variable & flag)?"T":"F")
#define REPORT_DISTUNITS(variable)          (((variable & CQ_UNITS_MASKDIST) == CQ_UNITS_BOHR)?CQ_TAG_BOHR\
                                           :(((variable & CQ_UNITS_MASKDIST) == CQ_UNITS_ANGSTROM)?CQ_TAG_ANGSTROM:"Unknown"))
#define REPORT_PARTITIONER(variable,flag)   ((variable & flag)?CQ_HILBERT_NAME:"Not "CQ_HILBERT_NAME)
#define REPORT_COORDINATE(index)            ((index == 0)?"x":((index == 1)?"y":((index == 2)?"z":"unknown")))

typedef struct
{
  int cores; /* Number of cores (or processes) where the job was run */
  int pt[3]; /* Grid points in the run (this could be larger than in density_t) */
  int blockpt[3]; /* Grid points per block */
  double lvm[3]; /* Module of lattice vectors (dimensions of the simulation box) */
  int numspecies; /* Number of atom types in the coordinates file */
  char **typename; /* List of atom types */
  int units; /* Flags about units */
  int have; /* Flags about available files */
  int read; /* Flags about loading conditions */
  char densityfl[MAXLIN]; /* Stub of the density files */
  char coordinatesfl[MAXLIN]; /* Name of the coordinates file */
  char inputfl[MAXLIN]; /* Name of the Conquest input file */
  char outputfl[MAXLIN]; /* Name of the Conquest output file */
  char blocksfl[MAXLIN]; /* Name of the list of blocks */
} cqruninfo_t;

/***** Information about "in memory" data *****/

typedef struct
{
  double realshft[3]; /* Shift in Angstroms, used to read/write a portion of a density */
  double realdims[3]; /* Dimensions of the (sub)box to be read */
  int gridshft[3]; /* Translation of the shift into grid points. This is what counts */
  int griddims[3]; /* Translation of the box dimensions */
  int strides[3]; /* Used to skip points of the input */
/* Note: gridshift and griddims are the final numbers (after striding) */
} limits_t;

// TODO Consider using enums for these
/* The following are used in makeLegend, where the actual scale (ticks, etc.) is created */
#define CQ_LEGEND_TICKS_MASK 0x0F /* Mask for the ticks method bits */
#define CQ_LEGEND_TICKS_LINEAR 0x01 /* Ticks from origin at constant linear steps (linear scale) */
#define CQ_LEGEND_TICKS_SIGNABSLOG 0x02 /* Ticks from origin at constant log10 steps (signed log of abs val) */
#define CQ_LEGEND_INTERPOL_MASK 0xF0 /* Mask for interpolation method bits */
#define CQ_LEGEND_INTERPOL_LINEAR 0x10 /* Interpolate colors linearly between color interval extremes */
#define CQ_LEGEND_INTERPOL_SIGNABSLOG 0x20 /* Interpolate colors using signed abs. logs. btw. color interval extremes */
#define CQ_LEGEND_INTERPOL_CONSTANT 0x30 /* Constant color (using lower end val) within each color interval */
/* The following are used in prepareLegend, where the coloring system (color intervals and values) is created */
/* Even when the color intervals are signabslog, the ticks could be linear, and vice versa */
#define CQ_LEGEND_INTERVAL_MASK 0x0F00 /* Mask for the interval method bits */
#define CQ_LEGEND_INTERVAL_LINEAR 0x0100 /* Create color intervals using a linear scale */
#define CQ_LEGEND_INTERVAL_SIGNABSLOG 0x0200 /* Create color intervals using orders of magnitude */
#define CQ_LEGEND_COLOR_MASK 0xF000 /* Mask for the color method bits */
#define CQ_LEGEND_COLOR_SATREDBLUE 0x1000 /* Create colors using saturated red for pos. and blue for neg. + proportional green */
#define CQ_LEGEND_COLOR_REDBLUE 0x2000 /* Create colors using proportional red for pos. and blue for neg.; green, not used */
#define CQ_LEGEND_COLOR_RAINBOW 0x3000 /* Create colors using a variety of colors */

/* Legend collective methods */

/* Implies CQ_LEGEND_TICKS_LINEAR, CQ_LEGEND_INTERPOL_LINEAR, CQ_LEGEND_INTERVAL_LINEAR */
#define CQ_LEGEND_COLLECTIVE_METHOD_LINEAR 0x01
/* Implies CQ_LEGEND_TICKS_SIGNABSLOG, CQ_LEGEND_INTERPOL_SIGNABSLOG, CQ_LEGEND_INTERVAL_SIGNABSLOG */
#define CQ_LEGEND_COLLECTIVE_METHOD_LOG 0x02

/* Default values, used to (re)set a ctrl structure */
#define CQ_LEGEND_DFLT_X 30
#define CQ_LEGEND_DFLT_Y 70
#define CQ_LEGEND_DFLT_W 20
#define CQ_LEGEND_DFLT_H 300
#define CQ_LEGEND_DFLT_TICKORIGIN 0.0

#define CQ_LEGEND_DFLT_METHOD_TICKS CQ_LEGEND_TICKS_LINEAR
//#define CQ_LEGEND_DFLT_METHOD_TICKS CQ_LEGEND_TICKS_SIGNABSLOG

//#define CQ_LEGEND_DFLT_METHOD_INTERPOL CQ_LEGEND_INTERPOL_CONSTANT
#define CQ_LEGEND_DFLT_METHOD_INTERPOL CQ_LEGEND_INTERPOL_LINEAR
//#define CQ_LEGEND_DFLT_METHOD_INTERPOL CQ_LEGEND_INTERPOL_SIGNABSLOG

#define CQ_LEGEND_DFLT_METHOD_INTERVAL CQ_LEGEND_INTERVAL_LINEAR
//#define CQ_LEGEND_DFLT_METHOD_INTERVAL CQ_LEGEND_INTERVAL_SIGNABSLOG

//#define CQ_LEGEND_DFLT_METHOD_COLOR CQ_LEGEND_COLOR_SATREDBLUE
//#define CQ_LEGEND_DFLT_METHOD_COLOR CQ_LEGEND_COLOR_REDBLUE
#define CQ_LEGEND_DFLT_METHOD_COLOR CQ_LEGEND_COLOR_RAINBOW

#define CQ_LEGEND_DFLT_COLLECTIVE_METHOD CQ_LEGEND_COLLECTIVE_METHOD_LINEAR

#define CQ_LEGEND_DFLT_BORDER_THICKNESS 2
#define CQ_LEGEND_DFLT_TICK_LENGTH 2
#define CQ_LEGEND_DFLT_APPROX_TICKS 10
#define CQ_LEGEND_DFLT_INTERVALS 10 /* This is not necessarily used; it depends on the interval method */
#define CQ_LEGEND_DFLT_FORMAT "%3.2e"
#define CQ_LEGEND_POINTS_BETWEEN_SIGNS 15 /* The number of points to be used to separate neg. from pos. in log. scales */

typedef struct
{
  char good; /* If NO, this is not a valid color, and should be represented in some special way */
  int r; /* Red */
  int g; /* Green */
  int b; /* Blue */
} rgbcolor_t;

typedef struct
{
  int pixels;
  double *val; /* Value attached to that pixel */
  char *tick; /* Which rows should be marked with ticks? */
  char *label; /* Which rows should be labeled? */
} scalemarks_t;

typedef struct
{
  int x; /* Position of the legend (including border) */
  int y;
  int w; /* Width of the legend (w/o borders) */
  int h; /* Height of the legend (w/o borders) */
  int bt; /* Border thickness */
  int tl; /* Tick length */
  double tickorigin; /* "Origin" of the ticks, used to make sure this tick is shown */
  double tickinterval; /* This is ignored or interpreted as a linear or logarithmic increment (depends on tick method) */
  int noticks; /* Number of ticks; this can be just a suggestion (in ctrl structs) or the actual number (in densities) */
  int nointer; /* Number of color intervals */
  int method; /* Interpolation method, log/linear scale, etc. */
  double min; /* Minimum datum */
  double max; /* Maximum datum */
  double *val; /* Values that define the intervals; val[0] <= min ; val[nointer] >= max */
  rgbcolor_t *col; /* Colors corresponding to the val array */
  char format[MAXLIN]; /* Format of the labels, in printf style */
  scalemarks_t marks; /* Information about marks (ticks and labels) */
  int warn; /* A flag to control the amount of output */
} legend_t;

typedef struct
{
  int pt[3]; /* Grid points */
  double lvm[3]; /* Module of lattice vectors */
  double angle[3]; /* Cell angles */
  double ***data; /* Density data */
  double pmin; /* Min value for positive data */
  double pmax; /* Max value for positive data */
  double nmin; /* Min (absolute) value for negative data */
  double nmax; /* Max (absolute) value for negative data */
  double del[3]; /* Spacing between grid points */
  legend_t legend; /* Legend */
} density_t;

typedef struct
{
  int numatoms;
  char **atname;
  double **xyz;
  char *move; /* Movement flags. Not strictly for XYZ files, but it's convenient here for Conquest files */
} xyzcoords_t;

typedef struct
{ /* This structure defines a box */
  double corner[8][3]; /* Vertices of the box */
  double refpoint[12][3]; /* Each edge is given by the initial point */
  double edge[12][3]; /* and a vector from that point to the final point */
} box_t;

typedef struct
{ /* This structure represents the intersection of a plane with a box */
  int nocorners;
  double corner[6][3]; /* A maximum of 6 vertices define the intersection */
  int refcorner[2]; /* Two corners sitting on the 1st edge of the rectangle, used for reference */
  GLdouble inplane[3]; /* Point in density slice; completes the plane */
  double xaxis[3]; /* Axes relative to the density cut, for convenience */
  double yaxis[3];
  GLdouble rectangle[4][3]; /* Minimum-area enclosing rectangle, in world coords. */
  double area;
  double dx; /* Spacing between grid points : Not necessarily the same as for density_t */
  double dy; /* Note that dz is common to all slices, and it is thus kept in the index */
  int pt[2]; /* Number of grid points in x and y directions */
  double **data; /* Interpolated data */
  double cturn[3]; /* Interpolation at the (potential) two turns (see below) */
  /* The third element is reserved for the "peak" in triangle slices */
  GLuint slice; /* Name for the precompiled list */
  int *lim[2]; /* Extremes of plottable data (min and max y coords.) */
  double *xbord[2]; /* Extreme x coordinates (in rectangle coords.) */
  int turn[2][2]; /* Potentially, there are two corners that need to be treated specially */
/* We store info about the scanning line (turns[LEFT/RIGHT][LINE]) */
/*  and corner (turns[LEFT/RIGHT][CORNER]) */
} densitycut_t;

typedef struct
{
  densitycut_t **dindex; /* Index of density cuts */
  GLuint cuts; /* Maximum number of cuts (slices) in current view */
  GLuint curr; /* Current (visible) slice */
  GLdouble surfnormal[3]; /* Normal to the density slice. Stored here because common to all cuts */
  GLdouble inplane[3]; /* This point is used as an origing for all cuts */
  GLuint orig; /* In this is to remember which slice was the first */
  double dz;
  FILE *fp; /* A pointer to a file where a representation will be written (usually a pymol script) */
} index_t;

/** Queue structures **/

/* Queue for names (used to store the py scripts requested in the command line) */
typedef struct _namequeuenode_t namequeuenode_t; /* Self-referencing: we need to pre-declare the type */

struct _namequeuenode_t
{
  char name[MAXLIN]; /* String to contain a name */
  namequeuenode_t *next; /* Pointer to next node */
};

typedef struct
{
  namequeuenode_t *head; /* Pointer to head node */
  namequeuenode_t *tail; /* Pointer to tail node */
} namequeue_t;

/* New type for the length of an array within a communicator */
typedef unsigned int commlen_t;

/* A "communicator" type, used to pass arbitrary parameters through a unified interface */
/* The main purpose is to have a way to stored commands in a queue using pointers to functions */
typedef struct
{
  int num; /* Number of p's */
  void **p; /* Array of pointers to data */
  commlen_t *length; /* Length of each array; zero for actual arrays */
  size_t *size; /* Size of each piece of data; if zero, keep the pointer */
  char *dealloc; /* Is the corresponding p element a true array that must be deallocated after use? */
} comm_t;

/* Communicator macros */
/*  These are provided for legibility: each call corresponds to a array of data/pointer of the communicator */
/*  It is easier to identify each piece if the length and size are visually separated */
#define COMM_DATA(a,b)  a,sizeof(b) /* a = length of data array; b = size per unit */
#define COMM_PTR NULL /* Placeholder to indicate that the next piece of data of a communicator is an actual pointer */

/* Queue for commands requested from python, but to be executed by the GL thread */
typedef struct _cmdqueuenode_t cmdqueuenode_t; /* Self-referencing: we need to pre-declare the type */

struct _cmdqueuenode_t
{
  comm_t *params; /* Pointer to communicator with parameter information */
  int (*cmd)(comm_t *); /* Pointer to a function to be executed */
  cmdqueuenode_t *next; /* Pointer to next node */
};

typedef struct
{
  cmdqueuenode_t *head; /* Pointer to head node */
  cmdqueuenode_t *tail; /* Pointer to tail node */
} cmdqueue_t;

#define CQ_CTRL_HAS_DISPLAY      0x01
#define CQ_CTRL_HAS_DENSITY      0x02
#define CQ_CTRL_HAS_COORDS       0x04
#define CQ_CTRL_HAS_COORDS_XYZ   0x08
#define CQ_CTRL_SHOW_MOLECULE    0x01
#define CQ_CTRL_SHOW_BOX         0x02
#define CQ_CTRL_SHOW_AXES        0x04
#define CQ_CTRL_SHOW_RECTANGLE   0x08
#define CQ_CTRL_SHOW_SHAPE       0x10
#define CQ_CTRL_SHOW_LEGEND      0x20

typedef struct
{
  pthread_mutex_t mutex; /* Mutex for protection between threads */
  pthread_cond_t cond; /* Condition */
  int redisplays; /* Number of requested redisplays */
  int x; /* A point in the window */
  int y;
  double z; /* Auxiliary coordinate on virtual trackball */
  int w; /* Width of the window */
  int h; /* Height */
  int verbose;
  GLdouble center[3];
  GLdouble up[3];
  GLdouble normal[3];
  GLdouble upnormal[3]; /* The third basis vector of the world coordinates = up x normal */
  GLdouble distance;
  GLdouble size;
  GLdouble field; /* Field of vision angle (degrees) */
  GLdouble near; /* Near clipping plane */
  GLdouble far; /* Far clipping plane */
  int has; /* Flags to indicate what objects we have */
  int display; /* Flags to indicate what objects must be displayed */
  int action;
  int shftcoords; /* Flag to indicate whether coordinates need to be shifted */
  // TODO Move spherescale, etc. to a molecule structure
  GLdouble spherescale;
  int selectedatoms; /* Used for selections */
  int selection[3]; /* Identity of selected atoms */
  pid_t pid; /* Our process ID, for signals */
  legend_t legend; /* Values used to create a legend for newly-loaded data */
} ctrl_t;

// TODO Only if these functions are to be called in several files should they be declared here
// TODO Move to implementation files if they must be "private"; make static, even
void init(int argc, char **argv, index_t *idx, density_t *d);
void processCommandLine(int argc, char **argv);
void checkArgumentLength(char *arg);
void showUsageAndExit(void);
void initializeQueue(namequeue_t *q);
void pushNameToQueue(namequeue_t *q, char *name);
char *popNameFromQueue(namequeue_t *q);
int findConquestAvailability(int verbose);
int getFlag(FILE *fp, char *flag, char *value);
int findBlockLabel(FILE *fp, char *label);
int isTrue(char *value);
int isEqualStr(char *value1, char *value2);
int scanPattern(FILE *fp, char *pattern, ...);
void iluminate(void);
//void color(density_t *d, int x, int y, int z);
void prepareBoundingBox(GLuint *bb);
void prepareAxes(GLuint *a);
void prepareMol(GLuint molecule, double scale);
void prepareLegend(ctrl_t *c, density_t *d);
void resetView(ctrl_t *c);
void resetConquestInfo(cqruninfo_t *info);
void resetLimits(limits_t *l);
void setView(ctrl_t *c);
void reorientRecenter(ctrl_t *c);
void recenter(ctrl_t *c);
void makeLegend(GLuint *legend, ctrl_t *c, density_t *d);
void drawRectangle(char *image, int sx, int sy, int x0, int y0, int x1, int y1, int color);
void allocateMarks(scalemarks_t *s, int p);
void deallocateMarks(scalemarks_t *s);
void deallocateIntervals(density_t *d);
void decideMarks(scalemarks_t *marks, ctrl_t *c, density_t *d);
void makeLinearMarks(scalemarks_t *m, ctrl_t *c, density_t *d);
void makeSignedAbsoluteLogMarks(scalemarks_t *m, ctrl_t *c, density_t *d);
void makeLinearColorIntervals(ctrl_t *c, density_t *d);
void makeSignedAbsoluteLogColorIntervals(ctrl_t *c, density_t *d);
void makeSaturatedRedBlueColors(ctrl_t *c, density_t *d);
void makeRedBlueColors(ctrl_t *c, density_t *d);
void makeRainbowColors(ctrl_t *c, density_t *d);
void getColorFromScale(density_t *d, double val, rgbcolor_t *c);
void display(void);
void reshape(int w, int h);
void keyboard(unsigned char key, int x, int y);
void specialKeys(int key, int x, int y);
void writePpmImage(char *file);
int printPymolScript(index_t *idx, density_t *d, char *file);
int printPymolScriptHead(index_t *idx, density_t *d);
int printPymolScriptTail(index_t *idx, density_t *d);
int writePymolVertex(index_t *idx, double x, double y, double z);
int drawPymolQuad(double **pt, int **c, index_t *idx, density_t *d);
int drawPymolTriangle(double **pt, int **c, index_t *idx, density_t *d);
int drawPymolCorner(double **pt, int x, int y, int side, index_t *idx, density_t *d);
int drawPymolPeak(double **pt, int y, index_t *idx, density_t *d);
int colorifyPymolGradient(index_t *idx, density_t *d, int x, int y);
int colorifyPymolGradient2(index_t *idx, density_t *d, double v);
void mouse(int button, int state, int x, int y);
void processHits(int h);
void calculateDrawSlice(index_t *idx, box_t *b, density_t *d);
void populateData(index_t *idx, box_t *b, density_t *d);
double trilinearInterpolation(double *pt, density_t *d);
double findBorder(double *pt, densitycut_t *spt, double *dir, int corner, int *flag);
int drawSlice(index_t *idx, density_t *d);
int drawQuad(double **pt, index_t *idx, density_t *d, int x, int y);
int drawTriangle(double **pt, index_t *idx, density_t *d, int x, int y, int ini, int sign);
int drawCorner(double **pt, index_t *idx, density_t *d, int x, int y, int side);
int drawPeak(double **pt, index_t *idx, density_t *d, int y);
void colorifyGradient(index_t *idx, density_t *d, int x, int y);
void colorifyGradient2(index_t *idx, density_t *d, double v);
double dotProduct(double *v1, double *v2);
void resetClickCount(int v);
void mouseMotion(int x, int y);
void crossProduct(double *a, double *b, double *result);
void rotateVector(GLdouble angle, GLdouble *axis, GLdouble *v);
void setBox(density_t *d, box_t *b);
void resetDensity(density_t *d, int *old);
void resetIndex(GLdouble *normal, GLdouble *inplane, index_t *idx, density_t *d, box_t *b);
int readXplorDensity(char *name, density_t *d);
int readCqDensity(density_t *d);
int readXyzCoordinates(char *name, xyzcoords_t *c);
int readCqCoordinates(cqruninfo_t *cq, xyzcoords_t *c);
int writeCqCoordinates(cqruninfo_t *cq, density_t *d, xyzcoords_t *c, char *name);
int writeXyzCoordinates(cqruninfo_t *cq, density_t *d, xyzcoords_t *c, char *name);
int makeFilename(char *name, cqruninfo_t *info, int coreno, char *ext);
double normalize(double *v);
void *pythonConsole(void *id);
void pythreadFinalize(void);
void requestRedisplay(void);
void awakeMainThread(int v);
void checkWindowVisibility(int has);
void runCommandQueue(void);
void waitForRedisplay(ctrl_t *c);
void initializeCommands(cmdqueue_t *q);
void pushCmdToQueue(cmdqueue_t *q, int(*c)(comm_t *), comm_t*p);
cmdqueuenode_t *popCmdFromQueue(cmdqueue_t *q);
void writeQueue(cmdqueue_t *q);
comm_t *allocateCommunicator(int num, ...);
void deallocateCommunicator(comm_t *c);
/* Service functions, to be scheduled using a command queue */
/* This functions are scheduled from the python thread and executed in the openGL (main) thread */
#ifdef DEBUG
int srvPrintNumbers(comm_t *c);
#endif
int srvResetIndex(comm_t *c);
int srvWritePpmImage(comm_t *c);
int srvWritePymolSlice(comm_t *c);
//int srvWriteCqCoordinates(comm_t *c);
//int srvWriteXplorDensity(comm_t *c);
int srvRemakeCell(comm_t *c);
int srvDeleteMoleculeList(comm_t *c);
int srvCreateMoleculeList(comm_t *c);
int srvCalculateDrawSlice(comm_t *c);
int srvResetRedraw(comm_t *c);
int srvRecreateLegend(comm_t *c);

//GLuint slice;
GLuint boundingbox;
GLuint axes;
GLuint molecule;
GLuint legend;
//GLuint currentcut;

int modifiers;
static int clicks = 0; /* RECONSIDER: Does 'static' make sense here? It means file scope */

// TODO Make as many of these as possible local
cqruninfo_t cqinfo;
limits_t limits;
density_t density;
//densitycut_t densitycut;
xyzcoords_t coords;
index_t idx;
box_t box;
ctrl_t ctrl;
GLuint pickbuffer[PICKBUFSIZE];
int rendermode;
namequeue_t scripts;
cmdqueue_t pycmds;

//int flagRedisplay;

#endif /* CQ_UT_CONQTOUR_H_ */
