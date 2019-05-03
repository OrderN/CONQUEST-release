
#if defined(DARWIN) || defined(FREEBSD)
#else
#include <malloc.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

static int ctrlCset = 0;
void catch22(int sig);
static int *iwesc;
static int *nproc;
static int *taskid;

#define NUMAT 2000
#define MAXBND  2*NUMAT
#define MAXANG  3*NUMAT
#define MAXTORS 4*NUMAT
#define MXCON 10
#define MAX13 MXCON*3
#define MAX14 MXCON*9
#define NUMADD 50000
#define NUMRES 10000
static int addat = NUMADD;
#define MSAVE 7

typedef struct { float coo[NUMAT][3];
		 float ctmp[NUMAT][3];
		 float q[NUMAT];
		 float fr[NUMAT][3];
		 float ftmp[NUMAT];
		 int iaton[NUMAT];
		 int iopt[NUMAT];
		 int n13[NUMAT];
		 int i13[NUMAT][MAX13];
		 int n14[NUMAT];
		 int i14[NUMAT][MAX14];
		 int iresid[NUMAT];
		 int iconn[NUMAT][MXCON+1];
		 short int ityp[NUMAT];
		 
               } ATOMSTRU;

static ATOMSTRU *atomptr;

typedef struct { 
		  float bl[MAXBND];
		  float bk[MAXBND];
		  float ango[MAXANG];
		  float ak[MAXANG];
		  float trs1[MAXTORS][4];
		  float trs2[MAXTORS][4];
		  float trs3[MAXTORS][4];
		  float trs4[MAXTORS][4];
		  float trsi1[MAXTORS][4];
		  float trsi2[MAXTORS][4];
		  int maxbnd;
		  int maxang;
		  int maxtors;
		  int nbnd;
		  int ibnd[MAXBND][2];
		  int nang;
		  int iang[MAXANG][3];
		  int nt;
		  int it[MAXTORS][4];
		  int nti;
		  int iti[MAXTORS][4];
               } FORSTRU;

static FORSTRU *forvarptr;

typedef struct { float coot[NUMAT][3];
		 float fort[NUMAT][3];
		 float zr[NUMAT][3];
		 float y[NUMAT][3];
		 float yt[NUMAT][3];
		 float pt[NUMAT][3];
		 float s[NUMAT][3];
               } OPTSCRSTRU;

static OPTSCRSTRU *optscr;

typedef struct { float *coo;
		 float *ctmp;
		 float *q;
		 float *fr;
		 float *ftmp;
		 float *coot;
		 float *frt;
		 float *zr;
		 float *y;
		 float *yt;
		 float *pt;
		 float *s;
		 int *iaton;
		 int *iopt;
		 int *n13;
		 int *i13;
		 int *n14;
		 int *i14;
		 int *iconn;
		 int *iresid;
		 short int *ityp;
		 int *nbnd;
		 int *nang;
		 int *nt;
		 int *nti;
		 int *ibnd;
		 int *iang;
		 int *it;
		 int *iti;
		 int *watprot;
		 float *bl;
		 float *bk;
		 float *ango;
		 float *ak;
		 float *trs1;
		 float *trs2;
		 float *trs3;
		 float *trs4;
		 float *trsi1;
		 float *trsi2;
		 float *work;
		 int *mxnat;
		 int *iatoms;
		 int mxorg;
		 int *maxbnd;
		 int *maxang;
		 int *maxtors;
               } COOSTRU;

static COOSTRU xyz;
static COOSTRU TMPxyz;

typedef struct { float *v;
		 float *a;
		 float *m;
               } MDSTRU;

static MDSTRU md;
	
#define MXNEIB 200
static int *nlst = NULL;
static int *lst = NULL;

static float water[648][3] = {
      3.009,    6.289,    1.536,
      2.346,    5.879,    2.126,
      2.541,    6.968,    1.000,
      2.147,    3.888,   -7.775,
      2.450,    3.249,   -7.100,
      1.818,    3.297,   -8.479,
     -0.200,    3.894,    7.128,
     -0.006,    3.280,    7.875,
     -0.614,    3.297,    6.466,
      5.038,   -5.416,   -6.189,
      4.421,   -5.875,   -6.806,
      5.509,   -4.776,   -6.766,
     -2.797,   -3.106,    7.212,
     -3.242,   -3.604,    7.924,
     -3.295,   -3.471,    6.443,
     -1.191,    6.280,    8.994,
     -0.626,    6.780,    8.379,
     -0.577,    5.706,    9.484,
     -7.041,    6.418,    7.260,
     -6.964,    6.968,    6.458,
     -7.539,    6.996,    7.871,
     -1.185,    6.280,    1.920,
     -1.048,    6.013,    0.967,
     -1.614,    5.510,    2.316,
      7.841,   -7.248,   -7.751,
      7.528,   -6.641,   -7.072,
      7.873,   -6.639,   -8.549,
      8.126,    8.499,   -0.316,
      8.639,    9.076,    0.293,
      8.789,    8.211,   -0.977,
      6.718,   -8.969,    7.092,
      7.277,   -8.669,    7.851,
      6.069,   -8.242,    6.961,
      6.637,   -7.256,   -0.316,
      7.087,   -6.641,    0.286,
      6.138,   -6.678,   -0.932,
      3.374,   -4.346,   -8.071,
      2.557,   -4.509,   -8.589,
      3.955,   -3.844,   -8.702,
     -4.259,   -3.031,   -5.785,
     -3.817,   -2.309,   -5.289,
     -4.916,   -3.372,   -5.128,
      4.833,    3.894,   -4.259,
      5.485,    3.351,   -4.793,
      4.429,    3.239,   -3.659,
     -3.024,    4.023,    7.219,
     -3.231,    4.544,    6.429,
     -3.175,    4.678,    7.908,
      1.969,   -0.811,    4.987,
      2.544,   -0.513,    5.719,
      2.542,   -1.352,    4.409,
     -6.207,   -5.387,   -1.870,
     -6.566,   -4.648,   -1.357,
     -6.767,   -6.124,   -1.549,
      2.034,   -4.164,    5.170,
      2.544,   -3.631,    5.784,
      1.456,   -4.609,    5.815,
      3.103,   -7.256,    9.173,
      2.745,   -8.058,    9.597,
      2.629,   -6.554,    9.649,
      6.513,   -0.807,   -4.193,
      6.069,   -1.529,   -3.707,
      6.069,   -0.028,   -3.792,
      5.208,   -8.757,    1.476,
      5.684,   -8.374,    0.696,
      5.684,   -8.284,    2.204,
     -2.684,    1.973,   -4.268,
     -3.242,    2.479,   -4.894,
     -3.292,    1.512,   -3.659,
      0.641,    6.289,    5.157,
      0.387,    5.519,    5.701,
      0.455,    7.002,    5.787,
     -0.219,   -0.811,    3.176,
      0.318,   -0.964,    3.993,
      0.277,   -0.123,    2.696,
      0.944,   -3.022,    1.825,
      0.397,   -3.805,    2.059,
      0.466,   -2.326,    2.323,
     -1.193,    9.784,    7.092,
     -0.614,    9.289,    6.474,
     -1.728,   10.288,    6.446,
     -1.494,   -7.248,    5.042,
     -1.612,   -6.678,    4.262,
     -1.881,   -6.678,    5.743,
      1.018,    7.300,   -2.181,
      1.770,    7.611,   -2.708,
      0.316,    7.411,   -2.844,
      6.513,   -5.358,    7.239,
      5.906,   -6.056,    7.581,
      6.069,   -5.005,    6.443,
     -5.868,   -5.591,   -4.383,
     -5.336,   -6.135,   -4.985,
     -5.588,   -5.940,   -3.504,
      8.119,   -0.753,    5.056,
      8.439,   -0.120,    4.389,
      7.582,   -0.184,    5.651,
      1.015,   -8.918,   -6.247,
      1.828,   -9.338,   -6.614,
      0.334,   -9.361,   -6.768,
     -2.931,   -8.757,    1.666,
     -3.242,   -9.488,    2.215,
     -3.178,   -7.998,    2.215,
     -6.310,   -0.786,   -2.263,
     -6.575,   -0.054,   -2.855,
     -6.575,   -1.533,   -2.843,
      3.332,    2.063,   -0.316,
      3.947,    2.756,   -0.683,
      2.736,    2.601,    0.253,
     -4.477,   -5.167,   -9.459,
     -5.033,   -4.539,   -8.944,
     -4.889,   -6.010,   -9.164,
     -2.581,   -0.811,    5.051,
     -2.214,   -1.362,    4.325,
     -2.259,   -1.311,    5.830,
      5.146,    7.419,    3.176,
      5.672,    7.841,    2.483,
      5.626,    7.654,    3.989,
     -5.993,    3.894,   -8.722,
     -6.767,    3.313,   -8.902,
     -5.289,    3.223,   -8.655,
     -2.566,   -2.739,    3.134,
     -1.899,   -2.343,    2.544,
     -3.324,   -2.824,    2.504,
      6.717,    6.280,   -4.324,
      6.974,    5.847,   -3.485,
      7.050,    5.659,   -4.985,
     -4.198,    5.067,   -2.051,
     -3.626,    5.828,   -1.882,
     -3.747,    4.452,   -1.460,
     -8.592,    8.506,   -2.263,
     -8.754,    7.762,   -2.865,
     -8.923,    9.233,   -2.836,
     -8.410,    3.946,   -6.152,
     -8.994,    4.454,   -5.562,
     -8.982,    3.300,   -6.612,
     -8.366,    6.289,    4.987,
     -8.933,    6.738,    5.651,
     -8.998,    5.815,    4.408,
      1.018,    8.499,    7.264,
      1.666,    8.044,    7.851,
      0.522,    9.084,    7.873,
      1.018,   -0.580,   -2.447,
      0.785,   -0.484,   -1.476,
      0.334,   -0.051,   -2.870,
      9.085,    3.894,    1.557,
      9.548,    3.278,    0.957,
      8.697,    3.280,    2.206,
     -1.191,    6.280,   -0.553,
     -1.659,    6.984,   -1.062,
     -1.249,    5.534,   -1.182,
     -4.347,   -7.302,    4.987,
     -3.821,   -7.868,    4.389,
     -4.916,   -6.752,    4.389,
      5.146,    2.054,    8.732,
      5.285,    2.399,    9.646,
      5.548,    2.758,    8.200,
     -8.295,   -4.318,    5.042,
     -8.853,   -3.844,    4.371,
     -8.913,   -4.499,    5.770,
      3.443,   -0.811,    9.276,
      4.051,   -1.457,    8.872,
      4.029,   -0.126,    9.644,
     -0.247,   -7.337,    8.994,
     -0.608,   -7.890,    8.246,
     -0.614,   -7.843,    9.747,
     -6.301,   -3.031,   -4.259,
     -6.287,   -3.649,   -3.511,
     -6.767,   -3.582,   -4.921,
      3.210,   -8.998,    3.089,
      4.016,   -9.166,    2.543,
      2.544,   -9.432,    2.543,
     -4.272,    7.635,    7.218,
     -4.916,    7.183,    7.805,
     -3.555,    7.787,    7.857,
      4.759,   -5.549,   -2.071,
      4.269,   -5.794,   -2.885,
      4.226,   -6.024,   -1.409,
     -7.114,   -0.895,    1.536,
     -6.860,   -1.431,    2.307,
     -6.769,   -1.431,    0.797,
      1.025,    2.061,   -3.850,
      0.277,    2.633,   -3.591,
      1.695,    2.451,   -3.252,
     -2.797,   -5.416,    1.278,
     -3.242,   -4.585,    1.051,
     -3.230,   -6.033,    0.670,
      5.078,   -0.811,    3.062,
      5.504,   -0.478,    3.896,
      5.716,   -1.431,    2.655,
     -4.289,    1.759,    4.869,
     -5.022,    1.340,    4.389,
     -3.525,    1.359,    4.417,
      2.036,   -0.811,   -6.158,
      2.491,   -1.060,   -5.317,
      2.478,   -1.431,   -6.768,
      9.139,   -0.855,    8.994,
      9.544,   -1.462,    9.649,
      8.697,   -1.494,    8.392,
      3.332,    6.289,    4.987,
      3.869,    5.725,    4.389,
      2.749,    6.800,    4.389,
     -2.652,    2.054,   -0.646,
     -2.740,    2.161,    0.345,
     -3.242,    1.317,   -0.855,
      6.697,    1.824,    1.239,
      7.516,    1.314,    1.088,
      6.069,    1.400,    0.633,
      8.194,    2.054,   -2.004,
      8.652,    1.414,   -1.422,
      7.696,    2.624,   -1.395,
      9.064,   -3.090,   -7.486,
      8.308,   -3.582,   -7.106,
      9.760,   -3.639,   -7.104,
     -5.865,    6.280,   -4.268,
     -5.281,    5.828,   -3.610,
     -5.282,    6.972,   -4.673,
      8.119,   -5.358,    9.159,
      8.703,   -4.655,    9.483,
      7.706,   -4.975,    8.341,
      9.064,    5.265,    9.068,
      8.777,    4.720,    8.288,
      9.626,    4.664,    9.576,
      4.970,    0.553,    5.157,
      4.394,    0.840,    5.902,
      4.600,    1.068,    4.389,
     -4.304,   -7.256,   -0.152,
     -3.809,   -7.919,    0.404,
     -4.916,   -6.787,    0.455,
      3.290,    0.341,    1.575,
      3.799,    0.884,    0.933,
      3.975,   -0.027,    2.204,
     -9.584,    2.063,    4.891,
    -10.319,    2.486,    4.417,
     -8.845,    2.540,    4.486,
     -8.401,    6.280,   -4.140,
     -8.939,    5.690,   -3.589,
     -7.492,    5.938,   -3.954,
     -6.033,   -2.920,    7.239,
     -6.566,   -2.343,    6.648,
     -5.505,   -2.314,    7.805,
      4.833,   -7.337,    7.205,
      4.277,   -7.666,    7.954,
      4.277,   -7.632,    6.443,
     -6.016,   -0.558,    5.215,
     -5.266,   -0.178,    5.719,
     -6.767,   -0.182,    5.721,
     -8.293,    5.344,   -0.243,
     -8.992,    5.154,    0.408,
     -7.577,    5.661,    0.336,
      1.024,    6.094,   -4.390,
      1.769,    5.828,   -4.957,
      1.270,    5.727,   -3.523,
      5.109,    2.054,   -7.383,
      5.547,    1.362,   -6.856,
      4.291,    2.248,   -6.856,
      8.120,    3.940,    7.050,
      7.571,    4.468,    6.446,
      8.559,    3.280,    6.469,
     -8.293,    0.312,    7.100,
     -8.920,    0.113,    6.390,
     -8.781,   -0.052,    7.889,
     -4.272,   -7.567,    8.767,
     -4.947,   -7.962,    8.166,
     -3.858,   -6.839,    8.239,
     -4.361,    6.200,   -9.292,
     -5.079,    5.657,   -8.902,
     -3.575,    5.766,   -8.906,
     -0.226,    8.499,    3.219,
      0.228,    7.879,    3.832,
     -0.591,    7.879,    2.543,
      4.950,   -5.295,    1.374,
      4.400,   -5.193,    2.204,
      4.364,   -4.885,    0.714,
     -1.470,   -5.416,    7.014,
     -0.719,   -5.791,    6.500,
     -1.728,   -4.632,    6.489,
     -4.321,    6.280,    1.476,
     -4.882,    5.729,    2.079,
     -3.679,    5.664,    1.084,
      4.915,    2.061,    3.104,
      5.688,    1.857,    2.523,
      4.394,    2.669,    2.543,
     -0.849,   -5.416,   -2.240,
     -0.614,   -6.164,   -2.810,
     -0.726,   -4.648,   -2.812,
     -8.286,   -5.416,   -5.845,
     -9.076,   -5.607,   -5.289,
     -7.615,   -5.912,   -5.320,
     -2.613,   -2.010,   -2.071,
     -2.260,   -1.548,   -2.865,
     -1.792,   -2.308,   -1.627,
     -7.324,   -7.252,   -7.809,
     -7.610,   -6.654,   -7.096,
     -7.485,   -6.678,   -8.596,
     -2.674,   -7.256,   -7.775,
     -3.026,   -6.641,   -7.096,
     -3.223,   -7.012,   -8.558,
      9.092,    0.553,   -4.153,
      9.626,    0.717,   -4.973,
      9.543,    1.162,   -3.540,
     -8.366,   -7.256,   -1.853,
     -8.923,   -6.552,   -1.453,
     -8.845,   -8.034,   -1.516,
      6.658,   -7.615,    3.465,
      5.937,   -7.910,    4.043,
      7.419,   -7.985,    3.951,
      3.424,    8.708,   -6.135,
      3.902,    9.127,   -5.382,
      3.701,    7.778,   -6.001,
     -5.868,    3.946,   -6.019,
     -5.454,    4.663,   -5.491,
     -6.649,    4.394,   -6.445,
      1.018,   -3.031,   -4.351,
      0.619,   -2.338,   -4.922,
      1.694,   -3.482,   -4.921,
      3.317,   -7.256,    5.004,
      3.750,   -6.643,    4.371,
      2.888,   -7.910,    4.389,
      5.192,    4.964,    8.998,
      5.311,    5.495,    8.164,
      5.745,    5.467,    9.644,
     -5.993,   -5.287,    3.841,
     -5.276,   -4.703,    4.167,
     -6.767,   -4.858,    4.286,
     -5.980,   -3.967,   -7.654,
     -5.282,   -3.655,   -7.028,
     -6.724,   -4.231,   -7.095,
     -2.686,    6.359,    5.042,
     -3.178,    6.968,    5.651,
     -2.237,    6.968,    4.414,
      3.443,   -3.032,   -5.950,
      3.593,   -3.587,   -6.768,
      3.857,   -3.620,   -5.281,
      3.573,    6.289,   -2.475,
      3.957,    7.109,   -2.844,
      4.198,    5.653,   -2.852,
     -2.652,    7.485,   -2.271,
     -3.402,    7.825,   -2.810,
     -1.917,    7.762,   -2.867,
     -8.286,    0.553,   -0.724,
     -7.755,    1.328,   -0.975,
     -7.686,   -0.146,   -1.041,
      4.896,    0.473,   -2.106,
      4.277,    0.151,   -2.810,
      4.290,    0.943,   -1.483,
      6.513,    8.730,   -7.834,
      6.118,    9.184,   -8.603,
      6.079,    9.184,   -7.099,
     -4.614,    8.539,   -8.033,
     -4.889,    9.266,   -8.621,
     -4.772,    7.762,   -8.593,
     -6.020,    6.026,    3.465,
     -5.295,    5.661,    4.006,
     -6.769,    5.671,    3.989,
      3.044,    2.054,   -5.845,
      2.544,    2.669,   -5.281,
      2.662,    1.208,   -5.535,
      0.718,   -5.364,   -5.711,
      0.525,   -6.145,   -5.134,
      0.387,   -4.632,   -5.170,
      6.839,    6.205,   -8.033,
      7.487,    5.729,   -8.606,
      6.988,    7.146,   -8.289,
     -1.186,   -8.931,   -0.331,
     -1.773,   -9.165,   -1.068,
     -1.728,   -9.168,    0.444,
      3.135,    0.553,    7.070,
      2.450,    1.035,    6.578,
      2.662,    0.144,    7.841,
     -2.613,    3.946,    3.491,
     -2.070,    3.311,    4.029,
     -2.760,    4.678,    4.125,
     -8.286,    8.212,   -9.462,
     -7.541,    7.883,   -8.933,
     -9.025,    7.762,   -9.024,
     -8.293,   -3.030,   -2.208,
     -8.920,   -2.573,   -2.836,
     -7.896,   -3.653,   -2.836,
     -0.178,   -3.030,   -2.071,
      0.347,   -3.297,   -2.865,
      0.157,   -3.649,   -1.400,
      3.374,    8.503,   -0.316,
      4.004,    8.964,    0.286,
      3.947,    7.905,   -0.843,
     -2.613,    2.054,    1.825,
     -3.243,    1.342,    2.094,
     -2.870,    2.800,    2.417,
     -5.940,    8.499,   -0.217,
     -5.363,    7.879,    0.286,
     -5.339,    8.826,   -0.937,
     -7.761,   -5.364,    8.895,
     -7.855,   -6.083,    8.223,
     -7.468,   -4.632,    8.329,
      3.443,    8.499,   -9.654,
      3.946,    8.038,  -10.363,
      3.640,    7.879,   -8.906,
      3.367,    6.250,   -8.448,
      2.555,    5.726,   -8.626,
      4.030,    5.671,   -8.902,
      5.125,    3.952,   -0.463,
      5.276,    4.426,    0.408,
      5.510,    4.603,   -1.074,
     -7.179,    3.952,   -2.071,
     -7.755,    4.452,   -1.440,
     -6.381,    4.465,   -1.890,
     -7.107,    2.058,   -4.254,
     -6.955,    2.669,   -3.507,
     -6.861,    2.633,   -5.008,
     -0.074,   -5.416,    1.536,
     -0.614,   -5.973,    0.934,
      0.236,   -6.063,    2.215,
      0.944,   -0.779,    9.069,
      0.337,   -1.548,    9.005,
      1.750,   -1.174,    9.468,
      5.082,   -3.030,    8.994,
      5.722,   -2.337,    9.294,
      5.626,   -3.585,    8.390,
     -8.348,    8.482,   -6.135,
     -7.856,    7.879,   -5.529,
     -8.902,    7.879,   -6.661,
      0.944,    6.280,   -7.485,
      1.771,    5.925,   -7.106,
      0.316,    5.703,   -7.046,
      2.132,    3.952,   -2.182,
      2.542,    4.607,   -2.778,
      1.825,    4.535,   -1.460,
     -4.566,   -9.021,   -2.071,
     -5.020,   -8.283,   -1.612,
     -5.084,   -9.225,   -2.867,
      8.119,    6.221,    3.123,
      8.581,    5.656,    2.483,
      8.060,    7.076,    2.650,
     -2.592,    6.195,   -6.129,
     -2.289,    5.729,   -5.321,
     -2.201,    5.661,   -6.827,
      3.103,   -3.019,   -1.943,
      2.583,   -2.316,   -1.506,
      2.736,   -3.806,   -1.506,
      6.697,    3.946,   -6.207,
      6.136,    3.351,   -6.768,
      6.988,    4.641,   -6.827,
      9.064,   -2.261,   -4.259,
      8.697,   -2.201,   -5.143,
      8.752,   -1.431,   -3.860,
      5.112,    4.993,    1.838,
      5.695,    5.471,    2.454,
      4.290,    5.552,    1.896,
      1.024,   -7.340,   -4.238,
      1.283,   -7.894,   -3.449,
      1.223,   -7.982,   -4.975,
     -5.937,   -8.757,   -6.135,
     -5.356,   -9.251,   -6.768,
     -6.652,   -8.426,   -6.740,
     -5.978,    8.506,    3.176,
     -6.498,    7.882,    3.732,
     -5.451,    7.879,    2.631,
     -0.842,   -1.078,   -6.281,
     -0.008,   -1.268,   -6.766,
     -1.516,   -1.457,   -6.885,
     -6.171,   -0.810,   -7.737,
     -5.531,   -1.227,   -8.378,
     -6.535,   -1.548,   -7.224,
      0.944,   -5.416,   -8.138,
      0.695,   -5.639,   -7.187,
      0.387,   -6.052,   -8.640,
      6.697,   -2.010,   -1.968,
      6.069,   -1.431,   -1.501,
      7.109,   -1.457,   -2.664,
      9.069,    0.379,    1.541,
      9.588,    0.172,    0.721,
      9.644,   -0.021,    2.214,
     -1.191,    1.974,    5.042,
     -0.730,    1.400,    4.391,
     -1.614,    1.348,    5.662,
     -4.361,   -3.022,    5.042,
     -3.793,   -2.510,    4.408,
     -4.909,   -2.343,    5.493,
      1.018,   -0.786,   -0.034,
      0.555,   -1.431,    0.529,
      1.693,   -0.370,    0.539,
     -5.868,    9.877,    6.945,
     -5.365,    9.168,    6.500,
     -5.497,   10.635,    6.443,
      8.119,   -5.397,   -4.323,
      7.651,   -4.915,   -3.589,
      7.429,   -5.959,   -4.718,
      6.697,   -3.022,    1.333,
      7.056,   -2.481,    0.616,
      6.080,   -3.649,    0.902,
     -2.593,   -0.807,   -4.419,
     -2.030,   -1.431,   -4.921,
     -2.317,    0.004,   -4.894,
     -0.054,    1.973,   -9.635,
      0.316,    1.325,   -9.013,
     -0.661,    1.416,  -10.159,
      3.127,    3.849,    7.142,
      2.542,    3.517,    6.446,
      2.830,    3.326,    7.904,
     -6.220,   -3.030,   -0.316,
     -6.767,   -3.609,    0.263,
     -6.833,   -2.628,   -0.970,
      8.438,   -8.757,   -4.351,
      8.710,   -9.471,   -4.939,
      8.697,   -8.008,   -4.903,
      9.078,    8.585,    3.176,
      9.648,    9.123,    2.570,
      9.673,    8.084,    3.761,
     -4.232,   -3.030,    1.277,
     -4.916,   -3.482,    0.717,
     -3.877,   -2.341,    0.672,
      3.289,   -5.461,   -4.323,
      4.028,   -5.744,   -4.921,
      2.544,   -5.959,   -4.731,
     -2.652,   -5.364,   -5.885,
     -2.224,   -5.996,   -5.277,
     -3.193,   -4.775,   -5.317,
     -8.293,   -7.337,    7.239,
     -8.872,   -7.923,    7.786,
     -7.852,   -7.950,    6.612,
     -8.293,    2.008,   -9.546,
     -7.927,    1.457,  -10.286,
     -8.855,    1.387,   -9.053,
      6.681,   -4.100,    4.959,
      7.280,   -4.569,    4.325,
      6.016,   -3.758,    4.344,
     -0.862,    3.946,   -4.460,
     -0.521,    4.722,   -4.921,
     -0.701,    3.210,   -5.100,
      3.136,   -5.194,    3.176,
      2.544,   -4.764,    2.536,
      2.715,   -4.842,    4.021,
      1.018,   -8.721,   -2.047,
      0.319,   -8.440,   -1.418,
      1.554,   -9.306,   -1.460,
      6.717,    1.990,   -4.194,
      7.393,    1.302,   -4.028,
      7.067,    2.669,   -3.591,
      1.018,    3.894,    1.227,
      0.387,    3.375,    0.684,
      1.258,    4.655,    0.677,
     -4.374,    8.524,   -4.430,
     -5.033,    9.086,   -4.921,
     -3.572,    8.689,   -4.985,
      1.018,    1.973,    3.126,
      0.628,    2.667,    2.543,
      1.632,    1.491,    2.543,
     -6.175,   -5.364,    1.278,
     -6.819,   -6.033,    1.029,
     -6.427,   -5.134,    2.214,
      5.195,    8.524,   -4.323,
      5.745,    9.258,   -4.638,
      5.745,    7.762,   -4.610,
     -4.165,   -0.811,   -9.464,
     -3.684,   -0.003,   -9.215,
     -3.613,   -1.485,   -9.024,
      6.743,    6.280,   -1.743,
      7.021,    7.158,   -1.406,
      7.330,    5.690,   -1.252,
     -2.690,   -7.249,   -2.218,
     -3.193,   -7.035,   -1.402,
     -3.173,   -6.707,   -2.852,
      3.331,   -0.811,   -3.907,
      2.475,   -1.004,   -3.450,
      3.863,   -1.548,   -3.557,
     -2.492,    8.752,   -6.212,
     -2.336,    7.854,   -6.594,
     -2.322,    9.324,   -6.976,
     -5.870,    3.876,    4.876,
     -5.306,    3.278,    4.352,
     -5.660,    3.520,    5.787,
      5.021,    6.571,    7.047,
      5.681,    7.001,    6.474,
      4.240,    6.537,    6.442,
     -1.470,    3.894,   -2.008,
     -1.584,    3.455,   -2.881,
     -1.794,    3.180,   -1.406,
     -1.064,    8.524,   -4.120,
     -0.615,    9.259,   -3.659,
     -1.479,    8.940,   -4.920,
      9.069,   -5.358,   -0.316,
      8.789,   -4.819,    0.441,
      8.697,   -4.855,   -1.070,
      6.513,   -0.725,   -9.626,
      6.069,   -0.105,  -10.242,
      7.096,   -0.144,   -9.100,
     -4.291,   -0.807,   -0.515,
     -5.020,   -0.622,   -1.152,
     -3.667,   -1.355,   -1.079,
      8.451,   -8.757,    8.994,
      9.249,   -9.285,    9.240,
      8.309,   -8.233,    9.822,
     -8.286,   -8.837,    1.105,
     -7.986,   -8.008,    0.717,
     -7.615,   -9.438,    0.717,
     -0.241,   -3.050,    8.994,
     -0.582,   -3.639,    8.291,
      0.234,   -3.658,    9.597,
      2.182,   -5.358,   -0.206,
      1.944,   -6.157,    0.286,
      1.734,   -4.681,    0.334,
      1.058,   -5.433,    7.220,
      0.557,   -6.013,    7.851,
      1.624,   -6.061,    6.713,
      7.031,   -4.531,   -2.071,
      6.069,   -4.746,   -1.936,
      7.079,   -3.621,   -1.730,
     -0.241,    1.973,   -6.143,
      0.234,    1.416,   -5.482,
     -0.608,    1.328,   -6.768,
      4.942,    3.894,    5.051,
      4.624,    3.278,    4.361,
      4.394,    3.634,    5.825,
     -2.480,   -3.031,   -7.834,
     -2.201,   -3.582,   -8.593,
     -2.322,   -3.646,   -7.095,
      3.443,    6.112,   -5.887,
      4.054,    5.661,   -5.279,
      3.804,    5.828,   -6.768,
      8.089,   -3.022,    7.047,
      8.061,   -2.249,    6.458,
      7.582,   -3.658,    6.491,
     -4.477,    2.054,   -7.485,
     -4.916,    1.268,   -7.116,
     -4.962,    2.757,   -6.987,
     -8.293,    0.392,   -6.345,
     -8.993,   -0.051,   -6.828,
     -7.490,    0.079,   -6.807,
      8.126,   -5.416,    3.176,
      7.582,   -6.033,    2.643,
      8.616,   -6.015,    3.776,
      6.658,   -2.921,   -6.135,
      6.118,   -2.342,   -5.555,
      7.050,   -2.291,   -6.768,
     -4.477,    0.379,    1.718,
     -4.661,   -0.029,    0.831,
     -5.111,   -0.077,    2.282,
     -5.993,    3.952,    7.282,
     -6.649,    4.585,    6.893,
     -5.506,    4.519,    7.904,
     -2.797,    2.054,    8.980,
     -3.242,    2.449,    9.752,
     -3.101,    2.633,    8.237,
      0.944,   -7.337,    3.176,
      1.451,   -7.951,    2.613,
      0.439,   -7.941,    3.774,
      9.085,   -2.745,    3.291,
      8.696,   -2.343,    2.504,
      8.697,   -2.223,    4.021
};

static int *listw = NULL;
static float *qpot = NULL;
static float *woo = NULL;
static float *whh = NULL;
static float *woh = NULL;
static float *dewoo = NULL;
static float *dewhh = NULL;
static float *dewoh = NULL;
static float *apprerf = NULL;
static float *appexp = NULL;
static float *ewoo = NULL;
static float *ewhh = NULL;
static float *ewoh = NULL;
static float *edwoo = NULL;
static float *edwhh = NULL;
static float *edwoh = NULL;
static float *drb = NULL;
static float *rsolv = NULL;
static float *rborn = NULL;
static float *shct = NULL;

void catch22(int sig)
{

	*iwesc = 1;
	(void) signal(SIGINT, catch22);

}


#if defined(VMS) || defined(UNDERSC)
void wexit()
#else
#ifdef CRAY
void WEXIT()
#else
void wexit_()
#endif
#endif
{
#ifndef MD

#if defined(VMS) || defined(UNDERSC)
	wrtesc();
#else
#ifdef CRAY
	WRTESC();
#else
	wrtesc_();
#endif
#endif

#else

#if defined(VMS) || defined(UNDERSC)
	wmdesc();
#else
#ifdef CRAY
	WMDESC();
#else
	wmdesc_();
#endif
#endif

#endif
	exit(1);
}

#if defined(VMS) || defined(UNDERSC)
parptr(nptr, fptr, ffptr, nitem)
#else
#ifdef CRAY
PARPTR(nptr, fptr, ffptr, nitem)
#else
parptr_(nptr, fptr, ffptr, nitem)
#endif
#endif

int *nptr;
float *fptr;
float *ffptr;
int *nitem;
{

    switch(*nptr) {
    case 0: atomptr   = (ATOMSTRU *) fptr; 
	    xyz.coo   = (float *) atomptr->coo;
	    xyz.ctmp   = (float *) atomptr->ctmp;
	    xyz.q     = (float *) atomptr->q;
	    xyz.fr   = (float *) atomptr->fr;
	    xyz.ftmp  = (float *) atomptr->ftmp;
	    xyz.iaton = (int *) atomptr->iaton;
	    xyz.iopt  = (int *) atomptr->iopt;
	    xyz.n13   = (int *) atomptr->n13;
	    xyz.i13   = (int *) atomptr->i13;
	    xyz.n14   = (int *) atomptr->n14;
	    xyz.i14   = (int *) atomptr->i14;
	    xyz.iresid= (int *) atomptr->iresid;
	    xyz.iconn = (int *) atomptr->iconn;
	    xyz.ityp  = (short int *) atomptr->ityp;
	    break;
    case 1: forvarptr = (FORSTRU *) fptr; 
	    xyz.maxbnd =  (int *) &forvarptr->maxbnd;
	    xyz.maxang =  (int *) &forvarptr->maxang;
	    xyz.maxtors = (int *) &forvarptr->maxtors;
	    xyz.nbnd  = (int *) &forvarptr->nbnd;
	    xyz.nang  = (int *) &forvarptr->nang;
	    xyz.nt    = (int *) &forvarptr->nt;
	    xyz.nti   = (int *) &forvarptr->nti;
	    xyz.ibnd  = (int *) forvarptr->ibnd;
	    xyz.iang  = (int *) forvarptr->iang;
	    xyz.it    = (int *) forvarptr->it;
	    xyz.iti   = (int *) forvarptr->iti;
	    xyz.bl    = (float *) forvarptr->bl;
	    xyz.bk    = (float *) forvarptr->bk;
	    xyz.ango  = (float *) forvarptr->ango;
	    xyz.ak    = (float *) forvarptr->ak;
	    xyz.trs1  = (float *) forvarptr->trs1;
	    xyz.trs2  = (float *) forvarptr->trs2;
	    xyz.trs3  = (float *) forvarptr->trs3;
	    xyz.trs4  = (float *) forvarptr->trs4;
	    xyz.trsi1 = (float *) forvarptr->trsi1;
	    xyz.trsi2 = (float *) forvarptr->trsi2;
	    break;
    case 2: 
	    xyz.iatoms = nitem;
	    xyz.mxorg = NUMAT;
	    break;
    case 3: optscr = (OPTSCRSTRU *) fptr; 
	    xyz.coot  = (float *) optscr->coot;
	    xyz.frt   = (float *) optscr->fort;
	    xyz.zr    = (float *) optscr->zr;
	    xyz.y     = (float *) optscr->y;
	    xyz.yt    = (float *) optscr->yt;
	    xyz.pt    = (float *) optscr->pt;
	    xyz.s     = (float *) optscr->s;
	    xyz.work  = NULL;
	    break;
    case 4: 
	    xyz.mxnat = nitem;
	    break;
    case 5: 
	    iwesc = nitem;
	    break;
    case 6: 
	    nproc = nitem;
	    break;
    case 7: 
	    taskid = nitem;
	    break;
    }
    if (!ctrlCset) {
	(void) signal(SIGINT, catch22);
	ctrlCset = 1;
        unlink("ambfor.wr");
    }
}

#if defined(VMS) || defined(UNDERSC)
allcoo(ZSizep)
#else
#ifdef CRAY
ALLCOO(ZSizep)
#else
allcoo_(ZSizep)
#endif
#endif
int *ZSizep;
{
   int memstat;
   float d;
   int i;
   short int j;
   int ZSize, wsize;

   ZSize = *xyz.mxnat + *ZSizep;
   memstat = 1;

   TMPxyz = xyz;

   if ((xyz.coo = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   if ((xyz.ctmp = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   if ((xyz.q = (float *) malloc((sizeof d)*ZSize)) == NULL) {
	memstat = 0;
   }

   if ((xyz.fr = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   if ((xyz.coot = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   if ((xyz.frt = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   if ((xyz.zr = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   if ((xyz.y = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   if ((xyz.yt = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   if ((xyz.pt = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   if ((xyz.s = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   wsize = ZSize*3*(2*MSAVE+1)+2*MSAVE;

   if ((xyz.work = (float *) malloc((sizeof d)*wsize)) == NULL) {
	memstat = 0;
   }

   if ((xyz.iaton = (int *) malloc((sizeof i)*ZSize)) == NULL) {
	memstat = 0;
   }

   if ((xyz.iopt = (int *) malloc((sizeof i)*ZSize)) == NULL) {
	memstat = 0;
   }

   if ((xyz.watprot = (int *) malloc((sizeof i)*ZSize)) == NULL) {
	memstat = 0;
   }

   if ((xyz.ftmp = (float *) malloc((sizeof d)*ZSize)) == NULL) {
	memstat = 0;
   }

   if ((xyz.n13 = (int *) malloc((sizeof i)*ZSize)) == NULL) {
	memstat = 0;
   }

   if ((xyz.i13 = (int *) malloc((sizeof i)*ZSize*(MAX13))) == NULL) {
	memstat = 0;
   }

   if ((xyz.n14 = (int *) malloc((sizeof i)*ZSize)) == NULL) {
	memstat = 0;
   }

   if ((xyz.i14 = (int *) malloc((sizeof i)*ZSize*(MAX14))) == NULL) {
	memstat = 0;
   }

   if ((xyz.iresid = (int *) malloc((sizeof i)*ZSize)) == NULL) {
	memstat = 0;
   }

   if ((xyz.iconn = (int *) malloc((sizeof i)*ZSize*(MXCON+1))) == NULL) {
	memstat = 0;
   }

   if ((xyz.ityp = (short int *) malloc((sizeof j)*ZSize)) == NULL) {
	memstat = 0;
   }

   if ((xyz.ibnd = (int *) malloc((sizeof i)*ZSize*2*2)) == NULL) {
	memstat = 0;
   }

   if ((xyz.iang = (int *) malloc((sizeof i)*ZSize*3*3)) == NULL) {
	memstat = 0;
   }

   if ((xyz.it = (int *) malloc((sizeof i)*ZSize*4*4)) == NULL) {
	memstat = 0;
   }

   if ((xyz.iti = (int *) malloc((sizeof i)*ZSize*4*4)) == NULL) {
	memstat = 0;
   }

   if ((xyz.bl = (float *) malloc((sizeof d)*ZSize*2)) == NULL) {
	memstat = 0;
   }

   if ((xyz.bk = (float *) malloc((sizeof d)*ZSize*2)) == NULL) {
	memstat = 0;
   }

   if ((xyz.ango = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   if ((xyz.ak = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   if ((xyz.trs1 = (float *) malloc((sizeof d)*ZSize*4*4)) == NULL) {
	memstat = 0;
   }

   if ((xyz.trs2 = (float *) malloc((sizeof d)*ZSize*4*4)) == NULL) {
	memstat = 0;
   }

   if ((xyz.trs3 = (float *) malloc((sizeof d)*ZSize*4*4)) == NULL) {
	memstat = 0;
   }

   if ((xyz.trs4 = (float *) malloc((sizeof d)*ZSize*4*4)) == NULL) {
	memstat = 0;
   }

   if ((xyz.trsi1 = (float *) malloc((sizeof d)*ZSize*4*4)) == NULL) {
	memstat = 0;
   }

   if ((xyz.trsi2 = (float *) malloc((sizeof d)*ZSize*4*4)) == NULL) {
	memstat = 0;
   }


   if (!memstat) {
	fprintf(stderr,"Out of memory\n");
	xyz = TMPxyz;
   } else {
	for (i=0; i < *xyz.mxnat; i++) {
	   for (j=0; j<3; j++) {
		xyz.coo[i*3+j] = TMPxyz.coo[i*3+j];
		xyz.ctmp[i*3+j] = TMPxyz.ctmp[i*3+j];
		xyz.fr[i*3+j] = TMPxyz.fr[i*3+j];
	   }
	   xyz.q[i] = TMPxyz.q[i];
	   xyz.iaton[i] = TMPxyz.iaton[i];
	   xyz.iopt[i] = TMPxyz.iopt[i];
	   xyz.n13[i] = TMPxyz.n13[i];
	   xyz.n14[i] = TMPxyz.n14[i];
	   xyz.iresid[i] = TMPxyz.iresid[i];
	   xyz.ityp[i] = TMPxyz.ityp[i];
	   for (j=0; j<MXCON; j++)
		xyz.iconn[i*(MXCON+1)+j] = TMPxyz.iconn[i*(MXCON+1)+j];
	   for (j=0; j<MAX13; j++)
		xyz.i13[i*MAX13+j] = TMPxyz.i13[i*MAX13+j];
	   for (j=0; j<MAX14; j++)
		xyz.i14[i*MAX14+j] = TMPxyz.i14[i*MAX14+j];
	}
	for (i=0; i < *xyz.nbnd; i++) {
	   for (j=0; j<2; j++)
		xyz.ibnd[i*2+j] = TMPxyz.ibnd[i*2+j];
	   xyz.bl[i] = TMPxyz.bl[i];
	   xyz.bk[i] = TMPxyz.bk[i];
	}
	for (i=0; i < *xyz.nang; i++) {
	   for (j=0; j<3; j++)
		xyz.iang[i*3+j] = TMPxyz.iang[i*3+j];
	   xyz.ango[i] = TMPxyz.ango[i];
	   xyz.ak[i] = TMPxyz.ak[i];
	}
	for (i=0; i < *xyz.nt; i++) {
	   for (j=0; j<4; j++)
		xyz.it[i*4+j] = TMPxyz.it[i*4+j];
	   for (j=0; j<4; j++) {
		xyz.trs1[i*4+j] = TMPxyz.trs1[i*4+j];
		xyz.trs2[i*4+j] = TMPxyz.trs2[i*4+j];
		xyz.trs3[i*4+j] = TMPxyz.trs3[i*4+j];
		xyz.trs4[i*4+j] = TMPxyz.trs4[i*4+j];
	   }
	}
	for (i=0; i < *xyz.nti; i++) {
	   for (j=0; j<4; j++)
		xyz.iti[i*4+j] = TMPxyz.iti[i*4+j];
	   for (j=0; j<4; j++) {
		xyz.trsi1[i*4+j] = TMPxyz.trsi1[i*4+j];
		xyz.trsi2[i*4+j] = TMPxyz.trsi2[i*4+j];
	   }
	}

	if (*TMPxyz.mxnat != TMPxyz.mxorg) {
	   free(TMPxyz.coo);
	   free(TMPxyz.ctmp);
	   free(TMPxyz.q);
	   free(TMPxyz.fr);
	   free(TMPxyz.coot);
	   free(TMPxyz.frt);
	   free(TMPxyz.zr);
	   free(TMPxyz.y);
	   free(TMPxyz.yt);
	   free(TMPxyz.pt);
	   free(TMPxyz.s);
	   free(TMPxyz.ftmp);
	   free(TMPxyz.iaton);
	   free(TMPxyz.iopt);
	   free(TMPxyz.watprot);
	   free(TMPxyz.n13);
	   free(TMPxyz.i13);
	   free(TMPxyz.n14);
	   free(TMPxyz.i14);
	   free(TMPxyz.iresid);
	   free(TMPxyz.iconn);
	   free(TMPxyz.ityp);
	   free(TMPxyz.ibnd);
	   free(TMPxyz.iang);
	   free(TMPxyz.it);
	   free(TMPxyz.iti);
	   free(TMPxyz.bl);
	   free(TMPxyz.bk);
	   free(TMPxyz.ango);
	   free(TMPxyz.ak);
	   free(TMPxyz.trs1);
	   free(TMPxyz.trs2);
	   free(TMPxyz.trs3);
	   free(TMPxyz.trs4);
	   free(TMPxyz.trsi1);
	   free(TMPxyz.trsi2);
	}
	*xyz.mxnat   = ZSize;
	*xyz.maxbnd  = ZSize*2;
	*xyz.maxang  = ZSize*3;
	*xyz.maxtors = ZSize*4;
   }
}

#if defined(VMS) || defined(UNDERSC)
allq()
#else
#ifdef CRAY
ALLQ()
#else
allq_()
#endif
#endif
{
   int memstat;
   float d;
   int i;
   short int j;
   int ZSize, wsize;

   ZSize = *xyz.mxnat;

   if ((qpot = (float *) malloc((sizeof d)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

}

#if defined(VMS) || defined(UNDERSC)
angarr(istat)
#else
#ifdef CRAY
ANGARR(istat)
#else
angarr_(istat)
#endif
#endif
int *istat;
{
#if defined(VMS) || defined(UNDERSC)
	angard(istat,
#else
#ifdef CRAY
	ANGARD(istat,
#else
	angard_(istat,
#endif
#endif
	xyz.nang,xyz.iang,xyz.ango,xyz.ak,xyz.maxang,xyz.iconn,xyz.ityp,xyz.iopt);
}

#if defined(VMS) || defined(UNDERSC)
bndarr(istat)
#else
#ifdef CRAY
BNDARR(istat)
#else
bndarr_(istat)
#endif
#endif
int *istat;
{
#if defined(VMS) || defined(UNDERSC)
	bndard(istat,
#else
#ifdef CRAY
	BNDARD(istat,
#else
	bndard_(istat,
#endif
#endif
	xyz.nbnd,xyz.ibnd,xyz.bl,xyz.bk,xyz.maxbnd,xyz.iconn,xyz.ityp,xyz.iopt);
}

#if defined(VMS) || defined(UNDERSC)
itrarr(istat)
#else
#ifdef CRAY
ITRARR(istat)
#else
itrarr_(istat)
#endif
#endif
int *istat;
{
#if defined(VMS) || defined(UNDERSC)
	itrard(istat,
#else
#ifdef CRAY
	ITRARD(istat,
#else
	itrard_(istat,
#endif
#endif
		xyz.maxtors,xyz.nti,xyz.iti,
		xyz.trsi1,xyz.trsi2,xyz.iconn,xyz.ityp,xyz.iopt);
}

#if defined(VMS) || defined(UNDERSC)
torarr(istat)
#else
#ifdef CRAY
TORARR(istat)
#else
torarr_(istat)
#endif
#endif
int *istat;
{
#if defined(VMS) || defined(UNDERSC)
	torard(istat,
#else
#ifdef CRAY
	TORARD(istat,
#else
	torard_(istat,
#endif
#endif
		xyz.maxtors,xyz.nbnd,xyz.ibnd,
		xyz.nt,xyz.it,xyz.trs1,xyz.trs2,xyz.trs3,xyz.trs4,
		xyz.iconn,xyz.ityp,xyz.iopt);
}

#if defined(VMS) || defined(UNDERSC)
prttor(istat)
#else
#ifdef CRAY
PRTTOR(istat)
#else
prttor_(istat)
#endif
#endif
int *istat;
{
#if defined(VMS) || defined(UNDERSC)
	prttod(istat,
#else
#ifdef CRAY
	PRTTOD(istat,
#else
	prttod_(istat,
#endif
#endif
		xyz.maxtors,
		xyz.nt,xyz.it,xyz.trs1,xyz.trs2,xyz.trs3,xyz.trs4,
		xyz.iconn,xyz.ityp,xyz.iopt);
}

#if defined(VMS) || defined(UNDERSC)
pritor(istat)
#else
#ifdef CRAY
PRITOR(istat)
#else
pritor_(istat)
#endif
#endif
int *istat;
{
#if defined(VMS) || defined(UNDERSC)
	pritod(istat,
#else
#ifdef CRAY
	PRITOD(istat,
#else
	pritod_(istat,
#endif
#endif
		xyz.maxtors,
		xyz.nti,xyz.iti,xyz.trsi1,xyz.trsi2,
		xyz.iconn,xyz.ityp,xyz.iopt);
}

#if defined(VMS) || defined(UNDERSC)
asschg(idebug)
#else
#ifdef CRAY
ASSCHG(idebug)
#else
asschg_(idebug)
#endif
#endif
int *idebug;
{
#if defined(VMS) || defined(UNDERSC)
	asschd(idebug,
#else
#ifdef CRAY
	ASSCHD(idebug,
#else
	asschd_(idebug,
#endif
#endif
		xyz.ityp,xyz.iconn,xyz.q);
}

#if defined(VMS) || defined(UNDERSC)
conn34(istat)
#else
#ifdef CRAY
CONN34(istat)
#else
conn34_(istat)
#endif
#endif
int *istat;
{
#if defined(VMS) || defined(UNDERSC)
	cond34(istat,
#else
#ifdef CRAY
	COND34(istat,
#else
	cond34_(istat,
#endif
#endif
		xyz.n13,xyz.i13,xyz.n14,xyz.i14,xyz.iconn);
}

#if defined(VMS) || defined(UNDERSC)
getinp(igtinp)
#else
#ifdef CRAY
GETINP(igtinp)
#else
getinp_(igtinp)
#endif
#endif
int *igtinp;
{
	int i;

#if defined(VMS) || defined(UNDERSC)
	getind(igtinp,
#else
#ifdef CRAY
	GETIND(igtinp,
#else
	getind_(igtinp,
#endif
#endif
	 xyz.iconn,xyz.iopt,xyz.iresid,xyz.ityp,xyz.coo,xyz.q);

	if (*igtinp == -1) {
	   addat = *xyz.iatoms+2;

/* +2 (2*3=6) for cell parameters */

#if defined(VMS) || defined(UNDERSC)
	   allcoo(&addat);
#else
#ifdef CRAY
	   ALLCOO(&addat);
#else
	   allcoo_(&addat);
#endif
#endif
#if defined(VMS) || defined(UNDERSC)
	   getind(igtinp,
#else
#ifdef CRAY
	   GETIND(igtinp,
#else
	   getind_(igtinp,
#endif
#endif
	   xyz.iconn,xyz.iopt,xyz.iresid,xyz.ityp,xyz.coo,xyz.q);
	}
	for (i=0; i < 3*(*xyz.iatoms); i++) {
	    xyz.ctmp[i] = xyz.coo[i];
	}
}

#if defined(VMS) || defined(UNDERSC)
wrmon(ncyc,emin)
#else
#ifdef CRAY
WRMON(ncyc,emin)
#else
wrmon_(ncyc,emin)
#endif
#endif
int *ncyc;
float *emin;
{
	FILE *out;
	char lckfile[80];
	int status;

	if (*nproc > 0 && *taskid != 0) return;

        sprintf(lckfile,"%d.ambforw",*ncyc);

	out = fopen(lckfile,"w");
	fprintf(out,"lock\n");
        fclose(out);

#if defined(VMS) || defined(UNDERSC)
	wrmod(ncyc,emin);
#else
#ifdef CRAY
	WRMOD(ncyc,emin);
#else
	wrmod_(ncyc,emin);
#endif
#endif

	status = remove(lckfile);
}

#if defined(VMS) || defined(UNDERSC)
wrtout(iun,emin)
#else
#ifdef CRAY
WRTOUT(iun,emin)
#else
wrtout_(iun,emin)
#endif
#endif
int *iun;
float *emin;
{

	if (*nproc > 0 && *taskid != 0) return;

#if defined(VMS) || defined(UNDERSC)
	wrtoud(iun,emin,
#else
#ifdef CRAY
	WRTOUD(iun,emin,
#else
	wrtoud_(iun,emin,
#endif
#endif
		xyz.iconn,xyz.ityp,xyz.coo,xyz.q,xyz.iopt);

}

#if defined(VMS) || defined(UNDERSC)
wrtbin(iun,emin)
#else
#ifdef CRAY
WRTBIN(iun,emin)
#else
wrtbin_(iun,emin)
#endif
#endif
int *iun;
float *emin;
{
#if defined(VMS) || defined(UNDERSC)
	wrtbid(iun,emin,
#else
#ifdef CRAY
	WRTBID(iun,emin,
#else
	wrtbid_(iun,emin,
#endif
#endif
		xyz.coo);
}

#if defined(VMS) || defined(UNDERSC)
restr(coo)
#else
#ifdef CRAY
RESTR(coo)
#else
restr_(coo)
#endif
#endif
float *coo;
{
   int i,j,k;

   for (i=0; i < *xyz.iatoms; i++ ) {
     if (!xyz.iopt[i]) {
	for (j=0; j < 3; j++) {
	    coo[i*3+j] = xyz.ctmp[i*3+j];
	}
     }
   }
}

static float *dftmp = NULL;
static float *fwat = NULL;
static float *fx = NULL;
static float *fy = NULL;
static float *fz = NULL;

#if defined(VMS) || defined(UNDERSC)
enegrd(energy,par,forces)
#else
#ifdef CRAY
ENEGRD(energy,par,forces)
#else
enegrd_(energy,par,forces)
#endif
#endif
float *energy;
float *par;
float *forces;
{
    int memstat,wsize;
    float f;
    float d;

    memstat = 1;
    wsize = (*xyz.iatoms);

    if (fx == NULL) {
	if ((fx = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((fy = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((fz = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
    }
    if (dftmp == NULL) {
	if ((dftmp = (float *) malloc((sizeof d)*3*wsize)) == NULL) 
		memstat = 0;
    }

    if (fwat == NULL) {
	if ((fwat = (float *) malloc((sizeof d)*3*wsize)) == NULL) 
		memstat = 0;
    }

    if (!memstat) fprintf(stderr,"enegrd: Out of memory\n");

#if defined(VMS) || defined(UNDERSC)
	enegdd(energy,par,forces,
#else
#ifdef CRAY
	ENEGDD(energy,par,forces,
#else
	enegdd_(energy,par,forces,
#endif
#endif
		dftmp,fwat,fx,fy,fz,
		xyz.ftmp,xyz.iopt,xyz.iresid,
		xyz.n13,xyz.i13,xyz.n14,xyz.i14,
		xyz.nbnd,xyz.ibnd,xyz.bl,xyz.bk,
		xyz.nang,xyz.iang,xyz.ango,xyz.ak,
		xyz.nt,xyz.it,xyz.trs1,xyz.trs2,xyz.trs3,xyz.trs4,
		xyz.nti,xyz.iti,xyz.trsi1,xyz.trsi2,
		xyz.q,xyz.iconn,xyz.ityp,xyz.watprot,nlst,lst);
}

#ifndef MD
#if defined(VMS) || defined(UNDERSC)
optimise(energy,gtol,nsd)
#else
#ifdef CRAY
OPTIMISE(energy,gtol,nsd)
#else
optimise_(energy,gtol,nsd)
#endif
#endif
float *energy;
float *gtol;
int *nsd;
{
#if defined(VMS) || defined(UNDERSC)
	optimisd(energy,gtol,nsd,
#else
#ifdef CRAY
	OPTIMISD(energy,gtol,nsd,
#else
	optimisd_(energy,gtol,nsd,
#endif
#endif
		xyz.coo,xyz.coot,xyz.fr,xyz.frt,xyz.zr,xyz.y,xyz.yt,
		xyz.pt,xyz.s,xyz.iopt,xyz.ityp);
}

#if defined(VMS) || defined(UNDERSC)
lbfgs(energy,gtol,nsd)
#else
#ifdef CRAY
LBFGS(energy,gtol,nsd)
#else
lbfgs_(energy,gtol,nsd)
#endif
#endif
float *energy;
float *gtol;
int *nsd;
{
   float d;
   int wsize, memstat;

   memstat = 1;

   wsize = (*xyz.iatoms)*3*(2*MSAVE+1)+2*MSAVE;
   if (xyz.work == NULL) {
	if ((xyz.work = (float *) malloc((sizeof d)*wsize)) == NULL) {
	   memstat = 0;
	   fprintf(stderr,"sdrive: Out of memory\n");
	}
   }

#if defined(VMS) || defined(UNDERSC)
	lbfgd(energy,gtol,nsd,
#else
#ifdef CRAY
	LBFGD(energy,gtol,nsd,
#else
	lbfgd_(energy,gtol,nsd,
#endif
#endif
		xyz.coo,xyz.fr,xyz.frt,xyz.work);
}
#endif

static float *cx = NULL;
static float *cy = NULL;
static float *cz = NULL;
static float *vdwe = NULL;
static float *vdwr = NULL;
static float *potq = NULL;
static float *potv = NULL;
static float *eps = NULL;
static float *pre6 = NULL;
static float *ftmp = NULL;
static int *iag = NULL;
static int *ikl = NULL;
static int *idoq = NULL;
static int *idov = NULL;

static float *dvdwe = NULL;
static float *dvdwr = NULL;
static float *dpotq = NULL;
static float *dpotv = NULL;
static float *deps = NULL;
static float *dpre6 = NULL;
static float *dfftmp = NULL;
static float *frac = NULL;
static float *exyz = NULL;

#if defined(VMS) || defined(UNDERSC)
qvdwr(
#else
#ifdef CRAY
QVDWR(
#else
qvdwr_(
#endif
#endif
	ec,ev,nac,iac,nad,iad,iresid,coo,iconn,q,iopt,ityp,fx,fy,fz)
float *ec;
float *ev;
int *nac;
int *iac;
int *nad;
int *iad;
int *iresid;
float *coo;
int *iconn;
float *q;
int *iopt;
short int *ityp;
float *fx;
float *fy;
float *fz;
{
   float f;
   int i,wsize, memstat;

   memstat = 1;

   wsize = (*xyz.iatoms);

   if (cx == NULL) {
	if ((cx = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((cy = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((cz = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((potq = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((potv = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((vdwe = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((vdwr = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((eps = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((pre6 = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((iag = (int *) malloc((sizeof i)*wsize)) == NULL) memstat = 0;
	if ((ikl = (int *) malloc((sizeof i)*wsize)) == NULL) memstat = 0;

	if (!memstat) fprintf(stderr,"qvdwr: Out of memory\n");
   }

#if defined(VMS) || defined(UNDERSC)
   qvdwd(
#else
#ifdef CRAY
   QVDWD(
#else
   qvdwd_(
#endif
#endif
	ec,ev,nac,iac,nad,iad,iresid,coo,iconn,q,iopt,ityp,fx,fy,fz,
	cx,cy,cz,vdwe,vdwr,potq,potv,eps,pre6,iag,ikl,ftmp);
}


#if defined(VMS) || defined(UNDERSC)
qvdwrd(
#else
#ifdef CRAY
QVDWRD(
#else
qvdwrd_(
#endif
#endif
	ec,ev,fintr,nac,iac,nad,iad,iresid,coo,iconn,q,iopt,ityp)
float *ec;
float *ev;
float *fintr;
int *nac;
int *iac;
int *nad;
int *iad;
int *iresid;
float *coo;
int *iconn;
float *q;
int *iopt;
short int *ityp;
{
   float f;
   int i,wsize, memstat;

   memstat = 1;

   wsize = (*xyz.iatoms);
   if (cx != NULL) {
	free(cx);
	free(cy);
	free(cz);
	free(fx);
	free(fy);
	free(fz);
	free(potq);
	free(potv);
	free(vdwe);
	free(vdwr);
	free(eps);
	free(pre6);
	free(ftmp);
	cx = NULL;
	cy = NULL;
	cz = NULL;
/*
	fx = NULL;
	fy = NULL;
	fz = NULL;
*/
	vdwe = NULL;
	vdwr = NULL;
	potq = NULL;
	potv = NULL;
	eps = NULL;
	pre6 = NULL;
	ftmp = NULL;
   }

   if (dpotq == NULL) {
	if ((dpotq = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dpotv = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dvdwe = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dvdwr = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((deps = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dpre6 = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((iag = (int *) malloc((sizeof i)*wsize)) == NULL) memstat = 0;
	if ((ikl = (int *) malloc((sizeof i)*wsize)) == NULL) memstat = 0;
	if ((dfftmp = (float *) malloc((sizeof f)*3*wsize)) == NULL) memstat = 0;

	if (!memstat) fprintf(stderr,"qvdwd: Out of memory\n");
   }

#if defined(VMS) || defined(UNDERSC)
   qvdwdd(
#else
#ifdef CRAY
   QVDWDD(
#else
   qvdwdd_(
#endif
#endif
	ec,ev,fintr,nac,iac,nad,iad,iresid,coo,iconn,q,iopt,ityp,
	dvdwe,dvdwr,dpotq,dpotv,deps,dpre6,iag,ikl,dfftmp);
}


#if defined(VMS) || defined(UNDERSC)
qvc(
#else
#ifdef CRAY
QVC(
#else
qvc_(
#endif
#endif
	ec,ev,nac,iac,nad,iad,iresid,coo,iconn,q,forces,cellder,iopt,ityp)
float *ec;
float *ev;
int *nac;
int *iac;
int *nad;
int *iad;
int *iresid;
float *coo;
int *iconn;
float *q;
float *forces;
float *cellder;
int *iopt;
short int *ityp;
{
   float f;
   int i,wsize, memstat;

   memstat = 1;

   wsize = (*xyz.iatoms);
   if (cx != NULL) {
	free(cx);
	free(cy);
	free(cz);
	free(fx);
	free(fy);
	free(fz);
	free(potq);
	free(potv);
	free(vdwe);
	free(vdwr);
	free(eps);
	free(pre6);
	free(ftmp);
	cx = NULL;
	cy = NULL;
	cz = NULL;
	fx = NULL;
	fy = NULL;
	fz = NULL;
	vdwe = NULL;
	vdwr = NULL;
	potq = NULL;
	potv = NULL;
	eps = NULL;
	pre6 = NULL;
	ftmp = NULL;
   }

   if (dpotq == NULL) {
	if ((dpotq = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dpotv = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dvdwe = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dvdwr = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((deps = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dpre6 = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((iag = (int *) malloc((sizeof i)*wsize)) == NULL) memstat = 0;
	if ((ikl = (int *) malloc((sizeof i)*wsize)) == NULL) memstat = 0;
	if ((dfftmp = (float *) malloc((sizeof f)*3*wsize)) == NULL) memstat = 0;
	if ((frac = (float *) malloc((sizeof f)*3*wsize)) == NULL) memstat = 0;

	if ((exyz = (float *) malloc((sizeof f)*3*wsize*125)) == NULL) memstat = 0;

	if (!memstat) fprintf(stderr,"qvdwd: Out of memory\n");
   }

#if defined(VMS) || defined(UNDERSC)
   qvd(
#else
#ifdef CRAY
   QVD(
#else
   qvd_(
#endif
#endif
	ec,ev,nac,iac,nad,iad,iresid,coo,iconn,q,forces,cellder,iopt,ityp,
	dvdwe,dvdwr,dpotq,dpotv,deps,dpre6,iag,ikl,dftmp,dfftmp,frac,exyz);
}


#if defined(VMS) || defined(UNDERSC)
qenvdw(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
      q,forces,iopt,nlst,lst,ityp)
#else
#ifdef CRAY
QENVDW(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
      q,forces,iopt,nlst,lst,ityp)
#else
qenvdw_(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
      q,forces,iopt,nlst,lst,ityp)
#endif
#endif
float *ec;
float *ev;
int *nac;
int *iac;
int *nad;
int *iad;
int *iresid;
float *coo;
int *iconn;
float *q;
float *forces;
int *iopt;
int *nlst;
int *lst;
short int *ityp;
{
   float f;
   int i,wsize, memstat;

   memstat = 1;

   wsize = (*xyz.iatoms);
   if (cx != NULL) {
	free(cx);
	free(cy);
	free(cz);
	free(fx);
	free(fy);
	free(fz);
	free(potq);
	free(potv);
	free(vdwe);
	free(vdwr);
	free(eps);
	free(pre6);
	free(ftmp);
	cx = NULL;
	cy = NULL;
	cz = NULL;
	fx = NULL;
	fy = NULL;
	fz = NULL;
	vdwe = NULL;
	vdwr = NULL;
	potq = NULL;
	potv = NULL;
	eps = NULL;
	pre6 = NULL;
	ftmp = NULL;
   }

   if (dpotq == NULL) {
	if ((dpotq = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dpotv = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dvdwe = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dvdwr = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((deps = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dpre6 = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((iag = (int *) malloc((sizeof i)*wsize)) == NULL) memstat = 0;
	if ((ikl = (int *) malloc((sizeof i)*wsize)) == NULL) memstat = 0;
	if ((idoq = (int *) malloc((sizeof i)*wsize)) == NULL) memstat = 0;
	if ((idov = (int *) malloc((sizeof i)*wsize)) == NULL) memstat = 0;
	if ((dftmp = (float *) malloc((sizeof f)*3*wsize)) == NULL) memstat = 0;
	if ((dfftmp = (float *) malloc((sizeof f)*3*wsize)) == NULL) memstat = 0;

	if (!memstat) fprintf(stderr,"qenvdw: Out of memory\n");
   }

#if defined(VMS) || defined(UNDERSC)
   qenvdd(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
         q,forces,iopt,nlst,lst,ityp,
#else
#ifdef CRAY
   QENVDD(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
         q,forces,iopt,nlst,lst,ityp,
#else
   qenvdd_(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
         q,forces,iopt,nlst,lst,ityp,
#endif
#endif
	 idoq,idov,dvdwe,dvdwr,dpotq,dpotv,deps,dpre6,iag,ikl,dftmp,dfftmp);
}

#if defined(VMS) || defined(UNDERSC)
qvnoe(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
      q,forces,iopt,ityp,iwtpr)
#else
#ifdef CRAY
QVNOE(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
      q,forces,iopt,ityp,iwtpr)
#else
qvnoe_(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
      q,forces,iopt,ityp,iwtpr)
#endif
#endif
float *ec;
float *ev;
int *nac;
int *iac;
int *nad;
int *iad;
int *iresid;
float *coo;
int *iconn;
float *q;
float *forces;
int *iopt;
int *ityp;
int *iwtpr;
{
   float f;
   int i,wsize, memstat;

   memstat = 1;

   wsize = (*xyz.iatoms);
   if (cx != NULL) {
	free(cx);
	free(cy);
	free(cz);
	free(fx);
	free(fy);
	free(fz);
	free(potq);
	free(potv);
	free(vdwe);
	free(vdwr);
	free(eps);
	free(pre6);
	free(ftmp);
	cx = NULL;
	cy = NULL;
	cz = NULL;
	fx = NULL;
	fy = NULL;
	fz = NULL;
	vdwe = NULL;
	vdwr = NULL;
	potq = NULL;
	potv = NULL;
	eps = NULL;
	pre6 = NULL;
	ftmp = NULL;
   }

   if (dpotq == NULL) {
	if ((dpotq = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dpotv = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dvdwe = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dvdwr = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((deps = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((dpre6 = (float *) malloc((sizeof f)*wsize)) == NULL) memstat = 0;
	if ((iag = (int *) malloc((sizeof i)*wsize)) == NULL) memstat = 0;
	if ((ikl = (int *) malloc((sizeof i)*wsize)) == NULL) memstat = 0;
	if ((dftmp = (float *) malloc((sizeof f)*3*wsize)) == NULL) memstat = 0;
	if ((dfftmp = (float *) malloc((sizeof f)*3*wsize)) == NULL) memstat = 0;

	if (!memstat) fprintf(stderr,"qvnoe: Out of memory\n");
   }

#if defined(VMS) || defined(UNDERSC)
   qvnod(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
         q,forces,iopt,ityp,iwtpr,
#else
#ifdef CRAY
   QVNOD(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
         q,forces,iopt,ityp,iwtpr,
#else
   qvnod_(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
         q,forces,iopt,ityp,iwtpr,
#endif
#endif
	 dvdwe,dvdwr,dpotq,dpotv,deps,dpre6,iag,ikl,dftmp,dfftmp);
}

#if defined(VMS) || defined(UNDERSC)
esolb()
#else
#ifdef CRAY
ESOLB()
#else
esolb_()
#endif
#endif
{
   int memstat,ZSize;
   float d;

   ZSize = *xyz.mxnat;
   memstat = 1;

   if (drb == NULL) {
   if ((drb = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }
   if ((rsolv = (float *) malloc((sizeof d)*ZSize)) == NULL) {
	memstat = 0;
   }
   if ((rborn = (float *) malloc((sizeof d)*ZSize)) == NULL) {
	memstat = 0;
   }
   if ((shct = (float *) malloc((sizeof d)*ZSize)) == NULL) {
	memstat = 0;
   }
   }

#if defined(VMS) || defined(UNDERSC)
   esolv(
#else
#ifdef CRAY
   ESOLV(
#else
   esolv_(
#endif
#endif
	xyz.coo,xyz.fr,drb,rsolv,rborn,shct,xyz.q,xyz.iconn,xyz.ityp);
}

#if defined(VMS) || defined(UNDERSC)
int mseed()
#else
#ifdef CRAY
int MSEED()
#else
int mseed_()
#endif
#endif
{
  return (int)(((long)time(NULL)+(long)getpid()) % (long)1000000);
}


#if defined(VMS) || defined(UNDERSC)
makbox()
#else
#ifdef CRAY
MAKBOX()
#else
makbox_()
#endif
#endif
{
#if defined(VMS) || defined(UNDERSC)
	makbod(
#else
#ifdef CRAY
	MAKBOD(
#else
	makbod_(
#endif
#endif
		xyz.coo);
}

#if defined(VMS) || defined(UNDERSC)
filbox()
#else
#ifdef CRAY
FILBOX()
#else
filbox_()
#endif
#endif
{
#if defined(VMS) || defined(UNDERSC)
	filbod(
#else
#ifdef CRAY
	FILBOD(
#else
	filbod_(
#endif
#endif
		water,xyz.coo,xyz.iconn,xyz.iresid,xyz.ityp,
		xyz.n13,xyz.i13,xyz.n14,
		xyz.nbnd,xyz.ibnd,xyz.bl,xyz.bk,
		xyz.nang,xyz.iang,xyz.ango,xyz.ak,
		xyz.q,xyz.iopt,xyz.watprot);
}

#if defined(VMS) || defined(UNDERSC)
cntwat()
#else
#ifdef CRAY
CNTWAT()
#else
cntwat_()
#endif
#endif
{
#if defined(VMS) || defined(UNDERSC)
	cntwad(
#else
#ifdef CRAY
	CNTWAD(
#else
	cntwad_(
#endif
#endif
		water,xyz.coo,xyz.iconn,xyz.iresid,xyz.ityp,
		xyz.q,xyz.iopt,qpot);
}

#if defined(VMS) || defined(UNDERSC)
watcnt()
#else
#ifdef CRAY
WATCNT()
#else
watcnt_()
#endif
#endif
{

#if defined(VMS) || defined(UNDERSC)
   watcnd(
#else
#ifdef CRAY
   WATCND(
#else
   watcnd_(
#endif
#endif
	xyz.coo);
}

#if defined(VMS) || defined(UNDERSC)
watlst(niwat)
#else
#ifdef CRAY
WATLST(niwat)
#else
watlst_(niwat)
#endif
#endif
int *niwat;
{
   int i,nelem;

   nelem = 2*(*niwat);

   if (listw == NULL) {
	if ((listw = (int *) malloc((sizeof i)*2*nelem)) == NULL) {
	   fprintf(stderr,"Out of memory\n");
	}
   }

#if defined(VMS) || defined(UNDERSC)
   watlsd(
#else
#ifdef CRAY
   WATLSD(
#else
   watlsd_(
#endif
#endif
	xyz.coo,xyz.watprot,listw);
}

#if defined(VMS) || defined(UNDERSC)
stutab()
#else
#ifdef CRAY
STUTAB()
#else
stutab_()
#endif
#endif
{

   float d;
   int ZSize;

   ZSize = 10000;

   if ((woo = (float *) malloc((sizeof d)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

   if ((dewoo = (float *) malloc((sizeof d)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

   if ((whh = (float *) malloc((sizeof d)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

   if ((dewhh = (float *) malloc((sizeof d)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

   if ((woh = (float *) malloc((sizeof d)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

   if ((dewoh = (float *) malloc((sizeof d)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

#if defined(VMS) || defined(UNDERSC)
   stutad(
#else
#ifdef CRAY
   STUTAD(
#else
   stutad_(
#endif
#endif
	woo,whh,woh,dewoo,dewhh,dewoh);
}

#if defined(VMS) || defined(UNDERSC)
etutab(nsize)
#else
#ifdef CRAY
ETUTAB(nsize)
#else
etutab_(nsize)
#endif
#endif
int *nsize;
{

   float r;
   int ZSize;

   ZSize = *nsize;

   if ((ewoo = (float *) malloc((sizeof r)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

   if ((edwoo = (float *) malloc((sizeof r)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

   if ((ewhh = (float *) malloc((sizeof r)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

   if ((edwhh = (float *) malloc((sizeof r)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

   if ((ewoh = (float *) malloc((sizeof r)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

   if ((edwoh = (float *) malloc((sizeof r)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

#if defined(VMS) || defined(UNDERSC)
   etutad(nsize,
#else
#ifdef CRAY
   ETUTAD(nsize,
#else
   etutad_(nsize,
#endif
#endif
	ewoo,ewhh,ewoh,edwoo,edwhh,edwoh);
}

#if defined(VMS) || defined(UNDERSC)
fstwat(ew,coo,forces,q)
#else
#ifdef CRAY
FSTWAT(ew,coo,forces,q)
#else
fstwat_(ew,coo,forces,q)
#endif
#endif
float *ew;
float *coo;
float *forces;
float *q;
{
   float f;
   int wsize, memstat;

   memstat = 1;

   wsize = (*xyz.iatoms);

   if (dfftmp == NULL) {
	if ((dfftmp = (float *) malloc((sizeof f)*3*wsize)) == NULL) memstat = 0;

	if (!memstat) fprintf(stderr,"fstwat: Out of memory\n");
   }

#if defined(VMS) || defined(UNDERSC)
	fstwad(ew,coo,forces,q,
#else
#ifdef CRAY
	FSTWAD(ew,coo,forces,q,
#else
	fstwad_(ew,coo,forces,q,
#endif
#endif
		dfftmp,listw,woo,whh,woh,dewoo,dewhh,dewoh);
}

#if defined(VMS) || defined(UNDERSC)
bldlst()
#else
#ifdef CRAY
BLDLST()
#else
bldlst_()
#endif
#endif
{
   int memstat,i,size;

   size = *xyz.iatoms;
   memstat = 1;

   if (nlst == NULL) {
	if ((nlst = (int *) malloc((sizeof i)*size)) == NULL) {
	   memstat = 0;
	   fprintf(stderr,"bldlst: Out of memory\n");
	}
   }

   size = MXNEIB*size;
   if (lst == NULL) {
	if ((lst = (int *) malloc((sizeof i)*size)) == NULL) {
	   memstat = 0;
	   fprintf(stderr,"bldlst: Out of memory\n");
	}
   }

#if defined(VMS) || defined(UNDERSC)
   bldlsd(
#else
#ifdef CRAY
   BLDLSD(
#else
   bldlsd_(
#endif
#endif
	xyz.coo,xyz.iresid,nlst,lst,&memstat);
}

#if defined(VMS) || defined(UNDERSC)
allerf()
#else
#ifdef CRAY
ALLERF()
#else
allerf_()
#endif
#endif
{

   float f;
   int ZSize;
/*
   RC = 15.0
   ZSize = 3000000;
*/

   ZSize = 2700000;
   if ((apprerf = (float *) malloc((sizeof f)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

/*
   RC = 15.0
   ZSize = 9000000;
*/
   ZSize = 7290000;

   if ((appexp = (float *) malloc((sizeof f)*ZSize)) == NULL) {
	fprintf(stderr,"Out of memory\n");
   }

#if defined(VMS) || defined(UNDERSC)
   allerd(
#else
#ifdef CRAY
   ALLERD(
#else
   allerd_(
#endif
#endif
	apprerf,appexp);
}

#if defined(VMS) || defined(UNDERSC)
apperfc(r,fc,fc2)
#else
#ifdef CRAY
APPERFC(r,fc,fc2)
#else
apperfc_(r,fc,fc2)
#endif
#endif
float *r;
float *fc;
float *fc2;
{
      float r1000,x;
      int ir;
      float rr,rr2;

      rr = *r;
      rr2 = rr*rr;

      if (*r <= 2.7) {
          r1000 = (float) (*r) * 1000000.0;
          ir = (int) r1000;
          x = (r1000 - (float) ir);
          *fc = ((1.0-(float) x)*apprerf[ir] + (float) x*apprerf[ir+1]);
      } else {
          *fc = apprerf[2699999];
      }

      if (rr2 <= 7.29) {
          r1000 = (float) (rr2) * 1000000.0;
          ir = (int) r1000;
          x = (r1000 - (float) ir);
          *fc2 = ((1.0-(float) x)*appexp[ir] + (float) x*appexp[ir+1]);
      } else {
          *fc2 = appexp[7289999];
      }

}

/*
#if defined(VMS) || defined(UNDERSC)
apperfc(r,fc,fc2)
#else
#ifdef CRAY
APPERFC(r,fc,fc2)
#else
apperfc_(r,fc,fc2)
#endif
#endif
float *r;
float *fc;
float *fc2;
{
      float r1000,x;
      int ir;
      float rr,rr2;

      rr = *r;
      rr2 = rr*rr;

      if (*r <= 3.0) {
          r1000 = (float) (*r) * 1000000.0;
          ir = (int) r1000;
          x = (r1000 - (float) ir);
          *fc = ((1.0-(float) x)*apprerf[ir] + (float) x*apprerf[ir+1]);
      } else {
          *fc = apprerf[2999999];
      }

      if (rr2 <= 9.0) {
          r1000 = (float) (rr2) * 1000000.0;
          ir = (int) r1000;
          x = (r1000 - (float) ir);
          *fc2 = ((1.0-(float) x)*appexp[ir] + (float) x*appexp[ir+1]);
      } else {
          *fc2 = appexp[8999999];
      }

}
*/

#if defined(VMS) || defined(UNDERSC)
appbnx()
#else
#ifdef CRAY
APPBNX()
#else
appbnx_()
#endif
#endif
{
#if defined(VMS) || defined(UNDERSC)
	appbnd(
#else
#ifdef CRAY
	APPBND(
#else
	appbnd_(
#endif
#endif
		xyz.coo,xyz.ityp);
}

#ifdef MD

#if defined(VMS) || defined(UNDERSC)
allmd(ZSizep)
#else
#ifdef CRAY
ALLMD(ZSizep)
#else
allmd_(ZSizep)
#endif
#endif
int *ZSizep;
{
   int memstat,i;
   float d;
   int ZSize;

   ZSize = *ZSizep;
   memstat = 1;

   if ((md.v = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   if ((md.a = (float *) malloc((sizeof d)*ZSize*3)) == NULL) {
	memstat = 0;
   }

   if ((md.m = (float *) malloc((sizeof d)*ZSize)) == NULL) {
	memstat = 0;
   }

   for (i=0; i < ZSize*3; i++) {
	md.a[i] = 0.0;
   }
}

#if defined(VMS) || defined(UNDERSC)
assmas()
#else
#ifdef CRAY
ASSMAS()
#else
assmas_()
#endif
#endif
{
#if defined(VMS) || defined(UNDERSC)
	assmad(
#else
#ifdef CRAY
	ASSMAD(
#else
	assmad_(
#endif
#endif
		xyz.ityp,md.m);
}

#if defined(VMS) || defined(UNDERSC)
velini()
#else
#ifdef CRAY
VELINI()
#else
velini_()
#endif
#endif
{
#if defined(VMS) || defined(UNDERSC)
	velind(
#else
#ifdef CRAY
	VELIND(
#else
	velind_(
#endif
#endif
		xyz.coo,xyz.fr,md.v,md.a,md.m);
}

#if defined(VMS) || defined(UNDERSC)
verlet(istep)
#else
#ifdef CRAY
VERLET(istep)
#else
verlet_(istep)
#endif
#endif
int *istep;
{
#if defined(VMS) || defined(UNDERSC)
	verled(istep,
#else
#ifdef CRAY
	VERLED(istep,
#else
	verled_(istep,
#endif
#endif
		xyz.coo,xyz.fr,md.v,md.a,md.m,xyz.ityp);
}
#endif

prtmas_(iopt)
int *iopt;
{
 int i;
    
	for (i=0; i < *xyz.nang; i++) {
	fprintf(stderr,"%d ang %d %d %d %f %f\n",i,xyz.iang[i*3],xyz.iang[i*3+1],xyz.iang[i*3+2],xyz.ango[i],xyz.ak[i]);
	}
}
