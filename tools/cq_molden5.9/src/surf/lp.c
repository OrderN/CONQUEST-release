/* Copyright amitabh varshney@cs.unc.edu. All Rights Reserved. March 16, 1994 */
/*--------------------------------------------------------------------------------
lp.c

This file contains the function for linear programming in three dimensions.
This is used to check if the feasible cell is non-null. The algorithm that
is being used has a linear expected time complexity and has been observed
to be quite fast in practice. The algorithm for this has been given in:
"Linear Programming and Convex Hulls Made Easy" by Raimund Seidel, Proceedings 
of the Sixth Annual Symposium on Computational Geometry, June 6-8, 1990, 
Berkeley, CA, pp 211--215. To understand this code, I suggest you should read 
that paper. The terminology in the code is basically the same as in that paper, 
with a minor difference -- here array B is maintained virtually in A itself;
B goes from 0 .. i-1, A is from i .. num_A-1. The pseudo-code with that paper 
had a couple of minor errors that I have corrected in this implementation.
----------------------------------------------------------------------------------*/
#include "surf.h"

/*------------------------------------------------------------------------ 
returns 0 if infeasible, 1 if feasible
--------------------------------------------------------------------------*/
int linear_prog(dim, obj, extents, A, num_A, sol, sol_cons)
int		dim, num_A;
Vector		obj;
Vector		*sol;
Vector		*A;
float		extents[3][2];
int		sol_cons[3];
{ int		i, j, k, m, n;
  float		high, low;
  double	z;
  int		high_cons, low_cons;
  float		tempf;
  Vector	*A_new;
  Vector	factor;
  Vector	obj_new, sol_new;
  int		num_A_new;
  float	        extents_new[3][2];
  int		feasible;

  low_cons = high_cons = sol_cons[dim-1];

  if (dim == 1)
  { high = extents[0][MAXC];
    low = extents[0][MINC];
    z = HUGE_VAL;
    sol_cons[0] = -1;
    for(i = 0; i < num_A; i++)
    { if ((A[i].coeff[0] < LP_EPS) && (A[i].coeff[0] > -LP_EPS))
      { if (z > A[i].coeff[1])  z = A[i].coeff[1]; }
      else
      { tempf = A[i].coeff[1]/A[i].coeff[0];
        if (A[i].coeff[0] >= LP_EPS)
        { if (high > tempf) 
	  { high = tempf; 
	    high_cons =  (A[i].atom_id >= 0)? i: -1;
	  }
	}
        else if (A[i].coeff[0] <= -LP_EPS)
	{ if (low < tempf) 
	  { low = tempf; 
	    low_cons = (A[i].atom_id >= 0)? i: -1;
	  }
	}
      }
    }
    if ((z < 0) || (high < low)) 
      return(0);			/* problem infeasible */
    else
    { sol->coeff[0] = (obj.coeff[0] < -LP_EPS)? low: high;
      sol_cons[0] = (obj.coeff[0] < -LP_EPS)? low_cons: high_cons;
      return(1);			 /* return a feasible solution */
    }
  }
  else		/* higher dimensional problem */
  {
    A_new = &temp_cons[dim-2][0];
  
    /* compute optimum solution for the constraints given by extents */
    for(i = 0; i < dim; i++)
      if (obj.coeff[i] >= 0) sol->coeff[i] = extents[i][MAXC];
      else	             sol->coeff[i] = extents[i][MINC];
    
    /* Add all constraints of A one by one */
    for(i = 0; i < num_A; i++)
    { tempf = 0;
      for(j = 0; j < dim; j++)
	tempf += A[i].coeff[j]*sol->coeff[j];
      if (tempf > A[i].coeff[dim])	/* sol violates this constraint */
      { k = -1; 
	sol_cons[dim-1] = (A[i].atom_id >= 0)? i: -1;

	for (j = dim - 1; j >= 0; j--)
	{ if ((A[i].coeff[j] >= LP_EPS) || (A[i].coeff[j] <= -LP_EPS))
	  { k = j; 
	    break;
          }
        }
	if (k == -1)			/* no such k found */
	  return(0);			/* problem infeasible */

	/* if such a k is found */
	/* find factor */
	tempf = 1/A[i].coeff[k];
	for(j = 0; j <= dim; j++)
	  factor.coeff[j] = A[i].coeff[j]*tempf;

	/* Eliminate x_{k} from constraints in B and from C */
	num_A_new = i;
	for(j = 0; j < num_A_new; j++)
	{ for (m = n = 0; m <= dim; m++)
	  { if (m != k)
	      A_new[j].coeff[n++] = A[j].coeff[m]-A[j].coeff[k]*factor.coeff[m];
          }
        }
	for(m = n = 0; m < dim; m++)
	{ if (m != k)
	    obj_new.coeff[n++] = obj.coeff[m] - obj.coeff[k]*factor.coeff[m];
        }
        /* Incorporate the constraints l_{k} <= x_{k} <= u_{k} into A_new*/
	for(m = n = 0; m <= dim; m++)
	{ if (m != k)
	  { A_new[num_A_new].coeff[n] = -factor.coeff[m];
	    A_new[num_A_new+1].coeff[n] = factor.coeff[m];
	    n++;
          }
        }
	A_new[num_A_new].coeff[dim-1] += extents[k][MAXC];
	A_new[num_A_new].atom_id = -1;
	A_new[num_A_new+1].coeff[dim-1] -= extents[k][MINC];
	A_new[num_A_new+1].atom_id = -1;
        num_A_new += 2;
        /* Generate new upper and lower bounds */
	for(m = n = 0; m < dim; m++)
	{ if (m != k)
	  { extents_new[n][MINC] = extents[m][MINC];
	    extents_new[n][MAXC] = extents[m][MAXC];
	    n++;
          }
        }
        /* solve the (dim-1) dimensional problem and "lift" the solution */
        feasible = linear_prog(dim-1,obj_new, extents_new, A_new, num_A_new, 
			       &sol_new, sol_cons);
        if (!feasible) 
	  return (0);

	/* if feasible solution is returned */
        /* Insert sol_new into sol with 0 as the kth component */
	for(m = n = 0; m < dim; m++)
	  sol->coeff[m] = (m == k)? 0 : sol_new.coeff[n++];

	tempf = 0;
	for(m = 0; m < dim; m++)
	  tempf += A[i].coeff[m]*sol->coeff[m];
	sol->coeff[k] = (A[i].coeff[dim] - tempf)/A[i].coeff[k];
      }
    }
    return(1);
  }
}
