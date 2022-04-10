/* Nsp
 * Copyright (C) 2014-2022 Jean-Philippe Chancelier ENPC/Cermics
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 */

#include <ilcplex/cplex.h>

#include <nsp/nsp.h>
#include <nsp/objects.h>
#include <nsp/imatrix.h>
#include "nspcplex.h"


static int nsp_cplex_redirect_channels(CPXENVptr env,int loglevel);
static void CPXPUBLIC nsp_cplex_msg(void *handle, const char *msg);
static void CPXPUBLIC nsp_cplex_error(void *handle, const char *msg);

double  nsp_cplex_dbl_max()
{
  return CPX_INFBOUND;
}

/* creates a cplex problem and solve the problem.
 * If a filname is given the problem is also saved 
 * if save_only if true the problem is just saved.
 */ 

int nsp_cplex_solve(const char* problemName, int sense, int ncols, int nrows, int neq,
		    NspIMatrix*Cmatbeg, NspIMatrix *Cmatcount, NspIMatrix *Cmatind, NspMatrix *Cmatval, 
		    NspMatrix *lower, NspMatrix *upper, NspMatrix *Objective,
		    NspIMatrix *Qmatbeg,NspIMatrix *Qmatcnt, NspIMatrix *Qmatind, NspMatrix *Qmatval, 
		    NspMatrix *Rhs,const char *columnType,  NspMatrix *X,NspMatrix *Lambda,
		    NspMatrix *RetCost,NspMatrix *Retcode,const char *rowType,
		    int semiCount, int *semiIndex,NspHash *Options,int loglevel,
		    const char *filename,int save_only)
{
  /* matrix A part */
  int *matrixBegin = (int *) Cmatbeg->Iv;
  int *matrixIndex = (int *) Cmatind->Iv;
  int *matrixCount = (int *) Cmatcount->Iv;
  double *matrixValues = Cmatval->R; 
  int colCount = ncols, rowCount = nrows; 
  int objectSense = (sense == 0 ) ?  CPX_MIN: CPX_MAX;
  double *objectCoeffs = Objective->R;
  double *lowerBounds = lower->R;
  double *upperBounds = upper->R;
  double *rhsValues = Rhs->R;

  CPXENVptr     env = NULL;
  CPXLPptr      lp = NULL;
  int status;
  int cur_numrows,cur_numcols;
  double *dj=NULL,*slack=NULL;
  int solstat=0;
  double objval;
  int ret = FAIL;

  /* Initialize the CPLEX environment */

  env = CPXopenCPLEX (&status);
   
  /* If an error occurs, the status value indicates the reason for
     failure.  A call to CPXgeterrorstring will produce the text of
     the error message.  Note that CPXopenCPLEX produces no output,
     so the only way to see the cause of the error is to use
     CPXgeterrorstring.  For other CPLEX routines, the errors will
     be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON.  */
  
  if ( env == NULL ) {
    char  errmsg[CPXMESSAGEBUFSIZE];
    Scierror("Could not open CPLEX environment.\n");
    CPXgeterrorstring (env, status, errmsg);
    Scierror("%s", errmsg);
    goto error;
  }

  if (nsp_cplex_redirect_channels(env,loglevel) == FAIL)
    {
      goto error;
    }

  /* Turn on/off output to the screen */

  status = CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_OFF);
  if ( status ) {
    Scierror("Failure to set screen indicator, error %d.\n", status);
    goto error;
  }
  
  /* Fill in the data for the problem.  */
     /* Create the problem. */

  lp = CPXcreateprob (env, &status, "nsplp");

  /* A returned pointer of NULL may mean that not enough memory
     was available or there was some other problem.  In the case of
     failure, an error message will have been written to the error
     channel from inside CPLEX.  In this example, the setting of
     the parameter CPXPARAM_ScreenOutput causes the error message to
     appear on stdout.  */
  
  if ( lp == NULL ) {
    Scierror("Failed to create a cplex problem.\n");
    goto error;
  }

  /* Now copy the LP part of the problem data into the lp */
  status = CPXcopylp (env, lp, colCount, rowCount , objectSense, objectCoeffs, rhsValues,
		      rowType, matrixBegin, matrixCount, matrixIndex, matrixValues,
		      lowerBounds, upperBounds, NULL);

  if ( status ) {
    Scierror("Failed to copy problem data.\n");
    goto error;
  }
  
  if ( Qmatbeg != NULL ) 
    {
      /* we have a quadratic cost */
      status = CPXcopyquad (env, lp,(const int*) Qmatbeg->Iv,(const int*) Qmatcnt->Iv, (const int*) Qmatind->Iv,
			    Qmatval->R);
      if ( status ) {
	Scierror("Failed to copy quadratic matrix.\n");
	goto error;
      }
    }

  /* type of variables */
  if (columnType) 
    {
      status = CPXcopyctype (env, lp, columnType);
      if ( status ) {
	Scierror("Failed to copy ctype\n");
	goto error;
      }
    }

  if ( filename != NULL) 
    {
      status = CPXwriteprob (env, lp, filename, NULL);
      if ( status ) 
	{
	  Scierror("CPXwriteprob failed.\n");
	  goto error;
	}
    }
  
  if ( save_only == TRUE ) 
    {
      /* do not solve the pb just return */
      goto ok;
    }
  
  /* solves the pb */

  if ( columnType ) 
    {
      /* MIP LP or QP */
      status = CPXmipopt (env, lp);
      if ( status ) {
	Scierror("Failed to optimize MIP or QMIP.\n");
	goto error;
      }
    }
  else
    {
      if ( Qmatbeg == NULL )
	{
	  status = CPXlpopt (env, lp);
	  if ( status ) {
	    Scierror("Failed to optimize LP.\n");
	    goto error;
	  }
	}
      else
	{
	  status = CPXqpopt (env, lp);
	  if ( status ) {
	    Scierror("Failed to optimize QP.\n");
	    goto error;
	  }
	}
    }
  
  cur_numrows = CPXgetnumrows (env, lp);
  cur_numcols = CPXgetnumcols (env, lp);

  if ( cur_numrows != rowCount || cur_numcols != colCount) 
    {
      Scierror("Error: cplex pb has row=%d and column=%d while expecting row=%d and column=%d\n",
	       cur_numrows,cur_numcols,rowCount,colCount);
      goto error;
    }
  
  solstat = CPXgetstat (env, lp);
  status = CPXgetobjval (env, lp, &objval);
  if ( status ) {
    Scierror("objective value not available\n");
    goto error;
  }
  
  status = CPXgetx (env, lp, X->R, 0, cur_numcols-1);
  if ( status ) {
    Scierror("Failed to get optimal integer x.\n");
    goto error;
  }
  
  if (  columnType ) 
    {
      int i;
      /* no multipliers for mip */
      for ( i = 0 ; i < cur_numrows; i++)
	Lambda->R[i]=0.0;
    }
  else
    {
      /* first the inequality constraints */
      status = CPXgetpi (env, lp, Lambda->R, neq, cur_numrows-1);
      if ( status ) {
	Scierror("Failed to get optimal pi values.\n");
	goto error;
      }
      /* then the equality constraints */
      status = CPXgetpi (env, lp, Lambda->R+(cur_numrows-neq), 0, neq-1);
      if ( status ) {
	Scierror("Failed to get optimal pi values.\n");
	goto error;
      }
    }

 ok: 
  Retcode->R[0]= solstat; 
  RetCost->R[0]= objval;
  ret = OK ;
  
 error:
   /* Free up the solution */
   if ( dj != NULL) free(dj);
   if ( slack != NULL) free(slack);
   if ( lp != NULL ) 
     {
       status = CPXfreeprob (env, &lp);
       if ( status ) {
	 Scierror("CPXfreeprob failed, error code %d.\n", status);
	 ret = FAIL;
       }
     }
   /* Free up the CPLEX environment, if necessary */
   if ( env != NULL ) 
     {
       status = CPXcloseCPLEX (&env);
       /* Note that CPXcloseCPLEX produces no output,
	  so the only way to see the cause of the error is to use
	  CPXgeterrorstring.  For other CPLEX routines, the errors will
	  be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON. 
       */
       if ( status ) {
	 char  errmsg[CPXMESSAGEBUFSIZE];
	 Scierror("Could not close CPLEX environment.\n");
	 CPXgeterrorstring (env, status, errmsg);
	 Scierror("%s", errmsg);
	 ret = FAIL;
       }
     }
   return ret;
}

/* 
   method   is the optimization method
   o          default
   p          primal simplex
   d          dual   simplex
   n          network with dual simplex cleanup
   h          barrier with crossover
   b          barrier without crossover
   s          sifting
   c          concurrent
*/

int nsp_cplex_set_method(CPXENVptr env, char c)
{
  int method, status;
  switch (c) {
  case 'o':
    method = CPX_ALG_AUTOMATIC;
    break;
  case 'p':
    method = CPX_ALG_PRIMAL;
    break;
  case 'd':
    method = CPX_ALG_DUAL;
    break;
  case 'n':
    method = CPX_ALG_NET;
    break;
  case 'h':
    method = CPX_ALG_BARRIER;
    break;
  case 'b':
    method = CPX_ALG_BARRIER;
    status = CPXsetintparam (env, CPXPARAM_Barrier_Crossover, CPX_ALG_NONE);
    if ( status ) {
      Scierror("Failed to set the crossover method, error %d.\n", status);
      return FAIL;
    }
    break;
  case 's':
    method = CPX_ALG_SIFTING;
    break;
  case 'c':
    method = CPX_ALG_CONCURRENT;
    break;
  default:
    method = CPX_ALG_NONE;
    break;
  }
  status = CPXsetintparam (env, CPXPARAM_LPMethod, method);
  if ( status ) {
    Scierror("Failed to set the optimization method, error %d.\n", status);
    return FAIL;
  }
  return OK;
}


static int nsp_cplex_redirect_channels(CPXENVptr env,int loglevel)
{
  char          errmsg[CPXMESSAGEBUFSIZE];
  CPXCHANNELptr  cpxerror   = NULL;
  CPXCHANNELptr  cpxwarning = NULL;
  CPXCHANNELptr  cpxresults = NULL;
  CPXCHANNELptr  cpxlog = NULL;

  static char errorlabel[] = "cpxerror";
  static char warnlabel[]  = "cpxwarning";
  static char reslabel[]   = "cpxresults";
  static char loglabel[]   = "cpxlog";

  int status;
  
   /* Now get the standard channels.  If an error, just call our
    * message function directly. 
    */
   
  status = CPXgetchannels (env, &cpxresults, &cpxwarning, &cpxerror, &cpxlog);
  if ( status ) {
    Scierror("Could not get standard channels.\n");
    CPXgeterrorstring (env, status, errmsg);
    Scierror("%s\n",errmsg);
    goto error;
  }

  status = CPXaddfuncdest (env, cpxerror, errorlabel, nsp_cplex_error);
  if ( status ) {
    Scierror("Could not set up error message handler.\n");
    CPXgeterrorstring (env, status, errmsg);
    Scierror("%s\n",errmsg);
    goto error;
  }

  /* Now that we have the error message handler set up, all CPLEX
     generated errors will go through nsp_cplex_msg.  So we don't have
     to use CPXgeterrorstring to determine the text of the message.
     We can also use CPXmsg to do any other printing.  */

  status = CPXaddfuncdest (env, cpxwarning, warnlabel, nsp_cplex_msg);
  if ( status ) {
    Scierror("Failed to set up handler for cpxwarning.\n");
    goto error;
  }

  if ( loglevel > 0 ) 
    {
      status = CPXaddfuncdest (env, cpxresults, reslabel, nsp_cplex_msg);
      if ( status ) {
	Scierror("Failed to set up handler for cpxresults.\n");
	goto error;
      }

      status = CPXaddfuncdest (env, cpxlog, loglabel, nsp_cplex_msg);
      if ( status ) {
	Scierror("Failed to set up handler for cpxresults.\n");
	goto error;
      }
    }

  /* Now turn on the iteration display. */
  status = CPXsetintparam (env, CPXPARAM_Simplex_Display, loglevel);
  if ( status ) {
    Scierror("Failed to turn on simplex display level.\n");
    goto error;
  }
   
  return OK;
   
 error:
   
  if ( cpxresults != NULL ) {
    int  chanstat;
    chanstat = CPXdelfuncdest (env, cpxresults, reslabel, nsp_cplex_msg);
    if ( chanstat && !status ) {
      Scierror("Failed to delete cpxresults function.\n");
    }
  }
   
  if ( cpxwarning != NULL ) {
    int  chanstat;
    chanstat = CPXdelfuncdest (env, cpxwarning, warnlabel, nsp_cplex_msg);
    if ( chanstat && !status ) {
      Scierror("Failed to delete cpxwarning function.\n");
    }
  }

  if ( cpxlog != NULL ) 
    {
      int  chanstat;
      chanstat = CPXdelfuncdest (env, cpxlog, loglabel, nsp_cplex_msg);
      if ( chanstat && !status ) {
	Scierror("Failed to delete cpxlog function.\n");
      }
    }

  if ( cpxerror != NULL ) {
    int  chanstat;
    chanstat = CPXdelfuncdest (env, cpxerror, errorlabel, nsp_cplex_error);
    if ( chanstat && !status ) {
      Scierror("Failed to delete cpxerror function.\n");
    }
  }
  return FAIL;
} 

static void CPXPUBLIC nsp_cplex_msg(void *handle, const char *msg)
{
  /* char const *label = (char const *) handle; */
  Sciprintf("%s",msg);
} 

static void CPXPUBLIC nsp_cplex_error(void *handle, const char *msg)
{
  /* char const *label = (char const *) handle;*/
  Scierror("Error:%s",msg);
} 



