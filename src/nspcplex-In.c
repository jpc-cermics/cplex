/* Nsp
 * Copyright (C) 2014-2022 J.-Ph. Chancelier Cermics/ENPC 
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
 * Interface for calling cplex solver 
 */

#include <nsp/nsp.h>
#include <nsp/objects.h>
#include <nsp/imatrix.h>
#include <nsp/spcolmatrix.h>
#include <nsp/interf.h>

#include "nspcplex.h"

int int_cplex_solve(Stack stack, int rhs, int opt, int lhs)
{
  int save_only=FALSE;
  char * filename = NULL;
  int rep;
  int loglevel = 0;
  const double cplex_dbl_max= nsp_cplex_dbl_max();
  char *sense_str = "min";
  NspMatrix *X,*Lambda, *Retcode, *RetCost;
  NspIMatrix *Cmatbeg=NULL,*Cmatind=NULL,*Cmatcount=NULL;
  NspMatrix *Cmatval=NULL;
  NspObject *ObjA,*ObjAe;
  NspMatrix *Objective, *Rhs, *Rhse, *B, *Lhs, *lb=NULL, *ub=NULL;
  NspMatrix *SemiCont = NULL;

  NspIMatrix *Qmatbeg=NULL,*Qmatcnt=NULL,*Qmatind=NULL;
  NspMatrix *Qmatval=NULL;
  NspObject *ObjQ=NULL;

  NspHash *Options = NULLHASH;
  NspSMatrix *var_type = NULLSMAT;
  int neq=0, ncols,nrows,i, sense=0;
  nsp_string columnType= NULL, rowType= NULL;

  int_types T[] = {realmat, obj , realmat, obj , realmat, new_opts, t_end} ;
  
  nsp_option opts[] ={
    {"Q", obj,NULLOBJ,-1},
    {"lb",realmatcopy,NULLOBJ,-1},
    {"ub",realmatcopy,NULLOBJ,-1},
    {"sense",string,NULLOBJ,-1},
    {"var_type",smat,NULLOBJ,-1},
    {"options", hash, NULLOBJ,-1},
    {"semi_cont",realmatcopy,NULLOBJ,-1},
    {"loglevel", s_int,NULLOBJ, -1},
    {"save_only",s_bool,NULLOBJ,-1},
    {"fname",string,NULLOBJ,-1},
    { NULL,t_end,NULLOBJ,-1}};
  
  if ( GetArgs(stack,rhs,opt,T,&Objective, &ObjA, &Rhs, &ObjAe, &Rhse, &opts,&ObjQ,
	       &lb,&ub,&sense_str,&var_type,&Options,&SemiCont,&loglevel,&save_only,&filename) == FAIL) 
    return RET_BUG;
  /* 
  if ( get_solver_options(stack, Options, &options) == FAIL )
    return RET_BUG;
  */

  /* optional argument sense : min or max default min */
  if ( strcmp(sense_str,"min") == 0 )
    sense = 0;
  else if ( strcmp(sense_str,"max") == 0 )
    sense = 1;
  else
    {
      Scierror("Error: sense should be 'min' or 'max'\n");
      return RET_BUG;
    }

  /* Matrix A and Ae : both sparses or both not sparses */
  
  if ( IsMat(ObjA) &&  IsMat(ObjAe))
    {
      NspMatrix *Ae;
      if ( ((NspMatrix *) ObjA)->m != Rhs->mn || ((NspMatrix *) ObjAe)->m != Rhse->mn )
	{
	  Scierror("Error: incompatible dimensions between matrices and Rhs\n",NspFname(stack));
	  return RET_BUG;
	}
      neq = ((NspMatrix *) ObjAe)->m;
      Ae = (neq == 0) ? NULL : (NspMatrix *) ObjAe;
      if ( nsp_matrix_to_sparse_triplet((NspMatrix *)ObjA,Ae, &Cmatbeg,&Cmatind,&Cmatval) == FAIL)
	return RET_BUG;
    }
  else if ( IsSpColMat(ObjA) && IsSpColMat(ObjAe) )
    {
      NspSpColMatrix *Ae;
      if ( ((NspSpColMatrix *) ObjA)->m != Rhs->mn || ((NspSpColMatrix *) ObjAe)->m != Rhse->mn )
	{
	  Scierror("Error: incompatible dimensions between matrices and Rhs\n",NspFname(stack));
	  return RET_BUG;
	}
      neq = ((NspSpColMatrix *) ObjAe)->m;
      Ae = (neq == 0) ? NULL : (NspSpColMatrix *) ObjAe;
      /* Take care that equality constraints are placed at the begining */
      if ( nsp_spcolmatrix_to_sparse_triplet((NspSpColMatrix *)ObjA,Ae, &Cmatbeg,&Cmatind,&Cmatval) == FAIL)
	return RET_BUG;
    }
  else
    {
      Scierror("Error: second and fourth arguments of function %s should be both full or both sparse matrices\n",
	       NspFname(stack));
      return RET_BUG;
    }

  ncols = Objective->mn; /* Length of c == number of columns*/
  nrows = Rhs->mn + Rhse->mn ; /* length of b == number of rows*/

  /* SemiCont : to be done */
  
  if ( SemiCont != NULL) 
    {
      Scierror("Error: not yet implemented in cplex interface \n");
      return RET_BUG;
      int *Sc = (int *) SemiCont->R;
      for ( i= 0 ; i < SemiCont->mn ; i++)
	{
	  if ( SemiCont->R[i] < 1 || SemiCont->R[i] > ncols) 
	    {
	      Scierror("Error: semi-cont index %d is not in the range [1,%d]\n",i,ncols);
	      return RET_BUG;
	    }
	}
      for ( i= 0 ; i < SemiCont->mn ; i++) Sc[i]= SemiCont->R[i] -1 ;
    }
  
  /* Cmatcount: extra matrix which counts number of non-null elements in each column 
   */

  if ( ( Cmatcount = nsp_imatrix_create(NVOID,1,ncols,nsp_gint32)) == NULLIMAT )
    {
      Scierror("Error: running out of memory\n");
      return RET_BUG;
    }

  for (i = 0; i < ncols ; i++) 
    {
      Cmatcount->Gint[i] = Cmatbeg->Gint[i+1] - Cmatbeg->Gint[i];
    }

  if ( var_type != NULL )
    {
      if ( var_type->mn != ncols ) 
	{
	  Scierror("Error: var_type should be of size %d\n", ncols ); 
	  return RET_BUG;
 	}
      if (( columnType =new_nsp_string_n(ncols+1)) == (nsp_string) 0)
	{
	  Scierror("Error: running out of memory\n");
	  return RET_BUG;
	}
      for (i = 0 ; i < ncols ; i++)
	{
	  if ( strlen(var_type->S[i]) > 0 ) 
	    columnType[i]= var_type->S[i][0];
	  else
	    columnType[i]= 'C';
	}
      columnType[ncols]='\0';
    }

  /* rowType should be of size rowcount and should contain
   * 'L', 'E', 'G', 'R', 'N' 
   * L: (Ax)_i <= b_i 
   * G: (Ax)_i >= b_i
   * E: (Ax)_i == b_i  
   * R: b_i -Abs(range_i) <= (Ax)_i <= b_i 
   * N: no constraint 
   */

  if (( rowType =new_nsp_string_n(nrows+1)) == (nsp_string) 0)
    {
      Scierror("Error: running out of memory\n");
      return RET_BUG;
    }
  for (i = 0 ; i < nrows ; i++)
    {
      if ( i < Rhse->mn) 
	rowType[i]= 'E';
      else
	rowType[i]= 'L';
    }
  rowType[nrows]='\0';

  /* Check that ObjA and ObjAe are compatible with Rhs and Rhse */

  /* agregates the Rhs */
  if ( ( B = nsp_matrix_create(NVOID,'r',1,nrows)) == NULLMAT) return RET_BUG;
  memcpy(B->R,Rhse->R,Rhse->mn*sizeof(double));
  memcpy(B->R+Rhse->mn,Rhs->R,Rhs->mn*sizeof(double));
  
  /* Create lower bounds if not available */
  
  if ( lb == NULL ) 
    {
      if (( lb =  nsp_matrix_create(NVOID,'r',1,ncols)) == NULL)
	return RET_BUG;
      for (i = 0; i < ncols; i++){
	lb->R[i] = 0; /* to fit with glpk - cplex_dbl_max; */
      }
    }
  else
    {
      if ( lb->mn != ncols) 
	{
	  Scierror("Error: lb should be of size %d\n",ncols);
	  return RET_BUG;
	}
      for (i = 0; i < lb->mn ; i++)
	{
	  if ( isinf(lb->R[i]) != 0 ) lb->R[i] = - cplex_dbl_max; 
	}
    }
  
  if ( ub == NULL ) 
    {
      if (( ub =  nsp_matrix_create(NVOID,'r',1,ncols)) == NULL)
	return RET_BUG;
      for (i = 0; i < ncols; i++){
	ub->R[i] = cplex_dbl_max;
     }        
    }
  else
    {
      if ( ub->mn != ncols) 
	{
	  Scierror("Error: ub should be of size %d\n",ncols);
	  return RET_BUG;
	}
      for (i = 0; i < ub->mn; i++)
	{
	  if ( isinf(ub->R[i]) != 0 ) ub->R[i] = cplex_dbl_max;
	}
    }

  if (( Lhs =  nsp_matrix_create(NVOID,'r',nrows,1)) == NULL)
    return RET_BUG;
  
  for (i = 0; i < neq; i++)
    {
      Lhs->R[i] = B->R[i];
    }
  for (i = neq; i < nrows; i++)
    {
      Lhs->R[i] = -cplex_dbl_max;
    }

  /* do we have a quadratic cost */
  
  if ( ObjQ != NULL )
    {
      if ( IsMat(ObjQ) ) 
	{
	  if ( ((NspMatrix *) ObjQ)->m != ncols && ((NspMatrix *) ObjQ)->n != ncols )
	    {
	      Scierror("Error: optional argument Q of function %s should be of size %dx%d\n",
		       NspFname(stack),ncols,ncols);
	      return RET_BUG;
	    }
	  if ( nsp_matrix_to_sparse_triplet((NspMatrix *)ObjQ,NULL, &Qmatbeg,&Qmatind,&Qmatval) == FAIL)
	    return RET_BUG;
	}
      else if ( IsSpColMat(ObjQ) ) 
	{
	  if ( ((NspSpColMatrix *) ObjQ)->m != ncols && ((NspSpColMatrix *) ObjQ)->n != ncols )
	    {
	      Scierror("Error: optional argument Q of function %s should be of size %dx%d\n",
		       NspFname(stack),ncols,ncols);
	      return RET_BUG;
	    }
	  if ( nsp_spcolmatrix_to_sparse_triplet((NspSpColMatrix *)ObjQ,NULL, &Qmatbeg,&Qmatind,&Qmatval) == FAIL)
	    return RET_BUG;
	}
      else
	{
	  Scierror("Error: optional argument Q of function %s should be a real full or sparse matrix\n",
		   NspFname(stack));
	  return RET_BUG;
	}

      /* extra matrix */
      if ( ( Qmatcnt = nsp_imatrix_create(NVOID,1,ncols,nsp_gint32)) == NULLIMAT )
	{
	  Scierror("Error: running out of memory\n");
	  return RET_BUG;
	}
      for (i = 0; i < ncols ; i++) 
	{
	  Qmatcnt->Gint[i] = Qmatbeg->Gint[i+1] - Qmatbeg->Gint[i];
	}
    }
  
  if (( X= nsp_matrix_create(NVOID,'r', ncols,1)) == NULL ) 
    return RET_BUG;
  if (( Lambda= nsp_matrix_create(NVOID,'r', nrows,1)) == NULL ) 
    return RET_BUG;
  if (( Retcode= nsp_matrix_create(NVOID,'r', 1,1)) == NULL ) 
    return RET_BUG;
  if (( RetCost= nsp_matrix_create(NVOID,'r', 1,1)) == NULL ) 
    return RET_BUG;

  rep = nsp_cplex_solve("Pb", sense, ncols, nrows, neq,Cmatbeg,Cmatcount,Cmatind, Cmatval,
			lb,ub,Objective, Qmatbeg, Qmatcnt, Qmatind, Qmatval,
			B, columnType,  X, Lambda,RetCost, Retcode,rowType, 
			(SemiCont != NULL) ? SemiCont->mn : 0, 
			(SemiCont != NULL) ? (int *) SemiCont->R: NULL, Options,loglevel, 
			filename, save_only);

  /* destroy allocated */

  if ( Cmatbeg !=NULL) nsp_imatrix_destroy(Cmatbeg);
  if ( Cmatcount !=NULL) nsp_imatrix_destroy(Cmatcount);
  if ( Cmatind !=NULL) nsp_imatrix_destroy(Cmatind);
  if ( Cmatval !=NULL) nsp_matrix_destroy(Cmatval);

  if ( Qmatbeg !=NULL) nsp_imatrix_destroy(Qmatbeg);
  if ( Qmatcnt !=NULL) nsp_imatrix_destroy(Qmatcnt);
  if ( Qmatind !=NULL) nsp_imatrix_destroy(Qmatind);
  if ( Qmatval !=NULL) nsp_matrix_destroy(Qmatval);

  nsp_matrix_destroy(B);
  nsp_matrix_destroy(Lhs);
  nsp_matrix_destroy(lb);
  nsp_matrix_destroy(ub);
  if ( SemiCont != NULL) nsp_matrix_destroy(SemiCont);

  if ( columnType != NULL) nsp_string_destroy(&columnType);
  if ( rowType != NULL) nsp_string_destroy(&columnType);

  if ( save_only == TRUE ) 
    {
      nsp_matrix_destroy(X);
      nsp_matrix_destroy(RetCost);
      nsp_matrix_destroy(Retcode);
      nsp_matrix_destroy(Lambda);
      if ( lhs > 0) 
	{
	  Scierror("Error: when save_only is set to %t then lhs should be 0\n");
	  return RET_BUG;
	}
      return 0;
    }

  if ( rep == FAIL ) 
    {
      nsp_matrix_destroy(X);
      nsp_matrix_destroy(RetCost);
      nsp_matrix_destroy(Retcode);
      nsp_matrix_destroy(Lambda);
      return RET_BUG;
    }
  else
    {
      MoveObj(stack,1,NSP_OBJECT(X));
      if ( lhs >= 2) 
	MoveObj(stack,2,NSP_OBJECT(RetCost));
      else
	nsp_matrix_destroy(RetCost);
      if ( lhs >= 3)
	MoveObj(stack,3,NSP_OBJECT(Retcode));
      else
	nsp_matrix_destroy(Retcode);
      if ( lhs >= 4)
	MoveObj(stack,4,NSP_OBJECT(Lambda)); 
      else
	nsp_matrix_destroy(Lambda);
    }
  return Max(lhs,1);
}

static OpTab libcplex_func[] = {
  {"linprog_cplex", int_cplex_solve},
  {(char *) 0, NULL},
};

int libcplex_Interf (int i, Stack stack, int rhs, int opt, int lhs)
{
  return (*(libcplex_func[i].fonc)) (stack, rhs, opt, lhs);
}

void libcplex_Interf_Info (int i, char **fname, function (**f))
{
 *fname = libcplex_func[i].name;
 *f = libcplex_func[i].fonc;
}
