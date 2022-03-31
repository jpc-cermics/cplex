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

#ifndef NSP_CPLEX_CPP_H 
#define NSP_CPLEX_CPP_H 

#ifdef __cplusplus
extern "C" {
#endif 

#include <nsp/nsp.h>
#include <nsp/objects.h>
#include <nsp/imatrix.h>
  
  extern double nsp_cplex_dbl_max();
  
  extern int nsp_cplex_solve(const char* problemName, int sense, int ncols, int nrows, int neq,
			     NspIMatrix*Cmatbeg, NspIMatrix *Cmatcount, NspIMatrix *Cmatind, NspMatrix *Cmatval, 
			     NspMatrix *lower, NspMatrix *upper, NspMatrix *Objective,
			     NspIMatrix *Qmatbeg,NspIMatrix *Qmatcnt, NspIMatrix *Qmatind, NspMatrix *Qmatval, 
			     NspMatrix *Rhs,const char *columnType,  NspMatrix *X,NspMatrix *Lambda,
			     NspMatrix *RetCost,NspMatrix *Retcode,const char *rowType,
			     int semiCount, int *semiIndex,NspHash *Options,int loglevel,
			     const char *filename,int save_only);


#ifdef __cplusplus
}
#endif 

#endif /*NSP_CPLEX_CPP_H */


