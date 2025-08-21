// -*- Mode: nsp -*-
// Copyright (C) 2014-2022 Bruno Pinçon ESIAL/IECN
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

load_toolbox('cplex');

/////////////////////////////////////////////////////////////////////////////////
// test 1-1:   max c'*x
//             Ax <= b
//              x >= 0
A = [ 2 4 8  6;
     10 8 6 10;
      1 1 2  2];
b = [100; 160;  20];
c= [50; 40; 70; 80];

[xopt,fopt,flag,lambda] = linprog_cplex(c,A,b,[],[],sense="max",loglevel=0);

// exact sol, function value and basis information
XM = [12;0;0;4]; FM = 920; final_basis = hash(aux = [1;3;3], str = [1;2;2;1]);

if norm(xopt-XM)/norm(XM) >= 4*%eps then, pause, end
if abs(fopt-FM)/abs(FM) >= 4*%eps then, pause, end
if abs((dot(lambda,b)-fopt)/fopt) >= 4*%eps then, pause, end

/////////////////////////////////////////////////////////////////////////////
// test 2-1  (unbounded solution)
//   max c'x
//   Ax <= b
//    x >= 0
A= [-1 -1; -1 2; -2 1];
b=[-3; -5; 5];
c=[1; 3];
ok=execstr('[xopt,fopt,flag] = linprog_cplex(c,A,b,[],[],sense=""max"");',errcatch=%t);
if ok then pause;else lasterror();end 

//////////////////////////////////////////////////////////////////////////////
// tests on easy mips (from netlib)

// bal8x12.mps is a mip with x >= 0 (no need to provide lb=0)
[c,A,b,Ae,be,sense,lb,ub,binprog,intprog,var_type] = readlp("NSP/tests/bal8x12.mps",verb=0);
Fe = 471.55; 
[xopt,fopt,flag] = linprog_cplex(c,A,b,Ae,be,ub=ub,var_type=var_type);
if abs((fopt-Fe)/Fe) >= 4*%eps then, pause, end

// gr4x6.mps is a mip with x >= 0 (no need to provide lb=0)
[c,A,b,Ae,be,sense,lb,ub,binprog,intprog,var_type] = readlp("NSP/tests/gr4x6.mps",verb=0);
Fe = 202.35; 
[xopt,fopt,flag] = linprog_cplex(c,A,b,Ae,be,ub=ub,var_type=var_type);
if abs((fopt-Fe)/Fe) >= 4*%eps then, pause, end

// bk4x3.mps is a mip with x >= 0 (no need to provide lb=0)
[c,A,b,Ae,be,sense,lb,ub,binprog,intprog,var_type] = readlp("NSP/tests/bk4x3.mps",verb=0);
Fe = 350.0; 
[xopt,fopt,flag] = linprog_cplex(c,A,b,Ae,be,ub=ub,var_type=var_type);
if abs((fopt-Fe)/Fe) >= 4*%eps then, pause, end


// basic test for quadratic programming 
//     Maximize
//       obj: x1 + 2 x2 + 3 x3
//              - 0.5 ( 33x1*x1 + 22*x2*x2 + 11*x3*x3
//                   -  12*x1*x2 - 23*x2*x3 )
//      Subject To
//       c1: - x1 + x2 + x3 <= 20
//       c2: x1 - 3 x2 + x3 <= 30
//      Bounds
//       0 <= x1 <= 40
//      End
// this is a test case of cplex 

c=[1,2,3];
A=sparse([-1,1,1;1,-3,1]);
b=[20;30];

lb=[0,-%inf,-%inf];
ub=[40,%inf,%inf];

Q= -[ 33,-06,00;
      -06,22,-23/2;
      00,-23/2,11];

sol= 2.015617;
xopt=[0.139115;0.598465;0.898396];
lambda=[18.642254;30.757886];

[xopt1,fopt1,flag1,lambda1] = linprog_cplex(c,A,b,sparse([]),[],ub=ub,lb=lb,sense="max",Q=Q);

if (fopt1- 2.015617) > 1.e-5 then pause;end
if norm(xopt1-xopt1) > 1.e-5 then pause;end
// if norm(lambda1-lambda) > 1.e-5 then pause;end ?

