// -*- Mode: nsp -*-
// Copyright (C) 2014-2022 Jean-Philippe Chancelier Enpc/Cermics
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
//
// a set of examples compared with linprog 

load_toolbox('nspcplex');

// example 1 
//-----------

c= [50; 40; 70; 80]; 
A = [ 2 4 8  6; 
      10 8 6 10; 
      1 1 2  2]; 
b = [100; 160;  20]; 

[xopt,fopt,flag,extra] = linprog(c,A,b,[],[],sense="max");
// Take care that default value for lb is 0 and ub = %inf 

// check optimality conditions 
L=extra.lambda;
T= and(-c + A'*L >= 0) && (-c+A'*L)'*xopt == 0 && L'*(A*xopt-b) == 0 ;
if ~T then pause;end 
// check optimal cost and dual cost 
if L'*b <> c'*xopt then pause;end
// check optimal cost
if fopt <> c'*xopt then pause;end 
if fopt <> 920 then pause;end 

[xopt1,fopt1,flag1,extra1] = linprog_cplex(c,A,b,[],[],sense="max");
// check optimality conditions 
L=extra1;
T= and(-c + A'*L >= 0) && (-c+A'*L)'*xopt1 == 0 && L'*(A*xopt1-b) == 0 ;
if ~T then pause;end 
// check optimal cost and dual cost 
if L'*b <> c'*xopt1 then pause;end
if norm(xopt-xopt1) >= 1.e-8 then pause;end 
if norm(extra1- extra.lambda) >= 1.e-8 then pause;end 

// example 2 
//---------- 
c = -[6;5];                
A = [1,4; 6,4; 2,-5];      
b = [16;28;6];    
lb = [0;0];                
ub = [10;10];

[xopt,fopt,flag,extra] = linprog(c,A,b,[],[],ub=ub,lb=lb,sense="min");
// check optimality conditions 
L=extra.lambda;
T= and(-c + A'*L >= 0) && abs((-c+A'*L)'*xopt) < 1.e-8 && abs(L'*(A*xopt-b)) < 1.e-8 ;
if ~T then pause;end 
// check optimal cost and dual cost 
if abs(L'*b - c'*xopt) > 1.e-8 then pause;end
[xopt1,fopt1,flag1,extra1] = linprog_cplex(c,A,b,[],[],ub=ub,lb=lb,sense="min");
if norm(xopt-xopt1) >= 1.e-8 then pause;end 
if norm(extra1- extra.lambda) >= 1.e-8 then pause;end 

// example 3
//------------
c = -[1 2 3];
A = [-1 , 1 , 1; 
     1 , -3 , 1];
Ae=[1  1  1];
b = [20;30];
be=40;
lb = [0;0;0];
ub = [40;%inf;%inf];

[xopt,fopt,flag,extra] = linprog(c,A,b,Ae,be,ub=ub,lb=lb,sense="min");
[xopt1,fopt1,flag1,extra1] = linprog_cplex(c,A,b,Ae,be,ub=ub,lb=lb,sense="min");

if norm(xopt-xopt1) >= 1.e-8 then pause;end 
if norm(extra1- extra.lambda) >= 1.e-8 then pause;end 

// We code here in nsp the problems used for coinmp test
//-----------------------------------------------------
// CoinTest 
//----------

objectConst = 0.0;
n=8;
m=5;
c = ones(1,8);
lb = 0*ones(1,n);
ub = 1000000*ones(1,n);
// char rowType[5] = [ 'L', 'L', 'L', 'L', 'L' ];
b = [14., 80., 50., 50., 50.];
beg=[0,2,4,6,8,10,11,12,14];
count=[2,2,2,2,2,1,1,2];
ind=[0,4,0,1,1,2,0,3,0,4,2,3,0,4];
val=[3., 5.6, 1., 2., 1.1, 1., -2., 2.8, -1., 1., 1., -1.2, -1., 1.9];
A=spfrommtlb(beg,ind,val,[m,n]);

optimalValue = 1428729.2857143;

[xopt,fopt,flag,extra] = linprog(c,A,b,[],[],ub=ub,lb=lb,sense="max");
if abs(fopt - optimalValue) > 1.e-7 then pause;end

[xopt1,fopt1,flag1,extra1] = linprog_cplex(c,A,b,sparse([]),[],ub=ub,lb=lb,sense="max");
if norm(xopt-xopt1) >= 1.e-8 then pause;end 
if abs(fopt1 - optimalValue) > 1.e-7 then pause;end
if norm(extra1- extra.lambda) >= 1.e-8 then pause;end 

// "Bakery"
// --------

n=2;
m=3;
cte=  - 4000.0 / 30.0;
c =[ 0.05 , 0.08];
lb =[ 0, 0 ];
ub =[ 1000000, 1000000 ];
b =[1400 , 8000 , 5000 ];

beg =  [ 0 , 2, 4 ];
ind =  [ 0, 1, 0, 2];
val =  [ 0.1, 1, 0.2, 1];
A=spfrommtlb(beg,ind,val,[m,n]);

optimalValue = 506.66666667 -cte ;

[xopt,fopt,flag,extra] = linprog(c,A,b,[],[],ub=ub,lb=lb,sense="max");
if abs(fopt - optimalValue) > 1.e-7 then pause;end

[xopt1,fopt1,flag1,extra1] = linprog_cplex(c,A,b,sparse([]),[],ub=ub,lb=lb,sense="max");
if norm(xopt-xopt1) >= 1.e-8 then pause;end 
if abs(fopt1 - optimalValue) > 1.e-7 then pause;end
if norm(extra1- extra.lambda) >= 1.e-8 then pause;end 

// Afiro
//-------
n = 32;
m = 27;
sense = "min";
c =[0, -0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.32, 0, 0, 0, -0.6, ... 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.48, 0, 0, 10];
lb = zeros(1,n);
ub = %inf*ones(1,n);

ct = ['E', 'E', 'L', 'L', 'E', 'E', 'L', 'L', 'L', 'L', 'E', 'E', 'L', ...
      'L', 'E', 'E', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L'];

b = [0, 0, 80, 0, 0, 0, 80, 0, 0, 0, 0, 0, 500, 0, 0, 44, 500, 0, ...
		0, 0, 0, 0, 0, 0, 0, 310, 300];
beg = [0, 4, 6, 8, 10, 14, 18, 22, 26, 28, 30, 32, 34, 36, 38, 40, ...
       44, 46, 48, 50, 52, 56, 60, 64, 68, 70, 72, 74, 76, 78, 80, 82, ...
       83];
ind=[0, 1, 2, 23, 0, 3, 0, 21, 1, 25, 4, 5, 6, 24, 4, 5, 7, 24, 4, 5, ...
     8, 24, 4, 5, 9, 24, 6, 20, 7, 20, 8, 20, 9, 20, 3, 4, 4, 22, 5, 26, 10, 11, ...
     12, 21, 10, 13, 10, 23, 10, 20, 11, 25, 14, 15, 16, 22, 14, 15, 17, 22, 14, ...
     15, 18, 22, 14, 15, 19, 22, 16, 20, 17, 20, 18, 20, 19, 20, 13, 15, 15, 24, ...
     14, 26, 15];

val =[-1, -1.06, 1, 0.301, 1, -1, 1, -1, 1, 1, -1, -1.06, 1, 0.301, ...
      -1, -1.06, 1, 0.313, -1, -0.96, 1, 0.313, -1, -0.86, 1, 0.326, -1, 2.364, -1, ...
      2.386, -1, 2.408, -1, 2.429, 1.4, 1, 1, -1, 1, 1, -1, -0.43, 1, 0.109, 1, -1, ...
      1, -1, 1, -1, 1, 1, -0.43, 1, 1, 0.109, -0.43, 1, 1, 0.108, -0.39, 1, 1, ...
      0.108, -0.37, 1, 1, 0.107, -1, 2.191, -1, 2.219, -1, 2.249, -1, 2.279, 1.4, ...
      -1, 1, -1, 1, 1, 1];

optimalValue = -464.753142857;

A=spfrommtlb(beg,ind,val,[m,n]);
Eq=find(ct == 'E');
Ae=A(Eq,:);be=b(Eq);
Lq=find(ct == 'L');
Al=A(Lq,:);bl=b(Lq);

[xopt,fopt,flag,extra] = linprog(c,A,b,[],[],ub=ub,lb=lb,sense=sense);
if abs(fopt - optimalValue) > 1.e-7 then pause;end

[xopt1,fopt1,flag1,extra1] = linprog_cplex(c,A,b,sparse([]),[],ub=ub,lb=lb,sense=sense);
if norm(xopt-xopt1) >= 1.e-8 then pause;end 
if abs(fopt1 - optimalValue) > 1.e-7 then pause;end
// XXX if norm(extra1- extra.lambda) >= 1.e-8 then pause;end 

// [xopt1,fopt1,flag1,extra2] = linprog_clp(c,A,b,sparse([]),[],ub=ub,lb=lb,sense=sense);
// if norm(xopt-xopt1) >= 1.e-8 then pause;end 
// if abs(fopt1 - optimalValue) > 1.e-7 then pause;end
// if norm(extra2- extra.lambda) >= 1.e-8 then pause;end 

[xopt,fopt,flag,extra] = linprog(c,Al,bl,Ae,be,ub=ub,lb=lb,sense=sense);
if abs(fopt - optimalValue) > 1.e-7 then pause;end

[xopt1,fopt1,flag1,extra1] = linprog_cplex(c,Al,bl,Ae,be,ub=ub,lb=lb,sense=sense);
if norm(xopt-xopt1) >= 1.e-8 then pause;end 
if abs(fopt1 - optimalValue) > 1.e-7 then pause;end
// if norm(extra1- extra.lambda) >= 1.e-8 then pause;end 

// P0033
//--------

n = 33;
m = 15;
sense = "min";
c = [171, 171, 171, 171, 163, 162, 163, 69, 69, 183, 183, 183, ...
     183, 49, 183, 258, 517, 250, 500, 250, 500, 159, 318, 159, 318, 159, 318, 159, ...
     318, 114, 228, 159, 318];
lb= zeros(1,n);
ub= ones(1,n);

b=[1, 1, 1, 1, -5, 2700, -2600, -100, -900, -1656, -335, -1026, -5, -500, -270];

beg=[0, 3, 6, 10, 14, 19, 24, 26, 31, 36, 38, 41, 45, 49, 53, 54, ...
     55, 56, 58, 60, 61, 62, 66, 70, 73, 76, 80, 84, 87, 90, 93, 96, 97, 98];
ind=[0, 8, 9, 0, 12, 13, 0, 5, 6, 9, 0, 5, 6, 7, 1, 5, 6, 10, 11, 1, ...
		5, 6, 8, 9, 1, 14, 2, 5, 6, 10, 11, 2, 5, 6, 8, 9, 3, 4, 3, 10, 11, 3, 5, 6, ...
		11, 3, 5, 6, 9, 5, 6, 8, 9, 3, 4, 4, 12, 13, 12, 13, 13, 13, 5, 6, 10, 11, 5, ...
		6, 10, 11, 5, 6, 11, 5, 6, 11, 5, 6, 8, 9, 5, 6, 8, 9, 5, 6, 9, 5, 6, 9, 5, 6, ...
		7, 5, 6, 7, 14, 14];
val=[1, -300, -300, 1, -300, -300, 1, 300, -300, -300, 1, 300, ...
     -300, -300, 1, 285, -285, -285, -285, 1, 285, -285, -285, -285, 1, -285, 1, ...
     265, -265, -265, -265, 1, 265, -265, -265, -265, 1, -230, 1, -230, -230, 1, ...
     230, -230, -230, 1, 230, -230, -230, 190, -190, -190, -190, 1, -200, -400, ...
     -200, -200, -400, -400, -200, -400, 200, -200, -200, -200, 400, -400, -400, ...
     -400, 200, -200, -200, 400, -400, -400, 200, -200, -200, -200, 400, -400, ...
     -400, -400, 200, -200, -200, 400, -400, -400, 200, -200, -200, 400, -400, ...
     -400, -200, -400];

ctyp = smat_create(1,n,"B");
A=spfrommtlb(beg,ind,val,[m,n]);

Ae=sparse([]);
be =[];
optimalValue = 3089.0;
ctyp = smat_create(1,n,"B"); 

[xopt,fopt,flag] = linprog(c,A,b,Ae,be,binprog=%t,sense=sense);
if abs(fopt - optimalValue) > 1.e-7 then pause;end

[xopt1,fopt1,flag1,extra1] = linprog_cplex(c,A,b,Ae,be,var_type=ctyp, sense=sense);
// if norm(xopt-xopt1) >= 1.e-8 then pause;end 
if abs(fopt1 - optimalValue) > 1.e-7 then pause;end

// Exmip1
//--------------------------

n = 8;
m= 5;
// objectname = "z";
sense = "min";
objconst = 0.0;
c=[1, 0, 0, 0, 2, 0, 0, -1];
lb=[2.5, 0, 0, 0, 0.5, 0, 0, 0];
ub=[%inf, 4.1, %inf, %inf, 4, %inf, %inf, 4.3];
// char rtyp[5]= ['G', 'L', 'E', 'G', 'L'];
b=[2.5, 2.1, 4, 1.8, 15];
//drng[5]=[0, 0, 0, -3.2, 12];

beg=[0, 2, 4, 6, 8, 10, 11, 12, 14];
ind=[0, 4, 0, 1, 1, 2, 0, 3, 0, 4, 2, 3, 0, 4];
val=[3, 5.6, 1, 2, 1.1, 1, -2, 2.8, -1, 1, 1, -1.2, -1, 1.9];

ctyp = [ 'C', 'C', 'B', 'B', 'C', 'C', 'C', 'C'];
A=spfrommtlb(beg,ind,val,[m,n]);
Ae=A(3,:);be=b(3);
A=[-A(1,:);A(2,:);-A(4,:);A(5,:)];
b=[-b(1);b(2);-b(4);b(5)];
optimalValue = 3.23684210526;

// no 'B' in linprog 
// we use I with extra bounds 

ctyp1 = [ 'C', 'C', 'I', 'I', 'C', 'C', 'C', 'C'];
ub1=ub;ub1(3:4)=1;
[xopt,fopt,flag] = linprog(c,A,b,Ae,be,var_type=ctyp1,ub=ub1,lb=lb,sense=sense);
if abs(fopt - optimalValue) > 1.e-7 then pause;end

[xopt1,fopt1,flag1,extra1] = linprog_cplex(c,A,b,Ae,be,ub=ub,lb=lb, ...
					   var_type=ctyp, sense=sense);
// here xopt and xopt1 are optimal but not the same
// if norm(xopt-xopt1) >= 1.e-8 then pause;end 
if abs(fopt1 - optimalValue) > 1.e-7 then pause;end

// GamsSos1a
// Sos variables
//----------------------------------

// GamsSos2a
// Sos variables XXXX
// ----------------------------------

// char* probname = "SemiCont";
// --------------------------------

