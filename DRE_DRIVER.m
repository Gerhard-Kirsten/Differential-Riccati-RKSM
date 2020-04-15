%% A typical example for running the RKSM-DRE code (Algorithm 4.1) %%
clear all
close all
clc

%We provide this code without any guarantee that the code is bulletproof as input 
%parameters are not checked for correctness.

%We ask any user who receives the code at this point to please not
%circulate the code as it is still under active development.

%Finally, if the user encounters any problems with the code, either of the authors can be
%contacted via e-mail.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% LIST OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Coefficient matrices

% coeff.A            coeff. matrix.  n x n (A < 0)
% coeff.C            rhs factor   p x n
% coeff.B            second order term   n x s
% coeff.Z            The low rank factor of the init. cond. n x q
% coeff.E            The mass matrix (E > 0)
% coeff.EL           Lower triangular (Cholesky or LU) factor of the mass matrix E
% coeff.EU           Upper triangular (Cholesky or LU) factor of the mass matrix E


%%% Algebraic reduction Paramters

% params.m           max number of iterations allowed
% params.tol         Algebraic stopping tolerance 

% params.smin,         estimates for real spectral interval
% params.smax        associated with field of values of A -- e.g., smax=norm(A,1); smin=smax/condest(A); 

% params.ch          ch=1  complex poles  ch=0 real poles
% params.period      how often check convergence (period=1 means at each iteration)
% params.Hritz       Hritz=1 uses adaptive spectral region, Hritz=0 uses spectral 
%                    region only based on projected A
% params.iterative   iterative = 1 solves the linear system at each iteration
%                    by block pcg if A symm. Otherwise, backslash is used when
%                    A is nonsym or iterative = 0 


%%% Numerical integration Parameters (Differential Part)

% diff.tim             The final time Tf
% diff.stpsiz          The stepsize h used in the reduction procedure
% diff.refinestep      The stepsize that should be employed in the refinement
%                      phase (If the user would prefer to use a higher order 
%                      (different) integrator in the refinement phase they could do so)
% diff.alltimes        alltimes = 0 only the final solution at Tf is given in factored
%                      form. alltimes = 1, the small dimensional solution Y(t)
%                      for all refined t is given by Yset = [Y1(:); 
%                                                            Y2(:);
%                                                              ...]


%% %%%%%%%%%%%%%%%%%% DEFINE THE PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%

nh=  400;
n = nh^2; % This will be the size of A

% For a symmetric matrix A (It is required that A<0)

Inh=speye(nh);   
T=toeplitz([-2,1,sparse(1,nh-2)]); 
A= sparse(kron(T,Inh)+kron(Inh,T));      %p=symamd(A); A=A(p,p); %IF NEEDED

coeff.A = A;

%For a nonsymmetric matrix A (It is required that A<0)
% 
% Inh=speye(nh); 
% T=toeplitz([-2,1,sparse(1,nh-2)]); 
% T1=toeplitz([-3,1,sparse(1,nh-2)],[-3,-1,sparse(1,nh-2)]); 
% A= sparse(kron(T,Inh)+kron(Inh,T1));    %  p=symamd(A); A=A(p,p); %IF NEEDED
% 
% coeff.A = A;


%%% Randomized low rank coefficient factors

rng(7)
coeff.B = randn(n,1);
rng(2)
coeff.C = randn(6,n);
rng(3)
coeff.Z = randn(n,2);

%%% Mass matrix E. IF The problem has identity mass matrix, just set E = EL = EU = speye(n);
%%% (It is required that E>0)

coeff.E = speye(n);
LE=chol(coeff.E,'lower');
coeff.EL = LE;
coeff.EU = LE';


%%%% Algebraic PARAMETERS %%%%%%%%
   


params.m = 40;
params.tol=1e-7;


params.smax= 1e5;
params.smin= 1e0; % For this specific example, otherwise:
%The value can be selected a-priori if a rough estimate of the spectral
%interval is available, or it can be approximated via an eigenvalue solver
%with loose accuracy, e.g.,
%opts.tol=1e-1;
%params.smax = eigs(-coeff.A,coeff.E,1,'largestreal',opts);
%params.smin = eigs(-coeff.A,coeff.E,1,'smallestreal',opts) ;

params.ch=0;
params.period=1;
params.Hritz=1;
params.iterative=1; % Only available for symmetric A


%%%% DIFFERENTIAL PARAMETERS %%%%%%%%

diff.tim=1; 
diff.stpsiz=0.1; % Final time and reduction timestep 
diff.refinestep = 0.01; 
diff.alltimes = 1;





%%%% RUN THE CODE %%%%%%%

[ZZ,nrmrestot,VV,K,s]=RKSM_DRE_FINAL(coeff,params,diff); 
