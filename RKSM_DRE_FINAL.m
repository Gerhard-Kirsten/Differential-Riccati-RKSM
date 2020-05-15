function [Z,nrmrestot,VV,K,s,Yset]=RKSM_DRE_FINAL(coeff,params,diff)
%%

% Approximately Solve  
%               E'dX/dtE = A' XE  +  E'X A - E'XBB'XE + C'C
%
% by the Rational Krylov subspace method 
% (Order reduction onto the Rational Krylov subspace)
% This code implements Algorithm 4.1 in the corresponding Article.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% % Input: LIST OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%               [Z,nrmrestot,VV,K,s,Yset]=RKSM_DRE_FINAL(coeff,params,diff)
%                v.1.2

%%% Coefficient matrices

% coeff.A            coeff. matrix.  n x n
% coeff.C            rhs factor   p x n
% coeff.B            second order term   n x s
% coeff.Z            The low rank factor of the init. cond. n x q
% coeff.E            The mass matrix
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
% diff.alltimes             alltimes = 0 only the final solution at Tf is given in factored
%                      form. alltimes = 1, the small dimensional solution Y(t)
%                      for all refined t is given by Yset = [Y1(:); 
%                                                            Y2(:);
%                                                              ...]                                                    ...]


% This code is implemented such that the reduction and refinement procedures are performed
% by BDF(1,10) and BDF(1,refinestep) respectfully. The user can apply any other
% procedures.

%% %%%%%%%%%%%%%%%%%%%%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Z              factor of approximate solution at final timestep  X(Tf) = ZZ'
% nrmrestot      residual norm history
% VV             orthonormal basis for the approx space
% K              The coefficient matrix of the reduced system
% s              sequence of generated poles

% Yset           Contains the set of solutions Y(t)of the reduced problem
%                when alltimes = 1.

% Hints:
% 2) Provide "comfortable" (loose bounds) estimates s1, emax

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We provide this code without any guarantee that the code is bulletproof as input 
%parameters are not checked for correctness.
%Finally, if the user encounters any problems with the code, either of the authors can be
%contacted via e-mail.
% If this code is used, please cite:

% Kirsten, G., & Simoncini, V. (2019). Order reduction methods for solving 
% large-scale differential matrix Riccati equations. 
% arXiv preprint arXiv:1905.12119.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


A = coeff.A ;        
C = coeff.C ;         
B = coeff.B ;          
Z1 = coeff.Z ;           
E = coeff.E ;            
EL = coeff.EL ;         
EU = coeff.EU ;

   m=params.m;
   tol=params.tol;
   s1=params.smin;
   emax=params.smax;
   ch=params.ch;
   period=params.period;
   Hritz=params.Hritz;
   iterative  = params.iterative;
   
   refinestep = diff.refinestep;
   tim = diff.tim; 
   stpsiz = diff.stpsiz; 
   alltimes = diff.alltimes;
   
   dim1 = [];
if norm(Z1,1)<1e-12
    Lres= EL\C';
else
    Lres= [Z1,EL\C'];
end

timer  = [];
tic
ttime=cputime;
[n,n]=size(A);
C=full(C);
p=size(Lres,2);
I=speye(p);O=0*I;
In=speye(n);

if (norm(E-In,'fro')>1e-14)
  condestE=condest(E);
  singE=condestE/norm(E,'fro');
else
  singE=1;
end

bsize = size(B,2);
csize = size(C,1);
zsize = size(Z1,2);

[V1,rr1]=qr([EL\C',EL\B],0); 
nrmc=norm((rr1(:,1:csize)),'fro')^2; 
nrmb=norm((rr1(:,csize+1:end)),'fro')^2; 

[V,rr]=qr(Lres,0); 
nrmb2=norm((rr),'fro')^2; 
beta=V'*Lres; beta2=beta*beta';
errtot=[];

VV=V;

H=sparse(p*(m+2),p*(m+1));
nrmrestotnew=[];


if (norm(A-A',1)<1e-14), symm=1;
fprintf('The matrix A of dimension %d is symmetric \n', n)
else symm=0;
fprintf('The matrix A of dimension %d is nonsymmetric \n', n)
end

fprintf('\n')
 
 newEv = (EU\V);
 EV = newEv;
 newAv=EL\(A'*newEv);
 K=full(V'*newAv); 
 s=s1(1);
 eH=eig(K);
 eHpoints = sort([s1(:)',emax]);
 snew=newpolei(eHpoints,eH,s1(1)*ones(p,1));
 if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end
 s=[s,snew];

% additional steps
cmplxflag=0;
itsinner=0;


 if iterative & symm == 1
      
      fprintf('Applying PCG to symmetric linear solves \n')
      
 elseif symm == 1
     
     fprintf('Applying backslash to symmetric linear solves \n')
     
 else
     
     fprintf('Applying backslash to nonsymmetric linear solves \n')
     
 end
fprintf(' \n')
fprintf(' \n')
fprintf('   space dim   residual norm   PCG Its.    PCG Res. \n')
i=0;

tt = toc;
timer = [timer, tt];
itcheck = i;
BB = (EL\B);
CC = (EL\C');
ZZ = Z1;
while i < m

tic

  i=i+1;

  paired=0;
  itp=1;
  cflag = 0;
  Vwrk = V;
  while (paired==0),

    i1=i+1; it = 0; t=0.;
    w=Vwrk;
 
    
    %%%%%%% Use iterative solver %%%%%%%

    AA = (A'-snew*E);
  

%     
      if iterative & symm == 1
      
      ww = (EL*w);
      opts.type='ict'; opts.droptol=1e-4;
      
                  
    

      L=ichol(-real(AA),opts);


      [wrk1,its_cg,res_cg]=PBCG(AA,ww,zeros(n,p),25,1e-10,-L,-L');    
      wrk1 = EU*wrk1;
      
      else
   
    %%%%%% Use direct solver %%%%%%%%%
    
     its_cg=0; res_cg=0;
     
     ww = (EL*w);
     wrk1 = (A'-snew*E)\ww;
     wrk1 = EU*(wrk1);
     
      end
     
     
     %%%% All real basis implementation for RKSM %%%%%
     
     if imag(wrk1) ~= 0 & cflag == 0
         wrk = real(wrk1);
         cflag = 1;
     elseif imag(wrk1) ~= 0 & cflag == 1
         wrk = imag(wrk1);
     else
         wrk = wrk1;
     end
     
     
% Gram-Schmidt step
    jms=(i-1)*p+1;j1s=(i+1)*p;js=i*p;js1=js+1;
    Bm(jms:js,1:bsize)= V'*BB;
    Cm(jms:js,1:csize)= V'*CC;
    Zm(jms:js,1:zsize)= V'*ZZ;
 

    for it=1:2,
      for kk=1:i
        k1=(kk-1)*p+1; k2=kk*p; 
        ww=VV(1:n,k1:k2);
        gamma=ww'*wrk;
        H(k1:k2,jms:js) = H(k1:k2,jms:js)+ gamma;
        wrk = wrk - ww*gamma;
      end
    end
    
    [V,hinv]=qr(wrk,0); H(js1:j1s,jms:js)=hinv; hinv = inv(hinv);
    
    if (cmplxflag), snew=conj(snew); s=[s,snew];cmplxflag=0;
    
    newEv = (EU\V);
    newAv=EL\(A'*newEv); 
    D = kron(spdiags(s(2:end)),I);
    g = VV'*newAv;
    g1 = g; 
    EV = [EV];
    g22 = (EL\(A'*EV));
    g2 = V'*g22;
    g3 = V'*(newAv);
    K = [K g1; g2, g3];
    VV=[VV,V];
    EV = [EV, newEv];
    i=i+1; itp=itp+1;
    else, paired=1; end
  end

    VVnew = [VV, V];
    ih1=i1; ih=i;
    newEv = (EU\V);
    newAv=EL\(A'*newEv); 
    D = kron(spdiags(s(2:end)),I);
    g = VV'*newAv;
    
    
   if (symm), K=(K+K')/2; end

if (rem(i,period)==0)
    
    
    g22 = (EL\(A'*EV));
    KK = VVnew'*g22;
    EV = [EV, newEv];

    

% Solve the projected problem

%%%%%%%%%%%%%%%%%%%%%%%% Internal Solver %%%%%%%%%%%%%%%%%%%%%%%


    
   [Y,Ysum,Ysum2,Yset] = bdf1FR(K',Bm,Cm*Cm',Zm,tim,stpsiz,0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrmrestot = [];


% computed residual   (exact, in exact arithmetic) cheaper computation possible
     u1=newAv-VV*g;
     d=-VV*(Ysum*(H(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
     U=[-V*s(end),  d u1 ];
     rr=qr(full(U),0); rr=triu(rr(1:size(rr,2),:));
     

%backward error
nrmres = norm(rr*sparse([O I O; I O I; O I O ])*rr','fro');

% Integral norms for backward error

    nrmcint = tim*nrmc;
    nrmxint = norm(KK*Ysum,'fro');
    nrmxint2 =  norm(Ysum2,'fro');
    nrmresnew = (nrmres)/(nrmcint + 2*nrmxint + nrmxint2);
 
    nrmrestotnew = [nrmrestotnew, nrmresnew];

     dim = size(VV,2);
     dim1 = [dim1,dim];

     disp([real(size(VV,2)),real(nrmresnew),real(its_cg),real(res_cg)])
     if (nrmresnew<tol), 
    tt = toc;
    
    if i - itcheck == 1
    timer = [timer, timer(itcheck+1) + tt];
    else
    timer = [timer, timer(itcheck+1) + tt];   
    end
    itcheck  = i;
    
    
               %%%% REFINEMENT PROCEDURE %%%%
   
   fprintf(' \n')
   fprintf(' \n')
    
   fprintf('Tolerance Reached \n')
   fprintf('Algebraic reduction converged with a space dimension of %d \n',size(VV,2))
   fprintf('Algebraic Reduction took %d seconds \n',timer(end))
   
   fprintf(' \n')
   fprintf(' \n')
   
   tic  
   

       
   [Y,Ysum,Ysum2,Yset] = bdf1FR(K',Bm,Cm*Cm',Zm,tim,refinestep,alltimes);    
 
   
   refinetime = toc;

   %%% Compute refinement residual
   
   
   u1=newAv-VV*g;
   d=-VV*(Ysum*(H(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
   U=[-V*s(end),  d u1 ];
   rr=qr(full(U),0); rr=triu(rr(1:size(rr,2),:));
     
   fprintf('Refinement begins \n')
   fprintf('Refinement ended \n')
   nrmres = norm(rr*sparse([O I O; I O I; O I O ])*rr','fro');

   nrmcint = tim*nrmc;
   nrmxint = norm(KK*Ysum,'fro');
   nrmxint2 =  norm(Ysum2,'fro');
 
   nrmresnew = (nrmres)/(nrmcint + 2*nrmxint + nrmxint2);
   fprintf(' \n')
    
   fprintf('Refinement took %d seconds',refinetime)
   fprintf(' \n')
     break,
     end
end


% New poles and zeros
  if Hritz
    eH=sort(eig(K'-Y*Bm*Bm'));
    [ih,hh]=find(real(eH)>0);
    eH(ih)=-abs(real(eH(ih)))+i*imag(eH(ih));
  else
    eH=sort(eig(K));
  end
  eHorig=eH;
  eK=sort(eig(full(K)));

 if (ch)                     % Complex poles. Compute set for next complex pole of r_m

    if (any(imag(eH)) ~=0 & max(abs(imag(eH)))>1e-5 & length(eH)>2) % Roots lambdas come from convex hull too
     eH=[eH;-emax];
      ij=convhull(real(eH),imag(eH)); eH=eH(ij);
      ieH=length(eH); missing=ih*p-ieH;
      while missing>0,                         % include enough points from the border
        neweH=(eH(1:ieH-1)+eH(2:ieH))/2;missing=ih*p-length(eH);
        eH=[eH;neweH];
      end
      eHpoints=-eH;
      eH=eHorig;
    else                                  % if all real eigs, no convex hull possible
      eHpoints = sort([s1; emax.';-real(eH)]);
    end


 else   % Real poles s from real set. Compute complex roots of r_m via Ritz convex hull
     
     if (any(imag(eH)) ~=0 & length(eH)>2)    % Roots lambdas come from convex hull too
       eH=[eH;-s1;-emax.'];
       ij=convhull(real(eH),imag(eH)); eH=eH(ij);
       ieH=length(eH); missing=ih*p-ieH;
       while missing>0, % include enough points from the border
         neweH=(eH(1:ieH-1)+eH(2:ieH))/2;
         eH=[eH;neweH];
         missing=ih*p-length(eH);
       end
       eH=eH(1:ih*p);
     end
      eHpoints = sort([s1; emax.';-real(eH)]);
      eH=eHorig;
 end


 gs=kron(s(2:end),ones(1,p))';

 snew = newpolei(eHpoints,eH,gs);
 if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end


% If pole is complex, include its conjugate

 if (imag(snew) ~=0), cmplxflag=1;end
 s=[s,snew];

 g1 = g; 
 
 g2 = V'*g22;
 g3 = V'*newAv;

 K = [K g1; g2, g3];
 VV=[VV,V];
    
    tt = toc;
    if i - itcheck == 1
    timer = [timer, timer(itcheck+1) + tt];
    else
    timer = [timer, timer(itcheck+1) + tt];   
    end
    itcheck  = itcheck + 1;
    
end


% factored solution at Tf  (A similar implementation can be done for the other timesteps)


[uY,sY]=eig(Y); [sY,id]=sort(diag(sY));
sY=flipud(sY); uY=uY(:,id(end:-1:1));
is=sum(abs(sY)>1e-8);
Y0 = uY(:,1:is)*diag(sqrt(sY(1:is))); 
Z = EL'\VV(:,1:size(Y0,1))*Y0; 
final_rank=is;
RKStotal_time=cputime-ttime;
 fprintf(' \n')
fprintf('Subspace dimension and solution rank (at Tf): \n')
fprintf('Space dim %d  Solution rank %d \n',size(VV,2),is);
fprintf('Done\n')
sze = size(nrmrestot,2);
 

%%% Plot space dimension and times

figure(1) 
semilogy(timer(2:end),nrmrestotnew,'r->')
ylabel('Backward Error')
xlabel('Computational time')
hold on

figure(2)
semilogy(dim1,nrmrestotnew,'r->')
ylabel('Backward Error')
xlabel('Space dimension')
hold on
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfun(x,eH,s)
 
for j=1:length(x)
r(j)=abs(prod( (x(j)-s)./(x(j)-eH) ));
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfuni(x,eH,s)

for j=1:length(x)
r(j)=1./abs(prod( (x(j)-s)./(x(j)-eH) ));
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function snew=newpolei(eHpoints,eH,s)

for j=1:length(eHpoints)-1

    sval=linspace(eHpoints(j),eHpoints(j+1),20);

    [sf,jx] = max (abs(ratfun(sval,eH,s)));
    snew(j)=sval(jx);
end
[sn,jx]=max(abs(ratfun(snew,eH,s)));
snew=snew(jx);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RA=compute_RA(A,eH,s);

ne=length(eH);
s=[1e20,s];
I=speye(size(A));
RA=I;
for k=1:ne,
   RA = (A-eH(k)*I)/(A-s(k)*I)*RA;
end
return
