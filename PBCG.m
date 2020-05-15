function [x,its_cg,res_cg]=PBCG(a,b,x,maxit,tol,L,U)
%function [x,its_cg,res_cg]=PBCG(a,b,x,maxit,tol,L,U)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We provide this code without any guarantee that the code is bulletproof as input 
%parameters are not checked for correctness.

%We ask any user who receives the code at this point to please not
%circulate the code as it is still under active development.

%Finally, if the user encounters any problems with the code, either of the authors can be
%contacted via e-mail.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[n,s]=size(b);
beta=0; 
r=b; ap=x;
gamma=r'*r;
res0=sqrt(sum(r.*r));
res=res0;
mem=ones(1,s);
k=0;
while (res./res0 > tol & k<maxit)
 
  z= U\(L\r);
  k=k+1;
  gamma=r'*z;
if (k==1), [p,ig]=qr(z,0); g=inv(ig);else, beta=g\(gamma0\gamma);[p,ig]=qr(z+p*beta,0);g=inv(ig);end,  ap=a*p;
  delta=p'*ap;
  alfa = delta\(g'*gamma);
  x = x + p*alfa;
  r = r - ap*alfa;
  gamma0=gamma;
  res=sqrt(sum(r.*r));
  mem=[mem;res./res0];

 if rem(k,6)==0

 end

 
end 

   res_cg=max(res./res0);
   its_cg=k;
 

