function [Y,Ysum,Ysum2,Yset] = bdf1FR(Tm,Bm,Q,Zm1,t,h,alltimes)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    [Y,Ysum,Ysum2,Yset] = bdf1FR(Tm,Bm,Q,Zm1,t,h,alltimes)  v. 1.2

%We provide this code without any guarantee that the code is bulletproof as input 
%parameters are not checked for correctness.
%Finally, if the user encounters any problems with the code, either of the authors can be
%contacted via e-mail.
% If this code is used, please cite:

% Kirsten, G., & Simoncini, V. (2019). Order reduction methods for solving 
% large-scale differential matrix Riccati equations. 
% arXiv preprint arXiv:1905.12119.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solves the projected  DRE Ym_dot = Tm'Ym + YmTm - YmBmBm'Ym + Q, Ym(0) = Zm1Zm1', using the backward Euler scheme. 
% for the final time t  with stepsize h.

Xm1 = Zm1*Zm1';

% Integration coefficients

m = t/h;
beta = 1;
a0 = 1;


% Initialise integration - Solution at time t = t0


    [n1,~] = size(Xm1);
if alltimes == 1
    
    Ysum =  h*Xm1;
    Ysum2 =  h*Xm1*Bm*Bm'*Xm1;;
    Yset(:,:,1) = Xm1;
    
else
    
    Ysum =h*Xm1;
    Ysum2 =  h*Xm1*Bm*Bm'*Xm1;
end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now iterate m-1 times using BDF(1)

Anew = h*beta*Tm - 0.5*eye(n1);
Bnew = sqrt(h*beta)*Bm;

j=1;

while j < m

    Qnew = h*beta*Q + a0*Xm1;

    [Xm1] = care(full(Anew),full(Bnew),full(Qnew));
    
 
        %Update integral solution

    
if alltimes == 1
    
    Ysum = Ysum + h*Xm1;
    Ysum2 = Ysum2 + h*Xm1*Bm*Bm'*Xm1;
    Yset(:,:,j) = Xm1;
    
else
    
    Ysum = Ysum + h*Xm1;
    Ysum2 = Ysum2 + h*Xm1*Bm*Bm'*Xm1;
    Yset = 0;
    
end
    %Update iteration
    
  j=j+1;  
end
 
 Y = Xm1;   
end
 
    
    
    