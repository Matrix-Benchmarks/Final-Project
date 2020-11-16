function V = accproj(Vp,A0,B0,D0inv,r,eta,beta,maxiter)
%ACCPROJ subroutine used in AFGD
%Inputs:
% Vp         = V prior to projection onto U0 space - { X : X'U0 is PSD}
% [A0,D0,B0] = r-SVD of U0
% r          = target rank of V
% maxiter    = # of iterations to run algorithm internal process for
%Outputs:
%V           = projection of Vp onto U0 space
% =========================================================================
T = A0'*Vp*B0;
Sig1 = 0; Sig2 = 0;

for t = 1:maxiter
    Sigp = Sig2 - eta*D0inv*(D0inv*Sig2 - T); %easy inverse cuz D0 is diag
    [A,D] = eig((Sigp + Sigp')/2);
    Sig1new = A * max(D,0) * A';
    
    Sig2 = Sig1 + beta*(Sig1new - Sig1);
    Sig1 = Sig1new;
end
V = (eye - A0*A0')*Vp + A0*D0inv*Sig1*B0';

end