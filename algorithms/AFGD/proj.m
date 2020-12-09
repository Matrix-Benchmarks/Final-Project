function [yy] = proj (V, S1,D1,S2, T)
%%%Y is the base, Y = S1 D1 S2'
[~,rr] = size(V);
sig = zeros(rr);
sig1 = S1'*V*S2;

%%%minimize |inv(D1)*sig - sig1|
for t=1:T

%%%calculate the gradient

sig = sig - min(diag(D1))^2* inv(D1)*(inv(D1)*sig - sig1)/3;

%%%%truncate

sig = (sig+sig')/2;
%disp(norm(sig^(-1)));


end

yy = V - S1*(S1'*V) + S1*(inv(D1)*sig*S2');





end