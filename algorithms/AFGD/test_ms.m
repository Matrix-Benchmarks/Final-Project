clc
clear 
close all
% %load underlying matrix X_star
 d1=1000;
 d2=1000;
 r=5;
 c=10;
n=c*r*d1;
p=c*r*d1/(d1*d2);

 load(['./Data/data_d1_' num2str(d1) '_d2_' num2str(d2) '_r_' num2str(r) '.mat']);
% Calculate the incoherence and spectral norm of X_star
[mu,sigma1,kappa] = calc_para(X_star,r);
%get linear operator A
%A=generateA(d1,d2,n); 

fraction_shown = 0.2;
num_observations = round(fraction_shown*d1*d2);  
% mask = zeros(d1*d2,1);
% mask(1:num_observations) = ones(num_observations,1);
% [~,ind] = sort(randn(matrix_size*matrix_size,1));
% mask = mask(ind);
% phi = reshape(mask,d1,d2);
% Omega = find(phi);
Omega = randsample(d1*d2,num_observations);
[I, J] = ind2sub([d1 d2],Omega);

%phi = masking operator (Omega is masking indices in vector form)
%y = non-zero entries of X_star in vector form


idx=1:1:n;
A_id = @(z) z(Omega);
At_id = @(z) sparse(I,J,z,d1,d2,num_observations);
y = A_id(X_star);
f_grad = @(x) (-1)*At_id(y - A_id(x));

f= @(u) (1/2) * norm(y - A_id(u*u'), 2)^2;
%%tuning parameter
T=10000;
TS = 70;
rd=1;
kk = 0.2;
nn=length(kk);
tau=2;%interation number for initialization
err=10^(-7);
alpha = 0.2;

%initial1 nalbaf(0)
%init=1;
[U0,~] = initialization_ms(y,r,d1,d2,Omega);


for k = kk
eta = k/sigma1;
gamma = alpha^2/eta;
[X_hat1,dist1,loss1,time1] =MS_FGD(U0,X_star,f_grad,eta,gamma,T,TS,20);

%[X_hat2,dist2,loss2,time2] =MS_AGD(U0,X_star,f_grad,eta,gamma,T,TS,20);

%[X_hat3,dist3,loss3,time3] =MS_AFGD(U0,X_star,f_grad,eta,gamma,T,TS,20);

end





semilogy(time1,dist1,'-.*','LineWidth',2,'MarkerSize',4,'color',[0,0,1]);
%hold on
%semilogy(time2,dist2,'-.^','LineWidth',2,'MarkerSize',4,'color',[0,0.5,0]);
%hold on
%semilogy(time3,dist3,'-.o','LineWidth',2,'MarkerSize',4,'color',[1,0,0]);
xlabel('Time (s)','FontSize', 18);
ylabel('Frobenius norm error','FontSize', 18);

legend({'FGD', 'AGD', 'AFGD'},'interpreter','latex','FontSize', 20,'Location','northeast')

