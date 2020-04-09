% This example computes spectral measure of a magnetic model of graphene.

clear
n=5000; % number of basis sites used
epsilon=0.05; % smoothing parameter
X=-3.1:0.05:3.1; % vector of points where we compute the smooth measure
PHI=0.25; % magnetic field

% construct the matrix from precomputation
load('Pre_computed_mats/graphene_lattice.mat')
N=(2*n1+1)*(2*n2+1)*2;

d1=spdiags(ones(2*n1,1),-1,2*n1+1,2*n1+1);
d2=spdiags(ones(2*n2,1),-1,2*n2+1,2*n2+1);
D2=spdiags(transpose(exp(2i*pi*PHI*(-n1:1:n1))),0,2*n1+1,2*n1+1);

T1=[sparse(round(N/2),round(N/2)), speye(round(N/2));
    speye(round(N/2)), sparse(round(N/2),round(N/2))];
T2=[sparse(round(N/2),round(N/2)), kron(speye(2*n2+1),d1');
    kron(speye(2*n2+1),d1),sparse(round(N/2),round(N/2)) ];
T3=[sparse(round(N/2),round(N/2)), kron(d2',D2);
    kron(d2,D2'),sparse(round(N/2),round(N/2)) ];

H0=T1+T2+T3;
H0=H0(ord,ord); % reorder for operator on l^2(N)
H0=H0(1:f(n),1:n); % take rectangular truncation
b=zeros(f(n),1);
b(1)=1;

% perform computation
mu = infmatMeas(H0,b,X,epsilon,'order',2);

%% Plot
plot(X,mu)
ax = gca; ax.FontSize = 14;
