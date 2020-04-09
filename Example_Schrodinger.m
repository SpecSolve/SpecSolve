% This example computes spectral measure of Schrodinger operator on real line.

clear
% close all
a=cell(3,1); % these specify the coefficients
a{1}=@(x) x.^2./(1+x.^6); % 0th order
a{2}=@(x) 0*x; % 1st order
a{3}=@(x) -1+0*x; % 2nd order

f=@(x) x.^2./(1+x.^6)*sqrt(9/pi); % function we compute measure wrt
X=0:0.05:6; % points where we compute smoothed measure
epsilon=0.1; % smoothing parameter
N=10^4; % truncation parameter

mu1=diffMeas(a,f,X,epsilon,N,'order',1,'parallel','on');
mu2=diffMeas(a,f,X,epsilon,N,'order',2,'parallel','on');
mu6=diffMeas(a,f,X,epsilon,N,'order',6,'parallel','on');

%% Plot
figure
semilogy(X,mu1,'linewidth',2)
hold on
semilogy(X,mu2,'linewidth',2)
semilogy(X,mu6,'linewidth',2)
axis([0,6,10^(-6),1])
ax = gca; ax.FontSize = 14;
