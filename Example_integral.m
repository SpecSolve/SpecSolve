%% Experiment
X=-2.5:0.001:2.5;                           %evaluation pts
f=@(x) sqrt(3/2)*x;                         %measure wrt f(x)
a={@(x) x, @(x,y) exp(-(x.^2+y.^2))};       %integral op. coeffs

mu1=intMeas(a,f,X,0.1,'Order',1);           %epsilon = 0.1
mu2=intMeas(a,f,X,0.01,'Order',1);          %epsilon = 0.01
mu3=intMeas(a,f,X,0.001,'Order',1);         %epsilon = 0.001

%% Plot
figure
semilogy(X,mu3,'LineWidth',2)
hold on
semilogy(X,mu2,'LineWidth',2)
semilogy(X,mu1,'LineWidth',2)
xlim([X(1) X(end)])
ylim([6e-5 1e2])
ax = gca; ax.FontSize = 14;