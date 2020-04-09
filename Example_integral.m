%% Experiment
X=-2.5:0.001:2.5;                           %evaluation pts
f=@(x) sqrt(3/2)*x;                         %measure wrt f(x)
a={@(x) x, @(x,y) exp(-(x.^2+y.^2))};       %integral op. coeffs
smooth_meas1=intMeas(a,f,X,0.1,100,'order',1);      %epsilon=0.1
smooth_meas2=intMeas(a,f,X,0.01,300,'order',1);     %epsilon=0.01
smooth_meas3=intMeas(a,f,X,0.001,2250,'order',1,'PoleType','equi');   %epsilon=0.001

%% Plot
figure
semilogy(X,smooth_meas3,'LineWidth',2)
hold on
semilogy(X,smooth_meas2,'LineWidth',2)
semilogy(X,smooth_meas1,'LineWidth',2)
xlim([X(1) X(end)])
ylim([6e-5 1e2])
ax = gca; ax.FontSize = 14;