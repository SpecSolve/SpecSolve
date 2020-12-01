% Experiment 1: spectral measure of Schrodinger operator on real line
X=0:0.05:6;                                     %Evaluation pts
f=@(x) x.^2./(1+x.^6)*sqrt(9/pi);               %Measure wrt f(r)
a={@(x) x.^2./(1+x.^6), @(x) 0*x, @(x) -1+0*x}; %Schrodinger op

tic
mu1=diffMeas(a,f,X,0.1,'order',1,'N',10^4);  %epsilon=0.1, m=1
t1=toc;
mu2=diffMeas(a,f,X,0.1,'order',2,'N',10^4); %epsilon=0.1, m=2
t2=toc;
mu6=diffMeas(a,f,X,0.1,'order',6,'N',10^4); %epsilon=0.1, m=6
t3=toc;
t_min=(t1+t2+t3)/60

%% Plot
figure
semilogy(X,mu1,'linewidth',2)
hold on
semilogy(X,mu2,'linewidth',2)
semilogy(X,mu6,'linewidth',2)
axis([0,6,10^(-6),1])
ax = gca; ax.FontSize = 14;
