%% Experiment 1: Calculate P(0.5 < E < 2) = mu_f([0.5,2])
[ccn,ccw]=chebpts(20,[0.5 2]);                      %quadrature rule
indx=2:4; normf=sqrt(pi/8)*(2-igamma(1/2,2*indx.^2)/gamma(1/2));
f1 = @(r) exp(-(r-indx(1)).^2)/sqrt(normf(1));
f2 = @(r) exp(-(r-indx(2)).^2)/sqrt(normf(2));
f3 = @(r) exp(-(r-indx(3)).^2)/sqrt(normf(3));
V={@(r) 0, @(r) exp(-r)-1, 1};                      %Schrodinger potential
tic
mu1=rseMeas(V,f1,ccn,0.1,'Order',4);                %measure wrt f_1(x)
mu2=rseMeas(V,f2,ccn,0.1,'Order',4);                %measure wrt f_2(x)
mu3=rseMeas(V,f3,ccn,0.1,'Order',4);                %measure wrt f_3(x)
t=toc
P1=ccw*mu1, P2=ccw*mu2, P3=ccw*mu3                  %compute probabilities

%% Experiment 2: Plot densities w/shaded area
%(~3.2 min runtime on laptop intel i7 proc.)
X=0:0.05:5;                                       %evaluation pts
indx=2:4; normf=sqrt(pi/8)*(2-igamma(1/2,2*indx.^2)/gamma(1/2));
f1 = @(r) exp(-(r-indx(1)).^2)/sqrt(normf(1));
f2 = @(r) exp(-(r-indx(2)).^2)/sqrt(normf(2));
f3 = @(r) exp(-(r-indx(3)).^2)/sqrt(normf(3));
V={@(r) 0, @(r) exp(-r)-1, 1};                      %Schrodinger potential
tic;
mu1=rseMeas(V,f1,X,0.1,'Order',4);                  %measure wrt f_1(x)
t1=toc;
mu2=rseMeas(V,f2,X,0.1,'Order',4);                  %measure wrt f_2(x)
t2=toc;
mu3=rseMeas(V,f3,X,0.1,'Order',4);                  %measure wrt f_3(x)
t3=toc;
tot_min=(t1+t2+t3)/60

%plot
figure
fill([X(11:41) flip(X(11:41))],[1.75*ones(size(11:41)) zeros(size(11:41))],'k','LineStyle','none')
hold on
plot(X,mu1,'LineWidth',2)
plot(X,mu2,'LineWidth',2)
plot(X,mu3,'LineWidth',2)
xlim([0 5])
ylim([0 1.75])
clear alpha
alpha(0.2)
ax = gca; ax.FontSize = 14;