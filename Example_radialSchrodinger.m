%% Experiment (~2 min runtime on dual core i7 proc., e.g., laptop)
X=0:0.05:5;                                       %evaluation pts
indx=2:4; normf=sqrt(pi/8)*(2-igamma(1/2,2*indx.^2)/gamma(1/2));
f1 = @(r) exp(-(r-indx(1)).^2)/sqrt(normf(1));
f2 = @(r) exp(-(r-indx(2)).^2)/sqrt(normf(2));
f3 = @(r) exp(-(r-indx(3)).^2)/sqrt(normf(3));
V={@(r) 0, @(r) exp(-r)-1, 1};                      %Potential coeffs
tic;
smooth_meas1=rseMeas(V,f1,X,0.1,1e4,'Scale',20,'Order',4);          %measure wrt f_1(x)
t1=toc;
smooth_meas2=rseMeas(V,f2,X,0.1,1e4,'Scale',20,'Order',4);          %measure wrt f_2(x)  
t2=toc;
smooth_meas3=rseMeas(V,f3,X,0.1,1e4,'Scale',20,'Order',4);          %measure wrt f_3(x)
t3=toc;
tot_min=(t1+t2+t3)/60

%% Plot
figure
fill([X(11:41) flip(X(11:41))],[1.75*ones(size(11:41)) zeros(size(11:41))],'k','LineStyle','none')
hold on
plot(X,smooth_meas3,'LineWidth',2)
plot(X,smooth_meas2,'LineWidth',2)
plot(X,smooth_meas1,'LineWidth',2)
xlim([0 5])
ylim([0 1.75])
clear alpha
alpha(0.2)
ax = gca; ax.FontSize = 14;