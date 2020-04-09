function [mu] = rseMeas(V,f,X,epsilon,N,varargin)
% This code computes smoothed spectral measures of a radial Schrodinger 
% operator, with form Lu(r)=-D2u(r) + V(r), with respect
% to a function f(r).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% V: Cell array with function handle for regular, singular and angular
%    parts of potential V(r), i.e.
%    V = {@(r) V1, @(r) V2, l}
%    so that V(r) = V1(r) + V2(r)/r + l(l+1)u(r)/r^2
%    Note: l is a nonnegative integer.
% f: function handle for function we compute measure wrt
% X: evaluation points
% order: order of kernel used
% epsilon: distance to real axis from where we evaluate the resolvent
% N: size of N x N discretization from Chebyshev collocation scheme
% L: scale for mapping to [-1,1]

% OPTIONAL LABELLED INPUTS
% Scale: scale for mapping to [-1,1], default is 10
% PoleType: where to place poles, default is equispaced
% Parallel: parfor (on) or normal for (off) loop, default is "off"
% Order: order of kernel used, default is 2
%
% OUTPUTS
% mu: smoothed measure wrt f at points X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect the optional inputs
p = inputParser;
addRequired(p,'V',@iscell);
addRequired(p,'f',@(x) isa(x,'function_handle'));
addRequired(p,'X',@isnumeric);
addRequired(p,'epsilon',@isnumeric);
addRequired(p,'N',@isnumeric);

validPole = {'cheb','roots','extrap','equi'};
checkPole = @(x) any(validatestring(x,validPole));
validPar = {'on','off'};
checkPar = @(x) any(validatestring(x,validPar));

addParameter(p,'Scale',10,@(x) x>0)
addParameter(p,'PoleType','equi',checkPole)
addParameter(p,'Parallel','off',checkPar)
addParameter(p,'Order',2,@(x) x==floor(x))

p.CaseSensitive = false;
parse(p,V,f,X,epsilon,N,varargin{:})

%%Unpack potential
%V1=V{1}; V2=V{2}; l=V{3};

%% Construct discretization
%Derivative and conversion matrices
D1=ultraS.diffmat(N,1); D2=ultraS.diffmat(N,2);
S12=ultraS.convertmat(N,1,1); S02=ultraS.convertmat(N,0,1);

%Variable coeffs and mult. matrices for mapped operator
x=chebfun('x'); V1_x=chebfun(@(x) V{1}(p.Results.Scale*(1+x)./(1-x))); V2_x=chebfun(@(x) V{2}(p.Results.Scale*(1+x)./(1-x)));
a_2=-((1-x^2)^2)*(1-x)^2; a_1=2*((1-x^2)^2)*(1-x); a_0=4*p.Results.Scale*(1+x)^2*V1_x + 4*p.Results.Scale*(1+x)*(1-x)*V2_x + 4*V{3}*(V{3}+1)*(1-x)^2;
M_2=ultraS.multmat(N,a_2,2); M_1=ultraS.multmat(N,a_1,1); M_0=ultraS.multmat(N,a_0,0);

%Radial Schrodinger operator after mapping to [-1,1]
H=M_2*D2+S12*M_1*D1+S02*M_0;

%Self-adjoint boundary condition u(0)=0 when origin is a regular point
if V{3}==0 && sum(abs(V2_x.coeffs))==0
    H=[[1 zeros(1,N-1)]; H(1:N-1,:)]; S02=[zeros(1,N); S02(1:N-1,:)];
    a=[1; zeros(N-1,1)];
    b=[0 (-1).^(1:N-1)];      %exploit (H+ab) for fast almost-banded solves
end

%Variable coeff for shifts
s=4*p.Results.Scale^2*(1+x)^2; S=ultraS.multmat(N,s,0); S=S02*S;

%Right hand side
[ccn,ccw]=chebpts(N,1);
f_x=chebfun(@(x) f(p.Results.Scale*(1+x)./(1-x))); F_vals=f_x(ccn);
rhs=chebtech1.vals2coeffs(F_vals); rhs=S*rhs;

%Inner product in values space to mitigate impact of singular weight at endpoint
f_vals=2*p.Results.Scale*ccw'.*f_x(ccn)./(1-ccn).^2; f_vals=f_vals';

%% Compute smoothed measure at evaluation points
[pts,alpha]=rational_kernel(p.Results.Order,p.Results.PoleType);

mu=zeros(length(X),1);
pf = parfor_progress(length(X));
pfcleanup = onCleanup(@() delete(pf));
warning('off','all')
if p.Results.Parallel=="off"
    if V{3}==0 && sum(abs(V2_x.coeffs))==0 %Woodbury for fast almost-banded solves
        for n=1:length(X)
            u_coeffs=zeros(size(rhs));
            for m=1:p.Results.Order
                z=X(n)+pts(m)*epsilon;
                u_coeffs1=(H-z*S)\rhs; c=(H-z*S)\a;
                u_coeffs=u_coeffs+alpha(m)*(u_coeffs1-(b*u_coeffs1)*c/(1+b*c));
            end
            u_vals=chebtech1.coeffs2vals(u_coeffs);     %values for inner prod.
            mu(n)=imag(f_vals*u_vals)/pi;              %inner product
            parfor_progress(pf);
        end
    else        %Sparse LU for banded solves.
        for n=1:length(X)
            u_coeffs=zeros(size(rhs));
            for m=1:p.Results.Order
                z=X(n)+pts(m)*epsilon;
                u_coeffs=u_coeffs+alpha(m)*((H-z*S)\rhs);
            end
            u_vals=chebtech1.coeffs2vals(u_coeffs);     %values for inner product
            mu(n)=imag(f_vals*u_vals)/pi;              %inner product
            parfor_progress(pf);
        end
    end
else
    if V{3}==0 && sum(abs(V2_x.coeffs))==0 %Woodbury for fast almost-banded solves
        parfor n=1:length(X)
            u_coeffs=zeros(size(rhs));
            for m=1:p.Results.Order
                u_coeffs1=(H-(X(n)+pts(m)*epsilon)*S)\rhs; c=(H-(X(n)+pts(m)*epsilon)*S)\a;
                u_coeffs=u_coeffs+alpha(m)*(u_coeffs1-(b*u_coeffs1)*c/(1+b*c));
            end
            u_vals=chebtech1.coeffs2vals(u_coeffs);     %values for inner prod.
            mu(n)=imag(f_vals*u_vals)/pi;              %inner product
            parfor_progress(pf);
        end
    else        %Sparse LU for banded solves.
        parfor n=1:length(X)
            u_coeffs=zeros(size(rhs));
            for m=1:p.Results.Order
                u_coeffs=u_coeffs+alpha(m)*((H-(X(n)+pts(m)*epsilon)*S)\rhs);
            end
            u_vals=chebtech1.coeffs2vals(u_coeffs);     %values for inner product
            mu(n)=imag(f_vals*u_vals)/pi;              %inner product
            parfor_progress(pf);
        end
    end
end
warning('on','all')
end