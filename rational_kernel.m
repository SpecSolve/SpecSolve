function [poles,res] = rational_kernel(m,type)
%%rational_kernel 
% Inputs: 
%         - order of kernel 'm'
%         - Pole location 'roots','cheb'
% Outputs: - poles of rational function 'poles' 
%         - residues of rational function 'res'

%Poles
if type=="cheb"
    z=1i+chebpts(m,1);  %Chebyshev points
    %z=1i+(2*(1:m)/(m+1)-1);%.^2.*sign((2*(1:m)/(m+1)-1));
elseif type=="roots"
    z=exp(1i*pi*((1:m)'-1/2)/m); %Roots of unity
elseif type=="extrap"
    z=1i*(0.5.^(0:m-1));
elseif type=="equi"
    z=1i+(2*(1:m)/(m+1)-1);
end

%Vandermonde matrix for poles
V=zeros(m);
for i=1:m 
    V(:,i)=z.^(i-1);
end

%Get residues directly
rhs=eye(m,1);
res=transpose(V)\rhs;
poles=z;

end

