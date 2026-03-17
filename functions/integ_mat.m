function f = integ_mat(f,x,dim,method)
% Integrates a field f along a given dimension, for input f(x) output is int_x1^x f(z) dz
% - f: field
% - x: domain, x in [x1, x2]
% - dim: dimension to integrate through, (default: 1)
% - method: 1 - trapezoidal method (default)
%
% ----------------------------------------------------------------------------
% Note: Additional methods may or may not be added later.
%
% Note: This script calculates the integral from x1 to x at each x, contrast
%       with integ.m which calculates integral from x1 to x2.
% ----------------------------------------------------------------------------


s = size(f);

if nargin < 2; x = 1:s(1); end
if nargin < 3; dim = 1; end
if nargin < 4; method = 1; end

s2=[prod(s(1:dim-1)) s(dim) prod(s(dim+1:end))];

f=reshape(f,s2);
g=zeros(s2);

if method==1
    for i1 = 2:s(dim)
        g(:,i1,:) = g(:,i1-1,:)+(x(i1)-x(i1-1))*(f(:,i1,:)+f(:,i1-1,:))/2;
    end
end

if method==2
    % include additional methods
end

f=reshape(g,s);

end

