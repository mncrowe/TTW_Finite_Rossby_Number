function [M, x] = grid_spectral(type,n,params,flip,c)
% Creates the spectral collocation differentiation matrix, M, and grid, x, for a chosen set of polynomials, p_n(x), and weight function, w(x)
% - type:   1 - Chebyshev   (first kind, roots)
%           2 - Chebyshev   (second kind, maxima plus end-points)
%           3 - Legendre
%           4 - Laguerre    (weight: w(x) = exp(-c*x/2))
%           5 - Hermite     (weight: w(x) = exp(-c^2*x^2/2))
%           0 - Custom      (enter x, w(x), w'(x) using 'params')
% - n: number of gridpoints, default: 16
% - params: grid parameters, depends on type, for some type value may be vector or scalar
%           type: 1,2,3 - L, for coordinate interval [0 L], may have L < 0
%                       - [a b], coordinate interval, may have b < a, default: [-1 1]
%           type: 4     - L, for coordinate interval [0 L], if L = 0, uses unscaled Laguerre points, default: 0
%                       - [a b], coordinate interval, if a = b, uses unscaled Laguerre points
%           type: 5     - L, for coordinate interval [-L L], if L = 0, uses unscaled Hermite points, default: 0
%                       - [a b], coordinate interval, if a = b, uses unscaled Hermite points
%           type: 0     - A, an n x 3 matrix given by A = [x w(x) w'(x)] with vectors x, w, w' as columns
% - flip:   1 - flips the basis and differentiation matrix (i.e. order of points reversed)
%           2 - flips the basis and differentiation matrix if x(end) < x(1)
%           0 - does not flip (default)
% - c: enter stretching parameter c directly, Laguerre and Hermite only (optional)
%
% ----------------------------------------------------------------------------
% Note: The method relies on the expansion f(x) ~ Sum_n [f(x_n) phi_n(x) w(x)]
%       therefore the weight function w(x) is the square root of the weight
%       function in the orthogonality relation for the given polynomials
%
% Note: The parameter 'c' in w(x) for Laguerre and Hermite polynomials is
%       given by c = x_n/L for largest root x_n. Therefore we can either
%       choose L for a given n such that the decay matches the expected
%       decay of the solution, or, we can choose L such that the function
%       can be taken to vanish for |x| > L. Entering L = 0 sets c = 1. The
%       weight functions given above for Laguerre and Hermite polynomials
%       correspond to the intervals [0 L] and [-L L] respectively; if using
%       other domains use [a b] input instead.
%
% Note: For Laguerre and Hermite we may set c directly using the 'c'
%       argument. If set directly, this will override any limits L or [a b],
%       however, the coordinates will still be centred on (a+b)/2;
%
% Note: Higher order derivative matrices can be calculated similarly
%       however using d^n/dx^n = M^n is generally correct to a high degree
%       of accuracy.
% ----------------------------------------------------------------------------

if nargin < 2; n = 16; end
if nargin < 4; flip = 0; end
if nargin < 5; c = NaN; end

% define grid params, w = w(x) is the weight function and wp = w'(x)

if type == 1    % Chebyshev, points are roots (points of first kind)
    if nargin < 3; params = [-1 1]; end
    if numel(params) == 1; params = [0 params]; end
    x = cos((1/2:(n-1/2))*pi/n)';
    a = params(1); b = params(2);
    x = (a+b)/2-(b-a)/2*x;
    w = ones(n,1);
    wp = zeros(n,1);
end

if type == 2    % Chebyshev, points are maxima (points of second kind), plus end-points
    if nargin < 3; params = [-1 1]; end
    if numel(params) == 1; params = [0 params]; end
    x = cos((0:(n-1))*pi/(n-1))';
    a = params(1); b = params(2);
    x = (a+b)/2-(b-a)/2*x;
    w = ones(n,1);
    wp = zeros(n,1);
end

if type == 3    % Legendre, points are roots
    if nargin < 3; params = [-1 1]; end
    if numel(params) == 1; params = [0 params]; end
    J = diag(0.5./sqrt(1-(2*(1:n-1)).^-2),1)+diag(0.5./sqrt(1-(2*(1:n-1)).^-2),-1);
    x = sort(eig(J));   % get Legendre points as eigenvalues of a Jacobi matrix
    a = params(1); b = params(2);
    x = (a+b)/2+(b-a)/2*x;
    w = ones(n,1);
    wp = zeros(n,1);
end

if type == 4    % Laguerre, points are roots, plus x = 0
    if nargin < 3; params = 0; end
    shf = 0;
    if numel(params) == 2; shf = params(1); params = params(2)-params(1); end
    J = diag(1:2:2*n-3)-diag(1:n-2,1)-diag(1:n-2,-1);
    x = sort(eig(J));   % get Laguerre points as eigenvalues of a Jacobi matrix
    if params == 0; lambda = 1; else; lambda = x(end)/params; end
    if ~isnan(c); lambda = c; end
    x = [0; x/lambda];
    w = exp(-lambda*x/2);
    wp = -lambda/2*exp(-lambda*x/2);
    x = x + shf;
end

if type == 5    % Hermite, points are roots
    if nargin < 3; params = 0; end
    shf = 0;
    if numel(params) == 2; shf = (params(1)+params(2))/2; params = (params(2)-params(1))/2; end
    J = diag(sqrt((1:n-1)/2),1)+diag(sqrt((1:n-1)/2),-1);
    x = sort(eig(J));   % get Hermite points as eigenvalues of a Jacobi matrix
    if params == 0; lambda = 1; else; lambda = x(end)/params; end
    if ~isnan(c); lambda = c; end
    x = x/lambda;
    w = exp(-lambda^2*x.^2/2);
    wp = -lambda^2*x.*exp(-lambda^2*x.^2/2);
    x = x + shf;
end

if type == 0    % Custom
    if nargin < 3; error('[x, w(x) ,w''(x)] required'); end
    x = params(:,1);
    w = params(:,2);
    wp = params(:,3);
end

% create matrix

D = x'-x;               % matrix of x_n - x_m
D1 = D; D1(~D) = 1;     % version of D with 0 -> 1 for performing products
a_n = ones(n,1)*(w'.*prod(D1));

M = a_n'./(a_n.*D');
M(1:(n+1):n^2) = sum(1./D1)-1+(wp./w)'; % extra -1 removes the additional 1 gained from using D1

if flip==1 || (flip==2 && x(end)<x(1)); x = x(end:-1:1); M = M(end:-1:1,end:-1:1); end

end