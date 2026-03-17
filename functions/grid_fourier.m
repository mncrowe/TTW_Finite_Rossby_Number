function [M,x] = grid_fourier(type,n,interval,flip)
% Creates the spectral differentiation matrix, M, and grid, x, for a discrete Fourier transform on a periodic domain
% - type:   1 - periodic    (periodic sinc) (default)
%           2 - unbounded   (sinc)
% - n: number of gridpoints, default: 16, (should be even for type = 1)
% - interval: [a b], coordinate interval, may have b < a, default: [-pi pi]
% - flip:   1 - flips the basis and differentiation matrix (i.e. order of points reversed)
%           2 - flips the basis and differentiation matrix if x(end) < x(1)
%           0 - does not flip (default)
%
% ----------------------------------------------------------------------------
% Note: If the grid is assumed to be periodic (type 1), it will only include
%       the first end-point, e.g. interval of [-1 1] with n = 4 gives the grid
%       x = [-1, -0.5, 0, 0.5].
%
% Note: The value of n should be taken as even with using type 1. Results are
%       not strictly valid for odd n and may give Gibbs-like oscillations.
%
% Note: If using type 2, the function is assumed to be 0 outside the interval.
%       Generally, type 1 can be used wherever type 2 is used if the function
%       could be taken to be periodic and remain sufficiently smooth. If this
%       doesn't hold then type 2 will often not work either.
%
% Note: Higher order derivative matrices can be calculated similarly
%       however using d^n/dx^n = M^n is generally correct to a high degree
%       of accuracy.
%
% Note: Differentiation by this method relies on O(n^2) matrix multiplication
%       so is slower than O(n log n) FFT methods. Use FFT if speed is required.
% ----------------------------------------------------------------------------

if nargin < 1; type = 1; end
if nargin < 2; n = 16; end
if nargin < 3; interval = [-pi pi]; end
if nargin < 4; flip = 0; end

if type == 1
    
    if rem(n,2) == 1; fprintf(2,'Warning: even values of n are recommended\n'); end
    x = linspace(interval(1),interval(2),n+1)'; x = x(1:n);
    h = x(2)-x(1); L = interval(2)-interval(1);
    M = pi/L*(-1).^((x-x')/h)./tan(pi*(x-x')/L);
    
end

if type == 2
    
    x = linspace(interval(1),interval(2),n)';
    h = x(2)-x(1);
    M = (-1).^((x-x')/h)./(x-x');
    
end

M(eye(n)==1) = 0;
M = real(M);

if flip==1 || (flip==2 && x(end)<x(1)); x = x(end:-1:1); M = M(end:-1:1,end:-1:1); end

end