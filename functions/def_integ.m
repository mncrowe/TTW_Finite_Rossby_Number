function I = def_integ(M,i)
% Calculates the integration matrix corresponding to the differentiation matrix M
% - M: differentiation matrix
% - i: index of position to replace with integration constant, default: 1
%
% ----------------------------------------------------------------------------
% Note: Inverts the system f = M*g for integration matrix I with g = I*f.
%       Since differentiation maps constant vectors to 0 the matrix M must
%       have a zero eigenvalue (with eigenvector 1) and hence a zero
%       determinant. This represents the need to impose an integration
%       constant.
%
% Note: The constant is set by replacing the i^th row of the matrix I with
%       [0 0 .. 1 .. 0] where the 1 is in position i. This gives f_i = g_i
%       so the i^th position of f should be replaced with the value of g at
%       that point.
% ----------------------------------------------------------------------------

if nargin < 2; i = 1; end

M(i,:) = 0; M(i,i) = 1;
I = M^-1;

end

