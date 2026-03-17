function D = create_operator(n,varargin)
% Creates a discretised operator in R^n using the differential operators and boundary conditions for each dimension
% - n: space dimension
% - M, M_bc: enter n pairs of inputs;
%       M - N x N matrix operator corresponding to each dimension
%       M_bc - boundary conditions; 
%               - 2 x N matrix; here row 1 is the top BC and row 2 the bottom BC
%               - 2 x 1 vector; value of 'a' for Dirichlet condition a*f = b
%               - m; any scalar, do not apply BC, e.g. periodic condition

s = zeros(n,1);

for in = 1:n                        % define inputs as cell array
    
    M{in} = varargin{2*in-1};
    s(in) = length(M{in});
    
    M_bc{in} = varargin{2*in};
    
    if isequal(size(M_bc{in}),[2 1]); M_bc{in} = [M_bc{in}(1) zeros(1,s(in)-1); zeros(1,s(in)-1) M_bc{in}(2)]; end
    if numel(M_bc{in}) == 1; M_bc{in} = 0; end
    
end

D = M{1};                           % Use kronecker tensor product to create D
for in = 1:n-1
    D = kron(M{in+1},D);
end

r = prod(s);                        % number of rows in D

for in = 1:n                        % apply BCs
    
    if ~isequal(M_bc{in},0)
    
        b = prod([1 s(1:in-1)']);   % block width
        d = prod([1 s(1:in)']);     % distance between blocks

        rows1 = [];                 % top BC row
        rows2 = [];                 % bottom BC row

        for i2 = 1:(r/prod(s(1:in)))
            rows1 = [rows1 (1:b)+(i2-1)*d];
            rows2 = [rows2 (d-b+1:d)+(i2-1)*d];
        end

        ii = [];
        jj = [];
        kk = [];
        
        for r1 = rows1              % set top BC
            ii = [ii; r1*ones(s(in),1)];
            jj = [jj; (r1:b:r1+(s(in)-1)*b)'];
            kk = [kk; (M_bc{in}(1,:))'];
        end
        for r2 = rows2              % set bottom BC
            ii = [ii; r2*ones(s(in),1)];
            jj = [jj; (r2-(s(in)-1)*b:b:r2)'];
            kk = [kk; (M_bc{in}(2,:))'];
        end

        D2 = sparse(ii,jj,kk,r,r,max(nnz(D),length(ii)));   % create sparse matrix of boundary condition rows
        
        i1 = setxor(setxor(1:r,rows1),rows2);
        D2(i1,:) = D(i1,:);                 % set non BC rows to the rows of D
        
        D = D2;
    
    end
    
end
    
end