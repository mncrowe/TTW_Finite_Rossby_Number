function [U,V,W,B,P,x,z,M2,N2,Psi,err,dx,dz] = TTW_State(n,E,Ro,Pr,delta,Nx,Nz,inc_B0,method,dbdx,err_tol,err_max,out)
% Calculates the TTW steady state for a given set of parameters and grids.
%
% Inputs:
% n - number of iterations, enter 0 to get TTW solution (default: 10)
% E - Ekman number (default: 0.1)
% Ro - Rossby number (default: 1)
% Pr - Prandtl number (default: 1)
% delta - depth ratio (default: 0)
% Nx - number of horizontal gridpoints (default: 64 point, Fourier)
% Nz - number of vertial grispoints, should be odd (default: 31 point, type 2 Chebyshev)
% inc_B0 - include B_0 in B field if equal to 1 (default: 1)
% method - 1: seperate iteration (fast), 2: combined iteration (slow, default)
% dbdx - anonymous function for background gradient db_0/dx (default: @(x) exp(-(x*sqrt(pi)/2).^2))
% err_tol - method will terminate once err < err_tol (default: 0)
% err_max - method will terminate if err > err_max (default: inf)
% out - display iteration information if set to 1 (default: 1)
%
% Outputs:
% (U,V,W,B,P) - velocity, buoyancy abd pressure fields as (x,z) arrays
% (x,z) - coordinate grids as vectors
% M2 - horizontal buoyancy gradient, M^2 = dB/dx
% N2 - vertical buoyancy gradient, N^2 = dB/dz
% Psi - streamfunction for (U,W), U = dPsi/dz, W = -dPsi/dx
%
% -------------------------------------------------------------------------
% Notes:
%
% Script prints the iteration difference for (U,V) and B. This is the sum
% of the squared difference between the current and previous iteration
% values for the fields. Use this to tune the correct number of iterations.
% Failure to converge corresponds to non-decreasing difference and is
% usually due to frontal collapse.
%
% By default, this script uses:
%
%                   b(x) = erf[x*sqrt(pi)/2],
%
% where the factor of sqrt(pi)/2 ensures that db/dx has a maximum value of
% Ro in the centre of the front. Only dBdx is used in the calculations and 
% hence all quantities are periodic in the x direction and vanish far from 
% the frontal region. b(x) should be chosen s.t. b -> +/-1 as x -> +/-inf.
%
% Two methods are implemented. Method 1 calculates (U,V) initially, then
% uses mass conservation to find W. B is then found by inverting the
% advection diffusion operator before finally calculating P by vertically
% integrating the vertical momentum equation. This is fairly fast though
% the mismatch between calculating the horizontal and vertical velocities
% means it can be unstable for strongly non-hysrostatic cases (delta > 0).
% Method 2 solves the full semi-linear system simultaneously with the old
% solution for (U,W) used only in the advection operator. This is more
% stable and accurate but takes longer as the linear system is much larger.
% -------------------------------------------------------------------------

% set undefined parameters
if nargin < 1; n = 2; end
if nargin < 2; E = 0.1; end
if nargin < 3; Ro = 1; end
if nargin < 4; Pr = 1; end
if nargin < 5; delta = 0; end
if nargin < 6; Nx = 64; end    % 3.5 is sufficient to ensure dBdx is continuous across boundaries
if nargin < 7; Nz = 31; end   % use points of the second kind as they include boundary values
if nargin < 8; inc_B0 = 1; end
if nargin < 9; method = 2; end
if nargin < 10; dbdx = @(x) exp(-(x*sqrt(pi)/2).^2); end
if nargin < 11; err_tol = 0; end
if nargin < 12; err_max = inf; end
if nargin < 13; out = 1; end
if rem(Nz,2) == 0; error('Odd number of z points required'); end    % require odd Nz such that the DA(B) = 0 condition can be applied on z = 0

[dx,x] = grid_fourier(1,Nx,[-3.5,3.5]);
[dz,z] = grid_spectral(2,Nz,[-1/2 1/2]);

% define background horizontal buoyancy gradient
dBdx = dbdx(x);

% define vertical dependence of mixing
nu = @(z) (1+0*z)*E;        % viscosity
kappa = @(z) (1+0*z)*E/Pr;  % diffusivity

% define fields for vertically dependent Ekman numbers
E_v = ones(length(x),1)*nu(z');
E_k = ones(length(x),1)*kappa(z');

% define TTW solution to use as initial guess for iteration, this assumes nu, kappa are constant functions
K_def       
U = -sqrt(E)*dBdx*K_pp(z'/sqrt(E),1/sqrt(4*E));
V = -sqrt(E)*dBdx*K(z'/sqrt(E),1/sqrt(4*E));
W = E*(dx*dBdx)*K_p(z'/sqrt(E),1/sqrt(4*E));
B = -Ro*Pr*sqrt(E)*(dBdx.^2)*K(z'/sqrt(E),1/sqrt(4*E));
P = Ro*Pr*(dBdx.^2)*((z'.^2-1/12)/2+E*(K_ppp(z'/sqrt(E),1/sqrt(4*E))-2*sqrt(E)*K_pp(1/sqrt(4*E),1/sqrt(4*E))));
dBdx0 = dBdx*ones(1,length(z));
dPdx0 = dBdx*(z');

% reshape fields into (Nx*Nz)^2 diagonal operators
mat2vec = @(F) reshape(F,[Nx*Nz 1]);    % operator to turn Nx x Nz matrix into Nx*Nz x 1 vector
vec2mat = @(F) reshape(F,[Nx, Nz]);     % operator to turn Nx*Nz x 1 vector into Nx x Nz matrix

% define operators on x or z
Ix = eye(Nx); Iz = eye(Nz);                 % identities
Dz = def_integ(dz)*diag([0; ones(Nz-1,1)]); % depth integral, starting from 0
Az = ones(Nz,1)*Dz(Nz,:);                   % depth average, constant in z
Lz = eye(Nz)-Az;                            % depth-indep operator, i.e. Lz(f) = f - Az(f)

% combine into operators on (x,z) space
I = eye(Nx*Nz);                             % identity
dx_vec = create_operator(2,dx,0,Iz,0);      % x derivative
dz_vec = create_operator(2,Ix,0,dz,0);      % z derivative
Dz_vec = create_operator(2,Ix,0,Dz,0);      % depth interal
Az_vec = create_operator(2,Ix,0,Az,0);      % depth interal
Lz_vec = create_operator(2,Ix,0,Lz,0);      % subtract depth-independence

% turn all fields to vectors for iteration
U = mat2vec(U); V = mat2vec(V); W = mat2vec(W); B = mat2vec(B); P = mat2vec(P);
E_v = mat2vec(E_v); E_k = mat2vec(E_k); dBdx0 = mat2vec(dBdx0); dPdx0 = mat2vec(dPdx0);

i = 0;
err = err_tol+1e-6;
dA = mean(diff(x))*mean(diff(z)); % area element for domain averaged error
while i < n && err > err_tol && err < err_max
    i = i+1;
    if out == 1; disp(['Iteration ' num2str(i) ' of ' num2str(n) ':']); end
    if method == 1  % solve system in parts
        
        % build operator for (U,V), system is M*[U,V]' = R
        L = Ro*(dx_vec*Lz_vec*diag(U) + dz_vec*diag(W)) - dz_vec*diag(E_v)*dz_vec;  % define advection-diffusion operator
        M = [L -I; I L];                                                            % define linear system
        R = [-Lz_vec*(dPdx0+dx_vec*P); zeros(Nx*Nz,1)];                             % linear system RHS

        % apply boundary conditions
        M([1:Nx (Nz-1)*Nx+1:Nx*Nz],:) = [dz_vec(1:Nx,:) zeros(Nx,Nx*Nz); dz_vec((Nz-1)*Nx+1:Nx*Nz,:) zeros(Nx,Nx*Nz)];                    % U BC
        M([Nx*Nz+1:(Nz+1)*Nx (2*Nz-1)*Nx+1:2*Nx*Nz],:) = [zeros(Nx,Nx*Nz) dz_vec(1:Nx,:); zeros(Nx,Nx*Nz) dz_vec((Nz-1)*Nx+1:Nx*Nz,:)];   % V BC
        R([1:Nx (Nz-1)*Nx+1:Nx*Nz Nx*Nz+1:(Nz+1)*Nx (2*Nz-1)*Nx+1:2*Nx*Nz]) = 0;                                                          % BC values

        % Solve for (U,V), calculate differences and update (U,V)
        X = M\R; dU = sqrt(dA*sum((U-X(1:Nx*Nz)).^2)); dV = sqrt(dA*sum((V-X(Nx*Nz+1:2*Nx*Nz)).^2)); U = X(1:Nx*Nz); V = X(Nx*Nz+1:2*Nx*Nz);

        % solve for W, subject to W = 0 on bottom (and hopefully top)
        W = -Dz_vec*dx_vec*U;

        % build operator for B, enforce DA(B) = 0 on line describing z = 0 (we use z = 0 so system is symmetric)
        L = Ro*(dx_vec*Lz_vec*diag(U) + dz_vec*diag(W)) - dz_vec*diag(E_k)*dz_vec;  % define advection-diffusion operator
        R = -Ro*diag(U)*dBdx0;                                                         % linear system RHS
        L([1:Nx (Nz-1)*Nx+1:Nx*Nz],:) = dz_vec([1:Nx (Nz-1)*Nx+1:Nx*Nz],:);         % B BC
        R([1:Nx (Nz-1)*Nx+1:Nx*Nz]) = 0;                                            % BC values
        L((ceil(Nz/2)-1)*Nx+1:ceil(Nz/2)*Nx,:) = Az_vec((ceil(Nz/2)-1)*Nx+1:ceil(Nz/2)*Nx,:); % DA(B) BC
        R((ceil(Nz/2)-1)*Nx+1:ceil(Nz/2)*Nx) = 0;                                             % BC value

        % solve for B, calculate difference and update B
        Z = L\R; dB = sqrt(dA*sum((Z-B).^2)); B = Z;

        % solve for P
        P = Lz_vec*Dz_vec*(B-delta^2*(diag(U)*dx_vec*W+diag(W)*dz_vec*W-dz_vec*diag(E_v)*dz_vec*W));
    end
    if method == 2  % solve system in single step
        
        % define operators
        Lu = Ro*(dx_vec*Lz_vec*diag(U) + dz_vec*diag(W)) - dz_vec*diag(E_v)*dz_vec;
        Lb = Ro*(dx_vec*Lz_vec*diag(U) + dz_vec*diag(W)) - dz_vec*diag(E_k)*dz_vec;
        O = zeros(Nx*Nz);
        ix = 1:Nx;
        
        % define linear system
        M = [Lu -I O O dx_vec; I Lu O O O; O O delta^2*Lu -I dz_vec; Ro*diag(dBdx0) O O Lb O; dx_vec O dz_vec O O]; % linear operator on (U,V,W,B,P)
        R = [-dPdx0; zeros(4*Nx*Nz,1)];                 % RHS vector
        
        % define boundary conditions
        M([ix (Nz-1)*Nx+ix],:) = [dz_vec(ix,:) zeros(Nx,4*Nx*Nz); dz_vec((Nz-1)*Nx+ix,:) zeros(Nx,4*Nx*Nz)];                                          % U BC
        M([Nx*Nz+ix (2*Nz-1)*Nx+ix],:) = [zeros(Nx,Nx*Nz) dz_vec(ix,:) zeros(Nx,3*Nx*Nz); zeros(Nx,Nx*Nz) dz_vec((Nz-1)*Nx+ix,:) zeros(Nx,3*Nx*Nz)];  % V BC
        M([3*Nx*Nz+ix (4*Nz-1)*Nx+ix],:) = [zeros(Nx,3*Nx*Nz) dz_vec(ix,:) zeros(Nx,Nx*Nz); zeros(Nx,3*Nx*Nz) dz_vec((Nz-1)*Nx+ix,:) zeros(Nx,Nx*Nz)];% B BC
        M(2*Nz*Nx+ix,:) = [zeros(Nx,2*Nx*Nz) I(ix,:) zeros(Nx,2*Nx*Nz)];                                                                              % W BC
        M((3*Nz+ceil(Nz/2)-1)*Nx+ix,:) = [zeros(Nx,3*Nx*Nz) Az_vec((ceil(Nz/2)-1)*Nx+ix,:) zeros(Nx,Nx*Nz)];                                          % DA(B)
        M((4*Nz+ceil(Nz/2)-1)*Nx+ix,:) = [zeros(Nx,4*Nx*Nz) Az_vec((ceil(Nz/2)-1)*Nx+ix,:)];                                                          % DA(P)
        R([1:Nx (Nz-1)*Nx+ix]) = 0;                                             % set BC values for only non-zero part of RHS vector
        
        % calculate new solution, calculate differences and update solution
        Q = M\R;
        dU = sqrt(dA*sum((U-Q(1:Nx*Nz)).^2)); dV = sqrt(dA*sum((V-Q(Nx*Nz+1:2*Nx*Nz)).^2)); dB = sqrt(dA*sum((B-Q(3*Nx*Nz+1:4*Nx*Nz)).^2));
        U = Q(1:Nx*Nz); V = Q(Nx*Nz+1:2*Nx*Nz); W = Q(2*Nx*Nz+1:3*Nx*Nz); B = Q(3*Nx*Nz+1:4*Nx*Nz); P = Q(4*Nx*Nz+1:5*Nx*Nz);
    end
    if out == 1; disp(['(U,V,B) iteration difference: (' num2str(dU) ',' num2str(dV) ',' num2str(dB) ')']); end   % display differences

    % calculate error
    err = sqrt(dU^2+dV^2+dB^2);
end

% reshape outputs back to Nx x Nz arrays, calculate derived fields
U = vec2mat(U); V = vec2mat(V); W = vec2mat(W); B = vec2mat(B); P = vec2mat(P);
M2 = vec2mat(dBdx0+dx_vec*mat2vec(B)); N2 = vec2mat(dz_vec*mat2vec(B)); Psi = vec2mat(Dz_vec*mat2vec(U));

% add background buoyancy to B if inc_B0 = 1
if inc_B0 == 1; B = B + (integ_mat(dBdx,x)-1)*ones(1,length(z)); end

end