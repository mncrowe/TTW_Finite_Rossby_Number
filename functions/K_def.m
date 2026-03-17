% define functions K, K', K'', K''' used in all TTW work, here z0 = 1/(2*sqrt(E))

ss=@(t) sin(t).*sinh(t);
sc=@(t) sin(t).*cosh(t);
cc=@(t) cos(t).*cosh(t);
cs=@(t) cos(t).*sinh(t);
det=@(t) ss(t).^2+cc(t).^2;
A_p=@(t) (ss(t)+cc(t))/det(t);
A_m=@(t) (-ss(t)+cc(t))/det(t);

K = @(z,z0) -z+1/sqrt(2)*A_p(z0/sqrt(2)).*sc(z/sqrt(2))+1/sqrt(2)*A_m(z0/sqrt(2)).*cs(z/sqrt(2));
K_p = @(z,z0) -1+1/2*(A_p(z0/sqrt(2))-A_m(z0/sqrt(2))).*ss(z/sqrt(2))+1/2*(A_p(z0/sqrt(2))+A_m(z0/sqrt(2))).*cc(z/sqrt(2));
K_pp = @(z,z0) 1/sqrt(2)*A_p(z0/sqrt(2)).*cs(z/sqrt(2))-1/sqrt(2)*A_m(z0/sqrt(2)).*sc(z/sqrt(2));
K_ppp = @(z,z0) 1/2*(A_p(z0/sqrt(2))-A_m(z0/sqrt(2))).*cc(z/sqrt(2))-1/2*(A_p(z0/sqrt(2))+A_m(z0/sqrt(2))).*ss(z/sqrt(2));