K_def

z0 = 1/(2*sqrt(E));

A_coeff = 1/(K_pp(z0,z0)^2+(z0+K(z0,z0))^2)*(K_pp(z0,z0)*((11/80-3*Pr/4)*K_pp(z0,z0)-1/16*z0^2*(K(z0,z0)+z0))+(K(z0,z0)+z0)*((23/80-Pr/4)*K(z0,z0)+1/16*z0^2*K_pp(z0,z0)+(1/10-3*Pr/2)*z0));
B_coeff = 1/(K_pp(z0,z0)^2+(z0+K(z0,z0))^2)*(-(K(z0,z0)+z0)*((11/80-3*Pr/4)*K_pp(z0,z0)-1/16*z0^2*(K(z0,z0)+z0))+K_pp(z0,z0)*((23/80-Pr/4)*K(z0,z0)+1/16*z0^2*K_pp(z0,z0)+(1/10-3*Pr/2)*z0));

P1_da = integral(@(z) K_p(z/sqrt(E),z0).^2,-0.5,0.5);
P2_da = integral(@(z) K_pp(z/sqrt(E),z0).^2,-0.5,0.5)+Pr*(2*sqrt(E)*K_pp(z0,z0)+1/(24*E));

P1 = @(z) A_coeff*(K_p(z,z0)+1)+B_coeff*K_ppp(z,z0)+Pr/4*(z.*K_pp(z,z0)+2*K_p(z,z0)+6)+2/5*K_p(z,z0).^2-1/2*K(z,z0).*K_pp(z,z0)-1/10*K_ppp(z,z0).^2+1/16*z.*K_pp(z,z0)-1/5*K_p(z,z0)-1/16*z.^2.*K_ppp(z,z0)-1/10 - P1_da;
P2 = @(z) B_coeff*(K_p(z,z0)+1)-A_coeff*K_ppp(z,z0)+Pr/4*(z.*K(z,z0)+3*z.^2)-3/10*K(z,z0).^2+1/5*K_pp(z,z0).^2-1/2*K_p(z,z0).*K_ppp(z,z0)-1/16*z.^2.*K_p(z,z0)-23/80*z.*K(z,z0)-1/20*z.^2 - P2_da;

Q1_pp = Pr*(2/3*K_pp(z/sqrt(E),z0).*K(z/sqrt(E),z0)-1/3*K_p(z/sqrt(E),z0).^2 + P1_da)+2/3*P1(z/sqrt(E));  % integrate twice with no-flux BCs to get Q_1
