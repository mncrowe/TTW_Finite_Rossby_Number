% Identify the critical Rossby number E = E(Ro), uses Matlab 'parfor' to parallelise computation

clear
addpath("functions")

% define numerical parameters:

n_iter = 10;
err_tol = 1e-5;
err_max = 10;
Nx = 64;
Nz = 51;
method = 2;
out = 0;
N_E = 161;
N_Ro = 201;

savedata = 1;
savename = ["data_Pr_1.mat", "data_Pr_2.mat", "data_Pr_4.mat", "data_Pr_8.mat", "data_Pr_05.mat", "data_Pr_025.mat", "data_Pr_0125.mat"];
dirname = "test2/";
Pr_vec = [1 2 4 8 0.5 0.25 0.125];

% define physical parameters:

delta = 0;
dbdx = @(x) exp(-(x*sqrt(pi)/2).^2);
Ro_vec = 10.^(linspace(-2,3,N_Ro));
E_vec = 10.^(linspace(-4,0,N_E));

p = gcp('nocreate');
if isempty(p)
    parpool(length(Pr_vec))
end
clear p

mkdir(dirname)

parfor iPr = 1:length(Pr_vec)

    Pr = Pr_vec(iPr);

    % loop through parameters determining which cases converge:
    
    conv = zeros(length(E_vec),length(Ro_vec));
    
    tic
    for iE = 1:length(E_vec)
        for iRo = 1:length(Ro_vec)
    
            calc = 1;
    
            if iRo > 1
                if conv(iE,iRo-1) == 0
                    conv(iE,iRo) = 0; calc = 0;
                end
            end
    
            if iE > 1
                if conv(iE-1,iRo) == 1
                    conv(iE,iRo) = 1; calc = 0;
                end
            end
            
            if calc == 1
                [~,~,~,~,~,~,~,~,~,~,err] = TTW_State(n_iter,E_vec(iE),Ro_vec(iRo),Pr,delta,Nx,Nz,0,method,dbdx,err_tol,err_max,out);
                if err < err_tol; conv(iE,iRo) = 1; end
            end
    
            disp(['(Pr, Ro, E) = (' num2str(Pr) ', ' num2str(Ro_vec(iRo)) ', ' num2str(E_vec(iE)) '), conv = ' num2str(conv(iE,iRo))])
        end
    end
    toc
    
    % determine critical Rossby number as a function of E:
    
    Ro_i = sum(conv,2);
    Ro_c = E_vec;
    for iE = 1:length(E_vec)
        i1 = max(1,Ro_i(iE));
        i2 = min(i1+1,length(Ro_vec));
        Ro_c(iE) = 10^(1/2*(log(Ro_vec(i1))+log(Ro_vec(i2)))/log(10));
        %Ro_c(iE) = 1/2*(Ro_vec(i1)+Ro_vec(i2));
    end
    
    % save data:
    
    if savedata == 1
        parsave(dirname + savename(iPr), E_vec, Ro_c)
    end

end

function parsave(fname, E_vec, Ro_c)
    save(fname, 'E_vec', 'Ro_c')
end
