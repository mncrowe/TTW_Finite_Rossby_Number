% Identify the critical Rossby number E = E(Ro)
% Saves progress and can be resumed in case of Matlab crashes

clear
addpath("functions")

% define numerical parameters:

n_iter = 10;
err_tol = 1e-5;
err_max = 10;
Nx = 64;       % 64
Nz = 61;       % 51            
method = 2;
out = 0;
N_E = 161; %161;
N_Ro = 201; %201;

savedata = 1;
savename = ["data_Pr_1", "data_Pr_2", "data_Pr_4", "data_Pr_8", "data_Pr_05", "data_Pr_025", "data_Pr_0125"];
dirname = "data_61/";
Pr_vec = [1 2 4 8 0.5 0.25 0.125];
resume_on = 1;

% define physical parameters:

delta = 0;
dbdx = @(x) exp(-(x*sqrt(pi)/2).^2);
Ro_vec = 10.^(linspace(-2,3,N_Ro));
E_vec = 10.^(linspace(-4,0,N_E));

if ~exist(dirname, "dir")
    mkdir(dirname)
end

for iPr = 1:length(Pr_vec)

    Pr = Pr_vec(iPr);

    if ~exist(dirname + savename(iPr), "dir")
        mkdir(dirname + savename(iPr))
    end

    % loop through parameters determining which cases converge:
    
    conv = zeros(length(E_vec),length(Ro_vec));
    
    tic
    for iE = 1:length(E_vec)

        dat_file = dirname + savename(iPr) + "/" + num2str(iPr) + "_" + num2str(iE) + ".mat";

        if exist(dat_file,"file") & resume_on == 1

            conv_temp = load(dat_file, 'conv_temp').conv_temp;
            conv(iE, :) = conv_temp;

        else

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

            conv_temp = conv(iE, :);
            save(dirname + savename(iPr) + "/" + num2str(iPr) + "_" + num2str(iE) + ".mat", 'conv_temp')

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
        save(dirname + savename(iPr) + ".mat", "conv", "E_vec", "Ro_vec", "Pr", "Ro_c", "Nx", "Nz")
    end

end

for iPr = 1:length(Pr_vec)
    delete(dirname + savename(iPr) + "/*")
    rmdir(dirname + savename(iPr))
end
