% load and process data files produced by Critical_Ro.m, the data files
% produced by 'Critical_Ro.m' or 'Critical_Ro_Parallel.m' should be
% accessible on the Path

savename = ["data_Pr_0125.mat", "data_Pr_025.mat", "data_Pr_05.mat", "data_Pr_1.mat", "data_Pr_2.mat", "data_Pr_4.mat", "data_Pr_8.mat"];

N_E = 161;
N_Ro = 201;
N_Pr = 7;

Ro_c_temp = zeros(N_Pr, N_E);

for i = 1:length(savename)

    load(savename(i), 'Ro_c', 'E_vec', 'Ro_vec')

    Ro_c_temp(i, :) = smoothdata(Ro_c, "sgolay");

end

Ro_c = Ro_c_temp;

clear Ro_c_temp

%Ro_s = 0 * E_vec;

%for i = 1:length(E_vec)

    %Ro_s(i) = 1 / (sqrt(E_vec(i)) * abs(K(1/2/sqrt(E_vec(i)), 1/2/sqrt(E_vec(i))))) * sqrt(exp(1) / 2 / pi);

%end

figure
loglog(E_vec, Ro_c, 'LineWidth', 2)
xlabel('E'); ylabel('Ro_c'); grid
legend('Pr = 0.125', 'Pr = 0.25', 'Pr = 0.5', 'Pr = 1', 'Pr = 2', 'Pr = 4', 'Pr = 8', 'Location', 'northwest')
set(gca, 'FontSize', 12, 'linewidth', 0.7);

figure
loglog(E_vec, Ro_s, 'LineWidth', 2)
xlabel('E'); ylabel('Ro_s'); grid; ylim([1e-2 1e2])