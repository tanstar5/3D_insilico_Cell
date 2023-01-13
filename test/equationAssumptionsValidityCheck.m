tot_lips = 10000;
PA = 0.5;
PB = 1 - PA;
nA = round(6*rand(tot_lips,1));
nB = 6-nA;

cA = -10;
cB = 10;

H_spont = (nA*cA + nB*cB)/6;

H_mean = 0;
del_hij = H_mean - mean(H_spont)
var_H_spont = var(H_spont)


disp('Analytical')

time_points = 1000;
%generating_environment for t th time points
nA_t = zeros(time_points,1);
nB_t = zeros(time_points,1);
H_spont_t = zeros(time_points,1);
for i = 1:time_points
    nA_t(i) = binornd(6,PA);
    nB_t(i) = 6 - nA_t(i);
    H_spont_t(i) = (cA*nA_t(i) + cB*nB_t(i))/6; 
end
histogram(H_spont_t)
var_H_spont = var(H_spont_t)  %% varisnce of time averaged neighborhood of lipid
mean_H_spont = mean(H_spont_t)
H_mean_noteq = 1;
del_hij_mean = mean(H_spont_t - H_mean_noteq)
del_hij_square_mean = mean((H_spont_t - H_mean_noteq).^2)
disp('Analytical')
var_H_spont_A = 1/6*(cA-cB)^2*PA*PB
mean_H_spont_A = cA*PA + cB*PB
del_hij_mean_A = mean_H_spont_A-H_mean_noteq
del_hij_square_mean = del_hij_mean_A^2 + var_H_spont_A


