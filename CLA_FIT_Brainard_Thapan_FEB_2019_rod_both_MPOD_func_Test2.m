function [maxrsq] = CLA_FIT_Brainard_Thapan_FEB_2019_rod_both_MPOD_func_Test2(rodY,ofY,ofB,rodB,mp,ma,ivbd,fileStruct)

%% Thapan Data
supp_A = fileStruct.old_thapan.Supp;
A = [fileStruct.old_thapan.Wavelengths,fileStruct.old_thapan.SPD];
CLA_A = CLA_rod_both_MPOD_calculation_Test2(A, rodY, ofY, ofB, rodB, mp, ma,ivbd,fileStruct);

%% Brainard et al. monochromatic suppressions and spectra
supp_B = fileStruct.old_brainard.Supp;
B = [fileStruct.old_brainard.Wavelengths,fileStruct.old_brainard.SPD];
CLA_B = CLA_rod_both_MPOD_calculation_Test2(B, rodY, ofY, ofB, rodB, mp, ma,ivbd,fileStruct);
%% Compute logistic
CLA = vertcat(CLA_A', CLA_B');
supp = vertcat(supp_A, supp_B);
    

%% search for min error using least squares method
% Constants
A = 0.7;
B = 355.7;
C = 1.1024;
av = mean(supp);

% Functional Fit
fit = A*(1 - (1./(1 + (CLA./B).^C)));

% Calc error
err = (supp - fit).^2;
err(CLA < 90) = 30*err(CLA < 90);

% Total error
Serr = (supp - fit).^2;
Stot = (supp - av).^2;

total_err = sum(err);
maxrsq = 1 - (sum(Serr)/sum(Stot));

