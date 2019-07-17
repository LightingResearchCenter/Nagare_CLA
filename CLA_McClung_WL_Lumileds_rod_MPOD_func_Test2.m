function [maxrsq] = CLA_McClung_WL_Lumileds_rod_MPOD_func_Test2(rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct)

    typeoffit = 'original'; % best or original **** CHECK crange for BEST CASE

    %% Load Suppressions
    
    %% Load SPDs
    white_light_data = [fileStruct.white_light_data.Wavelengths,fileStruct.white_light_data.SPD];
    
    %% Combine SPDs 

    CLA = CLA_rod_both_MPOD_calculation_Test2(white_light_data, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct);

    CLA = CLA';
    supp = fileStruct.white_light_data.Supp;



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

