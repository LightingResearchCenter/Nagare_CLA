function [maxrsq] = CLA_McClung_WL_Lumileds_rod_MPOD_func_Test2(rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct)

    typeoffit = 'original'; % best or original **** CHECK crange for BEST CASE

    %% Load Suppressions
    ZA = fileStruct.McClung_1h_2700K;%
    supp_A = ZA(:,1)';

    ZB = fileStruct.McClung_1h_6500K;%
    supp_B = ZB(:,1)';

    ZC = fileStruct.WL_1h_2700K_corr;%
    supp_C = ZC(:,1)';

    ZD = fileStruct.WL_1h_5600K_corr;%
    supp_D = ZD(:,1)';

    ZE = fileStruct.Lumileds_1h_3000K;%
    supp_E = ZE(:,1)';

    ZF = fileStruct.Lumileds_1h_CG;%
    supp_F = ZF(:,1)';

    ZG = fileStruct.Lumileds_1h_4000K;%
    supp_G = ZG(:,1)';
    
    %% Load SPDs
    load('white_light_data.mat');
    
    %% Combine SPDs 

    CLA = CLA_rod_both_MPOD_calculation_Test2(white_light_data, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct);

    CLA = CLA';
    supp = vertcat(supp_A',supp_B',supp_C',supp_D',supp_E',supp_F',supp_G');



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

