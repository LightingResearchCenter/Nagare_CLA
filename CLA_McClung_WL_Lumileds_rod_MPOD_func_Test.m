function [maxrsq] = CLA_McClung_WL_Lumileds_rod_MPOD_func_Test(rodY, ofY, ofB, rodB, mp, ma,fileStruct)

    typeoffit = 'original'; % best or original **** CHECK crange for BEST CASE

    ZA = fileStruct.McClung_1h_2700K;%
    supp_A = ZA(:,1)';
    % CLA_A_test = ZA(:,2)';
    Melanopic_A = ZA(:,3)';
    lux_A = ZA(:,4)';    
    A = fileStruct.spd_McClung_2700K_1lux;

    ZB = fileStruct.McClung_1h_6500K;%
    supp_B = ZB(:,1)';
    % CLA_B_test = ZB(:,2)';
    Melanopic_B = ZB(:,3)';
    lux_B = ZB(:,4)';
    B = fileStruct.spd_McClung_6500K_1lux;

    ZC = fileStruct.WL_1h_2700K_corr;%
    supp_C = ZC(:,1)';
    % CLA_C_test = ZC(:,2)';
    Melanopic_C = ZC(:,3)';
    lux_C = ZC(:,4)';
    C = fileStruct.spd_WL_2700K_1lux;

    ZD = fileStruct.WL_1h_5600K_corr;%
    supp_D = ZD(:,1)';
    % CLA_D_test = ZD(:,2)';
    Melanopic_D = ZD(:,3)';
    lux_D = ZD(:,4)';
    D = fileStruct.spd_WL_5600K_1lux;

    ZE = fileStruct.Lumileds_1h_3000K;%
    supp_E = ZE(:,1)';
    % CLA_E_test = ZE(:,2)';
    Melanopic_E = ZE(:,3)';
    lux_E = ZE(:,4)';    
    E = fileStruct.spd_Lumileds_3000K_1lux;

    ZF = fileStruct.Lumileds_1h_CG;%
    supp_F = ZF(:,1)';
    % CLA_F_test = ZF(:,2)';
    Melanopic_F = ZF(:,3)';
    lux_F = ZF(:,4)';    
    F = fileStruct.spd_Lumileds_CG_1lux;

    ZG = fileStruct.Lumileds_1h_4000K;%
    supp_G = ZG(:,1)';
    % CLA_G_test = ZG(:,2)';
    Melanopic_G = ZG(:,3)';
    lux_G = ZG(:,4)';    
    G = fileStruct.spd_Lumileds_4000K_1lux;
    
    %% Combine SPDs 

    CLA_A = CLA_rod_both_MPOD_calculation_Test(A, lux_A, rodY, ofY, ofB, rodB, mp, ma,fileStruct);
    CLA_B = CLA_rod_both_MPOD_calculation_Test(B, lux_B, rodY, ofY, ofB, rodB, mp, ma,fileStruct);
    CLA_C = CLA_rod_both_MPOD_calculation_Test(C, lux_C, rodY, ofY, ofB, rodB, mp, ma,fileStruct);
    CLA_D = CLA_rod_both_MPOD_calculation_Test(D, lux_D, rodY, ofY, ofB, rodB, mp, ma,fileStruct);
    CLA_E = CLA_rod_both_MPOD_calculation_Test(E, lux_E, rodY, ofY, ofB, rodB, mp, ma,fileStruct);
    CLA_F = CLA_rod_both_MPOD_calculation_Test(F, lux_F, rodY, ofY, ofB, rodB, mp, ma,fileStruct);
    CLA_G = CLA_rod_both_MPOD_calculation_Test(G, lux_G, rodY, ofY, ofB, rodB, mp, ma,fileStruct);


    CLA = vertcat(CLA_A',CLA_B',CLA_C',CLA_D',CLA_E',CLA_F',CLA_G');
    supp = vertcat(supp_A',supp_B',supp_C',supp_D',supp_E',supp_F',supp_G');
    Brange = 100:1:600;



switch typeoffit
    
    case 'best'
            Crange = .5:.001:1.5;
    case 'original'
            Brange = 355.7:.001:355.7;
            Crange = 1.1026:.0001:1.1026;  
    otherwise
        error('No matching case')
end

    %search for min error using least squares method
    j = 1;
    k = 1;
    l = 1;
    min = 1000;
    minA = 0;
    minB = 0;
    minC = 0;
    maxrsq = 0;
    av = mean(supp);

    for A = .7:.01:.7
        j = 1;

        %B is the term in the denominator/half-saturation value (was 215.75)
        for B = Brange  % B = 200:1:600 ---- B = Brange -----
            k = 1;

            %C is the exponent (was .864)
            for C = Crange % C = .5:.001:1.5 best / 1.1026:0.0001:1.1026 original
                fit = (A*(1 - (1./(1 + (CLA./B).^C))))';
                err = (supp' - fit).^2;
                Serr = (supp' - fit).^2;
                Stot = (supp' - av).^2;
                
                total_err(j, k) = sum(err);
                rsq(j, k) = 1 - (sum(Serr)/sum(Stot));
                if(sum(err) < min)
                    min = sum(err);
                    minA = A;
                    minB = B;
                    minC = C;
                    maxrsq = 1 - (sum(Serr)/sum(Stot));
                    minSumSerr = sum(Serr);
                    minsumStot = sum(Stot);
                end
                k = k + 1;
            end
            j = j + 1;
        end
        l = l + 1;
    end

