function [maxrsq] = CLA_McClung_WL_Lumileds_warm_4K_rod_MPOD_func(rodY, ofY, ofB, rodB, mp, ma,fileStruct)

    typeoffit = 'original'; % best or original **** CHECK crange for BEST CASE

    ZA = fileStruct.McClung_1h_2700K;%
    supp_A = ZA(:,1)';
    % CLA_A_test = ZA(:,2)';
    Melanopic_A = ZA(:,3)';
    lux_A = ZA(:,4)';
    A = fileStruct.spd_McClung_2700K_1lux;

    ZC = fileStruct.WL_1h_2700K_corr;%
    supp_C = ZC(:,1)';
    % CLA_C_test = ZC(:,2)';
    Melanopic_C = ZC(:,3)';
    lux_C = ZC(:,4)';
    C = fileStruct.spd_WL_2700K_1lux;

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


    CLA_A = CLA_rod_both_MPOD_optimization_1luxspd(A, lux_A, rodY, ofY, ofB, rodB, mp, ma,fileStruct);
    CLA_C = CLA_rod_both_MPOD_optimization_1luxspd(C, lux_C, rodY, ofY, ofB, rodB, mp, ma,fileStruct);
    CLA_E = CLA_rod_both_MPOD_optimization_1luxspd(E, lux_E, rodY, ofY, ofB, rodB, mp, ma,fileStruct);
    CLA_F = CLA_rod_both_MPOD_optimization_1luxspd(F, lux_F, rodY, ofY, ofB, rodB, mp, ma,fileStruct);
    CLA_G = CLA_rod_both_MPOD_optimization_1luxspd(G, lux_G, rodY, ofY, ofB, rodB, mp, ma,fileStruct);


    CLA = vertcat(CLA_A',CLA_C',CLA_E',CLA_F',CLA_G');
    supp = vertcat(supp_A',supp_C',supp_E',supp_F',supp_G');
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
                for i = 1:length(supp)
                    fit(i) = A*(1 - (1/(1 + (CLA(i)/B)^C)));
                    err(i) = (supp(i) - fit(i))^2;
                    Serr(i) = (supp(i) - fit(i))^2;
                    Stot(i) = (supp(i) - av)^2;

                end
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

