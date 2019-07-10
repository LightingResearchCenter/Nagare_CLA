function [maxrsq] = CLA_McClung_WL_Lumileds_Cool_rod_MPOD_func(rodY, ofY, ofB, rodB, mp, ma,fileStruct)

    typeoffit = 'original'; % best or original **** CHECK crange for BEST CASE

    ZB = fileStruct.McClung_1h_6500K;%
    supp_B = ZB(:,1)';
    % CLA_B_test = ZB(:,2)';
    Melanopic_B = ZB(:,3)';
    lux_B = ZB(:,4)';
    B = fileStruct.spd_McClung_6500K_1lux;

    ZD = fileStruct.WL_1h_5600K_corr;%
    supp_D = ZD(:,1)';
    % CLA_D_test = ZD(:,2)';
    Melanopic_D = ZD(:,3)';
    lux_D = ZD(:,4)';
    D = fileStruct.spd_WL_5600K_1lux;


    CLA_B = CLA_rod_both_MPOD_optimization_1luxspd(B, lux_B, rodY, ofY, ofB, rodB, mp, ma,fileStuct);
    CLA_D = CLA_rod_both_MPOD_optimization_1luxspd(D, lux_D, rodY, ofY, ofB, rodB, mp, ma,fileStruct);

    CLA = vertcat(CLA_B',CLA_D');
    supp = vertcat(supp_B',supp_D');
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

