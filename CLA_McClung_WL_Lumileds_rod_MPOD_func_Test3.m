function [maxrsq,minA,minB,minC] = CLA_McClung_WL_Lumileds_rod_MPOD_func_Test3(rodY, ofY, ofB, rodB, mp, ma,ivdb,g,fileStruct,varargin)

    typeoffit = 'original'; % best or original **** CHECK crange for BEST CASE
if numel(varargin) == 0 
    testA2 = 0.7;
    testA3 = 0;
elseif numel(varargin) == 1 
    testA2 = 1;
    testA3 = varargin{1};
elseif numel(varargin) == 2
    testA2 = varargin{1};
    testA3 = varargin{2};
else
    error('Too many inputs')
end 
    %% Load Suppressions
    
    %% Load SPDs
    white_light_data = [fileStruct.white_light_data.Wavelengths,fileStruct.white_light_data.SPD(:,~fileStruct.white_light_data.warm)];
%     white_light_data = [fileStruct.white_light_data.Wavelengths,fileStruct.white_light_data.SPD];
    
    %% Combine SPDs 

    CLA = CLA_rod_both_MPOD_calculation_Test2(white_light_data, rodY, ofY, ofB, rodB, mp, ma,ivdb,g,fileStruct,testA2,testA3);

    CLA = CLA';
    supp = fileStruct.white_light_data.Supp(~fileStruct.white_light_data.warm);
%     supp = fileStruct.white_light_data.Supp;
% 


%% search for min error using least squares method
% Constants
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
        for B = 355.7:1:355.7  % B = 200:1:600 ---- B = Brange -----
            k = 1;

            %C is the exponent (was .864)
            for C = 1.1026:0.0001:1.1026 % C = .5:.001:1.5 best / 1.1026:0.0001:1.1026 original
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
