function [maxrsq,minA,minB,minC] = CLA_FIT_Brainard_Thapan_FEB_2019_rod_both_MPOD_func_Test3(rodY,ofY,ofB,rodB,mp,ma,ivbd,g,fileStruct, varargin)

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

%% Thapan Data
% supp_A = fileStruct.old_thapan.Supp(fileStruct.old_thapan.warm);
supp_A = fileStruct.old_thapan.Supp;
% A = [fileStruct.old_thapan.Wavelengths,fileStruct.old_thapan.SPD(:,fileStruct.old_thapan.warm)];
A = [fileStruct.old_thapan.Wavelengths,fileStruct.old_thapan.SPD];
CLA_A = CLA_rod_both_MPOD_calculation_Test6(A, rodY, ofY, ofB, rodB, mp, ma,ivbd,g,fileStruct,testA2,testA3);

%% Brainard et al. monochromatic suppressions and spectra
% supp_B = fileStruct.old_brainard.Supp(fileStruct.old_brainard.warm);
supp_B = fileStruct.old_brainard.Supp;
% B = [fileStruct.old_brainard.Wavelengths,fileStruct.old_brainard.SPD(:,fileStruct.old_brainard.warm)];
B = [fileStruct.old_brainard.Wavelengths,fileStruct.old_brainard.SPD];
CLA_B = CLA_rod_both_MPOD_calculation_Test6(B, rodY, ofY, ofB, rodB, mp, ma,ivbd,g,fileStruct,testA2,testA3);
%% Compute logistic
CLA = vertcat(CLA_A', CLA_B');
supp = vertcat(supp_A, supp_B);
    

%% search for min error using least squares method
% Constants
% A = 0.7;
% B = 355.7;
% C = 1.1026;
% av = mean(supp);
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
    for B = 355.7% B = 200:1:600 || 355.7:1:355.7
        k = 1;
        
        %C is the exponent (was .864)
        for C = 1.1026:0.001:1.1026   %  1.1026:0.001:1.1026 || .5:.001:1.55
            for i = 1:length(supp)
                fit(i) = A*(1 - (1/(1 + (CLA(i)/B)^C)));
                
                %weight small values more heavily to compensate for lack of
                %density there
%                 if(CLA(i) < 90)
%                     err(i) = 30*(supp(i) - fit(i))^2;
%                 else
                    err(i) = (supp(i) - fit(i))^2;
%                 end
%                 err(i) = (supp(i) - fit(i))^2;
                
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
% Functional Fit
% fit = A*(1 - (1./(1 + (CLA./B).^C)));

% Calc error
% err = (supp - fit).^2;
% err(CLA < 90) = 30*err(CLA < 90);
% 
% % Total error
% Serr = (supp - fit).^2;
% Stot = (supp - av).^2;
% 
% total_err = sum(err);
% maxrsq = 1 - (sum(Serr)/sum(Stot));

