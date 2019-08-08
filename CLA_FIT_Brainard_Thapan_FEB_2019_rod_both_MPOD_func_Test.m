function [maxrsq] = CLA_FIT_Brainard_Thapan_FEB_2019_rod_both_MPOD_func_Test(rodY,ofY,ofB,rodB,mp,ma,fileStruct,varargin)
if numel(varargin) == 0 
    testA2 = 0.7;
    testA3 = 0;
elseif numel(varargin) == 1 
    testA2 = 0.7;
    testA3 = varargin{1};
elseif numel(varargin) == 2
    testA2 = varargin{1};
    testA3 = varargin{2};
else
    error('Too many inputs')
end 

wave = (380:1:780)';
% Melanopic weighting function
Melanopsin = fileStruct.MelanopsinWlensBy2nm_02Oct2012; % lens data from Wyszecki and Stiles Table 1(2.4.6) Norren and Vos(1974) data
Melanopic = interp1(Melanopsin(:,1),Melanopsin(:,2),wave,'linear',0.0);
%M = M/macularTi;
Melanopic = Melanopic/max(Melanopic);
%M555 = M(wave==555);
%Melanopic = M/M555;

% Thapan et al. monochromatic suppressions and spectra
wavelengths = fileStruct.old_thapan_suppressions.Wavelengths;
irr = fileStruct.old_thapan_suppressions.Irr;
supp = fileStruct.old_thapan_suppressions.Supp;

% wavelengths = fileStruct.old_thapan_suppressions_below500.Wavelengths;
% irr = fileStruct.old_thapan_suppressions_below500.Irr;
% supp = fileStruct.old_thapan_suppressions_below500.Supp;

% wavelengths = fileStruct.old_thapan_suppressions_above500.Wavelengths;
% irr = fileStruct.old_thapan_suppressions_above500.Irr;
% supp = fileStruct.old_thapan_suppressions_above500.Supp;

wavelengthsThapan = wavelengths; % save for later plotting
%correct for dilated pupils
irr = 8.6*irr; 
%correct for time if thapan
irr = irr/1.8;
for i1 = 1:length(wavelengths)
    A = load(['T', num2str(wavelengths(i1)), '-1.txt']);
    MelanopicInterp = interp1(wave,Melanopic,A(:,1),'linear',0.0);
    %A(:,2) = A(:,2)/max(A(:,2));
    A(:,2) = (A(:,2)*irr(i1))/trapz(A(:,1), A(:,2));
    supp_A(i1) = supp(i1);
     lux_A(i1) = Lxy23Sep05([A(:,1) A(:,2)]);
    irrad_A(i1) = trapz(A(:,1), A(:,2));
    CLA_A(i1) = CLA_rod_both_MPOD_calculation_Test([A(:,1) A(:,2)], lux_A(i1), rodY, ofY, ofB, rodB, mp, ma,fileStruct,testA2,testA3);
    Melanopic_A(i1) = 843*trapz(A(:,1), MelanopicInterp.*A(:,2));
    %CLA(i) = spdtolux([A(:,1) A(:,2)]);
    %temp(i) = trapz(A(:,1), A(:,2));
end

% Brainard et al. monochromatic suppressions and spectra
wavelengths = fileStruct.old_brainard_suppressions.Wavelengths;
irr = fileStruct.old_brainard_suppressions.Irr;
supp = fileStruct.old_brainard_suppressions.Supp;

% wavelengths = fileStruct.old_brainard_suppressions_below500.Wavelengths;
% irr = fileStruct.old_brainard_suppressions_below500.Irr;
% supp = fileStruct.old_brainard_suppressions_below500.Supp;

% wavelengths = fileStruct.old_brainard_suppressions_above500.Wavelengths;
% irr = fileStruct.old_brainard_suppressions_above500.Irr;
% supp = fileStruct.old_brainard_suppressions_above500.Supp;

wavelengthsBrainard = wavelengths; % Save for later plotting 
%correct for dilated pupils
irr = 8.6*irr; 
for i1 = 1:length(wavelengths)
    A = load(['B', num2str(wavelengths(i1)), '-1.txt']);
    MelanopicInterp = interp1(wave,Melanopic,A(:,1),'linear',0.0);
    %A(:,2) = A(:,2)/max(A(:,2));
    A(:,2) = (A(:,2)*irr(i1))/trapz(A(:,1), A(:,2));
    supp_B(i1) = supp(i1);
    lux_B(i1) = Lxy23Sep05([A(:,1) A(:,2)]);
    irrad_B(i1) = trapz(A(:,1), A(:,2));
    CLA_B(i1) = CLA_rod_both_MPOD_calculation_Test([A(:,1) A(:,2)], lux_B(i1),rodY, ofY, ofB, rodB, mp, ma, fileStruct,testA2,testA3);
    Melanopic_B(i1) = 843*trapz(A(:,1), MelanopicInterp.*A(:,2));
end


%********************************************************************
% Compute logistic
%********************************************************************
%combine sets

CLA = vertcat(CLA_A', CLA_B');
supp = vertcat(supp_A', supp_B');
    

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
    for B = 355.7:1:355.7   % B = 200:1:600 || 355.7:1:355.7 | 330:0.05:370
        k = 1;
        
        %C is the exponent (was .864)
        for C = 1.1026:0.001:1.1026    %  1.1026:0.001:1.1026 || .5:.001:1.55 || 1.0500:.0001:1.1500
            for i = 1:length(supp)
                fit(i) = A*(1 - (1/(1 + (CLA(i)/B)^C)));
                
                %weight small values more heavily to compensate for lack of
                %density there
                if(CLA(i) < 90)
                    err(i) = 30*(supp(i) - fit(i))^2;
                else
                    err(i) = (supp(i) - fit(i))^2;
                end
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
