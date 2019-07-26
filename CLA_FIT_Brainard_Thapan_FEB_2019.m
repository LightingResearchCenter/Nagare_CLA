%
clear
close('all')
fileStruct = loadAllTextFiles();
fileStruct2 = loadAllTextFiles2();
wave = (380:1:780)';
% Melanopic weighting function
Melanopsin = fileStruct.MelanopsinWlensBy2nm_02Oct2012; % lens data from Wyszecki and Stiles Table 1(2.4.6) Norren and Vos(1974) data
Melanopic = interp1(Melanopsin(:,1),Melanopsin(:,2),wave,'linear',0.0);
%M = M/macularTi;
Melanopic = Melanopic/max(Melanopic);
%M555 = M(wave==555);
%Melanopic = M/M555;

% Thapan et al. monochromatic suppressions and spectra
% wavelengths = fileStruct.old_thapan_suppressions.Wavelengths;
% irr = fileStruct.old_thapan_suppressions.Irr;
% supp = fileStruct.old_thapan_suppressions.Supp;

% wavelengths = fileStruct.old_thapan_suppressions_below500.Wavelengths;
% irr = fileStruct.old_thapan_suppressions_below500.Irr;
% supp = fileStruct.old_thapan_suppressions_below500.Supp;

wavelengths = fileStruct.old_thapan_suppressions_below500.Wavelengths;
irr = fileStruct.old_thapan_suppressions_below500.Irr;
supp = fileStruct.old_thapan_suppressions_below500.Supp;
rodY = 1.55;
ofY = 1.5;
ofB = 0.71;
rodB = 1.48;
mp = 0.2;
ma = 0.35;
ivbd = 4.85;
wavelengthsThapan = wavelengths; % save for later plotting
%correct for dilated pupils
irr = 8.6*irr; 
%correct for time if thapan
irr = irr/1.8;
ind = 0;
for i1 = 1:length(wavelengths)
    try
    A = load(['T', num2str(wavelengths(i1)), '-1.txt']);
    catch
        continue
    end
    ind = ind+1;
    MelanopicInterp = interp1(wave,Melanopic,A(:,1),'linear',0.0);
    %A(:,2) = A(:,2)/max(A(:,2));
    A(:,2) = (A(:,2)*irr(ind))/trapz(A(:,1), A(:,2));
    supp_A(ind) = supp(ind);
%     lux_A(ind) = Lxy23Sep05([A(:,1) A(:,2)]);
    irrad_A(ind) = trapz(A(:,1), A(:,2));
%     CLA_A(i1) = CLA_postBerlinCorrMelanopsin_06Feb2014([A(:,1) A(:,2)],fileStruct);
    CLA_A(ind) = CLA_rod_both_MPOD_calculation_Test3([A(:,1) A(:,2)], rodY, ofY, ofB, rodB, mp, ma,ivbd,fileStruct2);
Melanopic_A(ind) = 843*trapz(A(:,1), MelanopicInterp.*A(:,2));
    %CLA(i) = spdtolux([A(:,1) A(:,2)]);
    %temp(i) = trapz(A(:,1), A(:,2));
end

% Brainard et al. monochromatic suppressions and spectra

wavelengths = fileStruct.old_brainard_suppressions_below500.Wavelengths;
irr = fileStruct.old_brainard_suppressions_below500.Irr;
supp = fileStruct.old_brainard_suppressions_below500.Supp;

% wavelengths = fileStruct.old_brainard_suppressions_below500.Wavelengths;
% irr = fileStruct.old_brainard_suppressions_below500.Irr;
% supp = fileStruct.old_brainard_suppressions_below500.Supp;

% wavelengths = fileStruct.old_brainard_suppressions_above500.Wavelengths;
% irr = fileStruct.old_brainard_suppressions_above500.Irr;
% supp = fileStruct.old_brainard_suppressions_above500.Supp;

wavelengthsBrainard = wavelengths; % Save for later plotting 
%correct for dilated pupils
irr = 8.6*irr;
ind = 0;
for i1 = 1:length(wavelengths)
    try
    A = load(['B', num2str(wavelengths(i1)), '-1.txt']);
    catch
        continue
    end
    ind = ind+1;
    
    MelanopicInterp = interp1(wave,Melanopic,A(:,1),'linear',0.0);
    %A(:,2) = A(:,2)/max(A(:,2));
    A(:,2) = (A(:,2)*irr(ind))/trapz(A(:,1), A(:,2));
    supp_B(ind) = supp(ind);
%     lux_B(i1) = Lxy23Sep05([A(:,1) A(:,2)]);
    irrad_B(ind) = trapz(A(:,1), A(:,2));
    CLA_B(ind) = CLA_rod_both_MPOD_calculation_Test3([A(:,1) A(:,2)], rodY, ofY, ofB, rodB, mp, ma,ivbd,fileStruct2);
    Melanopic_B(ind) = 843*trapz(A(:,1), MelanopicInterp.*A(:,2));
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
    for B = 355.7:1:355.7 % B = 200:1:600 || 355.7:1:355.7
        k = 1;
        
        %C is the exponent (was .864)
        for C = 1.1024:0.001:1.1024   %  1.1026:0.001:1.1026 || .5:.001:1.55
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
        j = j + 1
    end
    l = l + 1
end

%plot results
        xAxis_A = CLA_A;
        xAxis_B = CLA_B;
        xLim = [-1,5];
        xLabel = 'log10(CL_A)';
%         xLabel = 'log10(CL_A New Model [k = 0.2883])';
    
figure(1)
hold on
scatter(log10(xAxis_A), supp_A, [], [0 0 1],'s', 'filled');
scatter(log10(xAxis_B), supp_B, [], [0 .5 0], 'd', 'filled');

x = xLim(1):0.01:xLim(2); %-2:.01:5;
%x = 10.^x;
fit = minA*(1 - (1./(1 + (10.^x/minB).^minC)));
plot(x, fit, 'k', 'linewidth', 2)

HL1 = legend('Thapan 2001','Brainard 2001', 'location', 'northwest');%, , 'curve fit')
% HL1 = legend('Thapan 2001 > 500nm','Brainard 2001 > 500nm', 'location', 'northwest');%, , 'curve fit')

xlabel(xLabel,'FontSize',16)
ylabel('Melatonin suppression','FontSize',16)
titleStr = sprintf('ofY=%0.2f,rodY=%0.2f,ofB=%0.2f,rodB=%0.2f, 0.2 MPOD @35%',ofY,rodY,ofB,rodB); 
title(titleStr,'FontSize',16)

temp = find(fit > .11);
text(xLim(1)+0.5, .35, ['r^2: ', num2str(round(maxrsq,4))],'FontSize',22);
text(xLim(1)+0.5, .25, ['Threshold: ', num2str(10^x(temp(1)))],'FontSize',12);
text(xLim(1)+0.5, .20, ['half saturation: ', num2str(minB)],'FontSize',12);
text(xLim(1)+0.5, .15, ['rate: ', num2str(minC)],'FontSize',12)
set(gca,'XLim',xLim);
set(gca,'FontSize',16);
set(HL1,'FontSize',16);
hold off


