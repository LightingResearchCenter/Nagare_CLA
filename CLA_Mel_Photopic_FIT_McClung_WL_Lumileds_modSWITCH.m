%
clear
close('all')
% fileStruct = loadAllTextFiles();
fileStruct = loadAllTextFiles2();
typeoffit = 'original'; % best or original

ZA = [fileStruct.white_light_data.Wavelengths,fileStruct.white_light_data.SPD(:,logical(fileStruct.white_light_data.A))];%
supp_A = fileStruct.white_light_data.Supp(logical(fileStruct.white_light_data.A));
% CLA_A_test = ZA(:,2)';
% Melanopic_A = ZA(:,3)';
% lux_A = ZA(:,4)';    
% A = fileStruct.spd_McClung_2700K_1lux;

ZB = [fileStruct.white_light_data.Wavelengths,fileStruct.white_light_data.SPD(:,logical(fileStruct.white_light_data.B))];%
supp_B = fileStruct.white_light_data.Supp(logical(fileStruct.white_light_data.B));

% ZB = fileStruct.McClung_1h_6500K;%
% supp_B = ZB(:,1)';
% % CLA_B_test = ZB(:,2)';
% Melanopic_B = ZB(:,3)';
% lux_B = ZB(:,4)';
% B = fileStruct.spd_McClung_6500K_1lux;

ZC = [fileStruct.white_light_data.Wavelengths,fileStruct.white_light_data.SPD(:,logical(fileStruct.white_light_data.C))];%
supp_C = fileStruct.white_light_data.Supp(logical(fileStruct.white_light_data.C));
% CLA_C_test = ZC(:,2)';
% Melanopic_C = ZC(:,3)';
% lux_C = ZC(:,4)';
% C = fileStruct.spd_WL_2700K_1lux;


ZD = [fileStruct.white_light_data.Wavelengths,fileStruct.white_light_data.SPD(:,logical(fileStruct.white_light_data.D))];%
supp_D = fileStruct.white_light_data.Supp(logical(fileStruct.white_light_data.D));

% ZD = fileStruct.WL_1h_5600K_corr;%
% supp_D = ZD(:,1)';
% % CLA_D_test = ZD(:,2)';
% Melanopic_D = ZD(:,3)';
% lux_D = ZD(:,4)';
% D = fileStruct.spd_WL_5600K_1lux;

ZE = [fileStruct.white_light_data.Wavelengths,fileStruct.white_light_data.SPD(:,logical(fileStruct.white_light_data.E))];%
supp_E = fileStruct.white_light_data.Supp(logical(fileStruct.white_light_data.E));
% CLA_E_test = ZE(:,2)';
% Melanopic_E = ZE(:,3)';
% lux_E = ZE(:,4)';    
% E = fileStruct.spd_Lumileds_3000K_1lux;

ZF = [fileStruct.white_light_data.Wavelengths,fileStruct.white_light_data.SPD(:,logical(fileStruct.white_light_data.F))];%
supp_F = fileStruct.white_light_data.Supp(logical(fileStruct.white_light_data.F));
% CLA_F_test = ZF(:,2)';
% Melanopic_F = ZF(:,3)';
% lux_F = ZF(:,4)';    
% F = fileStruct.spd_Lumileds_CG_1lux;

ZG = [fileStruct.white_light_data.Wavelengths,fileStruct.white_light_data.SPD(:,logical(fileStruct.white_light_data.G))];%
supp_G = fileStruct.white_light_data.Supp(logical(fileStruct.white_light_data.G));
% CLA_G_test = ZG(:,2)';
% Melanopic_G = ZG(:,3)';
% lux_G = ZG(:,4)';    
% G = fileStruct.spd_Lumileds_4000K_1lux;

rodY = 1.8;
ofY = 1.3;
ofB = 1.3;
rodB = 1.8;
mp = 0.3;
ma = 0.2;
ivdb = 4.3;
a2 = .45;
a3 = 2.5;

% CLA_A = CLA_postBerlinCorrMelanopsin_06Feb2014_1luxspd(A, lux_A,fileStruct);
% CLA_B = CLA_postBerlinCorrMelanopsin_06Feb2014_1luxspd(B, lux_B,fileStruct);
% CLA_C = CLA_postBerlinCorrMelanopsin_06Feb2014_1luxspd(C, lux_C,fileStruct);
% CLA_D = CLA_postBerlinCorrMelanopsin_06Feb2014_1luxspd(D, lux_D,fileStruct);
% CLA_E = CLA_postBerlinCorrMelanopsin_06Feb2014_1luxspd(E, lux_E,fileStruct);
% CLA_F = CLA_postBerlinCorrMelanopsin_06Feb2014_1luxspd(F, lux_F,fileStruct);
% CLA_G = CLA_postBerlinCorrMelanopsin_06Feb2014_1luxspd(G, lux_G,fileStruct);



CLA_A = CLA_rod_both_MPOD_calculation_Test5(ZA, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct,a2,a3);
CLA_B = CLA_rod_both_MPOD_calculation_Test5(ZB, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct,a2,a3);
CLA_C = CLA_rod_both_MPOD_calculation_Test5(ZC, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct,a2,a3);
CLA_D = CLA_rod_both_MPOD_calculation_Test5(ZD, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct,a2,a3);
CLA_E = CLA_rod_both_MPOD_calculation_Test5(ZE, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct,a2,a3);
CLA_F = CLA_rod_both_MPOD_calculation_Test5(ZF, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct,a2,a3);
CLA_G = CLA_rod_both_MPOD_calculation_Test5(ZG, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct,a2,a3);

% %**********
% Compute logistic
%************

CLA = vertcat(CLA_A',CLA_B',CLA_C',CLA_D',CLA_E',CLA_F',CLA_G');
supp = vertcat(supp_A,supp_B,supp_C,supp_D,supp_E,supp_F,supp_G);
Brange = 200:1:600;


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
        for B = 355.7  % B = 200:1:600 ---- B = Brange -----
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

    %plot results
    xAxis_A = CLA_A;
    xAxis_B = CLA_B;
    xAxis_C = CLA_C;
    xAxis_D = CLA_D;
    xAxis_E = CLA_E;
    xAxis_F = CLA_F;
    xAxis_G = CLA_G;
    xLim = [-1,5];
    xLabel = 'log10(CL_A)';
                
    figure(1)
    hold on
    scatter(log10(xAxis_A), supp_A, [130], [0.9290 0.6940 0.1250],'d', 'filled');
    scatter(log10(xAxis_B), supp_B, [130], [0 0.4470 0.7410],'d', 'filled');
    scatter(log10(xAxis_C), supp_C, [130], [0.9290 0.6940 0.1250],'+');
    scatter(log10(xAxis_D), supp_D, [130], [0 0.4470 0.7410],'+');
    scatter(log10(xAxis_E), supp_E, [130], [0.9290 0.6940 0.1250],'o', 'filled');
    scatter(log10(xAxis_F), supp_F, [130], [0.9290 0.6940 0.1250],'*');
    scatter(log10(xAxis_G), supp_G, [130], [0.9290 0.6940 0.1250],'s');

    x = xLim(1):0.01:xLim(2); %-2:.01:5;
    %x = 10.^x;
    fit = minA*(1 - (1./(1 + (10.^x/minB).^minC)));
    plot(x, fit, 'k', 'linewidth', 2)
    ylim([-0.1 0.8])
    HL1 = legend('2700K (M)', '6500K (M)', '2700K (WL)', '5600K (WL)', '3000K (L)', 'Cyan-gap (L)', '4000K (L)','location', 'northwest');%, , 'curve fit')
    xlabel(xLabel,'FontSize',14)
    ylabel('Melatonin suppression','FontSize',14)
    titleStr = '1-h data; Polychromatic All, 0.2MPOD35'; % 

%     title({'Original fit'; titleStr},'FontSize',14)% *********************************
%     title({'Eratio c2 1: ofB 0.75,rodB 0.75, ofY 4.4,rodY 3.4, a3 0, 0.2@35'; titleStr},'FontSize',13)% *********************************
%     title({'Parameters: ofB 0.7,rodB 0.65, ofY 3.2,rodY 3.5'; titleStr},'FontSize',13)% *********************************
    title({sprintf('Parameters: of %0.2f,rod %0.2f, mp %0.2f, ma %0.2f',ofB,rodB,mp,ma); titleStr},'FontSize',13)% *********************************

    temp = find(fit > .11);
    % rsquare = num2str(maxrsq);
    % rsquare = round(rsquare,2);
    text(xLim(1)+4, .32, ['Polychromatic'],'FontSize',18);
    text(xLim(1)+4, .22, ['r^2: ', num2str(round(maxrsq,2))],'FontSize',18);
    text(xLim(1)+4, .12, ['Threshold: ', num2str(round(10^x(temp(1)),1))],'FontSize',16);
    % text(xLim(1)+0.25, .35, ['max. resp. (a): ', num2str(round(minA,2))],'FontSize',16);
    text(xLim(1)+4, .06, ['half-sat (b): ', num2str(round(minB,1))],'FontSize',16);
    text(xLim(1)+4, .0, ['rate (c): ', num2str(round(minC,4))],'FontSize',16)
    % text(xLim(1)+0.5, .15, ['max. response: ', num2str(minA)],'FontSize',12)
    set(gca,'XLim',xLim);
    set(gca,'FontSize',14);
    set(HL1,'FontSize',16);     % does not matter as all set by gca
    %         max_rsq = maxrsq
    hold off
