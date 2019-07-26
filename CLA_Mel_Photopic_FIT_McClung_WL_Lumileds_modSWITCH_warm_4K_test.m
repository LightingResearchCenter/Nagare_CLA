%
clear
close('all')
fileStruct = loadAllTextFiles2();
typeoffit = 'original'; % best or original

ZA = [fileStruct.white_light_data.Wavelengths,fileStruct.white_light_data.SPD(:,logical(fileStruct.white_light_data.A))];%
supp_A = fileStruct.white_light_data.Supp(logical(fileStruct.white_light_data.A));
% CLA_A_test = ZA(:,2)';
% Melanopic_A = ZA(:,3)';
% lux_A = ZA(:,4)';    
% A = fileStruct.spd_McClung_2700K_1lux;

ZC = [fileStruct.white_light_data.Wavelengths,fileStruct.white_light_data.SPD(:,logical(fileStruct.white_light_data.C))];%
supp_C = fileStruct.white_light_data.Supp(logical(fileStruct.white_light_data.C));
% CLA_C_test = ZC(:,2)';
% Melanopic_C = ZC(:,3)';
% lux_C = ZC(:,4)';
% C = fileStruct.spd_WL_2700K_1lux;

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

rodY=1.2;
ofY=1.5;
ofB=0.9500;
rodB=0.6000;
mp=0.2;
ma=0.35;
ivdb=.85;

CLA_A = CLA_rod_both_MPOD_calculation_Test2(ZA, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct);
CLA_C = CLA_rod_both_MPOD_calculation_Test2(ZC, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct);
CLA_E = CLA_rod_both_MPOD_calculation_Test2(ZE, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct);
CLA_F = CLA_rod_both_MPOD_calculation_Test2(ZF, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct);
CLA_G = CLA_rod_both_MPOD_calculation_Test2(ZG, rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct);


% %**********
% Compute logistic
%************

CLA = vertcat(CLA_A',CLA_C',CLA_E',CLA_F',CLA_G');
supp = vertcat(supp_A,supp_C,supp_E,supp_F,supp_G);
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

    %plot results
    xAxis_A = CLA_A;
    xAxis_C = CLA_C;
    xAxis_E = CLA_E;
    xAxis_F = CLA_F;
    xAxis_G = CLA_G;
    xLim = [-1,5];
    xLabel = 'log10(CL_A)';
                
    figure(1)
    hold on
    scatter(log10(xAxis_A), supp_A, [130], [0.9290 0.6940 0.1250],'d', 'filled');
    scatter(log10(xAxis_C), supp_C, [130], [0.9290 0.6940 0.1250],'+');
    scatter(log10(xAxis_E), supp_E, [130], [0.9290 0.6940 0.1250],'o', 'filled');
    scatter(log10(xAxis_F), supp_F, [130], [0.9290 0.6940 0.1250],'*');
    scatter(log10(xAxis_G), supp_G, [130], [0.9290 0.6940 0.1250],'s');

    x = xLim(1):0.01:xLim(2); %-2:.01:5;
    %x = 10.^x;
    fit = minA*(1 - (1./(1 + (10.^x/minB).^minC)));
    plot(x, fit, 'k', 'linewidth', 2)
    ylim([-0.1 0.8])
    HL1 = legend('2700K (M)', '2700K (WL)', '3000K (L)', 'Cyan-gap (L)', '4000K (L)','location', 'northwest');%, , 'curve fit')
    xlabel(xLabel,'FontSize',14)
    ylabel('Melatonin suppression','FontSize',14)
    titleStr = '1-h data; Polychromatic warm with 4K, 0.2MPOD35'; % 

%     title({'Original fit'; titleStr},'FontSize',14)% *********************************
%     title({'Eratio c2 1: ofY 4.4,rodY 3.4, a3 0, 0.2@35'; titleStr},'FontSize',13)% *********************************
    title({'Parameters: ofY 3.2,rodY 3.5'; titleStr},'FontSize',13)% *********************************
%     title({'Parameters: ofY 1.4,rodY 1.1'; titleStr},'FontSize',13)% *********************************

    temp = find(fit > .11);
    % rsquare = num2str(maxrsq);
    % rsquare = round(rsquare,2);
    text(xLim(1)+4, .22, ['r^2: ', num2str(round(maxrsq,4))],'FontSize',24);
    text(xLim(1)+4, .12, ['Threshold: ', num2str(round(10^x(temp(1)),1))],'FontSize',16);
    % text(xLim(1)+0.25, .35, ['max. resp. (a): ', num2str(round(minA,2))],'FontSize',16);
    text(xLim(1)+4, .06, ['half-sat (b): ', num2str(round(minB,1))],'FontSize',16);
    text(xLim(1)+4, .0, ['rate (c): ', num2str(round(minC,4))],'FontSize',16)
    % text(xLim(1)+0.5, .15, ['max. response: ', num2str(minA)],'FontSize',12)
    set(gca,'XLim',xLim);
    set(gca,'FontSize',16);
    set(HL1,'FontSize',16);     % does not matter as all set by gca
    %         max_rsq = maxrsq
    hold off
