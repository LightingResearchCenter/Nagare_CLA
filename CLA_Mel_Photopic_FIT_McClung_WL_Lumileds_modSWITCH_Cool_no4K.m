%
clear
close('all')
fileStruct = loadAllTextFiles();

typeoffit = 'original'; % best or original

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

CLA_B = CLA_postBerlinCorrMelanopsin_06Feb2014_1luxspd(B, lux_B);
CLA_D = CLA_postBerlinCorrMelanopsin_06Feb2014_1luxspd(D, lux_D);


% %**********
% Compute logistic
%************

CLA = vertcat(CLA_B',CLA_D');
supp = vertcat(supp_B',supp_D');
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
    xAxis_B = CLA_B;
    xAxis_D = CLA_D;
    xLim = [-1,5];
    xLabel = 'log10(CL_A)';
                
    figure(1)
    hold on
    scatter(log10(xAxis_B), supp_B, [130], [0 0.4470 0.7410],'d', 'filled');
    scatter(log10(xAxis_D), supp_D, [130], [0 0.4470 0.7410],'+');

    x = xLim(1):0.01:xLim(2); %-2:.01:5;
    %x = 10.^x;
    fit = minA*(1 - (1./(1 + (10.^x/minB).^minC)));
    plot(x, fit, 'k', 'linewidth', 2)
    ylim([-0.1 0.8])
    HL1 = legend('6500K (M)', '5600K (WL)','location', 'northwest');%, , 'curve fit')
    xlabel(xLabel,'FontSize',14)
    ylabel('Melatonin suppression','FontSize',14)
    titleStr = '1-h data; Polychromatic Cool (No 4K), 0.2MPOD35'; % 

%     title({'Original fit'; titleStr},'FontSize',14)% *********************************
%     title({'Eratio c2 1: ofB 0.75,rodB 0.75, a3 0, 0.2@35'; titleStr},'FontSize',13)% *********************************
    title({'Parameters: ofB 0.7,rodB 0.65'; titleStr},'FontSize',13)% *********************************
%     title({'Parameters: ofB 0.8,rodB 1.3'; titleStr},'FontSize',13)% *********************************

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
