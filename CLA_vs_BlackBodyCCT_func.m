
function CLA_vs_BlackBodyCCT_func(rodY, ofY, ofB, rodB, mp, ma,ivdb,g,fileStruct,TFplot,a2,a3)
% Lightning bolt for TRB
Tc = 1/7000:0.00001:1/2500;
Tc = 1./Tc;
wave = (380:1:780)';
lightLevel = 300; % lux
% fileStruct = loadAllTextFiles2;
% rodY = 1.8;
% ofY = 1.4;
% ofB = 0.6;
% rodB = 1.8;
% mp = 0.2;
% ma = 0.35;
% ivdb = 4.3;
% a2 = .2;
% a3 = 0.6;

for loop = 1:length(Tc)
    spd = blackBodySpectra23Sep05(Tc(loop),wave);
    spd = spd*lightLevel/Lxy23Sep05([wave spd]);
%     CLA(loop) = CLA_postBerlinCorrMelanopsin_06Feb2014([wave spd]);
    CLA(loop) =CLA_rod_both_MPOD_calculation_Test6([wave spd], rodY, ofY, ofB, rodB, mp, ma,ivdb,g,fileStruct,a2,a3);
end
if TFplot
    figure(2)
    plot(Tc,CLA,'b-','LineWidth',2)
    xlabel('Correlated color temperature (CCT), K')
    ylabel('CL_A')
    set(gca,'YLim',[0,600],'XLim',[0,7000])
    title(sprintf('Parameters: ofB %0.2f, rodB %0.2f, ofY %0.2f, rodY %0.2f',ofB,rodB,ofY,rodY));
end