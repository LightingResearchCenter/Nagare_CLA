% Lightning bolt for TRB
Tc = 1/7000:0.00001:1/500;
Tc = 1./Tc;
wave = (380:1:780)';
lightLevel = 300; % lux
fileStruct = loadAllTextFiles2;
rodY = 1.8*4.3;
ofY = 1.3*4.3;
ofB = 1.3*0.7;
rodB = 1.3*1.39*0.7;
mp = 0.2;
ma = 0.35;
ivdb = 1;
a2 = .45;
a3 = 1;
g = 0.2;

for loop = 1:length(Tc)
    spd = blackBodySpectra23Sep05(Tc(loop),wave);
    spd = spd*lightLevel/Lxy23Sep05([wave spd]);
%     CLA(loop) = CLA_postBerlinCorrMelanopsin_06Feb2014([wave spd]);
    CLA(loop) =CLA_rod_both_MPOD_calculation_Test2([wave spd], rodY, ofY, ofB, rodB, mp, ma,ivdb,g,fileStruct,a2,a3);
end

figure(1)
plot(Tc,CLA,'b-','LineWidth',2)
xlabel('Correlated color temperature (CCT), K')
ylabel('CL_A')
set(gca,'YLim',[0,600],'XLim',[0,7000])
