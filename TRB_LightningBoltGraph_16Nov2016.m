% Lightning bolt for TRB
Tc = 1/7000:0.00001:1/500;
Tc = 1./Tc;
wave = (380:1:780)';

Melanopsin = load('MelanopsinWlensBy2nm_02Oct2012.txt'); % lens data from Wyszecki and Stiles Table 1(2.4.6) Norren and Vos(1974) data
M = interp1(Melanopsin(:,1),Melanopsin(:,2),wave,'linear',0.0);
%M = M/macularTi;
M = M/max(M);
%M555 = M(wave==555);
%Melanopic = M/M555;
MelanopicFcn = M;
k = 1.5483e+03;

for loop = 1:length(Tc)
    spd = blackBodySpectra23Sep05(Tc(loop),wave);
    spd = spd*300/Lxy23Sep05([wave spd]);
    CLA(loop) = CLA_postBerlinCorrMelanopsin_06Feb2014([wave spd]);
    Melanopic(loop) = k*trapz(wave,MelanopicFcn.*spd);
end

%A = load('C:\AndyOct2012\CircadianModeling\F32T841K.txt');
A = load('C:\AndyOct2012\CircadianModeling\WhiteLEDsCLaCalculation\LumiledsLXML_PW51_4000K.txt');
A(:,2) = A(:,2)*300/Lxy23Sep05(A);
Tc1(1) = CCT23Sep05(A);
spd1 = interp1(A(:,1),A(:,2),wave,'linear',0.0);
CLA1(1) = CLA_postBerlinCorrMelanopsin_06Feb2014([wave spd1]);

A = load('C:\AndyOct2012\CircadianModeling\WhiteLEDsCLaCalculation\LumiledsLXM7_PW40_4000K.txt'); % This one
A(:,2) = A(:,2)*300/Lxy23Sep05(A);
Tc1(2) = CCT23Sep05(A);
spd2 = interp1(A(:,1),A(:,2),wave,'linear',0.0);
CLA1(2) = CLA_postBerlinCorrMelanopsin_06Feb2014([wave spd2]);

A = load('C:\AndyOct2012\CircadianModeling\WhiteLEDsCLaCalculation\LumiledsLXM8_PW30_3000K.txt'); % LumiledsLXM8_PW30_3000K
%A = load('C:\AndyNewStuff\MarkStuff\SEA_EfficacyRegulationPlot_Feb2016\whiteLED_Spectra\SeoulSTW8C2SA_3700To4700K.txt'); 
A(:,2) = A(:,2)*300/Lxy23Sep05(A);
Tc1(3) = CCT23Sep05(A);
spd3 = interp1(A(:,1),A(:,2),wave,'linear',0.0);
CLA1(3) = CLA_postBerlinCorrMelanopsin_06Feb2014([wave spd3]);

A = load('C:\AndyNewStuff\MarkStuff\SEA_EfficacyRegulationPlot_Feb2016\whiteLED_Spectra\Samsung362A_3000.txt'); %This One
A(:,2) = A(:,2)*300/Lxy23Sep05(A);
Tc1(4) = CCT23Sep05(A);
spd4 = interp1(A(:,1),A(:,2),wave,'linear',0.0);
CLA1(4) = CLA_postBerlinCorrMelanopsin_06Feb2014([wave spd4]);

A = load('SPD4000KHighCLA.txt');
A(:,2) = A(:,2)*300/Lxy23Sep05(A);
%Tc1(5) = CCT23Sep05(A);
spd5 = interp1(A(:,1),A(:,2),wave,'linear',0.0);
%spd5(wave>475 & wave<485) = 15*spd5(wave==480);
%q = (wave>470 & wave<=475) | (wave>=485&wave<490);
%spd5(q) = 5*spd5(wave==475);
%spd5(wave<450) = 0.0;
CLA1(5) = CLA_postBerlinCorrMelanopsin_06Feb2014([wave spd5]);
Tc1(5) = CCT23Sep05([wave spd5]);

A = load('SPD3000KLowCLA.txt');
A(:,2) = A(:,2)*300/Lxy23Sep05(A);
spd6 = interp1(A(:,1),A(:,2),wave,'linear',0.0);
%spd6(wave>455 & wave<500) = 0.0;
%spd6(wave<450) = 2*spd6(wave<450);
CLA1(6) = CLA_postBerlinCorrMelanopsin_06Feb2014([wave spd6]);
Tc1(6) = CCT23Sep05([wave spd6]);

figure(1)
plot(Tc,CLA,'b-','LineWidth',2)
hold on
plot(Tc1(1:4),CLA1(1:4),'bd','LineWidth',2)
plot(Tc1(5:6),CLA1(5:6),'rd','LineWidth',2)
hold off
xlabel('Correlated color temperature (CCT), K')
ylabel('CL_A')
set(gca,'YLim',[0,600],'XLim',[0,7000])

figure(3)
plot(wave,spd5,'b-')

[Tc1;CLA1]'

% Suppression at 18 and 27 lux
LED4000K = spd2;
LED3000K = spd4;
LED3000K = LED3000K*18/Lxy23Sep05([wave,LED3000K]);
LED3000KCLA = CLA_postBerlinCorrMelanopsin_06Feb2014([wave LED3000K])
LED3000KCS = CSCalc_postBerlin_12Aug2011(LED3000KCLA)
LED4000K = LED4000K*18/Lxy23Sep05([wave,LED4000K]);
LED4000KCLA = CLA_postBerlinCorrMelanopsin_06Feb2014([wave LED4000K])
LED4000KCS = CSCalc_postBerlin_12Aug2011(LED4000KCLA)

LED3000K = LED3000K*27/Lxy23Sep05([wave,LED3000K]);
LED3000KCLA = CLA_postBerlinCorrMelanopsin_06Feb2014([wave LED3000K])
LED3000KCS = CSCalc_postBerlin_12Aug2011(LED3000KCLA)
LED4000K = LED4000K*27/Lxy23Sep05([wave,LED4000K]);
LED4000KCLA = CLA_postBerlinCorrMelanopsin_06Feb2014([wave LED4000K])
LED4000KCS = CSCalc_postBerlin_12Aug2011(LED4000KCLA)

% CLA and lux at 30 scotopic lux
Vprime = load('C:\AndyOct2012\MatLabWork\PhotometryCalc\scotopicBy2.txt');
Vprime = interp1(Vprime(:,1),Vprime(:,2),wave,'linear',0.0);
LED3000K = LED3000K*30/(1700*trapz(wave,Vprime.*LED3000K));
LED3000KCLA = CLA_postBerlinCorrMelanopsin_06Feb2014([wave LED3000K])
LED3000Lux = Lxy23Sep05([wave,LED3000K])
LED4000K = LED4000K*30/(1700*trapz(wave,Vprime.*LED4000K));
LED4000KCLA = CLA_postBerlinCorrMelanopsin_06Feb2014([wave LED4000K])
LED4000Lux = Lxy23Sep05([wave,LED4000K])

figure(4)
plot(wave,LED3000K/max(LED3000K),'r-','LineWidth',2);
hold on
plot(wave,LED4000K/max(LED4000K),'b-','LineWidth',2)
hold off
xlabel('Wavelength (nm)','FontSize',14)
ylabel('Relative spectral power','FontSize',14)
set(gca,'XLim',[380,750])
set(gca,'FontSize',14)
legend('3000 K','4000 K');

% HPS
HPS = load('C:\AndyOct2012\CircadianModeling\HPS400W.txt');
HPS = interp1(HPS(:,1),HPS(:,2),wave,'linear',0.0);
HPS = HPS*27/Lxy23Sep05([wave,HPS]);
HPSCLA = CLA_postBerlinCorrMelanopsin_06Feb2014([wave HPS])
HPSCS = CSCalc_postBerlin_12Aug2011(HPSCLA)

% Plot SPD of commercially available LEDs
figure(5)
plot(wave,spd1/max(spd1),'b-','LineWidth',2)
set(gca,'XLim',[380,750])
xlabel('Wavelength (nm)','FontSize',12)
ylabel('Relative spectral power','FontSize',12)
title('LumiledsLXML-PW51-4000K')

figure(6)
plot(wave,spd2/max(spd2),'b-','LineWidth',2)
set(gca,'XLim',[380,750])
xlabel('Wavelength (nm)','FontSize',12)
ylabel('Relative spectral power','FontSize',12)
title('LumiledsLXM7-PW40-4000K')

figure(7)
plot(wave,spd3/max(spd3),'b-','LineWidth',2)
set(gca,'XLim',[380,750])
xlabel('Wavelength (nm)','FontSize',12)
ylabel('Relative spectral power','FontSize',12)
title('LumiledsLXM8-PW30-3000K')

figure(8)
plot(wave,spd4/max(spd4),'b-','LineWidth',2)
set(gca,'XLim',[380,750])
xlabel('Wavelength (nm)','FontSize',12)
ylabel('Relative spectral power','FontSize',12)
title('Samsung362A-3000')
