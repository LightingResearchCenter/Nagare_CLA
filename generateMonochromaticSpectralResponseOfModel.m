%Brainard and Thappan Data
% [wavelength  response]
BandT = [
%     420 0.256; ...
%     424 0.815; ...
%     440 0.953; ...
%     456 1.000; ...
%     460 1.000; ...
%     472 0.8560; ...
%     480 0.7259; ...
%     496 0.5869; ...
    505 0.7916; ...
    520 0.5202; ...
    530 0.3958; ...
    548 0.1428; ...
    555 0.1089; ...
    575 0.0554; ...
    600 0.0282
    ];
waveBT = BandT(:,1);
BT = BandT(:,2);
fileStruct = loadAllTextFiles2();
wave = (420:5:600)'; % Wavelength range for calculations

criterion = 300; % Criterion response in units of CLA ("Circadian Light"). CLA is numerically equal to lux for CIE Illuminant A

rodY=1.2;
ofY=1.5;
ofB=0.9500;
rodB=0.6000;
mp=0.2;
ma=0.35;
ivdb=0.85;

criterionIrrad = zeros(size(wave)); % initialize array
for j = 1:length(wave)
    spd = zeros(length(wave),1);
    spd(j) = 1;
    steps = -4:0.2:3; % steps of log10(irradiance) (W/m^2)
    CLA = zeros(length(steps),1);
    irrad = zeros(length(steps),1);
    for i = 1:length(steps)
        irrad(i) = 10^steps(i); % irradiance (W/m^2)
        P = spd*irrad(i); % scaled spd
 
        % CLA is "Circadian Light" as given by the phototranduction model
        CLA(i) = CLA_rod_both_MPOD_calculation_Test2([wave P], rodY, ofY, ofB, rodB, mp, ma,ivdb,fileStruct);
        
    end
    criterionIrrad(j) = interp1(CLA,irrad,criterion,'linear'); % find the irradiance that gives the criterion response
%     display(['Wavelength = ' num2str(wave(j),'%d') ' nm; Stimulus for criterion response = ' num2str(criterionIrrad(j)) ' W/m^2']);
end
for j = 1:length(wave)
    display(['Wavelength = ' num2str(wave(j),'%d') ' nm; Stimulus for criterion response = ' num2str(criterionIrrad(j)) ' W/m^2']);
end
Efficacy = 1./criterionIrrad; % efficacy  = 1/(criterion response)
Efficiency = Efficacy/max(Efficacy); % normalize to a maximum of 1

%figure(1)
%hold off
figure(2)
plot(wave,Efficiency,'k--')
axis([400 650 0 1.2])
hold on
plot(waveBT,BT,'rd')
% HL1 = legend('Original Model', 'Brainard and Thapan');
% HL1 = legend('Brainard and Thapan: Rod-Brightness model', 'rodB-ofB 0.1-1.0 a3-1.0, rodY-ofY 0.1-1.2, g3');
%         HL1 = legend('E-ratio; g= s/mel; ofB 0.85, rodB 0.75, ofY 1.5, rodY 1.1', 'Brainard and Thapan 0.2MPOD35');
HL1 = legend(sprintf('Parameters: ofB %0.2f, rodB %0.2f, ofY %0.2f, rodY %0.2f',ofB,rodB,ofY,rodY),'Brainard and Thapan 0.2MPOD35');
% HL1 = legend('Rod-brightness: ofY = ofB = 1.4, rodY = 0.2, rodB = 0.08', 'Brainard and Thapan: 0.2 MPOD at 35%');
xlabel('Wavelength (nm)','FontSize',16);
ylabel('Efficiency','FontSize',16);
set(gca,'FontSize',16);
set(HL1,'FontSize',13);
hold off
