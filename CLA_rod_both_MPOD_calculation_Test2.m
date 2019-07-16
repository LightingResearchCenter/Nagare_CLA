function CLA = CLA_rod_both_MPOD_calculation_Test2(spd, rodY, ofY, ofB, rodB, mp, ma,fileStruct)

% RN July 2019
wavelength_spd = spd(:,1);
spd = spd(:,2:end);

GAI_test = GamutArea23Sep05_test([wavelength_spd, spd])' * 13600;
vd = exp(1.1-(1.1./(1+GAI_test)));%1-(1./(1+GAI_test));%

rodY = rodY * vd;
ofY = ofY * vd;
ofB = ofB * vd;
rodB = rodB * vd;


Vlamda = fileStruct.Vlamda;
Vlambda = interp1(Vlamda(:,1),Vlamda(:,2),wavelength_spd,'linear',0.0);

Vprime = fileStruct.Vprime;
Vprime = interp1(Vprime(:,1),Vprime(:,2),wavelength_spd,'linear',0.0);
Vprime = Vprime/max(Vprime);

Scone = fileStruct.Scone;
Scone = interp1(Scone(:,1),Scone(:,2),wavelength_spd,'linear',0.0);

Macula = fileStruct.MacularPigmentODfromSnodderly;
thickness = 1.0; % macular thickness factor
macularT = 10.^(-Macula(:,2)*thickness);
macularTi = interp1(Macula(:,1),macularT,wavelength_spd,'linear',1.0);

Scone = Scone./macularTi;
Scone = Scone/(max(Scone));

Vlambda = Vlambda./macularTi;
Vlambda = Vlambda/max(Vlambda);

%Melanopsin = load('Melanopsin with corrected lens.txt');
Melanopsin = fileStruct.MelanopsinWlensBy2nm_02Oct2012; % lens data from Wyszecki and Stiles Table 1(2.4.6) Norren and Vos(1974) data
M = interp1(Melanopsin(:,1),Melanopsin(:,2),wavelength_spd,'linear',0.0);
%M = M/macularTi;
M = M/max(M);

%----------------CHANGE HERE for MPOD---------------------------------------
% p = 0.35; % Percent corneal stimulus passing through macula **** ma
% MPOD = 0.0; % Estimated MPOD of the subject to  put in calculations *** mp
thickness_exp = 2*mp; %  MPOD THICKNESS  -------------------------
macularT_exp = 10.^(-Macula(:,2)*thickness_exp);
macularTi_exp = interp1(Macula(:,1),macularT_exp,wavelength_spd,'linear',1.0);
%---------------------------------------------------------------------------
spd = ma*spd.*macularTi_exp + (1-ma)*spd;
%-------------------------------------------------------------------------

% weighted responses
vl_response = trapz(wavelength_spd,Vlambda.*spd);
scone_response = trapz(wavelength_spd,Scone.*spd);
rod_response = trapz(wavelength_spd,Vprime.*spd);
mel_response = trapz(wavelength_spd,M.*spd);
scone_over_mel = scone_response/mel_response;


BF_eff_func = fileStruct.CIE31by1;
wave = BF_eff_func(:,1);
BF_Vlambda = interp1(wave,BF_eff_func(:,3),wavelength_spd,'linear',0.0);
 g = 3; 
%g = scone_over_mel; 
BrightnessFunction = BF_Vlambda + g*Scone;
brightness = BrightnessFunction/max(BrightnessFunction); %  normalize to max=1 (luminous efficiency)

brightness_response = trapz(wavelength_spd,brightness.*spd);
rod_over_brightness = rod_response./brightness_response;
c1 = 0.81;
c2 = 0.3;
rod_over_brightness_E = c1*exp(1-c2./rod_over_brightness);

%--------------------------------------------------------------------------------------------------------------------

rodSat = 35000; % Scotopic Trolands
retinalE = [1 3 10 30 100 300 1000 3000 10000 30000 100000];
pupilDiam = [7.1 7 6.9 6.8 6.7 6.5 6.3 5.65 5 3.65 2.3];
diam = interp1(retinalE,pupilDiam,rodSat,'linear');
rodSat = rodSat/(diam^2/4*pi)*pi/1700;

a1 = 1.0;      % 1.0 originally        %0.285 % Melanopsin correction factor
b1 = 0.0;                             %0.01  
a2 = 0.7000;  % a_(b-y)  was 0.6201 prior to 06Feb2014     %0.2
b2 = 0.0;                    %0.001
k =  0.2616;                  %0.31 -- 0.2616 original / 0.2883 for 4K switch / 0.2436 for CG switch                           %0.001
a3 = 0*3.300;  % a_rod  was 3.3 originally - made 0 now as a new rod threshold term added 

%-------------------------------------------------------------------------
P = spd;
            %% Set B - Y
            byIdx = (trapz(wavelength_spd,Scone.*spd)-k*trapz(wavelength_spd,Vlambda.*spd)) >= 0;
            CS = nan(size(byIdx));
            
            %% B - Y Response
            CS1 = (a1*trapz(wavelength_spd,M.*spd)-b1).*byIdx;

            CS1(CS1 < 0) = 0; % remove negative values that are below threshold set by constant b1.

            CS2 = (a2*(trapz(wavelength_spd,Scone.*spd)-k*trapz(wavelength_spd,Vlambda.*spd))-b2).*byIdx;
            CS2(CS2 < 0) = 0; % This is the important diode operator, the (b-y) term cannot be less than zero

            Rod = (a3*(1-exp(-trapz(wavelength_spd,Vprime.*spd)./rodSat))).*byIdx; %*(1 - exp(-20*(trapz(wavelength_spd,Scone.*spd)-k*trapz(wavelength_spd,V10.*spd))));
            %disp(Rod)
            
    %         CS = (CS1 + CS2 - Rod);
            CSa = (ofB.*(CS1 + CS2 - Rod - rodB.*(rod_over_brightness).*(1-exp(-trapz(wavelength_spd,Vprime.*spd)/rodSat)))).*byIdx;
            %CS = ofB*(CS1 + CS2 - Rod - rodB*(rod_over_brightness_E)*(1-exp(-trapz(wavelength_spd,Vprime.*spd)/rodSat)));
            
            CSa(CSa < 0) = 0; % Rod inhibition cannot make the CS less than zero
            %disp('(B-Y) > 0')
        %% Melaniopsin Response
    %             CS = a1*trapz(wavelength_spd,M.*spd)-b1;
    %         CS = ofY*a1*trapz(wavelength_spd,M.*P)-b1 - rodY*(a3*(1-exp(-trapz(wavelength_spd,Vprime.*spd)/rodSat)));
            CSb = (ofY.*(a1*trapz(wavelength_spd,M.*P)-b1 - rodY.*(rod_over_brightness).*(1-exp(-trapz(wavelength_spd,Vprime.*spd)/rodSat)))).*(~byIdx);
            %CS = ofY*(a1*trapz(wavelength_spd,M.*P)-b1 - rodY*(rod_over_brightness_E)*(1-exp(-trapz(wavelength_spd,Vprime.*spd)/rodSat)));

            CSb(CSb < 0) = 0; % Negative values mean stimulus is below threshold set by constant b1
            %disp('(B-Y) < 0')
        %end
CS(byIdx) = CSa(byIdx);
CS(~byIdx) = CSb(~byIdx);
CLA = CS.*1547.9; % originally used to CLA equal to photopic value for 1000 lux of 2856 K. Was 1622.5 prior to 27-Jun-2014 
% CSe = 0.7*(1-(1/(1+(CLA/355.7)^1.1026)));
