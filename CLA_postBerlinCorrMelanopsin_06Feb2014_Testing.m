function CLA = CLA_postBerlinCorrMelanopsin_06Feb2014_Testing(spd,fileStruct, varargin)

% July 2019 RN 

%% Split imported SPD
[rows columns] = size(spd);
if columns > 2
    error('Not column oriented data. Try transposing spd');
end
wavelength_spd = spd(:,1);
spd = spd(:,2);

%% Calculate GAI
GAI = GamutArea23Sep05([wavelength_spd, spd]) * 13600;
vd = exp(1.1-(1.1/(1+GAI)));

%% Calculate Duv
[~,DC,~] = CRI23Sep05([wavelength_spd, spd]);

%% Display Values
tVal = table(GAI,vd,DC);
%disp(tVal);

%% Generate Physical Values
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

%% MPOD Calculations
%----------------CHANGE HERE for MPOD---------------------------------------
p = 0.35; % Percent corneal stimulus passing through macula
MPOD = 0.2; % Estimated MPOD of the subject to  put in calculations
thickness_exp = 2*MPOD; %  MPOD THICKNESS  -------------------------
macularT_exp = 10.^(-Macula(:,2)*thickness_exp);
macularTi_exp = interp1(Macula(:,1),macularT_exp,wavelength_spd,'linear',1.0);
%---------------------------------------------------------------------------
spd = p*spd.*macularTi_exp + (1-p)*spd;
%-------------------------------------------------------------------------

%% Weighted Photoreceptor responses
% weighted responses for photoreceptors
vl_response = trapz(wavelength_spd,Vlambda.*spd);
scone_response = trapz(wavelength_spd,Scone.*spd);
rod_response = trapz(wavelength_spd,Vprime.*spd);
mel_response = trapz(wavelength_spd,M.*spd);

%% Brightness Response Caluculation
scone_over_mel = scone_response/mel_response;

BF_eff_func = fileStruct.CIE31by1;
wave = BF_eff_func(:,1);
BF_Vlambda = interp1(wave,BF_eff_func(:,3),wavelength_spd,'linear',0.0);
g = 3; % original g
%g = scone_over_mel; % 
BrightnessFunction = BF_Vlambda + g*Scone;
brightness = BrightnessFunction/max(BrightnessFunction); % normalize to max=1 (luminous efficiency)
%--------------------------------------------------------------------------------------------------------------------

brightness_response = trapz(wavelength_spd,brightness.*spd);
rod_over_brightness = rod_response/brightness_response;
c1 = 0.81;
c2 = 0.3;
rod_over_brightness_E = c1*exp(1-c2/rod_over_brightness);

%--------------------------------------------------------------------------------------------------------------------
%% Define rodSat
rodSat = 35000; % Scotopic Trolands
retinalE = [1 3 10 30 100 300 1000 3000 10000 30000 100000];
pupilDiam = [7.1 7 6.9 6.8 6.7 6.5 6.3 5.65 5 3.65 2.3];
diam = interp1(retinalE,pupilDiam,rodSat,'linear');
rodSat = rodSat/(diam^2/4*pi)*pi/1700;

%% Define CLA Constants
% original
a1 = 1.0;      
b1 = 0.0;                           
a2 = 0.7000;
b2 = 0.0;                    
k =  0.2616;                                      
a3 = 0*3.300;  % a_rod  was 3.3 originally - made 0 in the new model and new rod threshold added

% new
ofB = 0.85;
rodB = 0.45; % rod inhibition multiplier on the cool side - is 1 originally
ofY = 3.1;  % of_last_used = 1.5; of = 0.55 for McClung optimized; 0.8 for MC+WL+Lum; 0.75 for for MC+WL+Lum (0.2 MPOD 35%)
rodY = 0.95; % rod inhibition multiplier on the warm side - is 0 originally

%% B - Y Switch
 if (trapz(wavelength_spd,Scone.*spd)-k*trapz(wavelength_spd,Vlambda.*spd)) >= 0
            %% B - Y Response
            %% CS1 Calculation
            CS1 = a1*trapz(wavelength_spd,M.*spd)-b1;

            if CS1 < 0
                CS1(CS1 < 0) = 0; % remove negative values that are below threshold set by constant b1.
            end
            
            %% CS2 Calculation
            CS2 = a2*(trapz(wavelength_spd,Scone.*spd)-k*trapz(wavelength_spd,Vlambda.*spd))-b2;

            if CS2 < 0
                CS2(CS2 < 0) = 0; % This is the important diode operator, the (b-y) term cannot be less than zero
            end
            
            %% Rod Calculation
            Rod = a3*(1-exp(-trapz(wavelength_spd,Vprime.*spd)/rodSat)); %*(1 - exp(-20*(trapz(wavelength_spd,Scone.*spd)-k*trapz(wavelength_spd,V10.*spd))));
%             Rod = rodB*rod_over_brightness*(1-exp(-trapz(wavelength_spd,Vprime.*spd)/rodSat)); %*(1 - exp(-20*(trapz(wavelength_spd,Scone.*spd)-k*trapz(wavelength_spd,V10.*spd))));
            %disp(Rod)
            
            
            %% Rod Calculation
%             CS = (CS1 + CS2 - Rod);
%             CS = ofB*(CS1 + CS2 - Rod - rodB*(rod_over_brightness)*(1-exp(-trapz(wavelength_spd,Vprime.*spd)/rodSat)));
            CS = ofB*(CS1 + CS2 - Rod - rodB*(rod_over_brightness_E)*(1-exp(-trapz(wavelength_spd,Vprime.*spd)/rodSat)));

            if CS < 0
                CS(CS < 0) = 0; % Rod inhibition cannot make the CS less than zero
            end
            %disp('(B-Y) > 0')
 else
            %% Melanopic Response
%             CS = a1*trapz(wavelength_spd,M.*spd)-b1;
%             CS = ofY*(a1*trapz(wavelength_spd,M.*spd)-b1 - rodY*(rod_over_brightness)*(1-exp(-trapz(wavelength_spd,Vprime.*spd)/rodSat)));
            CS = ofY*(a1*trapz(wavelength_spd,M.*spd)-b1 - rodY*(rod_over_brightness_E)*(1-exp(-trapz(wavelength_spd,Vprime.*spd)/rodSat)));

            if CS < 0
                CS(CS < 0) = 0; % Negative values mean stimulus is below threshold set by constant b1
            end
            %disp('(B-Y) < 0')
        end


CLA = CS*1547.9; % orignally used to sets CLA equal to photopic value for 1000 lux of 2856 K. Was 1622.5 prior to 27-Jun-2014
% CSe = 0.7*(1-(1/(1+(CLA/355.7)^1.1026)));
