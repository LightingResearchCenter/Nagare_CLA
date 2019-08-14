function [Duv] = calcDuv2(spd, fileStruct)

wavelength_spd = spd(:,1);
spd = spd(:,2:end);
%% Calc chromaticty 
Table = fileStruct.CIE31by1;
wavelength = Table(:,1);
xbar = Table(:,2);
ybar = Table(:,3);
zbar = Table(:,4);

xbark = interp1(wavelength,xbar,wavelength_spd);
xbark(isnan(xbark)) = 0.0;
ybark = interp1(wavelength,ybar,wavelength_spd);
ybark(isnan(ybark)) = 0.0;
zbark = interp1(wavelength,zbar,wavelength_spd);
zbark(isnan(zbark)) = 0.0;

Xk = trapz(wavelength_spd,spd .* xbark);
Yk = trapz(wavelength_spd,spd .* ybark);
Zk = trapz(wavelength_spd,spd .* zbark);

uk = (4*Xk) ./ (Xk + 15*Yk + 3*Zk);
vk = (6*Yk) ./ (Xk + 15*Yk + 3*Zk);

%% Calculate Lfp
Lfp = sqrt(((uk-0.292).^2) + ((vk-0.24).^2));

%% Calculate Lbb
a = acos((uk-0.292)/Lfp);
k = [-0.471106,1.925865,-2.4243787,1.5317403,-0.5179722,0.0893944,-0.00616793];

Lbb = (k(7).*a.^6)+(k(6).*a.^5)+(k(5).*a.^4)+(k(4).*a.^3)+(k(3).*a.^2)+(k(2).*a.^1)+(k(1));

%% Calculate Duv
Duv = Lfp - Lbb;

% %% Generate BlackBody
% spdref = blackBodySpectra23Sep05(3000,wavelength);
% 
% xbarr = interp1(wavelength,xbar,wavelength);
% xbarr(isnan(xbarr)) = 0.0;
% ybarr = interp1(wavelength,ybar,wavelength);
% ybarr(isnan(ybarr)) = 0.0;
% zbarr = interp1(wavelength,zbar,wavelength);
% zbarr(isnan(zbarr)) = 0.0;
% 
% Xr = trapz(wavelength,spdref .* xbarr);
% Yr = trapz(wavelength,spdref .* ybarr);
% Zr = trapz(wavelength,spdref .* zbarr);
% 
% ur = (4*Xr) ./ (Xr + 15*Yr + 3*Zr);
% vr = (9*Yr) ./ (Xr + 15*Yr + 3*Zr);
% 
% Duv = sqrt((uk-ur).^2 + (vk-vr).^2);

% %% Minimum Tint
% % Xr = 0.4431;%2700minTint
% % Yr = 0.3806;%2700minTint
% 
% Xr = 0.4212; %3000kMinTint
% Yr = 0.3716; %3000kMinTint
% % Xr = 0.3980; %3500kMinTint
% % Yr = 0.3710; %3500kMinTint
% ur = (4*Xr) ./ (-2*Xr + 12*Yr + 3);
% vr = (9*Yr) ./ (-2*Xr + 12*Yr + 3);
% 
% Duv = sqrt((uk-ur).^2 + (vk-vr).^2);


end