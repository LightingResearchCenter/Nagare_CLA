function [DC,sstar] = CRI23Sep05_PreLoaded(spd,fileStruct,varargin)
% CRI Color Rendering Indices (1-8)
% Calculates the Color Rendering Indices (CRI) for the first 8 CIE
% test color samples according to CIE 13.3 (1995).
% Interpolates xbar, ybar, and zbar values to same increment as spd data. 
% The following files are needed:
%	CIE31_1.mat, CIEDaySn.mat, CCT_1.m, isoTempLines.mat, CIEtcsa.mat
% Function arguments are:
%	spd - spectral power distribution vector,
%	startw - starting wavelength of spd in nanometers,
%	endw - ending wavelength of spd in nanometers,
%	incrementw - increment of wavelength data in nanometers.
%  OR
%   spd - spectral power distribution 2-column matrix
%         column 1 = wavelength in nm (can be any arbitrary spacing)
%         column 2 = intensity values
%   put in place holders for startw, endw and incrementw, e.g.  [R,DC,Tc] = CRI_1(spd_matrix,0,0,0)
% Function returns:
%   R = vector of the color rendering indicies for each of the eight test color samples
%   DC = scalar tolorance value having to do with the appropriateness of applying CRI to colored light sources (see CIE documentation)
%   Tc = The correlated color temperature of the light source in Kelvin units

% Calculate Correlated Color Temperature, Tc.
wavelength_spd = spd(:,1);
spd = spd(:,2:end);
Tc = CCT23Sep05_PreLoaded([wavelength_spd,spd],fileStruct);

%load('CIE31_1', 'wavelength','xbar','ybar','zbar');
Table = fileStruct.CIE31by1;
wavelength = Table(:,1);
xbar = Table(:,2);
ybar = Table(:,3);
zbar = Table(:,4);

xbar = interp1(wavelength,xbar,wavelength_spd);
xbar(isnan(xbar)) = 0.0;
ybar = interp1(wavelength,ybar,wavelength_spd);
ybar(isnan(ybar)) = 0.0;
zbar = interp1(wavelength,zbar,wavelength_spd);
zbar(isnan(zbar)) = 0.0;

% Calculate Reference Source Spectrum, spdref.
for i = 1:length(Tc)
    if (Tc < 5000)
        spdref(:,i) = blackBodySpectra23Sep05(Tc(i),wavelength_spd);
    else
        if (Tc <= 25000)
            spdref(:,i) = CieDaySpectra23Sep05(Tc(i),wavelength_spd);
        else
            %error('CCT above 25,000 K');
            spdref(:,i) = nan(size(wavelength_spd));
        end
    end
end
% TCS = load('Tcs14_23Sep09.txt'); % 'Tcs14.txt'
TCS = fileStruct.TCS;
%Interpolate TCS values from 5 nm to spd nm increments
TCS_1 = zeros(length(wavelength_spd),14);
wavelength_5 = TCS(:,1);
TCS = TCS(:,2:end); % Remove wavelength column
TCS = TCS/1000;
for i = 1:14
	TCS_1(:,i) = interp1(wavelength_5,TCS(:,i),wavelength_spd,'linear',0);
	%TCS_1(isnan(TCS_1(:,i)),i) = 0.0; % remove NaN from vector, but it
	%should't have any since a zero value was specified for extrapolated
	%values in above interp1 function.
end

% Calculate u, v chromaticity coordinates of samples under test illuminant, uk, vk and
% reference illuminant, ur, vr.
uki = zeros(1,14);
vki = zeros(1,14);
uri = zeros(1,14);
vri = zeros(1,14);
Xk = trapz(wavelength_spd,spd .* xbar);
Yk = trapz(wavelength_spd,spd .* ybar);
Zk = trapz(wavelength_spd,spd .* zbar);
Yknormal = 100 ./ Yk;
Yk = Yk.*Yknormal;
uk = 4.*Xk./(Xk+15.*Yk+3.*Zk);
vk = 9.*Yk./(Xk+15.*Yk+3.*Zk);

Xr = trapz(wavelength_spd,spdref .* xbar);
Yr = trapz(wavelength_spd,spdref .* ybar);
Zr = trapz(wavelength_spd,spdref .* zbar);
Yrnormal = 100 ./ Yr;
Yr = Yr.*Yrnormal;
ur = 4.*Xr./(Xr+15.*Yr+3.*Zr);
vr = 9.*Yr./(Xr+15.*Yr+3.*Zr);

for i = 1:length(Yr)
    if Yr(i)./Yk(i) <=((6/26).^3)
        LStar(i) = ((29/3).^3).*(Yr(i)./Yk(i));
    else
        LStar(i) = 116.*(Yr(i)./Yk(i)).^(1/3)-16;
    end
end
Ustar = 13.*LStar.*(uk-ur);
Vstar = 13.*LStar.*(vk-vr);
Cstar = sqrt(Ustar.^2+Vstar.^2);
sstar = Cstar./LStar;
% Check tolorence for reference illuminant
DC = sqrt((uk-ur).^2 + ((vk)-(vr)).^2);

% % Apply adaptive (perceived) color shift.
% ck = (4 - uk - 10.*vk) ./ vk;
% dk = (1.708.*vk + 0.404 - 1.481.*uk) ./ vk;
% cr = (4 - ur - 10.*vr) ./ vr;
% dr = (1.708.*vr + 0.404 - 1.481.*ur) ./ vr;
% 
% for i = 1:14
% 	cki = (4 - uki(i) - 10.*vki(i)) ./ vki(i);
% 	dki = (1.708.*vki(i) + 0.404 - 1.481.*uki(i)) ./ vki(i);
% 	ukip(:,i) = (10.872 + 0.404.*cr./ck.*cki - 4.*dr./dk.*dki) ./ (16.518 + 1.481.*cr./ck.*cki - dr./dk.*dki);
% 	vkip(:,i) = 5.520 ./ (16.518 + 1.481.*cr./ck.*cki - dr./dk.*dki);
% end
% 
% %  Transformation into 1964 Uniform space coordinates.
% for i = 1:14
% 	Wstarr(i,:) = 25.*Yri(:,i).^.333333 - 17;
% 	Ustarr(i,:) = 13.*Wstarr(:,i).*(uri(i) - ur);
% 	Vstarr(i,:) = 13.*Wstarr(:,i).*(vri(i) - vr);
% 	
% 	Wstark(i,:) = 25.*Yki(:,i).^.333333 - 17;
% 	Ustark(i,:) = 13.*Wstark(:,i).*(ukip(:,i) - ur); % after applying the adaptive color shift, u'k = ur
% 	Vstark(i,:) = 13.*Wstark(:,i).*(vkip(i,:) - vr); % after applying the adaptive color shift, v'k = vr
% end
% 
% % Determination of resultant color shift, delta E.
% deltaE = zeros(1,14);
% R = zeros(1,14);
% for i = 1:14
% 	deltaE(i) = sqrt((Ustarr(i) - Ustark(i)).^2 + (Vstarr(i) - Vstark(i)).^2 + (Wstarr(i) - Wstark(i)).^2);
% 	R(i) = 100 - 4.6.*deltaE(i);
% end
% Ra = sum(R(1:8))/8;
% end