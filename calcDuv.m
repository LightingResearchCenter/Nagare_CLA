function [Duv] = calcDuv(spd, fileStruct)

wavelength_spd = spd(:,1);
spd = spd(:,2:end);
%% Calc chromaticty 
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

Xk = trapz(wavelength_spd,spd .* xbar);
Yk = trapz(wavelength_spd,spd .* ybar);
Zk = trapz(wavelength_spd,spd .* zbar);

uk = (4*Xk) ./ (Xk + 15*Yk + 3*Zk);
vk = (6*Yk) ./ (Xk + 15*Yk + 3*Zk);
%% Calc CCT
Tc = CCT23Sep05_PreLoaded([wavelength_spd,spd],fileStruct);
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

Xr = trapz(wavelength_spd,spdref .* xbar);
Yr = trapz(wavelength_spd,spdref .* ybar);
Zr = trapz(wavelength_spd,spdref .* zbar);

ur = (4*Xr) ./ (Xr + 15*Yr + 3*Zr);
vr = (6*Yr) ./ (Xr + 15*Yr + 3*Zr);

Duv = sqrt((uk-ur).^2 +((2/3)*vk-(2/3)*vr).^2);
end