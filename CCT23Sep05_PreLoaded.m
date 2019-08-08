function Tc = CCT23Sep05_PreLoaded(spd,fileStruct)
% CCT (Correlated Color Temperature)
% Calculates the correlated color temperature (CCT) of spd.
% Uses method by Allen Robertson as references in W&S (1982).
% Needs file 'CIE31_1.mat' for CIE 1931 Std.Observer xbar,ybar and zbar data.
% Needs file 'isoTempLines.mat' for isotemperature lines from W&S (1982) or similarly calculated.
% Interpolates xbar, ybar, and zbar values to same increment as spd data.
% Function arguments are:
%	spd - spectral power distribution vector,
%	startw - starting wavelength of spd in nanometers,
%	endw - ending wavelength of spd in nanometers,
%	incrementw - increment of wavelength data in nanometers.
% OR needs spd file with 2 columns [wavelength value] in any increment, regular or not.
% When used with a 2-column spd arguement, startw, endw and increment arguements ar not used, but
% place holder values, such as zeros, must be specified. e.g., Tc = CCT_1(spd,0,0,0)

wavelength_spd = spd(:,1);
spd = spd(:,2:end);

%load('CIE31_1', 'wavelength','xbar','ybar','zbar');
Table = fileStruct.CIE31by1;
wavelength = Table(:,1);
xbar = Table(:,2);
ybar = Table(:,3);
zbar = Table(:,4);
%load('isoTempLinesNewestFine.mat','T','ut','vt','tt');
% Table = load('isoTempLinesNewestFine23Sep05.txt');
Table = fileStruct.isoTempLines;
T = Table(:,1);
ut = Table(:,2);
vt = Table(:,3);
tt = Table(:,4);

xbar = interp1(wavelength,xbar,wavelength_spd);
xbar(isnan(xbar)) = 0.0;
ybar = interp1(wavelength,ybar,wavelength_spd);
ybar(isnan(ybar)) = 0.0;
zbar = interp1(wavelength,zbar,wavelength_spd);
zbar(isnan(zbar)) = 0.0;

% Calculate Chromaticity Coordinates
%diffwave = diff(wavelength_spd);
%deltaWave = [diffwave(1)/2;(diffwave(1:end-1)+diffwave(2:end))/2;diffwave(end)/2];
%X = sum(spd .* xbar.*deltaWave);
%Y = sum(spd .* ybar.*deltaWave);
%Z = sum(spd .* zbar.*deltaWave);

X = trapz(wavelength_spd,spd.*xbar);
Y = trapz(wavelength_spd,spd.*ybar);
Z = trapz(wavelength_spd,spd.*zbar);
x = X./(X+Y+Z);
y = Y./(X+Y+Z);
u = 4.*x./(-2.*x+12.*y+3);
v = 6.*y./(-2.*x+12.*y+3);

% Find adjacent lines to (us,vs)
[j,k] = size(spd);
Tc = zeros(k,1);
for i1= 1:k
    n = length(T);
    index = 0;
    d1 = ((v(i1)-vt(1)) - tt(1).*(u(i1)-ut(1)))./sqrt(1+tt(1)*tt(1));
    for i=2:n
        d2 = ((v(i1)-vt(i)) - tt(i).*(u(i1)-ut(i)))/sqrt(1+tt(i).*tt(i));
        if (d1/d2 < 0)
            index = i;
            break;
        else
            d1 = d2;
        end
    end
    if index == 0
        Tc(i1) = -1; % Not able to calculate CCT, u,v oordinates outside range.
    else
        % Calculate CCT by interpolation between isotemperature lines
        Tc = 1/(1/T(index-1)+d1/(d1-d2)*(1/T(index)-1/T(index-1)));
    end
end

% Calculate CCT by interpolation between isotemperature lines

