function [RGSignal,BYSignal] = DKLchromaticity(SPD,fileStruct)

wave = fileStruct.reaChromaticity(:,1);
L10 = fileStruct.reaChromaticity(:,2);
M10 = fileStruct.reaChromaticity(:,3);
S10 = fileStruct.reaChromaticity(:,4);
vLambda  = fileStruct.reaChromaticity(:,5);


wavelengths = SPD(:,1);
SPD = SPD(:,2:end);

L10 = interp1(wave,L10,wavelengths);
M10 = interp1(wave,M10,wavelengths);
S10 = interp1(wave,S10,wavelengths);
vLambda = interp1(wave,vLambda,wavelengths);

N = 1./(trapz(wavelengths,vLambda.*SPD));

RGWeights = (0.28.*L10) - (0.37.*M10);
RGSignal = N.* trapz(wavelengths,RGWeights.*SPD);

BYWeights = (4.00.*S10)-((0.57.*L10)+(0.8.*M10));
BYSignal = N.* trapz(wavelengths,BYWeights.*SPD);
end