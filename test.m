i=1;
fileStruct = loadAllTextFiles2;
duv = [];
wavelengths = [
    424,440,456,460,472,480,496,...
%     505,520,530,548,575,600
    ];
for wave = wavelengths
    SPD = [wave-1,0;wave,1;wave+1,0];
    [,duv1(i)] = DKLchromaticity(SPD,fileStruct);
    i=i+1;
end
% vd = .6 + ((4.3 - 0.6)./(1 + (duv./0.06).^2.5));
% duv = calcDuv2([fileStruct.white_light_data.Wavelengths,fileStruct.white_light_data.SPD],fileStruct);

wave = fileStruct2.white_light_data.Wavelengths;
spd = fileStruct2.white_light_data.SPD(:,~fileStruct2.white_light_data.warm);
[~,duv2] = DKLchromaticity([wave,spd],fileStruct);