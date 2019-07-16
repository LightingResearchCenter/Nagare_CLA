function [CLA] = luxspd2claRohan(spd,lux,fileStruct)
%LUXSPD2CLAROHAN Summary of this function goes here
%   Detailed explanation goes here

%% Create scaled SPD
spd_scaled_value = (spd(:,2) .* lux)/Lxy23Sep05(spd);
scaled_spd = [spd(:,1),spd_scaled_value];

%% Calculate CLA
CLA = CLA_postBerlinCorrMelanopsin_06Feb2014(scaled_spd,fileStruct, varargin);

end

