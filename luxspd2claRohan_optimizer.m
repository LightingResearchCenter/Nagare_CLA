function [CLA] = luxspd2claRohan_optimizer(spd, lux, rodY, ofY, ofB, rodB, mp, ma,fileStruct)
%LUXSPD2CLAROHAN Summary of this function goes here
%   Detailed explanation goes here

%% Lux loop
for i1 = 1:numel(lux)
    %% Create scaled SPD
    spd_scaled_value = (spd(:,2) .* lux(i1))/Lxy23Sep05(spd);
    scaled_spd = [spd(:,1),spd_scaled_value];

    %% Calculate CLA
    CLA(i1) = CLA_rod_both_MPOD_calculation_Test(scaled_spd, rodY, ofY, ofB, rodB, mp, ma,fileStruct);
end

end

