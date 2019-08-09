function fileStruct = loadAllTextFiles2()
load('old_brainard.mat','old_brainard')
[~,y] = size(old_brainard);
brainard_step = old_brainard(end,1)-old_brainard(end-1,1);
new_wavlengths = old_brainard(end,1)+brainard_step:brainard_step:780;
[~,y2] = size(new_wavlengths);
newArr = zeros(y2,y);
newArr(:,1) =new_wavlengths'; 
old_brainard(end+1:end+y2,:) = newArr;


load('old_thapan.mat','old_thapan')
[~,y] = size(old_thapan);
brainard_step = old_thapan(end,1)-old_thapan(end-1,1);
new_wavlengths = old_thapan(end,1)+brainard_step:brainard_step:780;
[~,y2] = size(new_wavlengths);
newArr = zeros(y2,y);
newArr(:,1) =new_wavlengths'; 
old_thapan(end+1:end+y2,:) = newArr;

load('white_light_data.mat','white_light_data')
[~,y] = size(white_light_data);
brainard_step = white_light_data(end,1)-white_light_data(end-1,1);
new_wavlengths = white_light_data(end,1)+brainard_step:brainard_step:780;
[~,y2] = size(new_wavlengths);
newArr = zeros(y2,y);
newArr(:,1) =new_wavlengths'; 
white_light_data(end+1:end+y2,:) = newArr;

[~, ~, old_thapan_suppressionssupp] = textread('old thapan suppressions.txt', '%f %f %f', 'headerlines', 1);
[~, ~, old_brainard_suppressionssupp] = textread('old brainard suppressions.txt', '%f %f %f', 'headerlines', 1);
white_Light_wave = 380:780;
McClung_1h_2700K = load('McClung_1h_2700K.txt');
spd_McClung_2700K_1lux = load('spd_McClung_2700K_1lux.txt');
supp1 = McClung_1h_2700K(:,1)';
SPD1 = interp1(spd_McClung_2700K_1lux(:,1),spd_McClung_2700K_1lux(:,2),white_Light_wave,'linear',0)';
SPD1 = (SPD1 .* McClung_1h_2700K(:,4)')/Lxy23Sep05([white_Light_wave',SPD1]);

McClung_1h_6500K = load('McClung_1h_6500K.txt');
spd_McClung_6500K_1lux = load('spd_McClung_6500K_1lux.txt');
supp2 = McClung_1h_6500K(:,1)';
SPD2 = interp1(spd_McClung_6500K_1lux(:,1),spd_McClung_6500K_1lux(:,2),white_Light_wave,'linear',0)';
SPD2 = (SPD2 .* McClung_1h_6500K(:,4)')/Lxy23Sep05([white_Light_wave',SPD2]);

WL_1h_2700K_corr = load('WL_1h_2700K_corr.txt');
spd_WL_2700K_1lux = load('spd_WL_2700K_1lux.txt');
supp3 = WL_1h_2700K_corr(:,1)';
SPD3 = interp1(spd_WL_2700K_1lux(:,1),spd_WL_2700K_1lux(:,2),white_Light_wave,'linear',0)';
SPD3 = (SPD3 .* WL_1h_2700K_corr(:,4)')/Lxy23Sep05([white_Light_wave',SPD3]);

WL_1h_5600K_corr = load('WL_1h_5600K_corr.txt');
spd_WL_5600K_1lux = load('spd_WL_5600K_1lux.txt');
supp4 = WL_1h_5600K_corr(:,1)';
SPD4 = interp1(spd_WL_5600K_1lux(:,1),spd_WL_5600K_1lux(:,2),white_Light_wave,'linear',0)';
SPD4 = (SPD4 .* WL_1h_5600K_corr(:,4)')/Lxy23Sep05([white_Light_wave',SPD4]);

Lumileds_1h_3000K = load('Lumileds_1h_3000K.txt');
spd_Lumileds_3000K_1lux = load('spd_Lumileds_3000K_1lux.txt');
supp5 = Lumileds_1h_3000K(:,1)';
SPD5 = interp1(spd_Lumileds_3000K_1lux(:,1),spd_Lumileds_3000K_1lux(:,2),white_Light_wave,'linear',0)';
SPD5 = (SPD5 .* Lumileds_1h_3000K(:,4)')/Lxy23Sep05([white_Light_wave',SPD5]);

Lumileds_1h_CG = load('Lumileds_1h_CG.txt');
spd_Lumileds_CG_1lux = load('spd_Lumileds_CG_1lux.txt');
supp6 = Lumileds_1h_CG(:,1)';
SPD6 = interp1(spd_Lumileds_CG_1lux(:,1),spd_Lumileds_CG_1lux(:,2),white_Light_wave,'linear',0)';
SPD6 = (SPD6 .* Lumileds_1h_CG(:,4)')/Lxy23Sep05([white_Light_wave',SPD6]);

Lumileds_1h_4000K = load('Lumileds_1h_4000K.txt');
spd_Lumileds_4000K_1lux = load('spd_Lumileds_4000K_1lux.txt');
supp7 = Lumileds_1h_4000K(:,1)';
SPD7 = interp1(spd_Lumileds_4000K_1lux(:,1),spd_Lumileds_4000K_1lux(:,2),white_Light_wave,'linear',0)';
SPD7 = (SPD7 .* Lumileds_1h_4000K(:,4)')/Lxy23Sep05([white_Light_wave',SPD7]);
z = zeros(1,18);
A = z;
A(1:4)=true;
B=z;
B(5:8)=true;
C =z;
C(9:10)=true;
D=z;
D(11:12)=true;
E=z;
E(13:14)=true;
F=z;
F(15:16)=true;
G=z;
G(17:18)=true;
supp = horzcat(supp1,supp2,supp3,supp4,supp5,supp6,supp7);
SPD = horzcat(SPD1,SPD2,SPD3,SPD4,SPD5,SPD6,SPD7);

opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [1, 201];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "E05", "E05_1", "E04", "E05_2", "VarName6"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
reaChromaticity = readtable("C:\Users\plummt\Documents\Nagare_CLA\reaChromaticity.csv", opts);
% Convert to output type
reaChromaticity = table2array(reaChromaticity);

fileStruct = struct('Vlamda',load('Vlamda.txt'),...
                    'Vprime',load('Vprime.txt'),...
                    'Scone',load('Scone.txt'),...
                    'CIE31by1',load('CIE31by1.txt'),...
                    'isoTempLines',load('isoTempLinesNewestFine23Sep05.txt'),...
                    'MacularPigmentODfromSnodderly',load('MacularPigmentODfromSnodderly.txt'),...
                    'MelanopsinWlensBy2nm_02Oct2012',load('MelanopsinWlensBy2nm_02Oct2012.txt'),...
                    'TCS', load('Tcs14_23Sep09.txt'),...
                    'reaChromaticity',reaChromaticity,...
                    'old_thapan', struct('Wavelengths',old_thapan(:,1),...
                                         'SPD',old_thapan(:,2:end),...
                                         'Supp',old_thapan_suppressionssupp,...
                                         'warm',~[true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,false,false,false,false,false,false,false,false,false]),...
                    'old_brainard', struct('Wavelengths',old_brainard(:,1),...
                                           'SPD',old_brainard(:,2:end),...
                                           'Supp',old_brainard_suppressionssupp,...
                                           'warm',~[true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,false,false,false,false,false]),...
                    'white_light_data',struct('Wavelengths',white_Light_wave',...
                                           'SPD',SPD,...
                                           'Supp',supp',...
                                           'warm',~[false,false,false,false,true,true,true,true,false,false,true,true,false,false,false,false,false,false],...
                                           'A',A,...
                                           'B',B,...
                                           'C',C,...
                                           'D',D,...
                                           'E',E,...
                                           'F',F,...
                                           'G',G));

end