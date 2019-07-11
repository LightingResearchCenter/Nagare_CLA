%% Set Up
clear;
clc;
fileStruct = loadAllTextFiles();
%% Initialize Loop Variables

ofYtest = 1.40;
ofYRange =ofYtest:0.01:ofYtest;  
% ofYRange = 2.5:0.1:3.5;   % ofY = 1 for original model ON WARM SIDE

ofBtest = 0.81;
ofBRange =ofBtest:0.01:ofBtest;  
% ofBRange = 0.5:0.05:1.0;   % ofB = 1 for original model ON COOL SIDE

rodYtest = 1.1;
rodYRange =rodYtest:0.01:rodYtest;   
% rodYRange = 3.5:0.1:4.5;   % rodY = 0 for original model

rodBtest = 1.28;
rodBRange =rodBtest:0.01:rodBtest;   
% rodBRange = 0.5:0.05:1.0;   % rodB = 0 for original model

mptest = 0.2;       % MPOD
% mpRange =mptest:0.01:mptest;   
mpRange = 0:0.1:1;   % mp = 0 for original model

matest = 0.35;       % MPOD attenuation
% maRange =matest:0.01:matest;   
maRange = 0:0.1:1;   % ma = 0 for original model

rodYBest = 0;
rodBBest = 0;
ofYBest = 0;
ofBBest = 0;
mpBest = 0;
maBest = 0;
maxrsq = 0;

%% Loops
for irodY = rodYRange
    for iOFY = ofYRange
        for iOFB = ofBRange
            for irodB = rodBRange
                for imp = mpRange
                    for ima = maRange
            
                    rsq = CLA_McClung_WL_Lumileds_rod_MPOD_func(irodY,iOFY,iOFB,irodB,imp,ima,fileStruct);
%                     rsq = CLA_McClung_WL_Lumileds_warm_4K_rod_MPOD_func(irodY,iOFY,iOFB,irodB,imp,ima,fileStruct);                    
%                     rsq = CLA_McClung_WL_Lumileds_Cool_rod_MPOD_func(irodY,iOFY,iOFB,irodB,imp,ima,fileStruct); % no ofb yet

%                     rsq = CLA_FIT_Brainard_Thapan_FEB_2019_rod_both_MPOD_func(irodY,iOFY,iOFB,irodB,imp,ima,fileStruct); %
% 
                        if rsq > maxrsq
                            maxrsq = rsq
                            rodYBest = irodY;
                            ofYBest = iOFY;
                            ofBBest = iOFB;
                            rodBBest = irodB;
                            mpBest = imp;
                            maBest = ima;
                        end
                    end
                end
            end
        end
    end
end
max_R2 = round(maxrsq,4)
% ofY___rodY = [ofYBest rodYBest]
% ofB___rodB = [ofBBest rodBBest]
% ofY___rodY___ofB___rodB = [ofYBest rodYBest ofBBest rodBBest] % display optimized coefficients
