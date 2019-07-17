%% Set Up
clear;
clc;
fileStruct = loadAllTextFiles2();
%% Initialize Loop Variables

ofYtest = 3.2;%1.40;%
%ofYRange =ofYtest-.5:0.5:ofYtest+.5;  
ofYRange = 0.5:0.1:5.0;   % ofY = 1 for original model ON WARM SIDE

ofBtest = 0.7;%0.81;%
%ofBRange =ofBtest:0.5:ofBtest;  
ofBRange = 0.1:0.05:1.5;   % ofB = 1 for original model ON COOL SIDE

rodYtest = 3.5;%1.1;%
%rodYRange =rodYtest-.5:0.5:rodYtest+.5;   
rodYRange = 1.0:0.1:5.0;   % rodY = 0 for original model

rodBtest = 0.65;%1.28;%
%rodBRange =rodBtest:0.5:rodBtest;   
rodBRange = 1.0:0.05:1.5;   % rodB = 0 for original model

mptest = 0.2;       % MPOD
 mpRange =mptest:0.01:mptest;   
%mpRange = 0:0.5:1;   % mp = 0 for original model

matest = 0.35;       % MPOD attenuation
 maRange =matest:0.01:matest;   
%maRange = 0:0.5:1;   % ma = 0 for original model

vdBasetest = 0;
vdBaseRange = 1.0:0.5:5.0;

rsq1Best = 0;
rsq2Best = 0;
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
                        for ivdb = vdBaseRange
            
                          rsqs(1) = CLA_McClung_WL_Lumileds_rod_MPOD_func_Test2(irodY,iOFY,iOFB,irodB,imp,ima,ivdb,fileStruct);
    %                     rsq = CLA_McClung_WL_Lumileds_warm_4K_rod_MPOD_func(irodY,iOFY,iOFB,irodB,imp,ima,fileStruct);                    
    %                     rsq = CLA_McClung_WL_Lumileds_Cool_rod_MPOD_func(irodY,iOFY,iOFB,irodB,imp,ima,fileStruct); % no ofb yet

                           %rsq2 = CLA_FIT_Brainard_Thapan_FEB_2019_rod_both_MPOD_func_Test(irodY,iOFY,iOFB,irodB,imp,ima,fileStruct);
                           rsqs(2) = CLA_FIT_Brainard_Thapan_FEB_2019_rod_both_MPOD_func_Test2(irodY,iOFY,iOFB,irodB,imp,ima,ivdb,fileStruct);

                           rsqs(rsqs < 0) = 0;

                           rsq = sum(rsqs);

                           if rsq > maxrsq
                                maxrsq = rsq
                                rsq1Best = rsqs(1);
                                rsq2Best = rsqs(2);
                                rodYBest = irodY;
                                ofYBest = iOFY;
                                ofBBest = iOFB;
                                rodBBest = irodB;
                                mpBest = imp;
                                maBest = ima;
                                vdBaseBest = ivdb;
                           end
                        end
                    end
                end
            end
        end
    end
end
%max_R2 = round(maxrsq,4)
% ofY___rodY = [ofYBest rodYBest]
% ofB___rodB = [ofBBest rodBBest]
rqsBest = [rsq1Best rsq2Best]
ofY___rodY___ofB___rodB = [ofYBest rodYBest ofBBest rodBBest] % display optimized coefficients
