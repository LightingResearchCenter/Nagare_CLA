%% Set Up
clear;
clc;
fileStruct = loadAllTextFiles2();
%% Initialize Loop Variables


ofYtest = 1.5;%1.40;%
ofYRange =ofYtest:0.5:ofYtest;
%ofYRange = 1:.1:3;   % ofY = 1 for original model ON WARM SIDE

ofBtest = 0.81;%0.81;%
ofBRange =ofBtest:0.5:ofBtest;
%ofBRange = 0.75:0.5:1.25;   % ofB = 1 for original model ON COOL SIDE

rodYtest = 1.55;%1.1;%
rodYRange =rodYtest:0.5:rodYtest;
%rodYRange = 1:.1:3;   % rodY = 0 for original model

rodBtest = 1.28;%1.28;%
rodBRange =rodBtest:0.5:rodBtest;
%rodBRange = 0.25:0.5:1;   % rodB = 0 for original model

mptest = 0.2;       % MPOD
mpRange =mptest:0.01:mptest;
%mpRange = 0:0.5:1;   % mp = 0 for original model

matest = 0.35;       % MPOD attenuation
maRange =matest:0.01:matest;
%maRange = 0:0.5:1;   % ma = 0 for original model


vdBasetest = 4.85;
vdBaseRange = vdBasetest:0.1:vdBasetest;
%vdBaseRange = 3:0.01:5;

rsq1Best = 0;
rsq2Best = 0;
rsq3Best = 0;
rodYBest = 0;
rodBBest = 0;
ofYBest = 0;
ofBBest = 0;
mpBest = 0;
maBest = 0;
maxrsq = 0;

%% Loops
tic
for irodY = rodYRange
    irodY
    for iOFY = ofYRange
        iOFY
        for iOFB = ofBRange
            iOFB
            for irodB = rodBRange
                for imp = mpRange
                    for ima = maRange
                        for ivdb = vdBaseRange
                            t= toc;
                            rsqs(1) = CLA_McClung_WL_Lumileds_rod_MPOD_func_Test3(irodY,iOFY,iOFB,irodB,imp,ima,ivdb,fileStruct);
                            %rsq = CLA_McClung_WL_Lumileds_warm_4K_rod_MPOD_func(irodY,iOFY,iOFB,irodB,imp,ima,fileStruct);
                            %rsq = CLA_McClung_WL_Lumileds_Cool_rod_MPOD_func(irodY,iOFY,iOFB,irodB,imp,ima,fileStruct); % no ofb yet
                            
                            %rsq2 = CLA_FIT_Brainard_Thapan_FEB_2019_rod_both_MPOD_func_Test(irodY,iOFY,iOFB,irodB,imp,ima,fileStruct);
                            rsqs(2) = CLA_FIT_Brainard_Thapan_FEB_2019_rod_both_MPOD_func_Test3(irodY,iOFY,iOFB,irodB,imp,ima,ivdb,fileStruct);
                            %rsqs(3) = generateMonochromaticSpectralResponseOfModel_efficacy_Func(irodY, iOFY, iOFB, irodB, imp, ima,ivdb,fileStruct);
                            rsqs(3) = 0;%generateMonochromaticSpectralResponseOfModel_Func_Test3(irodY, iOFY, iOFB, irodB, imp, ima,ivdb,fileStruct,false);
                            
                            rsqs(isnan(rsqs) | rsqs < 0) = 0;
                            
                            rsq = sum(rsqs);
                            
                            if rsq > maxrsq
                                maxrsq = rsq
                                generateMonochromaticSpectralResponseOfModel_Func_Test3(irodY, iOFY, iOFB, irodB, imp, ima,ivdb,fileStruct,true);
                                rsq1Best = rsqs(1);
                                rsq2Best = rsqs(2);
                                rsq3Best = rsqs(3);
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
rqsBest = [rsq1Best rsq2Best rsq3Best]
ofY___rodY___ofB___rodB___VDBase = [ofYBest rodYBest ofBBest rodBBest vdBaseBest] % display optimized coefficients
