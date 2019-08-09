%% Set Up
clear;
clc;
fileStruct = loadAllTextFiles2();
%% Initialize Loop Variables


ofYtest = 1.3;%1.40;%
% ofYRange =ofYtest:0.5:ofYtest;
ofYRange = 0.1:.1:2;   % ofY = 1 for original model ON WARM SIDE

ofBtest = 1.4*.35;%0.81;%
% ofBRange =ofBtest:0.5:ofBtest;
ofBRange =0.5:0.1:2;   % ofB = 1 for original model ON COOL SIDE

rodYtest = 1.8;%1.1;%
% rodYRange =rodYtest:0.5:rodYtest;
rodYRange = 0.1:.1:2;   % rodY = 0 for original model

rodBtest = 1.8;%1.28;%
% rodBRange =rodBtest:0.5:rodBtest;
rodBRange = 0.1:0.1:2;   % rodB = 0 for original model

mptest = 0.2;       % MPOD
mpRange =mptest:0.01:mptest;
% mpRange = 0.1:0.05:.6;   % mp = 0 for original model

matest = 0.35;       % MPOD attenuation
maRange =matest:0.01:matest;
% maRange = 0.2:0.05:.4;   % ma = 0 for original model


vdBasetest = 4.3;
vdBaseRange = vdBasetest:0.1:vdBasetest;
% vdBaseRange = 0:0.05:1;

a2test = 0.45;
% a2Range = a2test:0.1:a2test;
a2Range = 0.1:0.05:1;

a3test = 2.55;
% a3Range = a3test:0.1:a3test;
a3Range = 1:0.1:3;


rsq1Best = 0;
rsq2Best = 0;
rsq3Best = 0;
rodYBest = 0;
rodBBest = 0;
ofYBest = 0;
ofBBest = 0;
vdBaseBest = 0;
mpBest = 0;
maBest = 0;
maxrsq = 0;
minABest = 0;
minBBest = 0;
minCBest = 0;

%% Loops
tic
indexRodY = 1;
msgRodY = sprintf('Loop %d of %d',0,length(rodYRange));
f1 = waitbar(0,msgRodY);
indexOfY = 1;
msgOfY = sprintf('Loop %d of %d',0,length(ofYRange));
f2 = waitbar(0,msgOfY);
indexOfB = 1;
msgOfB = sprintf('Loop %d of %d',0,length(ofBRange));
f3 = waitbar(0,msgOfB);
for irodY = rodYRange
    irodY
    toc
    indexRodYRatio = indexRodY/length(rodYRange);
    msgRodY = sprintf('irodY Loop %d of %d',indexRodY,length(rodYRange));
    waitbar(indexRodYRatio,f1,msgRodY);
    indexOfY = 1;
    for iOFY = ofYRange
        iOFY
        indexOfYRatio = indexOfY/length(ofYRange);
        msgOfY = sprintf('iOFY Loop %d of %d',indexOfY,length(ofYRange));
        waitbar(indexOfYRatio,f2,msgOfY);
        indexOfB = 1;
        for iOFB = ofBRange
            %iOFB
            indexOfBRatio = indexOfB/length(ofBRange);
            msgOfB = sprintf('iOFB Loop %d of %d',indexOfB,length(ofBRange));
            waitbar(indexOfBRatio,f3,msgOfB);
            for irodB = rodBRange
                for imp = mpRange
                    for ima = maRange
                        for ivdb = vdBaseRange
                            for ia2 = a2Range
                                for ia3 = a3Range
                                    t= toc;
                                    [rsqs(1),minA(1),minB(1),minC(1)] = CLA_McClung_WL_Lumileds_rod_MPOD_func_Test3(irodY,iOFY,iOFB,irodB,imp,ima,ivdb,fileStruct,ia2,ia3);
                                    %rsq = CLA_McClung_WL_Lumileds_warm_4K_rod_MPOD_func(irodY,iOFY,iOFB,irodB,imp,ima,fileStruct);
                                    %rsq = CLA_McClung_WL_Lumileds_Cool_rod_MPOD_func(irodY,iOFY,iOFB,irodB,imp,ima,fileStruct); % no ofb yet
                                    
%                                     rsq2 = CLA_FIT_Brainard_Thapan_FEB_2019_rod_both_MPOD_func_Test(irodY,iOFY,iOFB,irodB,imp,ima,fileStruct);
                                    [rsqs(2),minA(2),minB(2),minC(2)] = CLA_FIT_Brainard_Thapan_FEB_2019_rod_both_MPOD_func_Test3(irodY,iOFY,iOFB,irodB,imp,ima,ivdb,fileStruct,ia2,ia3);
                                    %rsqs(3) = generateMonochromaticSpectralResponseOfModel_efficacy_Func(irodY, iOFY, iOFB, irodB, imp, ima,ivdb,fileStruct);
%                                     rsqs(3) = generateMonochromaticSpectralResponseOfModel_Func_Test3(irodY, iOFY, iOFB, irodB, imp, ima,ivdb,fileStruct,true,ia2,ia3);
                                    
                                    rsqs(isnan(rsqs) | rsqs < 0) = 0;
                                    
                                    rsq = sum(rsqs);
                                    
                                    if rsq > maxrsq
                                        maxrsq = rsq
%                                         generateMonochromaticSpectralResponseOfModel_Func_Test3(irodY, iOFY, iOFB, irodB, imp, ima,ivdb,fileStruct,true,ia2,ia3);
%                                         CLA_vs_BlackBodyCCT_func(irodY, iOFY, iOFB, irodB, imp, ima,ivdb,fileStruct,true,ia2,ia3);
                                        rsq1Best = rsqs(1);
%                                         rsq2Best = rsqs(2);
%                                         rsq3Best = rsqs(3);
                                        rodYBest = irodY;
                                        ofYBest = iOFY;
                                        ofBBest = iOFB;
                                        rodBBest = irodB;
                                        mpBest = imp;
                                        maBest = ima;
                                        vdBaseBest = ivdb;
                                        a2Best = ia2;
                                        a3Best = ia3;
                                        minABest = minA;
                                        minBBest = minB;
                                        minCBest = minC;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            indexOfB = indexOfB +1;
        end
        indexOfY = indexOfY+1;
    end
    indexRodY = indexRodY+1;
end
%max_R2 = round(maxrsq,4)
% ofY___rodY = [ofYBest rodYBest]
% ofB___rodB = [ofBBest rodBBest]
rqsBest = [rsq1Best rsq2Best rsq3Best]
ofY___rodY___ofB___rodB___VDBase = [ofYBest rodYBest ofBBest rodBBest vdBaseBest] % display optimized coefficients
vd___a2___a3 = [vdBaseBest,a2Best,a3Best]
% mp___ma = [mpBest,maBest]