function [T1map] = EstimateT1fromMP2RAGE1(filename, jsonfileinv1,jsonfileinv2, TRFLASH, NzSlices1, NzSlices2, eff,outputfilename)    
% Example Usage:
% EstimateT1WithJson('/fileserver/fastscratch/prabh/datasets/thresherepilepsy/flywheel20230224/ForEpilepsyMynames/SUBJECTS/4066815/SESSIONS/3T/ACQUISITIONS/t1w/t1w_MP2RAGE_UNI/t1w_MP2RAGE_UNI.nii.gz','/fileserver/fastscratch/prabh/datasets/thresherepilepsy/flywheel20230224/ForEpilepsyMynames/SUBJECTS/4066815/SESSIONS/3T/ACQUISITIONS/t1w/t1w_MP2RAGE_INV1/t1w_MP2RAGE_INV1.json','/fileserver/fastscratch/prabh/datasets/thresherepilepsy/flywheel20230224/ForEpilepsyMynames/SUBJECTS/4066815/SESSIONS/3T/ACQUISITIONS/t1w/t1w_MP2RAGE_INV2/t1w_MP2RAGE_INV2.json',...
% 0.0064,[30 120],1.0)

    jsonString = fileread(jsonfileinv1);
    jsonString2 = fileread(jsonfileinv2);
    jsonData = jsondecode(jsonString);
    jsonData2 = jsondecode(jsonString2);
    MP2RAGE.B0=jsonData.MagneticFieldStrength;
    MP2RAGE.TR=jsonData.RepetitionTime;

    MP2RAGE.TRFLASH=str2num(TRFLASH);%0.0064; % TR of the GRE readout
    MP2RAGE.TIs=[jsonData.InversionTime jsonData2.InversionTime];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
    MP2RAGE.NZslices= [str2num(NzSlices1), str2num(NzSlices2)]; %[ 30 120];% Slices Per Slab * [PartialFourierInSlice-0.5  0.5]
    MP2RAGE.FlipDegrees=[jsonData.FlipAngle jsonData2.FlipAngle];% Flip
    MP2RAGE.filename=filename;%'MP2RAGE_UNI.nii' % file

    % check the properties of this MP2RAGE protocol... this happens to be a
    % very B1 insensitive protocol

    %plotMP2RAGEproperties(MP2RAGE)
%%
    % load the MP2RAGE data - it can be either the SIEMENS one scaled from
    % 0 4095 or the standard -0.5 to 0.5
    MP2RAGEimg=load_untouch_nii(MP2RAGE.filename);
    disp(MP2RAGE)
    disp(size(MP2RAGEimg.img))
    [T1map, R1map]=T1estimateMP2RAGE(MP2RAGEimg,MP2RAGE,str2num(eff));
    T1map.img=1000.*T1map.img;
    save_untouch_nii(T1map,outputfilename);
end
