addpath(genpath('.'))

subj=['5193338','5721998','5834888','5077519','4066815'];
for sub_no=1:5
addpath(genpath('/fileserver/fastscratch/prabh/work/Utilities/'))
subjectpath=strcat('/fileserver/fastscratch/prabh/datasets/thresherepilepsy/flywheel20230224/ForEpilepsyMynames/SUBJECTS/',num2str(subj(sub_no)),'/SESSIONS/3T/ANALYSIS/CRLPIPELINE/');
fs = '3T';
filename=strcat(subjectpath,'common-processed/anatomical/01-t1w-ref/t1w_ref.nii.gz');


t1 = load_untouch_nii(filename); 
MP2RAGEINV2img=load_untouch_nii('/fileserver/fastscratch/prabh/datasets/thresherepilepsy/flywheel20230224/ForEpilepsyMynames/SUBJECTS/4066815/SESSIONS/3T/ACQUISITIONS/t1w/t1w_MP2RAGE_INV2/orientedINV2.nii.gz');

% WHITE MATTER
parcwm     =   load_untouch_nii(strcat(subjectpath,'common-processed/modules/Parcellation/NVM/ParcellationNVM_WM.nii.gz')); parcwm.img=double(parcwm ...
            .img);

% GRAY MATTER
parcgm     =   load_untouch_nii(strcat(subjectpath,'common-processed/modules/Parcellation/NVM/ParcellationNVM_GM.nii.gz')); parcgm.img=double(parcgm ...
    .img);

if strcmp(fs,'7T')
    MP2RAGE.B0=7;           % in Tesla
    MP2RAGE.TR=4.57;           % MP2RAGE TR in seconds
    MP2RAGE.TRFLASH=0.0064; % TR of the GRE readout
    MP2RAGE.TIs=[840e-3 2370e-3];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
    MP2RAGE.NZslices=[48 96];% Slices Per Slab * [PartialFourierInSlice-0.5  0.5]
    MP2RAGE.FlipDegrees=[5 6];% Flip
    % angle of the two readouts in degrees
    MP2RAGE.filename=filename;%'MP2RAGE_UNI.nii' % file
else
    if strcmp(fs,'3T')
        MP2RAGE.B0=3;           % in Tesla
        MP2RAGE.TR=5;           % MP2RAGE TR in seconds
        MP2RAGE.TRFLASH=7.14e-3; % TR of the GRE readout
        MP2RAGE.TIs=[700e-3 2500e-3];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
        MP2RAGE.NZslices=[88 88];% Slices Per Slab * [PartialFourierInSlice-0.5  0.5]
        MP2RAGE.FlipDegrees=[4 5];% Flip
        % angle of the two readouts in degrees
        MP2RAGE.filename=filename;%'MP2RAGE_UNI.nii' % file
    end
end

T1map=EstimateT1func('3T',filename);
% [T1map , M0map , R1map]=T1M0estimateMP2RAGE(t1,MP2RAGEINV2img,MP2RAGE,0.96);

%%
wm= load_untouch_nii(strcat(subjectpath,'common-processed/modules/Parcellation/NVM/ParcellationNVM_WM.nii.gz'));
gm= load_untouch_nii(strcat(subjectpath,'common-processed/modules/Parcellation/NVM/ParcellationNVM_GM.nii.gz'));

%% Draw the expected signal intensities of tissues with T1= 800ms  and 1200ms, for different TIs. 

Intensitis = zeros(100,1);
cnt=1;
for tis=200:10:1000
    MP2RAGE.TIs= [tis*0.001 2.5] ; % first TI in MP2RAGE would be the our chosen TI like in 3D-EDGE acquisition design
    [Intensity, T1vector, ibc] = MP2RAGE_lookuptable(2,MP2RAGE.TR,MP2RAGE.TIs,MP2RAGE.FlipDegrees,MP2RAGE.NZslices,MP2RAGE.TRFLASH,'normal',0.96);
    % ibc is intensity before combining. :)  ibc(:,1) is for first TI, and
    % ibc(:,2) is for second TI.
    if not(any(ibc))
        tis % sometimes the ibc is null for given ti, like for very initial tis. 
    else
        cnt
        wmi(cnt)=ibc(17,1); % assuming t1 of wm = T1vector(17) which is 900
        gmi(cnt)=ibc(25,1); % assuming t1 of gm = T1vector(25) which is 1300
        tiss(cnt)=tis;
        % Generate the TI images
        TIimages(:,cnt)=interp1(T1vector,ibc(:,1),T1map.img(:));
        wmavg(cnt)=mean(wm.img(:).*TIimages(:,cnt));
        gmavg(cnt)=mean(gm.img(:).*TIimages(:,cnt));
        cnt=cnt+1;
    end
end

%% 
figure; plot(tiss,wmavg)
hold on
plot(tiss,gmavg)
grid on
tis
figure; plot(tiss,wmi)
hold on
plot(tiss,gmi)
grid on
figure; plot(tiss,abs(wmavg)-abs(gmavg))
grid on
[~,optimalTI]=min(abs(abs(wmavg)-abs(gmavg)));
optimalTI=tiss(optimalTI);
% optimal TI = 650 ms

%%
Intensitis = zeros(100,1);
cnt=1;
for tis=optimalTI-5:optimalTI+5%630:1:670
    MP2RAGE.TIs= [tis*0.001 2.5] ; % first TI in MP2RAGE would be the our chosen TI like in 3D-EDGE acquisition design
    [Intensity, T1vector, ibc] = MP2RAGE_lookuptable(2,MP2RAGE.TR,MP2RAGE.TIs,MP2RAGE.FlipDegrees,MP2RAGE.NZslices,MP2RAGE.TRFLASH,'normal',0.96);
    % ibc is intensity before combining. :)  ibc(:,1) is for first TI, and
    % ibc(:,2) is for second TI.
    if not(any(ibc))
        tis % sometimes the ibc is null for given ti, like for very initial tis. 
    else
        cnt
        wmi(cnt)=ibc(17,1); % assuming t1 of wm = T1vector(17) which is 900
        gmi(cnt)=ibc(25,1); % assuming t1 of gm = T1vector(25) which is 1300
        tiss(cnt)=tis;
        % Generate the TI images
        TIimages(:,cnt)=interp1(T1vector,ibc(:,1),T1map.img(:)); % SEe this is my assumption okay? If iterp can be used for combined intensity like this, can it be used in similar wa for uncombined single ti signal intensity too?
        TIimagesCombined(:,cnt)=interp1(T1vector,Intensity,T1map.img(:));
        wmavg(cnt)=mean(wm.img(:).*TIimages(:,cnt));
        gmavg(cnt)=mean(gm.img(:).*TIimages(:,cnt));
        cnt=cnt+1;
    end
end
mask=wm.img | gm.img;

%%  Middle like image
TIimagesCombined=mean(TIimagesCombined,2);
ao=(reshape(TIimagesCombined(:,1),[176,240,256]));
aaa=T1map; aaa.img=ao;
save_untouch_nii(aaa,strcat('edge/EDGE_',num2str(subj(sub_no)),'.nii.gz'))
end
% figure; imshow3Dfull(abs(ao),[0 0.35])
% figure; imshow3Dfull((ao),[-0.4 0.5])
% %% histogram of t1 values for edge'ed pixels
% th=0.05.*max(ao(:))
% bb=T1map.img(abs(ao)<th);
% figure; hist(bb)
% %% edge like image obtained via thresholding t1s
% tmin=min(bb(:));
% tmax=max(bb(:));
% ao2=ao;
% ao2(abs(ao2)<th)=0;
% figure; imshow(abs(ao2(:,:,144)),[0 0.02])
% figure; imshow3Dfull(abs(ao2(:,:,:)),[0 0.02])