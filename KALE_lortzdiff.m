%% the data prepare
% create by sunyaozong
% 2023,2,20
clear all;close all;clc
addpath('tools')
datapath = fullfile(pwd,'\data');
Savepath = fullfile(datapath,'test');
mkdir(Savepath)
scanNum = 8;
load(fullfile(datapath,['V_exp Scan',num2str(scanNum),'.mat']));
load(fullfile(datapath,['S0_Scan',num2str(scanNum),'E1_power0.7.mat']));
load(fullfile(datapath,['w_offset Scan',num2str(scanNum),'.mat']));
load(fullfile(datapath,['brainMask.mat']));
[V_exp_mask] = prepare(S0,V_exp,brainMask);
slice = 0;
Zspec_AddNoise = V_exp_mask;
w_offset_inter=[min(w_offset):0.1:max(w_offset)]'; 
clear V_norm_B0_0p7uT M0_stack Thmask_0p7uT V_exp
%%
K_value = 3;
iteration = 0;
method_fitting = 1; % 0 : class fitting     1 : after kmeans and voxels fitting
All_numClass = 1; 
medfilter = 1;
std_range = 76;    % Ts     
num_range = 12;    % Yn
noise_std = 0;
method = 'kmeanspp'; 
offset_choice_name = 'All_AddVoxels';
offset_choice_kmeans = w_offset;
filename = [offset_choice_name,'Slice',num2str(slice),'Kvalue',num2str(K_value),'interpolation'];
newSavepath = fullfile(Savepath,filename);
if exist(newSavepath) ~=7       
   mkdir(newSavepath) 
end
[Row,Column,numOffset] = size(Zspec_AddNoise);
oriindex = [1:1:Row*Column]';
Zspec_loli_vec = reshape(Zspec_AddNoise,Row*Column,numOffset);
Zspec_loli_vec_ori = Zspec_loli_vec;

cls_loli = zeros(Row*Column,1);
KALE_Ites = zeros(Row*Column,1);
if method_fitting == 0
    final_lortz_LDgroupk = zeros(Row*Column,length(w_offset_inter));
    final_lortzdiff_LDgroupk = zeros(Row*Column,length(w_offset_inter));
else
    final_lortz_LDk = zeros(Row*Column,length(w_offset_inter));
    final_lortzdiff_LDk = zeros(Row*Column,length(w_offset_inter));
end
%%
numCluster = 0;
ori_std = std2(Zspec_loli_vec)/std_range;   
iteration_data{1,1} = {Zspec_loli_vec};
iteration_data{1,2} = {oriindex};
par_0 = {[0,  50, 20, 1]};
tic
while 1
    num = 1;
    for numClass = 1 : size(iteration_data,1)
        [cls,~,ture_index,index_choice] = Clustering_KALE_fitting(iteration_data{numClass,1},...
            K_value,method,w_offset,offset_choice_kmeans,iteration_data{numClass,2});

        [ZsData_exp,w_offset1,mSigEStd] = ROIprocess_Zspecforkmeans(ture_index,w_offset,Zspec_AddNoise);

        [Lorentz_fit_Kmeans,Lor_fit,Lorentzian_diff_Kmeans,~,par_fit,~] =...
            LF_Poly(w_offset1,ZsData_exp,ture_index,-2,-6.25,-10,2,6.25,10,par_0{numClass});
        
        Mean_mSigEStd =  mean(mSigEStd(index_choice,:));
        for ii = 1 : length(ture_index)
            if (Mean_mSigEStd(ii) < ori_std)  || (length(ture_index{ii}) < (floor(Row/num_range))^2) 
                index = ture_index{ii};
                numCluster = numCluster + 1;
%%%%%%%%%%%%%%%%%%%%%% KALE-LD group-wise fitting %%%%%%%%%%%%%%%%%%%%%%%%
                if method_fitting == 0
                    for j = 1 : length(index)
                        final_lortz_LDgroupk(index(j),:) = Lorentz_fit_Kmeans(:,ii)';
                        final_lortzdiff_LDgroupk(index(j),:) = Lorentzian_diff_Kmeans(:,ii)';
                    end
%%%%%%%%%%%%%%%%%%    KALE-LD voxels-wise fitting   %%%%%%%%%%%%%%%%%%%
                else
                    if mSigEStd(ii) == 0
                        for j = 1 : length(index)
                            final_lortz_LDk(index(j),:) = Lorentz_fit_Kmeans(:,ii)';
                        end
                    else
                        for i = 1 : length(index)
                            ind = index(i);
                            x = Zspec_loli_vec(ind,:);
                            [Lorentzian_fit, Lorentzian_difference, ~,~,Iters]...
                                    = spectrum_LD_MTRrex_interpolation(x', w_offset,par_fit(ii,:),iteration);
                            final_lortz_LDk(ind,:) = Lorentzian_fit;
                            final_lortzdiff_LDk(ind,:) = Lorentzian_difference;
                            KALE_Ites(ind) = Iters;
                        end
                    end
                end
            else
                index = ture_index{ii};
                Zsdata_class_lolivec{num,1} = Zspec_loli_vec(index,:);
                par_0{num} = par_fit(ii,:);
                Zsdata_class_lolivec{num,2} = index;
                num = num + 1;
            end
        end
        All_numClass = All_numClass + 1;
    end
if exist('Zsdata_class_lolivec')
    iteration_data = Zsdata_class_lolivec;
else
    break
end
    clear Zsdata_class_lolivec
    iteration = iteration + 1
end
toc

%% adjust the LD maps and display it.
APT_DisplayRange_LD = [0.04,0.1];  
NOE_DisplayRange_LD = [0.06,0.16];
[~,index_APT] = min(abs(w_offset_inter-3.5));
[~,index_NOE] = min(abs(w_offset_inter+3.5));
Display_Range_LF = [0.6,1];
Zspec_AddNoise = reshape(Zspec_loli_vec,Row,Column,numOffset);
if method_fitting == 0
    final_lortz_LDgroupk = reshape(final_lortz_LDgroupk,Row,Column,length(w_offset_inter));
    final_lortzdiff_LDgroupk = reshape(final_lortzdiff_LDgroupk,Row,Column,length(w_offset_inter));
    if medfilter
        for numoffset = 1 : length(w_offset_inter)
            final_lortzdiff_LDgroupk(:,:,numoffset) = medfilt2(final_lortzdiff_LDgroupk(:,:,numoffset),[2,2]);
        end
    end
    display_image(final_lortzdiff_LDgroupk(:,:,index_APT),brainMask,['K class LD APT K=',num2str(K_value)],APT_DisplayRange_LD)
    display_image(final_lortzdiff_LDgroupk(:,:,index_NOE),brainMask,['K class LD NOE K=',num2str(K_value)],NOE_DisplayRange_LD)

    save(fullfile(newSavepath,['LDgroupk','.mat']),'final_lortz_LDgroupk','final_lortzdiff_LDgroupk')
else
    final_lortz_LDk = reshape(final_lortz_LDk,Row,Column,length(w_offset_inter));
    final_lortzdiff_LDk = reshape(final_lortzdiff_LDk,Row,Column,length(w_offset_inter));
    if medfilter
        for numoffset = 1 : length(w_offset_inter)
            final_lortzdiff_LDk(:,:,numoffset) = medfilt2(final_lortzdiff_LDk(:,:,numoffset),[2,2]);
        end
    end
    display_image(final_lortzdiff_LDk(:,:,index_APT),brainMask,'K Voxels LD APT',APT_DisplayRange_LD)
    display_image(final_lortzdiff_LDk(:,:,index_NOE),brainMask,'K Voxels LD NOE',NOE_DisplayRange_LD)
    save(fullfile(newSavepath,['LDk','.mat']),'final_lortz_LDk','final_lortzdiff_LDk')
end
%% traditional method for voxles by voxels fitting
Zspec_AddNoise_reloli = reshape(Zspec_loli_vec,Row,Column,numOffset);
Zspec_AddNoise_reloli = Zspec_AddNoise_reloli.*brainMask;
Zspec_AddNoise = Zspec_AddNoise.*brainMask;
tic
[~,Voxel_Lorentzian_difference,Voxel_Lorentzian_fit,freq_offsets,numIters]...
    =Voxel_LD_useful_pixels(Zspec_AddNoise,w_offset);
time_voxels = toc
%% display the traditional method LD mmps and save it.
medfilter = 1;
if medfilter
    for numoffset = 1 : length(freq_offsets)
        Voxel_Lorentzian_difference(:,:,numoffset) = medfilt2(Voxel_Lorentzian_difference(:,:,numoffset),[2,2]);
    end
end
[~,index_APTinterpolation] = min(abs(w_offset_inter-3.5));
[~,index_NOEinterpolation] = min(abs(w_offset_inter+3.5));
save(fullfile(newSavepath,['Slice#',num2str(slice),'Voxel_Lorentzian.mat']),'Voxel_Lorentzian_difference','Voxel_Lorentzian_fit',...
    'freq_offsets','time_voxels','numIters','meanIters')
display_image(Voxel_Lorentzian_difference(:,:,index_APTinterpolation),brainMask,'Traditional method LD APT',APT_DisplayRange_LD)
display_image(Voxel_Lorentzian_difference(:,:,index_NOEinterpolation),brainMask,'Traditional method LD NOE',NOE_DisplayRange_LD)
