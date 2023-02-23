%% the data prepare
clear all;clc;close all
ROIname = ['ROI 1','ROI 2'];
datapath = 'G:\实验数据\肿瘤数据\panlaiwang\apt_165\MATfile';
Savepath = fullfile(datapath,'kmeans');
%%%%%%%%%%%%%%% GBM DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(datapath,['test_0p7uT_B0_jiang.mat']));
slice = 3;
% load(fullfile(Savepath,['brainMaskS_4.mat']));
brainMask = Thmask_0p7uT(:,:,slice);
S0 = M0_stack(:,:,slice);
% brainMask = GeneratebrainMask(S0,Savepath);
V_norm = squeeze(V_norm_B0_0p7uT(:,:,slice,:));
% [~,index_APT] = min(abs(w_offset-3.5));
% disp("starting registration to S0 ........")
% [Vo_ori_reg] = Registration(V_norm,V_norm(:,:,index_APT),2,Savepath);
% disp("ending registration to S0 ........")
[V_exp_mask] = prepare(S0,V_norm,brainMask);
V_exp_mask(isnan(V_exp_mask)) = 0;
% load(fullfile(datapath,['test_0p7uT_B0_jiang.mat']));
% w_offset_inter=[min(w_offset):0.1:max(w_offset)]'; 
%%%%%%%%%%%%%%% human data %%%%%%%%%%%%%%%%%
% load(fullfile(datapath,['w_offset0.7uT.mat']));
% load(fullfile(datapath,['brainMask_s3.mat']));
% load(fullfile(datapath,['V_norm_B0_0.7uT.mat']));
% % load(fullfile(datapath,['S0.mat']));
% slice = 3;
% % S0 = M0_stack(:,:,slice);
% % [V_exp_mask] = prepare(S0,V_exp,brainMask);
% V_exp_mask = squeeze(V_norm(:,:,slice,:));
% V_exp_mask(isnan(V_exp_mask)) = 0;
% % w_offset_inter=[min(w_offset):0.1:max(w_offset)]'; 
% % [V_exp_mask] = Registration(V_exp_mask,S0,1,datapath);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load(fullfile(datapath,['w_offset.mat']));
% load(fullfile(datapath,['brainMask.mat']));
% load(fullfile(datapath,['V_exp_mask.mat']));
% % % load(fullfile(datapath,['S0.mat']));
% slice = 1;
% % % [V_exp_mask] = prepare(S0,V_exp,brainMask);
% % V_exp_mask = squeeze(V_norm(:,:,slice,:));
% V_exp_mask(isnan(V_exp_mask)) = 0;
%%%%%%%%%%% rat data %%%%%%%%%%%%%%%%%%%%%%%%%%%
% brainMask = Thmask;
% load(fullfile(datapath,['Vo_ori8.mat']));
% 
% load(fullfile(datapath,['S0_Scan8E1_power0.7.mat']));
% load(fullfile(datapath,['w_offset Scan8.mat']));
% load(fullfile(datapath,['Thmask_Scan8S1.mat']));
% % brainMask = Thmask;
% brainMask = GeneratebrainMask(S0,Savepath);
% [V_exp_mask] = prepare(S0,Vo_ori,brainMask);
% slice = 0;
% % brainMask = Thmask;
% % V_exp_mask = V_exp_mask(:,:,1:end-1);
% % w_offset = w_offset(1:end-1);
% % w_offset = w/300.33;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V_exp_mask = V_exp_mask.*brainMask;
% add the gaussian noise
Addnoise = 0;
if Addnoise
    for i = 1 : length(w_offset)
        Zdata_Noise = imnoise(V_exp_mask(:,:,i),'gaussian',0,0.08^2);
        Zspec_AddNoise(:,:,i) = Zdata_Noise;
    end
else
%     for i = 1 : length(w_offset)
%         temp = squeeze(V_exp_mask(:,:,i));
%         temp_filt = medfilt2(temp,[2,2]);
%         Zspec_AddNoise(:,:,i) = temp_filt;
%     end
    Zspec_AddNoise = V_exp_mask;
end
Zspec_AddNoise = Zspec_AddNoise.*brainMask;
[~,index_APT] = min(abs(w_offset-3.5));
display_image(Zspec_AddNoise(:,:,index_APT),brainMask,'ori',[0.6,1])
w_offset_inter=[min(w_offset):0.1:max(w_offset)]'; 

% figure;   imagesc(Zspec_AddNoise(:,:,index_APT),[0.1,0.9])    % [0.6,1]
% colorbar;axis off;colormap(jet);altered_colormap();title('Zspec\_AddNoise APT','FontSize',18) 
% clear V_exp_mask
%%
K_value = 3;
iteration = 0;
method_fitting = 0; % 0 : class fitting     1 : after kmeans and voxels fitting
All_numClass = 1; 
medfilter = 1;
std_range = 120;       % 50
num_range = 50;
noise_std = 0;
method = 'kmeanspp'; 
offset_choice_name = 'All_AddVoxels';
% offset_choice_kmeans  = [3,3.1,3.4,3.5,3.4,3.7];
offset_choice_kmeans = w_offset;
filename = [offset_choice_name,'Slice',num2str(slice),'Kvalue',num2str(K_value),'Noise',num2str(noise_std),'interpolation'];
newSavepath = fullfile(Savepath,filename);
if exist(newSavepath) ~=7       
   mkdir(newSavepath) 
end
[Row,Column,numOffset] = size(Zspec_AddNoise);
oriindex = [1:1:Row*Column]';
Zspec_loli_vec = reshape(Zspec_AddNoise,Row*Column,numOffset);
Zspec_loli_vec_ori = Zspec_loli_vec;
% % % add noise
if noise_std ~= 0
    Zspec_loli_vec = imnoise(Zspec_loli_vec,'gaussian',0,noise_std^2);
    temp = reshape(Zspec_loli_vec,Row,Column,numOffset);
    temp = temp.*brainMask;
    Zspec_loli_vec = reshape(temp,Row*Column,numOffset);
    clear 'temp'
%     Zspec_AddNoise = reshape(Zspec_loli_vec,Row,Column,numOffset);
end

cls_loli = zeros(Row*Column,1);
final_lortz = zeros(Row*Column,numOffset);
final_lortz_as = zeros(Row*Column,length(w_offset_inter));
final_lortz_as_Kvoxels = zeros(Row*Column,length(w_offset_inter));
final_mean = zeros(Row*Column,numOffset);
final_std = zeros(Row*Column,numOffset);
final_lortz_LD  = zeros(Row*Column,numOffset);
final_lortz_LD_as = zeros(Row*Column,length(w_offset_inter));
final_lortz_LD_as_AREX = zeros(Row*Column,length(w_offset_inter));


ori_std = std2(Zspec_loli_vec)/std_range;    % 40
iteration_data{1,1} = {Zspec_loli_vec};
iteration_data{1,2} = {oriindex};
par_0 = {[0,  50, 20, 1]};
tic
while 1
    num = 1;
    for numClass = 1 : size(iteration_data,1)
        [cls,~,ture_index,index_choice] = Clustering_IDEAL_fitting(iteration_data{numClass,1},K_value,method,w_offset,offset_choice_kmeans,iteration_data{numClass,2});
%%%%%%%%%%%%%%%%%%%%%%%%%%%为了画出分裂的不同的类%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         unique_cls = unique(cls);
%         for iii = 1 : length(ture_index)
%             if iteration == 0
%                 ind0 = ture_index{iii};
%                 clsvalue = unique_cls(iii);
%                 cls_loli(ind0) = clsvalue;
%             else
%                 ind0 = ture_index{iii};
%                 clsvalue = (cls_loli(ind0(1))/3)*unique_cls(iii);
%                 cls_loli(ind0) = clsvalue;
%             end
%         end
%         cls_loli = reshape(cls_loli,64,64);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [ZsData_exp,w_offset1,mSigEStd] = ROIprocess_Zspecforkmeans(ture_index,w_offset,Zspec_AddNoise);
        [Lorentz_fit_Kmeans,Lor_fit,Lorentzian_diff_Kmeans,~,par_fit,AREX_LD] = LF_Poly(w_offset1,ZsData_exp,ture_index,-2,-6.25,-10,2,6.25,10,par_0{numClass},iteration);
        Mean_mSigEStd =  mean(mSigEStd(index_choice,:));
%         All_mSigEStd{iteration+1,numClass} = Mean_mSigEStd;
        for ii = 1 : length(ture_index)
            if (Mean_mSigEStd(ii) < ori_std)  || (length(ture_index{ii}) < (floor(Row/num_range))^2) % 相应类别不用进行下采样 20
%             if All_mSigEStd{iteration+1,numClass}(ii)/Mean_mSigEStd(ii)
    %             final_index = Voxel_Index{ii};
%                 Lor_fit_loilvec = reshape(Lor_fit,x*y,numOffset);
                index = ture_index{ii};
%%%%%%%%%%%%%%%%%%%%%% 一次拟合作为最终结果 %%%%%%%%%%%%%%%%%%%%%%%%
                if method_fitting == 0
                    for j = 1 : length(index)
                        final_lortz(index(j),:) = Lor_fit(:,ii)';
                        final_lortz_as(index(j),:) = Lorentz_fit_Kmeans(:,ii)';
                        final_lortz_LD_as(index(j),:) = Lorentzian_diff_Kmeans(:,ii)';
                        final_lortz_LD_as_AREX(index(j),:) = AREX_LD(:,ii)';
                    end
%%%%%%%%%%%%%%%%%%使用原始像素点进行二次拟合数据%%%%%%%%%%%%%%%%%%%
                else
                    if mSigEStd(ii) == 0
                        for j = 1 : length(index)
                            final_lortz(index(j),:) = Lor_fit(:,ii)';
                            final_lortz_as(index(j),:) = Lorentz_fit_Kmeans(:,ii)';
                        end
                    else
                        for i = 1 : length(index)
                            ind = index(i);
                            x = Zspec_loli_vec(ind,:);
        %                     final_lortz(ind,:) = lorentz_N(par_fit(ii,:),x);
                            [Lorentzian_fit, Lorentzian_difference, AREX_Lorentzian_difference,~]...
                                    = spectrum_LD_MTRrex_interpolation(x', w_offset,par_fit(ii,:),iteration);
%                             final_lortz(ind,:) =  Lorentzian_fit;
                            final_lortz_as_Kvoxels(ind,:) = Lorentzian_fit;
                            final_lortz_LD_as(ind,:) = Lorentzian_difference;
                            final_lortz_LD_as_AREX(ind,:) = AREX_Lorentzian_difference;
                        end
                    end
                end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
            else   % 类别需要再分类
                % 把该类别的原始像素点取出来
                index = ture_index{ii};
                Zsdata_class_lolivec{num,1} = Zspec_loli_vec(index,:);
                par_0{num} = par_fit(ii,:);
                Zsdata_class_lolivec{num,2} = index;
%                 for j = 1 : length(index)
%                    final_lortz(index(j),:) = Lor_fit(:,ii)';
%                 end 
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
%     final_lortz_reshape = reshape(final_lortz,Row,Column,numOffset);
%     h1 = figure;   imagesc(final_lortz_reshape(:,:,index_APT),[0.6,1])
%     str = ['iteration',num2str(iteration),'KmeansLF'];
% %     colorbar;axis off;colormap(jet);altered_colormap();title('KmeansAPTLF','FontSize',18)    
%     display_image(final_lortz_reshape(:,:,index_APT),brainMask,str,[0.6,1])

end

toc

%% 
Zspec_AddNoise = reshape(Zspec_loli_vec,Row,Column,numOffset);
Z_spec_interpolation = interpolation(Zspec_AddNoise,w_offset);
[~,index_APT] = min(abs(w_offset_inter-3.5));
[~,index_NOE] = min(abs(w_offset_inter+3.5));
%%%%%%%%%%%%%%% 使用两次洛伦兹拟合时候使用 %%%%%%%%%%%%%%%%%%%%
APT_DisplayRange_LD = [0,0.03];   %human
NOE_DisplayRange_LD = [0.01,0.055];
% APT_DisplayRange_LD = [0.035,0.12];   %stork rat
% NOE_DisplayRange_LD = [0.05,0.2];
% APT_DisplayRange_LD = [0,0.04];   %gbm
% NOE_DisplayRange_LD = [0,0.06];
APT_DisplayRange_LDrex = [0,0.05];   %gbm
NOE_DisplayRange_LDrex = [0.01,0.07];

Display_Range_LF = [0.6,1];
final_lortz_as_im = reshape(final_lortz_as,Row,Column,length(w_offset_inter));
final_lortz_im_VOXELS = reshape(final_lortz_as_Kvoxels,Row,Column,length(w_offset_inter));
final_lortz_LD_as = reshape(final_lortz_LD_as,Row,Column,length(w_offset_inter));
final_lortz_LD_as_AREX = reshape(final_lortz_LD_as_AREX,Row,Column,length(w_offset_inter));

%%%%%%%%%%%  使用一次洛伦兹拟合时候使用  %%%%%%%%%%%%%%

if method_fitting == 0
    for numoffset = 1 : length(w_offset_inter)
        final_lortz_LDas_AREX_voxels(:,:,numoffset) = 1./(Z_spec_interpolation(:,:,numoffset) + 1e-5) - 1./(final_lortz_as_im(:,:,numoffset) + 1e-5);
        final_lortz_LDas_voxels(:,:,numoffset) = final_lortz_as_im(:,:,numoffset) - Z_spec_interpolation(:,:,numoffset);
    end
    if medfilter
        for numoffset = 1 : length(w_offset_inter)
            final_lortz_LD_as(:,:,numoffset) = medfilt2(final_lortz_LD_as(:,:,numoffset),[2,2]);
            final_lortz_LD_as_AREX(:,:,numoffset) = medfilt2(final_lortz_LD_as_AREX(:,:,numoffset),[2,2]);

            final_lortz_LDas_voxels(:,:,numoffset) = medfilt2(final_lortz_LDas_voxels(:,:,numoffset),[2,2]);
            final_lortz_LDas_AREX_voxels(:,:,numoffset) = medfilt2(final_lortz_LDas_AREX_voxels(:,:,numoffset),[2,2]);

        end
    end
    display_image(final_lortz_LD_as(:,:,index_APT),brainMask,'K class LD APT',APT_DisplayRange_LD)
    display_image(final_lortz_LD_as(:,:,index_NOE),brainMask,'K class LD NOE',NOE_DisplayRange_LD)
    display_image(final_lortz_LD_as_AREX(:,:,index_APT),brainMask,'K class LDrex APT',APT_DisplayRange_LDrex)
    display_image(final_lortz_LD_as_AREX(:,:,index_NOE),brainMask,'K class LDrex NOE',NOE_DisplayRange_LDrex)
    
    display_image(final_lortz_LDas_voxels(:,:,index_APT),brainMask,'K class(voxels) LD APT',APT_DisplayRange_LD)
    display_image(final_lortz_LDas_voxels(:,:,index_NOE),brainMask,'K class(voxels) LD NOE',NOE_DisplayRange_LD)
    display_image(final_lortz_LDas_AREX_voxels(:,:,index_APT),brainMask,'K class(voxels) LDrex APT',APT_DisplayRange_LDrex)
    display_image(final_lortz_LDas_AREX_voxels(:,:,index_NOE),brainMask,'K class(voxels) LDrex NOE',NOE_DisplayRange_LDrex)
    
    save(fullfile(newSavepath,['Kclass','.mat']),'final_lortz_LD_as','final_lortz_LDas_voxels',...
        'final_lortz_as_im','Zspec_AddNoise','final_lortz_LD_as_AREX','final_lortz_LDas_AREX_voxels')
else
    for numoffset = 1 : length(w_offset_inter)
        Kvoxels(:,:,numoffset) = final_lortz_im_VOXELS(:,:,numoffset) - Z_spec_interpolation(:,:,numoffset);
        Kvoxels_AREX(:,:,numoffset) =  1./(Z_spec_interpolation(:,:,numoffset) + 1e-5) - 1./(final_lortz_im_VOXELS(:,:,numoffset) + 1e-5) ;
    end
    if medfilter
        for numoffset = 1 : length(w_offset_inter)
            Kvoxels(:,:,numoffset) = medfilt2(Kvoxels(:,:,numoffset),[2,2]);
            Kvoxels_AREX(:,:,numoffset) = medfilt2(Kvoxels_AREX(:,:,numoffset),[2,2]);
        end
    end
    display_image(Kvoxels(:,:,index_APT),brainMask,'K Voxels LD APT',APT_DisplayRange_LD)
    display_image(Kvoxels(:,:,index_NOE),brainMask,'K Voxels LD NOE',NOE_DisplayRange_LD)
    display_image(Kvoxels_AREX(:,:,index_APT),brainMask,'K Voxels LDrex APT',APT_DisplayRange_LDrex)
    display_image(Kvoxels_AREX(:,:,index_NOE),brainMask,'K Voxels LDrex NOE',NOE_DisplayRange_LDrex)
    save(fullfile(newSavepath,['Kvoxels','.mat']),'Kvoxels','final_lortz_im_VOXELS','Kvoxels_AREX')
end

%%
% tic
% [Voxel_Lorentzian_difference,Voxel_Lorentzian_fit,freq_offsets]=Voxel_LD(Zspec_AddNoise,w_offset);
% time_voxels = toc
% Zspec_AddNoise_reloli = reshape(Zspec_loli_vec,Row,Column,numOffset);
% Zspec_AddNoise_reloli = Zspec_AddNoise_reloli.*brainMask;
% Zspec_AddNoise = Zspec_AddNoise.*brainMask;
tic
[Voxel_Lorentzian_differenceRex,Voxel_Lorentzian_difference,Voxel_Lorentzian_fit,freq_offsets]...
    =Voxel_LD_useful_pixels(Zspec_AddNoise,w_offset);
time_voxels = toc
if medfilter
    for numoffset = 1 : length(freq_offsets)
        Voxel_Lorentzian_difference(:,:,numoffset) = medfilt2(Voxel_Lorentzian_difference(:,:,numoffset),[2,2]);
        Voxel_Lorentzian_differenceRex(:,:,numoffset) = medfilt2(Voxel_Lorentzian_differenceRex(:,:,numoffset),[2,2]);
    end
end

[~,index_APTinterpolation] = min(abs(freq_offsets-3.5));
[~,index_NOEinterpolation] = min(abs(freq_offsets+3.5));
% [~,index_n10_inter] = min(abs(freq_offsets+10));
% [~,index_n10] = min(abs(w_offset+10));
% save(fullfile(newSavepath,['Slice#',num2str(slice),'Voxel_Lorentzian.mat']),'Voxel_Lorentzian_difference','Voxel_Lorentzian_fit',...
%     'Voxel_Lorentzian_differenceRex','freq_offsets','time_voxels')
% ssim_n10_voxelsfit = ssim(Voxel_Lorentzian_fit(:,:,index_n10_inter),Zspec_AddNoise(:,:,index_n10));
% ssim_n10_Kclassfit = ssim(final_lortz_im(:,:,index_n10),Zspec_AddNoise(:,:,index_n10));
display_image(Voxel_Lorentzian_difference(:,:,index_APTinterpolation),brainMask,'Traditional method LD APT',APT_DisplayRange_LD)
display_image(Voxel_Lorentzian_difference(:,:,index_NOEinterpolation),brainMask,'Traditional method LD NOE',NOE_DisplayRange_LD)
display_image(Voxel_Lorentzian_differenceRex(:,:,index_APTinterpolation),brainMask,'Traditional method LDrex APT',APT_DisplayRange_LDrex)
display_image(Voxel_Lorentzian_differenceRex(:,:,index_NOEinterpolation),brainMask,'Traditional method LDrex NOE',NOE_DisplayRange_LDrex)

% display_image(aa,brainMask,'diff NOE',[-0.00000001,0.00000001])
%% Multiple Pool lorentz Fitting 

[ampMaps, areaMaps, fwhmMaps, offsetMaps, VpoolNamesCellArr] = ...
          Voxel_MultiplePool_lF_MTRcontrastMap(Zspec_AddNoise,w_offset);


%% coefficient of variation
newSavepath = 'G:\实验数据\肿瘤数据\panlaiwang\apt_165\MATfile\kmeans\All_AddVoxelsSlice3Kvalue3Noise0interpolation';
ROInum = 2;
% % ROIname = ['ROI in CBM','ROI in contralateral WM'];
ROIname = ['ROI 1','ROI 2','ROI 3','ROI 4'];
% newSavepath = 'D:\dataset\20210521_brain\Kmeans\All_AddVoxelsSlice3Kvalue3Noise0.03';
% load(fullfile(newSavepath,'ROI_total_4X4.mat'));
load(fullfile(datapath,'resliced_DWI_img.mat'))
[ROI_total]=ROIplotSong(ROInum,ROIname,brainMask,resliced_DWI_img(:,:,3),newSavepath);
% ROI_total = ROI_square(ROInum, resliced_DWI_img(:,:,slice).*brainMask, 2);
% save(fullfile(newSavepath,['ROI_total_4X4','.mat']),'ROI_total')
save(fullfile(newSavepath,['ROI_total','.mat']),'ROI_total')
% load(fullfile(newSavepath,'ROI_total.mat'));
% load(fullfile(newSavepath,'ROI_total_4X4.mat'));

% ROI_total = {ROI_total{4};ROI_total{8}};
load(fullfile(newSavepath,['Slice#',num2str(3),'Voxel_Lorentzian.mat']));
load(fullfile(newSavepath,['Kvoxels','.mat']));
load(fullfile(newSavepath,['Kclass','.mat']));

% end
% ROInum = length(ROI_total);
% ROI_total = {ROI_total{1},ROI_total{2}};
ROIpath = fullfile(newSavepath,['ROIanaysis','num',num2str(ROInum)]);
if exist(ROIpath) ~= 7
    mkdir(ROIpath)
end
[Spec_Kclass_class,~,SpecStd_Kclass_class] = ROIprocess_Zspecforkmeans(ROI_total,w_offset_inter,final_lortz_LD_as);
[Spec_Kclass_voxels,~,SpecStd_Kclass_voxels] = ROIprocess_Zspecforkmeans(ROI_total,w_offset_inter,final_lortz_LDas_voxels);
[Spec_Kvoxels,~,SpecStd_Kvoxels] = ROIprocess_Zspecforkmeans(ROI_total,w_offset_inter,Kvoxels);
[Spec_Trad,~,SpecStd_Trad] = ROIprocess_Zspecforkmeans(ROI_total,w_offset_inter,Voxel_Lorentzian_difference);
temp1 = [Spec_Trad,Spec_Kclass_class,Spec_Kclass_voxels,Spec_Kvoxels];
%plot lf spectrum
[Z_Spec,~,Z_SpecStd] = ROIprocess_Zspecforkmeans(ROI_total,w_offset,Zspec_AddNoise);
[LF_Spec_Kclass,~,LF_SpecStd_Kclass] = ROIprocess_Zspecforkmeans(ROI_total,w_offset_inter,final_lortz_as_im);
[LF_Spec_Kvoxels,~,LF_SpecStd_Kvoxels] = ROIprocess_Zspecforkmeans(ROI_total,w_offset_inter,final_lortz_im_VOXELS);
[LF_Spec_Trad,~,LF_SpecStd_Trad] = ROIprocess_Zspecforkmeans(ROI_total,w_offset_inter,Voxel_Lorentzian_fit);
temp2 = [LF_Spec_Trad,LF_Spec_Kclass,LF_Spec_Kvoxels];

[figh, Lorentzian_fit_ref, Lorentzian_difference_ref, AREX_Lorentzian_difference_ref,freq_offsets]...
                            = perform_Lorentzian_fitting_LD_MTRrex(Z_Spec(:,1), w_offset);
% muitiple pool luotrz fit
[Z_spectrums, spectrumLorFit,spectrumLorFitPerPool,poolNamesCellArr] = ...
          Specturm_FourPool_lF_MTRspec(Z_Spec(:,1),w_offset);
                        
                        
                        
plot_z(newSavepath,ROInum,w_offset,Z_Spec,LF_Spec_Trad,Spec_Trad,2,1)
plot_z(newSavepath,ROInum,w_offset,Z_Spec,LF_Spec_Kclass,Spec_Kclass_class,2,2)
plot_z(newSavepath,ROInum,w_offset,Z_Spec,LF_Spec_Kclass,Spec_Kclass_voxels,2,3)
plot_z(newSavepath,ROInum,w_offset,Z_Spec,LF_Spec_Kvoxels,Spec_Kvoxels,2,4)
%%
[~,index_singleoft] = min(abs(w_offset_inter+1.4));
[~,index_singleoftNO] = min(abs(w_offset+1.4));
display_image(Zspec_AddNoise(:,:,index_singleoftNO),brainMask,'ORI',NOE_DisplayRange_LD)
display_image(Voxel_Lorentzian_fit(:,:,index_singleoft),brainMask,'Traditional method LF',NOE_DisplayRange_LD)
display_image(final_lortz_im_VOXELS(:,:,index_singleoft),brainMask,'Kvoxels method LF',NOE_DisplayRange_LD)
display_image(final_lortz_as_im(:,:,index_singleoft),brainMask,'Kclass method LF',NOE_DisplayRange_LD)

figure
imagesc(Zspec_AddNoise(:,:,index_singleoftNO) - Z_spec_interpolation(:,:,index_singleoft));

aa = Voxel_Lorentzian_fit(:,:,index_singleoft) - Zspec_AddNoise(:,:,index_singleoftNO);
display_image(aa,brainMask,'Traditional method LD',NOE_DisplayRange_LD)
%%
offset_choice_SNR  = [3.5,-3.5];
for i = 1:length(offset_choice_SNR)
   [~,index] = min(abs(w_offset_inter -offset_choice_SNR(i)));
   index2(i,1) = index;
end
mean_gbm_trad = Spec_Trad(index2,1);
mean_gbm_kclass_class = Spec_Kclass_class(index2,1);
mean_gbm_kclass_voxels = Spec_Kclass_voxels(index2,1);
mean_gbm_kvoxels = Spec_Kvoxels(index2,1);
mean_gbm = [mean_gbm_trad';mean_gbm_kclass_class';mean_gbm_kclass_voxels';mean_gbm_kvoxels'];

mean_wm_trad = Spec_Trad(index2,2);
mean_wm_kclass_class = Spec_Kclass_class(index2,2);
mean_wm_kclass_voxels = Spec_Kclass_voxels(index2,2);
mean_wm_kvoxels = SpecStd_Kvoxels(index2,2);
mean_wm = [mean_wm_trad';mean_wm_trad';mean_wm_trad';mean_wm_trad'];

std_gbm_trad = SpecStd_Trad(index2,1);
std_gbm_kclass_class = SpecStd_Kclass_class(index2,1);
std_gbm_kclass_voxels = SpecStd_Kclass_voxels(index2,1);
std_gbm_Kvoxels = SpecStd_Kvoxels(index2,1);
std_gbm = [std_gbm_trad';std_gbm_kclass_class';std_gbm_kclass_voxels';std_gbm_Kvoxels'];

std_wm_trad = SpecStd_Trad(index2,2);
std_wm_kclass_class = SpecStd_Kclass_class(index2,2);
std_wm_kclass_voxels = SpecStd_Kclass_voxels(index2,2);
std_wm_Kvoxels = SpecStd_Kvoxels(index2,2);
std_wm = [std_wm_trad';std_wm_kclass_class';std_wm_kclass_voxels';std_wm_Kvoxels'];

%%
% close all
% cnr(n) = abs([mean(inflamed) - mean(normal)] / sqrt((std(inflamed))^2 + (std(normal))^2));
offset_choice_SNR  = [3.5,-3.5];
for i = 1:length(offset_choice_SNR)
   [~,index] = min(abs(w_offset_inter -offset_choice_SNR(i)));
   index2(i,1) = index;
end
Trad_mean = mean(Spec_Trad(index2,:),2);
Trad_std = mean(SpecStd_Trad(index2,:),2);

Kclass_class_mean = mean(Spec_Kclass_class(index2,:),2);
Kclass_class_std = mean(SpecStd_Kclass_class(index2,:),2);

Kclass_voxels_mean = mean(Spec_Kclass_voxels(index2,:),2);
Kclass_voxels_std = mean(SpecStd_Kclass_voxels(index2,:),2);

Kvoxels_mean = mean(Spec_Kvoxels(index2,:),2);
Kvoxels_std = mean(SpecStd_Kvoxels(index2,:),2);

mean_data = [Trad_mean,Kclass_class_mean,Kclass_voxels_mean,Kvoxels_mean];
std_data = [Trad_std,Kclass_class_std,Kclass_voxels_std,Kvoxels_std];
% ori_snr_anaysis = [snr_ori_ld(index2,:), LDori(index2,:),mSigEStdori(index2,:)];
% kmeans_snr_anaysis = [snr_kmeans_ld(index2,:),LD_averagespecKmeans(index2,:),mSigEStd_LD(index2,:)];
save(fullfile(ROIpath,['snr_anaysis']),'mean_data','std_data')
%% compute CNR
R1 = 3
R2 = 4
[cnr_APT_Kclass, cnr_NOE_Kclass] = CNR_compute(ROI_total,final_lortz_LD,w_offset,R1,R2);
[cnr_APT_KclassV, cnr_NOE_KclassV] = CNR_compute(ROI_total,final_lortz_LD_voxels,w_offset,R1,R2);
[cnr_APT_Kvoxels, cnr_NOE_Kvoxels] = CNR_compute(ROI_total,Kvoxels,w_offset,R1,R2);
[cnr_APT_Trad, cnr_NOE_Trad] = CNR_compute(ROI_total,Voxel_Lorentzian_difference,w_offset,R1,R2);
all_CNRapt = [cnr_APT_Trad,cnr_APT_Kclass,cnr_APT_KclassV,cnr_APT_Kvoxels];
all_CNRnoe = [cnr_NOE_Trad,cnr_NOE_Kclass,cnr_NOE_KclassV,cnr_NOE_Kvoxels];

%%  画APT和NOE的不同方法的箱型图

% save(fullfile(newSavepath,['Slice#',num2str(slice),'all_LD_Nonoise.mat']),'final_lortz_LD','final_lortz_LD_voxels',...
%     'Voxel_Lorentzian_difference')
%  save(fullfile(newSavepath,['Slice#',num2str(slice),'Kvoxels_no_noise.mat']),'KvoxelsNonoise')
K = 2;
load(fullfile(Savepath,['Kvoxels',num2str(noise_std),'.mat']))
Display_ROI(ROI_total,Zspec_AddNoise(:,:,index_APT),K)
APT_data_Kclaaa = squeeze(final_lortz_LD(:,:,index_APT));
APT_data_Kvoxels = squeeze(Kvoxels(:,:,index_APT));
APT_data_KclassVoxels = squeeze(final_lortz_LD_voxels(:,:,index_APT));
APT_data_Tradition = squeeze(Voxel_Lorentzian_difference(:,:,index_APTinterpolation));

NOE_data_Kclaaa = squeeze(final_lortz_LD(:,:,index_NOE));
NOE_data_Kvoxels = squeeze(Kvoxels(:,:,index_NOE));
NOE_data_KclassVoxels = squeeze(final_lortz_LD_voxels(:,:,index_NOE));
NOE_data_Tradition = squeeze(Voxel_Lorentzian_difference(:,:,index_NOEinterpolation));

APT_box_data = [APT_data_Tradition(ROI_total{K}),APT_data_Kclaaa(ROI_total{K}),APT_data_KclassVoxels(ROI_total{K}),...
    APT_data_Kvoxels(ROI_total{K})];
NOE_box_data = [NOE_data_Tradition(ROI_total{K}),NOE_data_Kclaaa(ROI_total{K}),NOE_data_KclassVoxels(ROI_total{K}),...
    NOE_data_Kvoxels(ROI_total{K})];
box_data = [APT_box_data,NOE_box_data];
figure
boxplot(box_data)
title(['ROI',num2str(K)])
set(gca, 'FontWeight','bold','FontSize',16)
