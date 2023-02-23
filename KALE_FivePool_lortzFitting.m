%% load data and set the pars
% create by sunyaozong
% 2023,2,20
clear all;clc;close all
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
All_numClass = 1; 

K_value = 3;
method_fitting = 1; % 0 : class fitting     1 : after kmeans and voxels fitting
std_range = 76;     % 76
num_range = 12;     % 12

method = 'kmeanspp'; 
offset_choice_name = 'All_dynamic_maxstd';
offset_choice_kmeans = w_offset;
filename = [offset_choice_name,'Slice',num2str(slice),'Kvalue',num2str(K_value)];
newSavepath = fullfile(Savepath,filename);
mkdir(newSavepath) 
w_offset_inter=[min(w_offset):0.1:max(w_offset)]'; 
%% KALE
Addnoise = 0;
iteration = 0;
if Addnoise
    for i = 1 : length(w_offset)
        Zdata_Noise = imnoise(V_exp_mask(:,:,i),'gaussian',0,0.08^2);
        Zspec_AddNoise(:,:,i) = Zdata_Noise;
    end
else
    Zspec_AddNoise = V_exp_mask;
end
Zspec_AddNoise(isnan(Zspec_AddNoise)) = 0;
[Row,Column,numOffset] = size(Zspec_AddNoise);
oriindex = [1:1:Row*Column]';
Zspec_loli_vec = reshape(Zspec_AddNoise,Row*Column,numOffset);
Zspec_loli_vec_ori = Zspec_loli_vec;
% multiple pool loratzian fitting
Npools = 5;
ampMaps    =zeros(Row*Column,Npools);
areaMaps   =zeros(Row*Column,Npools);
fwhmMaps   =zeros(Row*Column,Npools);
KALE_fit   =zeros(Row*Column,length(w_offset_inter));
offsetMaps =zeros(Row*Column,Npools);
ori_std = std2(Zspec_loli_vec)/std_range;    % 40
iteration_data{1,1} = {Zspec_loli_vec};
iteration_data{1,2} = {oriindex};

par1=      [0.0621    1.49    3.63];   % initial guess of in NMR
par2=      [0.1505    4.28   -3.25];   % initial guess of in NMR
par3=      [0.71      1.35    0.02];   % initial guess of in NMR
par4=      [0.203     27.43  -1.48];   % initial guess of in NMR
par5=      [0.0951    2.23    2.04];   % initial guess of in NMR


par_0 ={[par1 par2 par3 par4 par5]};
tic
while 1
    num = 1;
    for numClass = 1 : size(iteration_data,1)
        [cls,~,ture_index,index_choice] = Clustering_KALE_fitting(iteration_data{numClass,1},...
            K_value,method,w_offset,offset_choice_kmeans,iteration_data{numClass,2});
       
        unique_cls = unique(cls);
        [ZsData_exp,w_offset1,mSigEStd] = ROIprocess_Zspecforkmeans(ture_index,w_offset,Zspec_AddNoise);
        Max_mSigEStd =  max(mSigEStd(index_choice,:));
        Mean_mSigEStd = mean(mSigEStd(index_choice,:));
        [Z_spectrums, spectrumLorFit,spectrumLorFitPerPool,poolNamesCellArr,par_fit] = ...
            Kmeans_FivePool_lF_MTRspec(ZsData_exp,w_offset,par_0{numClass},iteration,Max_mSigEStd);

        for ii = 1 : length(ture_index)
            if (Mean_mSigEStd(ii) < ori_std)  || (length(ture_index{ii}) < (floor(Row/num_range))^2) % 相应类别不用进行下采样 20
                index = ture_index{ii};
%%%%%%%%%%%%%%%%%%%%%% KALE group-wise fitting %%%%%%%%%%%%%%%%%%%%%%%%
                if method_fitting == 0
%                     for j = 1 : length(index)
%                         temp = 1;
%                         for curPool = 1 : Npools
%                             ampMaps(index(j),curPool) = par_fit(ii,temp);
%                             fwhmMaps(index(j),curPool)=par_fit(ii,temp+1);
%                             offsetMaps(index(j),curPool)=par_fit(ii,temp+2);
%                             temp=temp+3;
%                         end
%                     end
%%%%%%%%%%%%%%%%%%    KALE voxels-wise fitting      %%%%%%%%%%%%%%%%%%%
                else
                    if mean(mSigEStd(:,ii)) == 0
                        for j = 1 : length(index)
                            temp = 1;
                            for curPool = 1 : Npools
                                ampMaps(index(j),curPool) = par_fit(ii,temp);
                                fwhmMaps(index(j),curPool)=par_fit(ii,temp+1);
                                offsetMaps(index(j),curPool)=par_fit(ii,temp+2);
                                temp=temp+3;
                            end
                        end
                    else
                        for i = 1 : length(index)
                            ind = index(i);
                            x = Zspec_loli_vec(ind,:);
                            [Z_spectrums, spectrumLorFit,spectrumLorFitPerPool,~,par_fit_voxels] = ...
                                Kmeans_FivePool_lF_MTRspec(x',w_offset,par_fit(ii,:),iteration,Mean_mSigEStd);
                            temp = 1;
                            for curPool = 1 : Npools
                                ampMaps(ind,curPool) = par_fit_voxels(temp);
                                fwhmMaps(ind,curPool)=par_fit_voxels(temp+1);
                                offsetMaps(ind,curPool)=par_fit_voxels(temp+2);
                                temp=temp+3;
                            end
                            KALE_fit(ind,:) =  spectrumLorFit;
                        end
                    end
                end
            else
                % Categories need to be reclassified
                % Take the original pixels out of that category
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

%% adjust the size of ampMaps,fwhmMaps and offsetMaps.Compute the R2 and CC
ampMaps    =reshape(ampMaps,Row,Column,Npools);
fwhmMaps   =reshape(fwhmMaps,Row,Column,Npools);
offsetMaps =reshape(offsetMaps,Row,Column,Npools);
KALE_fit   =reshape(KALE_fit,Row,Column,length(w_offset_inter));
[R2_KALE,R2_map_KALE,~] = Goodness_Fitting(KALE_fit,V_exp_mask,w_offset,brainMask);
[CC_mean_KALE,CC_vec_KALE] = CCcalculate(KALE_fit,V_exp_mask,w_offset,brainMask);
save(fullfile(newSavepath,['KALE','.mat']),'ampMaps','fwhmMaps','offsetMaps','KALE_fit','R2_KALE');
%% display the KALE fitting maps
poolNamesCellArr{1}='Amide pool';
poolNamesCellArr{2}='NOE pool';
poolNamesCellArr{3}='Water pool';
poolNamesCellArr{4}='MT pool';
poolNamesCellArr{5}='Amine pool';
[px,py]=find(brainMask);
temp = 6;
for i = 1:length(poolNamesCellArr)
    figure
    imtemp = ampMaps(:,:,i).*brainMask;
    idx = find(brainMask == 0);
    imtemp(idx) = -10;
    imagesc(imtemp(min(px)-temp:max(px)+temp,min(py)-temp:max(py)+temp),[0.01,0.1])
    colorbar;axis off
    colormap(jet(256));
    altered_colormap();
    str = poolNamesCellArr{i};
    set(gca, 'FontWeight','bold','FontSize',20)
    title(str,'FontWeight','bold','FontSize',18);
    savefig(fullfile(newSavepath,['KALE_offsetMaps',poolNamesCellArr{i},'.fig']))
end

%% traditional method for voxels-voxels fitting
tic
[Voxels_ampMaps, Voxels_areaMaps, Voxels_fwhmMaps, Voxels_offsetMaps, VoxelsLorFit,poolNamesCellArr] = ...
  Voxel_FivePool_lF_MTRcontrastMap(V_exp_mask,w_offset);
times = toc
[R2_Trad,R2_map_Trad,~] = Goodness_Fitting(VoxelsLorFit,V_exp_mask,w_offset,brainMask);
save(fullfile(newSavepath,['Voxels_multLF_NMR10and2and0.1','.mat']),'Voxels_ampMaps','Voxels_fwhmMaps',...
    'Voxels_offsetMaps','VoxelsLorFit','times','R2_Trad');

%% display the CEST maps
[px,py]=find(brainMask);
temp = 6;
for i = 1:length(poolNamesCellArr)
    figure
    imtemp = Voxels_ampMaps(:,:,i).*brainMask;
    idx = find(brainMask == 0);
    imtemp(idx) = -10;
    imagesc(imtemp(min(px)-temp:max(px)+temp,min(py)-temp:max(py)+temp),[0.02,0.1])
    colorbar;axis off
    colormap(jet(256));
    altered_colormap();
    str = poolNamesCellArr{i};
    set(gca, 'FontWeight','bold','FontSize',20)
    title(str,'FontWeight','bold','FontSize',18);
    savefig(fullfile(newSavepath,['KALE_offsetMaps',poolNamesCellArr{i},'.fig']))
end
