%%
rng('default'); % For reproducibility
% X = [(randn(20,2)*0.75)+1;
%     (randn(20,2)*0.25)-1];
X = Zspec_loli_vec;
clear all;clc;close all
ROIname = ['ROI 1','ROI 2'];
Savepath = 'D:\CESTcode\CEST_20210115\savedatapath\kmeans';
datapath = 'D:\CESTcode\CEST_20210115\savedatapath';
% brainMask = Thmask;
% load(fullfile(datapath,['V_exp_mask.mat']));
load(fullfile(datapath,['w_offset0.7uT.mat']));
% load(fullfile(datapath,['S0_Scan8E1_power0.7.mat']));
load(fullfile(datapath,['brainMask_s5.mat']));
load(fullfile(datapath,['V_norm_B0_0.7uT.mat']));
slice = 5;
V_exp_mask = squeeze(V_norm(:,:,slice,:));
V_exp_mask(isnan(V_exp_mask)) = 0;
V_exp_mask = V_exp_mask.*brainMask;
% add the gaussian noise
Addnoise = 0;
if Addnoise
    for i = 1 : length(w_offset)
        Zdata_Noise = imnoise(V_exp_mask(:,:,i),'gaussian',0,0.020^2);
        Zspec_AddNoise(:,:,i) = Zdata_Noise;
    end
else
    Zspec_AddNoise = V_exp_mask;
end
clear V_exp_mask
%%
% [x,y,numOffset] = size(Zspec_AddNoise);
% Zspec_loli_vec = reshape(Zspec_AddNoise,x*y,numOffset);
% offset_choice  = [3.5,3.6];
% for i = 1:length(offset_choice)
%    [~,index] = min(abs(w_offset -offset_choice(i)));
%    index1(i,1) = index;
% end
% X = Zspec_loli_vec(:,index1);% Create a scatter plot of the data.
% X = X(find(X));
X = iteration_data{numClass,1};
X = X(:,index1);
scatter(X(:,1),X(:,2));
title('Randomly Generated Data');
% Create a hierarchical cluster tree using the ward linkage method.
Z = linkage(X,'ward');
% Create a dendrogram plot of the data.
dendrogram(Z)
% Cluster the data using a threshold of 3 for the inconsistency coefficient and looking to a depth of 4 below 
% each node. Plot the resulting clusters.
T = cluster(Z,'cutoff',2,'Depth',3);
gscatter(X(:,1),X(:,2),T)
index = iteration_data{numClass,2};
image_brain = zeros(x*y,1);
image_brain1 =  zeros(x*y,1);
image_brain(index) = X(:,1);
image_brain1(index) = T;
image_brain = reshape(image_brain,x,y);
h1 = figure;   imagesc(image_brain,[0,1.25])
colorbar;axis off;colormap(jet);altered_colormap();title('KmeansAPTLF','FontSize',18)    
image_brain1 = reshape(image_brain1,x,y);
h1 = figure;   imagesc(image_brain1)
colorbar;axis off;colormap(jet);altered_colormap();title('KmeansAPTLF','FontSize',18)    
