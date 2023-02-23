function [R2,R2_map,R2_MaskMean,R2_Voxels] = Goodness_Fitting(Z_spenturmfit,cestNormData,freq_offsets,brainMask)

%%

indexMask = find(brainMask == 1);

%% interpolation
do_interpolation = 0;
if do_interpolation
    ori_freq_offsets = freq_offsets;
    ori_cestNormData = cestNormData;
    [V_wxp_ROIori_row,V_wxp_ROIori_column,~] = size(ori_cestNormData);
    
    freq_stepsize = 0.1;
    freq_offsets = [min(freq_offsets) : freq_stepsize : max(freq_offsets)]';

    cestNormData = zeros(V_wxp_ROIori_row, V_wxp_ROIori_column, length(freq_offsets));
    for i = 1 : V_wxp_ROIori_row
        for j = 1 : V_wxp_ROIori_column
            cestNormData(i,j,:) = spline(ori_freq_offsets,ori_cestNormData(i,j,:),freq_offsets);
        end
    end
    
else
    ori_freq_offsets = freq_offsets;
    freq_stepsize = 0.1;
    freq_offsets = [min(freq_offsets) : freq_stepsize : max(freq_offsets)]';
    
    for i = 1:length(ori_freq_offsets)
       [~,index] = min(abs(freq_offsets -ori_freq_offsets(i)));
       index1(i,1) = index;
    end
    Z_spenturmfit = Z_spenturmfit(:,:,index1);
    freq_offsets = ori_freq_offsets;
end
[cestNormData_row,cestNormData_column,cestNormData_numoffset] = size(cestNormData);
%%
cestNormData = reshape(cestNormData,[cestNormData_row*cestNormData_column,cestNormData_numoffset]);
Z_spenturmfit = reshape(Z_spenturmfit,[cestNormData_row*cestNormData_column,cestNormData_numoffset]);
R2_Voxels = zeros(length(indexMask),1);


% [~,idx] = min(abs(freq_offsets-0));

[~,idxp1] = min(abs(freq_offsets-1));
[~,idxp6] = min(abs(freq_offsets-6));

for i = 1 : length(indexMask)
    Z_specturm = squeeze(cestNormData(indexMask(i),:));
    fit_specturm = squeeze(Z_spenturmfit(indexMask(i),:));
    % [-5,-1]   [1,5]
%     Z_specturm = Z_specturm(1:idx);
%     fit_specturm = fit_specturm(1:idx);
%     
%     Z_specturm = Z_specturm(idxp1:idxp6);
%     fit_specturm = fit_specturm(idxp1:idxp6);
    
    R2_Voxels(i) = 1 - sum((Z_specturm-fit_specturm).^2)/sum((Z_specturm-mean(Z_specturm)).^2);
end
R2_map = zeros(cestNormData_row*cestNormData_column,1);
R2_map(indexMask) = R2_Voxels;
R2_map = reshape(R2_map,[cestNormData_row,cestNormData_column]);
figure
imagesc(R2_map,[0.9,1])
axis off
colorbar
colormap(jet)
altered_colormap();
set(gca, 'FontWeight','bold','FontSize',20)
title('Goodness of fit (R2)','FontWeight','bold','FontSize',14)
R2 = mean(R2_Voxels);
%%
fitMask = Z_spenturmfit(indexMask,:);
ZspecMask = cestNormData(indexMask,:);
fitMaskMean = mean(fitMask,1);
ZspecMaskMean = mean(ZspecMask,1);
R2_MaskMean = 1 - sum((ZspecMaskMean-fitMaskMean).^2)/sum((ZspecMaskMean-mean(ZspecMaskMean)).^2);
end