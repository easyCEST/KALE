%CC calculate
function [CC_mean,CC_vec] = CCcalculate(DataFit,cestNormData,freq_offsets,brainMask)
%%
do_interpolation = 1;
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
end
% [~,~,~] = size(cestNormData);
%%
CC_vec = zeros(length(freq_offsets),1);
[r,c] = find(brainMask);
for cur_offset = 1 : length(freq_offsets)
    
%     [~,index] = min(abs(freq_offsets-freq_offsets(i)));
    L = (squeeze(DataFit(:,:,cur_offset)));
    Ori = squeeze(cestNormData(:,:,cur_offset));

    L_aver = mean2(L(r,c));
    Ori_aver = mean2(Ori(r,c));
    sum1 = 0;
    sum2 = 0;
    for indx=drange(1: length(r))
            CC1 = [L(r(indx),c(indx))-L_aver]*[Ori(r(indx),c(indx))-Ori_aver];
            sum1 = sum1+CC1;
            CC2 =sqrt([L(r(indx),c(indx))-L_aver]^2*[Ori(r(indx),c(indx))-Ori_aver]^2);
            sum2 = sum2+CC2;
    end
    CC = sum1/sum2;
    CC_vec(cur_offset) = CC;
end
CC_mean = mean(CC_vec);
end