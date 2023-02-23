function [Lorentzian_differenceRex,Lorentzian_difference,Lorentzian_fit,freq_offsets,numIters]=...
    Voxel_LD_useful_pixels(V_exp_ROI,freq_offsets)

%% 
do_interpolation = 1;
if do_interpolation
    ori_freq_offsets = freq_offsets;
    ori_V_exp_ROI = V_exp_ROI;
    [V_wxp_ROIori_row,V_wxp_ROIori_column,~] = size(ori_V_exp_ROI);
    
    freq_stepsize = 0.1;
    freq_offsets = [min(freq_offsets) : freq_stepsize : max(freq_offsets)]';

    V_exp_ROI = zeros(V_wxp_ROIori_row, V_wxp_ROIori_column, length(freq_offsets));
    for i = 1 : V_wxp_ROIori_row
        for j = 1 : V_wxp_ROIori_column
            V_exp_ROI(i,j,:) = spline(ori_freq_offsets,ori_V_exp_ROI(i,j,:),freq_offsets);
%             v_exp_mask_inter(i,j,:) = spline(ori_freq_offsets,v_exp_mask(i,j,:),freq_offsets);
        end
    end
end
[V_wxp_ROI_row,V_wxp_ROI_column,V_wxp_ROI_numoffset] = size(V_exp_ROI);
%
offset_l1 = -2;
offset_l2 = -6.25;
offset_l3 = -10;
offset_h1 = 2;
offset_h2 = 6.25;
offset_h3 = 10;

% offset_l1 = -2;
% offset_l2 = -4.5;
% offset_l3 = -6.7;
% offset_h1 = 2;
% offset_h2 = 4.5;
% offset_h3 = 5.1;

[~, ind_neg_10] = min(abs(freq_offsets - offset_l3));  
[~, ind_neg_625] = min(abs(freq_offsets - offset_l2));  
[~, ind_neg_05] = min(abs(freq_offsets - offset_l1));  
[~, ind_pos_10] = min(abs(freq_offsets - offset_h3));  
[~, ind_pos_625] = min(abs(freq_offsets - offset_h2));  
[~, ind_pos_05] = min(abs(freq_offsets - offset_h1));
%% perform Lorentzian fitting  
options = optimset('MaxFunEvals', 1000000, 'TolFun', 1e-10, 'TolX', 1e-10, 'Display',  'off' );
% options.Algorithm = 'trust-region-reflective';%'levenberg-marquardt';
options.Algorithm = 'levenberg-marquardt';
x_data = freq_offsets([ind_neg_10:ind_neg_625, ind_neg_05:ind_pos_05, ind_pos_625:ind_pos_10]);
Lorentzian_fit = zeros(V_wxp_ROI_row,V_wxp_ROI_column,V_wxp_ROI_numoffset);
Lorentzian_difference = zeros(V_wxp_ROI_row,V_wxp_ROI_column,V_wxp_ROI_numoffset);
% [r, c]= find(reshape(V_exp_ROI(:,:,end),[V_wxp_ROI_row,V_wxp_ROI_column]));  
numIters = zeros(V_wxp_ROI_row,V_wxp_ROI_column);
for i = 1 : V_wxp_ROI_row
% for indx=drange(1: length(r))
    for j = 1 : V_wxp_ROI_column
        single_voxel_Zspec = squeeze(V_exp_ROI(i,j,:));
        y_data = single_voxel_Zspec([ind_neg_10:ind_neg_625, ind_neg_05:ind_pos_05, ind_pos_625:ind_pos_10]);
        xall_da = freq_offsets;
        x0 = 0;
        lb = [];
        ub = [];
        par0 = [x0,  50, 20, 1]; 
        [par,~,~,~,output] = lsqcurvefit(@lorentz_N, par0, x_data, y_data, lb, ub, options);
        numIters(i,j) = output.iterations;
        Lorentzian_fit(i,j,:) = lorentz_N(par, xall_da);
%         Lorentzian_difference(i,j,:) = Lorentzian_fit(r(indx),c(indx),:) -  V_exp_ROI(r(indx),c(indx),:); 
       Lorentzian_difference(i,j,:) = Lorentzian_fit(i,j,:) -  V_exp_ROI(i,j,:); 
       Lorentzian_differenceRex(i,j,:) =1./(V_exp_ROI(i,j,:)+1e-5)-1./(Lorentzian_fit(i,j,:)+1e-5);
%       Lorentzian_differenceRex(i,j,:) = 1./(v_exp_mask_inter(i,j,:)) - 1./(Lorentzian_fit(i,j,:));
    end
end
end
function y_fit = lorentz_N(par, x)
    denum = 1+4*((x-par(1))./par(3)).^2;
    y_fit=par(4)+par(2)./denum;
end

