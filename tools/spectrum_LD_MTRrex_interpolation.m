function [Lorentzian_fit, Lorentzian_difference, AREX_Lorentzian_difference,freq_offsets,Iters]...
                            = spectrum_LD_MTRrex_interpolation(Z_spectrums, freq_offsets,par0,iter)
%%
%% 
do_interpolation = 1;

if do_interpolation
    ori_freq_offsets = freq_offsets;
    ori_Z_spectrums = Z_spectrums;
    
    freq_stepsize = 0.1;
    freq_offsets = [min(freq_offsets) : freq_stepsize : max(freq_offsets)]';
    
    Z_spectrums = zeros(length(freq_offsets), size(ori_Z_spectrums, 2));
    for i = 1:size(ori_Z_spectrums, 2)
        Z_spectrums(:, i) = spline(ori_freq_offsets, ori_Z_spectrums(:, i), freq_offsets);
    end
end

[nr_freq_offsets, nr_spectrum] = size(Z_spectrums);

%%

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
[~, ind_pos_05] = min(abs(freq_offsets - offset_h1));
[~, ind_neg_05] = min(abs(freq_offsets - offset_l1));  
[~, ind_pos_10] = min(abs(freq_offsets - offset_h3));  
[~, ind_pos_625] = min(abs(freq_offsets - offset_h2));  


%% perform Lorentzian fitting  
% cond_x = 1/(10^iter);
options = optimset('MaxFunEvals', 1000000, 'TolFun', 1e-10, 'TolX',1e-10 , 'Display',  'off' );
% options.Algorithm = 'trust-region-reflective';%'levenberg-marquardt';
options.Algorithm = 'levenberg-marquardt';

xdata = freq_offsets([ind_neg_10:ind_neg_625, ind_neg_05:ind_pos_05, ind_pos_625:ind_pos_10]);

Lorentzian_fit = zeros(nr_freq_offsets, nr_spectrum);
Lorentzian_difference = zeros(nr_freq_offsets, nr_spectrum);
Iters = 0;
for i = 1:nr_spectrum
    ydata = Z_spectrums([ind_neg_10:ind_neg_625, ind_neg_05:ind_pos_05, ind_pos_625:ind_pos_10], i);

    xall_da = freq_offsets;
%     x0 = 0;
    lb = [];   
    ub = [];
%     par0 = [x0,  50, 20, 1]; 
    [par,~,~,~,output] = lsqcurvefit(@lorentz_N, par0, xdata, ydata, lb, ub, options);
    Iters = output.iterations;
    Lorentzian_fit(:, i) = lorentz_N(par, xall_da);

    %求的拟合后和原始数据的差值信号
    Lorentzian_difference(:, i) = Lorentzian_fit(:, i) -  Z_spectrums(:, i); 
    AREX_Lorentzian_difference(:,i) = 1./Z_spectrums(:, i) - 1./Lorentzian_fit(:, i);
end

end

function y_fit = lorentz_N(par, x)
%     denum = 1 + (par(3) ./ (x - par(1))).^2;
%     y_fit = par(4) + par(2) ./ denum;
    
    denum = 1+4*((x-par(1))./par(3)).^2;
    y_fit=par(4)+par(2)./denum;
end
