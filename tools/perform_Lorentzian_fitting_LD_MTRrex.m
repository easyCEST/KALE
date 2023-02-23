function [R2,figh, Lorentzian_fit, Lorentzian_difference, AREX_Lorentzian_difference,freq_offsets]...
                            = perform_Lorentzian_fitting_LD_MTRrex(Z_spectrums, freq_offsets,cond_DS)

%%
do_adjustment = 0;

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

offset_l1 = -cond_DS;
offset_l2 = -8;
offset_l3 = -10;
offset_h1 = cond_DS;
offset_h2 = 8;
offset_h3 = 10;

% offset_l1 = -2;
% offset_l2 = -4.25;
% offset_l3 = -6;
% offset_h1 = 2;
% offset_h2 = 4.25;
% offset_h3 = 6;

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

xdata = freq_offsets([ind_neg_10:ind_neg_625, ind_neg_05:ind_pos_05, ind_pos_625:ind_pos_10]);

Lorentzian_fit = zeros(nr_freq_offsets, nr_spectrum);
Lorentzian_difference = zeros(nr_freq_offsets, nr_spectrum);
for i = 1:nr_spectrum
    ydata = Z_spectrums([ind_neg_10:ind_neg_625, ind_neg_05:ind_pos_05, ind_pos_625:ind_pos_10], i);

    xall_da = freq_offsets;
    x0 = 0;
    lb = [];   
    ub = [];
    par0 = [x0,  50, 20, 1]; 
    par = lsqcurvefit(@lorentz_N, par0, xdata, ydata, lb, ub, options);
    
    Lorentzian_fit(:, i) = lorentz_N(par, xall_da);
    if do_adjustment
         [~, ind_neg_01] = min(abs(freq_offsets + 0.1));  
         [~, ind_pos_01] = min(abs(freq_offsets - 0.1));  
         Lorentzian_fit([ind_neg_10:ind_neg_625], i) = Z_spectrums([ind_neg_10:ind_neg_625], i);
         Lorentzian_fit([ind_neg_01:ind_pos_01], i) = Z_spectrums([ind_neg_01:ind_pos_01], i);
         Lorentzian_fit([ind_pos_625:ind_pos_10], i) = Z_spectrums([ind_pos_625:ind_pos_10], i);
         
        
        ad_x = freq_offsets;
        ad_y = Lorentzian_fit;
        lb = [];   
        ub = [];
        par0 = par; 
        ad_par = lsqcurvefit(@lorentz_N, par0, ad_x, ad_y, lb, ub, options);
        
        Lorentzian_fit(:, i) = lorentz_N(ad_par, ad_x);
    end
    
    w_offset_n1 = -1.5;
    w_offset_p1 = 1.5;
    [~, ind_neg_1] = min(abs(freq_offsets - w_offset_n1));  
    [~, ind_pos_1] = min(abs(freq_offsets - w_offset_p1));
    Z_sub = Z_spectrums([ind_neg_1:ind_pos_1],i);
    LF_sub = Lorentzian_fit([ind_neg_1:ind_pos_1],i);
    R2(i) = 1 - sum((Z_sub(:,i)-LF_sub(:,i)).^2)/sum((Z_sub(:,i)-mean(Z_sub(:,i))).^2);
    %求的拟合后和原始数据的差值信号
    Lorentzian_difference(:, i) = Lorentzian_fit(:, i) -  Z_spectrums(:, i); 
    AREX_Lorentzian_difference(:,i) = 1./Z_spectrums(:, i) - 1./Lorentzian_fit(:, i);
end

figure
cla
hold on
% plot(freq_offsets, [Z_spectrums, Lorentzian_fit, Lorentzian_difference,AREX_Lorentzian_difference], '.-','linewidth',1.5);
plot(freq_offsets, [Z_spectrums, Lorentzian_fit], '.-','linewidth',1.5);
grid on
set(gca, 'Xdir', 'reverse', 'FontWeight','bold','FontSize',14);
legend('Z\_spectrum', 'Lorentzian fitting', 'Fitting difference','AREX Fitting difference')
axis([-10 10 0 1])
title(num2str(offset_h1))
figh = gcf;

end

function y_fit = lorentz_N(par, x)
    denum = 1+4*((x-par(1))./par(3)).^2;
    y_fit=par(4)+par(2)./denum;
end

function y_fit = Gaussian_N(par, x)

    y_fit=par(3)+par(2)*(exp(-(x-par(1)).^2).*(pi*(par(2)).^2));

end

function do_adjustment_fit(y_fit,z_data)

    


end