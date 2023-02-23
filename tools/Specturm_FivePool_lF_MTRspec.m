function [Z_spectrums, spectrumLorFit,spectrumLorFitPerPool,poolNamesCellArr] = ...
          Specturm_FivePool_lF_MTRspec(Z_spectrums,freq_offsets)

      
%%

colors = distinguishable_colors(50);
colors = colors(:,:);

%%
do_interpolation = 1;
if do_interpolation
    ori_freq_offsets = freq_offsets;
    ori_Z_spectrums = Z_spectrums;
    
    freq_stepsize = 0.1;
    freq_offsets = [min(freq_offsets) : freq_stepsize : max(freq_offsets)]';

    Z_spectrums = zeros(length(freq_offsets), size(ori_Z_spectrums, 2));
    for i = 1 : size(Z_spectrums, 2)
         Z_spectrums(:, i) = spline(ori_freq_offsets, ori_Z_spectrums(:, i), freq_offsets);
    end
end
[N_offsets, N_spectrums] = size(Z_spectrums);
            
%% Fit parameters: [p1 p2 p3], where p1 - [0 - 1] amplitude, p2 - [ppm] FWHM,p3 - [ppm] freq offset with respect to water at 0
% Amide
poolNamesCellArr{1}='Amide pool';
par1=      [0.025    0.5   3.5]; % initial guess
lb1=       [0        0.4    3]; % lower limits
ub1=       [0.2      3  4]; % upper limits

% NOE
poolNamesCellArr{2}='NOE pool';
par2=       [0.02    3       -3.5];
lb2=        [0     1     -4.5];
ub2=        [1   5     -2];

% Water
poolNamesCellArr{3}='Water pool';
par3=       [0.9     1.4     0];
lb3=        [0.02  0.3     -1];
ub3=        [1     10      1];

% MT
poolNamesCellArr{4}='MT pool';
par4=       [0.1       25     0];
lb4=        [0        10       -4];
ub4=        [1       100        4];

% Amine
poolNamesCellArr{5}='Amine pool';
par5=       [0.050     1     2];
lb5=        [0       0.3   1.7];
ub5=        [0.3     2.5   2.2];

% putting all parameters together
par0=[par1 par2 par3 par4,par5];
lb=[lb1 lb2 lb3 lb4,lb5];
ub=[ub1 ub2 ub3 ub4,ub5];
%%
Npools=numel(poolNamesCellArr);
spectrumLorFit=zeros(N_offsets,N_spectrums);
spectrumLorFitPerPool=zeros(N_offsets,N_spectrums,Npools);

for SpectrumsNo=1:N_spectrums
    spectrum=1-Z_spectrums(:,SpectrumsNo);
    options = optimset('MaxFunEvals', 10000, 'TolFun', 1e-10, 'TolX', 1e-10, 'Display',  'off' );
    options.Algorithm = 'levenberg-marquardt';
    fitPar=lsqcurvefit(@NpoolLorFitting,par0, freq_offsets(:),spectrum(:),lb,ub,options);
    spectrumLorFit(:,SpectrumsNo)=1-NpoolLorFitting(fitPar, freq_offsets);
    
    temp=1;
    for poolNo=1:Npools
        parOnlyCurrentPool=zeros(size(fitPar));
        parOnlyCurrentPool(temp:temp+2)=fitPar(temp:temp+2);
        spectrumLorFitPerPool(:,SpectrumsNo,poolNo)=NpoolLorFitting(parOnlyCurrentPool, freq_offsets);
        temp=temp+3;
    end
end

for SpectrumsNo=1:N_spectrums
    
    figName=strcat('CESTandLorFit-ROI-',num2str(SpectrumsNo));
    h=figure;
    set(gcf,'name',sprintf('%s',figName),'numbertitle','off')
    lgnd(1)=plot(freq_offsets, Z_spectrums(:, SpectrumsNo), 'b.','MarkerSize',10, 'DisplayName','Data');
    hold on, grid on
%     errorbar(freq_offsets, Z_spectrums(:, SpectrumsNo),'.b','LineWidth',1)
    lgnd(2)=plot(freq_offsets, spectrumLorFit(:, SpectrumsNo), 'k','LineWidth',1, 'DisplayName','full Lor fit');
    ylabel('M_z/M_0')
    xlabel('Frequency offset [ppm]')
%     axis(XYscaleSpectrum)
    set(gca, 'Xdir', 'reverse')
    set(gca,'FontSize',14)
    
    for poolNo=1:Npools
        lgnd(poolNo+2)=plot(freq_offsets, spectrumLorFitPerPool(:, SpectrumsNo,poolNo),'-','Color',colors(poolNo,:),'LineWidth',2,'DisplayName',poolNamesCellArr{poolNo});
    end
    legend(lgnd);
    
    %saveas(h,figName,'fig')
end

end

function fittedSpectrum = NpoolLorFitting(fitPar, offsets)
% fittedSpectrum = NpoolLorFitting(par, offsets)
%
% Implementation of an N-pool Lorentzian-lineshape model

count=1;
Npools=numel(fitPar)/3;
tempVar=zeros(numel(offsets),Npools);
for ii=1:Npools
    tempVar(:,ii)=fitPar(count)./(1+((offsets-fitPar(count+2))./fitPar(count+1)).^2+1e-06);
    count=count+3;
    tempVar(isnan(tempVar)) = 0;
end
fittedSpectrum=sum(tempVar,2);
end
