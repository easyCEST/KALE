function [Z_spectrums, spectrumLorFit,spectrumLorFitPerPool,poolNamesCellArr,par_fit] = ...
          Kmeans_FivePool_lF_MTRspec(Z_spectrums,freq_offsets,par0,iteration,Mean_mSigEStd)

      
%%
Z_spectrums(isnan(Z_spectrums)) = 0;

colors = distinguishable_colors(50);
colors = colors(:,:);
alpha = 0.1;

poolNamesCellArr{1}='Amide pool';
poolNamesCellArr{2}='NOE pool';
poolNamesCellArr{3}='Water pool';
poolNamesCellArr{4}='MT pool';
poolNamesCellArr{5}='Amine pool';


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
% % lamba = 0.1;
% if iteration == 0
% % Amide
% poolNamesCellArr{1}='Amide pool';
% % par1=      [0.025    0.5   3.5]; % initial guess
% lb1=       [0        0.4    3]; % lower limits
% ub1=       [0.2      3  4]; % upper limits
% 
% % NOE
% poolNamesCellArr{2}='NOE pool';
% % par2=      [0.02    3       -3.5];
% lb2=        [0     1     -4.5];
% ub2=        [1   5         -2];
% 
% % Water
% poolNamesCellArr{3}='Water pool';
% % par3=      [0.9     1.4     0];
% lb3=        [0.02  0.3     -1];
% ub3=        [1     10       1];
% 
% % MT
% poolNamesCellArr{4}='MT pool';
% % par4=      [0.1       25     0];
% lb4=        [0        10       -4];
% ub4=        [1       100        4];
% 
% % Amine
% poolNamesCellArr{5}='Amine pool';
% % par5=      [0.050     1     2];
% lb5=        [0       0.3   1.7];
% ub5=        [0.3     2.5   2.2];

%     
% else
% % % % Amide
% % poolNamesCellArr{1}='Amide pool';
% % % par1=      [0.025    0.5   3.5]; % initial guess
% % lb1=       [par0(1)*(1-lamba)    par0(2)*(1-lamba)     par0(3)*0.9 ]; % lower limits
% % ub1=       [par0(1)*(1+lamba)    par0(2)*(1+lamba)     par0(3)*1.1 ]; % upper limits
% % 
% % % NOE
% % poolNamesCellArr{2}='NOE pool';
% % % par2=       [0.02    3       -3.5];
% % lb2=       [par0(4)*(1-lamba)    par0(5)*(1-lamba)     par0(6)*0.9 ]; % lower limits
% % ub2=       [par0(4)*(1+lamba)    par0(5)*(1+lamba)     par0(6)*1.1 ]; % upper limits
% % 
% % % Water
% % poolNamesCellArr{3}='Water pool';
% % % par3=       [0.9     1.4     0];
% % lb3=       [par0(7)*(1-lamba)    par0(8)*(1-lamba)     par0(9)*0.9 ]; % lower limits
% % ub3=       [par0(7)*(1+lamba)    par0(8)*(1+lamba)     par0(9)*1.1 ]; % upper limits
% % 
% % % MT
% % poolNamesCellArr{4}='MT pool';
% % % par4=       [0.1       25     0];
% % lb4=       [par0(10)*(1-lamba)    par0(11)*(1-lamba)     par0(12)*0.9 ]; % lower limits
% % ub4=       [par0(10)*(1+lamba)    par0(11)*(1+lamba)     par0(12)*1.1 ]; % upper limits
% % 
% % % Amine
% % poolNamesCellArr{5}='Amine pool';
% % % par5=       [0.050     1     2];
% % lb5=       [par0(13)*(1-lamba)    par0(14)*(1-lamba)     par0(15)*0.9 ]; % lower limits
% % ub5=       [par0(13)*(1+lamba)    par0(14)*(1+lamba)     par0(15)*1.1 ]; % upper limits
% 
% % 
% % Amide
% % % lamba = 1 / iteration + 0.01;
% poolNamesCellArr{1}='Amide pool';
% % par1=      [0.025    0.5   3.5]; % initial guess
% lb1=       [par0(1)-lamba    par0(2)-lamba    par0(3)-lamba*0.5 ]; % lower limits
% ub1=       [par0(1)+lamba    par0(2)+lamba    par0(3)+lamba*0.5 ]; % upper limits
% 
% % NOE
% poolNamesCellArr{2}='NOE pool';
% % par2=       [0.02    3       -3.5];
% lb2=       [par0(4)-lamba    par0(5)-lamba     par0(6)-lamba*0.5 ]; % lower limits
% ub2=       [par0(4)+lamba    par0(5)+lamba     par0(6)+lamba*0.5 ]; % upper limits
% 
% % Water
% poolNamesCellArr{3}='Water pool';
% % par3=       [0.9     1.4     0];
% lb3=       [par0(7)-lamba   par0(8)-lamba     par0(9)-lamba*0.5 ]; % lower limits
% ub3=       [par0(7)+lamba   par0(8)+lamba     par0(9)+lamba*0.5 ]; % upper limits
% 
% % MT
% poolNamesCellArr{4}='MT pool';
% % par4=       [0.1       25     0];
% lb4=       [par0(10)-lamba    par0(11)-lamba     par0(12)-lamba*0.5 ]; % lower limits
% ub4=       [par0(10)+lamba    par0(11)+lamba     par0(12)+lamba*0.5 ]; % upper limits
% 
% % Amine
% poolNamesCellArr{5}='Amine pool';
% % par5=       [0.050     1     2];
% lb5=       [par0(13)-lamba    par0(14)-lamba     par0(15)-lamba*0.5 ]; % lower limits
% ub5=       [par0(13)+lamba    par0(14)+lamba     par0(15)+lamba*0.5 ]; % upper limits
% 
% % % % Amide
% % % lamba = 1 / iteration + 0.01;
% % poolNamesCellArr{1}='Amide pool';
% % % par1=      [0.025    0.5   3.5]; % initial guess
% % lb1=       [par0(1)-lamba    par0(2)*(1-alpha)    par0(3)-lamba*0.5 ]; % lower limits
% % ub1=       [par0(1)+lamba    par0(2)*(1+alpha)    par0(3)+lamba*0.5 ]; % upper limits
% % 
% % % NOE
% % poolNamesCellArr{2}='NOE pool';
% % % par2=       [0.02    3       -3.5];
% % lb2=       [par0(4)-lamba    par0(5)*(1-alpha)     par0(6)-lamba*0.5 ]; % lower limits
% % ub2=       [par0(4)+lamba    par0(5)*(1+alpha)     par0(6)+lamba*0.5 ]; % upper limits
% % 
% % % Water
% % poolNamesCellArr{3}='Water pool';
% % % par3=       [0.9     1.4     0];
% % lb3=       [par0(7)-lamba   par0(8)*(1-alpha)     par0(9)-lamba*0.5 ]; % lower limits
% % ub3=       [par0(7)+lamba   par0(8)*(1+alpha)     par0(9)+lamba*0.5 ]; % upper limits
% % 
% % % MT
% % poolNamesCellArr{4}='MT pool';
% % % par4=       [0.1       25     0];
% % lb4=       [par0(10)-lamba    par0(11)*(1-alpha)     par0(12)-lamba*0.5 ]; % lower limits
% % ub4=       [par0(10)+lamba    par0(11)*(1+alpha)     par0(12)+lamba*0.5 ]; % upper limits
% % 
% % % Amine
% % poolNamesCellArr{5}='Amine pool';
% % % par5=       [0.050     1     2];
% % lb5=       [par0(13)-lamba    par0(14)*(1-alpha)     par0(15)-lamba*0.5 ]; % lower limits
% % ub5=       [par0(13)+lamba    par0(14)*(1+alpha)     par0(15)+lamba*0.5 ]; % upper limits
% 
% 
% end
% % putting all parameters together
% % par0=[par1 par2 par3 par4,par5];
% lb=[lb1 lb2 lb3 lb4,lb5];
% ub=[ub1 ub2 ub3 ub4,ub5];
%%
Npools=numel(poolNamesCellArr);
spectrumLorFit=zeros(N_offsets,N_spectrums);
spectrumLorFitPerPool=zeros(N_offsets,N_spectrums,Npools);

for SpectrumsNo=1:N_spectrums
    spectrum=1-Z_spectrums(:,SpectrumsNo);
    [lb,ub] = setParameters(par0,Mean_mSigEStd(SpectrumsNo),iteration);
    options=optimset('MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-10,'TolX',1e-10,  'Display',  'off' );
    fitPar=lsqcurvefit(@NpoolLorFitting,par0, freq_offsets(:),spectrum(:),lb,ub,options);
    spectrumLorFit(:,SpectrumsNo)=1-NpoolLorFitting(fitPar, freq_offsets);
    par_fit(SpectrumsNo,:) = fitPar;
    temp=1;
    for poolNo=1:Npools
        parOnlyCurrentPool=zeros(size(fitPar));
        parOnlyCurrentPool(temp:temp+2)=fitPar(temp:temp+2);
        spectrumLorFitPerPool(:,SpectrumsNo,poolNo)=NpoolLorFitting(parOnlyCurrentPool, freq_offsets);
        temp=temp+3;
    end
end

end
function [lb,ub] = setParameters(par0,cluterStd,iteration)
% lamba = 0.02;
lamba = cluterStd * 10;
Bata = cluterStd * 100;
if iteration == 0
% Amide
lb1=       [0        0.4    3]; % lower limits
ub1=       [0.2      3      4]; % upper limits

% NOE
lb2=        [0     1     -4.5];
ub2=        [1   5         -2];

% Water
lb3=        [0.02  0.3     -1];
ub3=        [1     10       1];

% MT
lb4=        [0        10       -4];
ub4=        [1       100        4];

% Amine
lb5=        [0       0.3   1.7];
ub5=        [0.3     2.5   2.2];

    
else
% Amide
lb1=       [par0(1)-lamba    par0(2)-Bata    par0(3)-cluterStd*0.5 ]; % lower limits
ub1=       [par0(1)+lamba    par0(2)+Bata    par0(3)+cluterStd*0.5 ]; % upper limits

% NOE
lb2=       [par0(4)-lamba    par0(5)-Bata     par0(6)-cluterStd*0.5 ]; % lower limits
ub2=       [par0(4)+lamba    par0(5)+Bata     par0(6)+cluterStd*0.5 ]; % upper limits

% Water
lb3=       [par0(7)-lamba   par0(8)-Bata     par0(9)-cluterStd*0.5 ]; % lower limits
ub3=       [par0(7)+lamba   par0(8)+Bata     par0(9)+cluterStd*0.5 ]; % upper limits

% MT
lb4=       [par0(10)-lamba    par0(11)-Bata     par0(12)-cluterStd*0.5 ]; % lower limits
ub4=       [par0(10)+lamba    par0(11)+Bata     par0(12)+cluterStd*0.5 ]; % upper limits

% Amine
lb5=       [par0(13)-lamba    par0(14)-Bata     par0(15)-cluterStd*0.5 ]; % lower limits
ub5=       [par0(13)+lamba    par0(14)+Bata     par0(15)+cluterStd*0.5 ]; % upper limits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
% putting all parameters together
lb=[lb1 lb2 lb3 lb4,lb5];
ub=[ub1 ub2 ub3 ub4,ub5];
end
%%
function fittedSpectrum = NpoolLorFitting(fitPar, offsets)
% fittedSpectrum = NpoolLorFitting(par, offsets)
%
% Implementation of an N-pool Lorentzian-lineshape model

count=1;
Npools=numel(fitPar)/3;
tempVar=zeros(numel(offsets),Npools);
for ii=1:Npools
    tempVar(:,ii)=fitPar(count)./(1+ 4 * (((offsets-fitPar(count+2))./fitPar(count+1)).^2+1e-06));
    count=count+3;
    tempVar(isnan(tempVar)) = 0;
end
fittedSpectrum=sum(tempVar,2);
end
