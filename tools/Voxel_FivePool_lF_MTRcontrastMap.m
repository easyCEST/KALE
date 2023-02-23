function [ampMaps, areaMaps, fwhmMaps, offsetMaps, spectrumLorFit,poolNamesCellArr] = ...
          Voxel_FivePool_lF_MTRcontrastMap(cestNormData,freq_offsets)
%

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
[cestNormData_row,cestNormData_column,cestNormData_numoffset] = size(cestNormData);
            
%% Fit parameters: [p1 p2 p3], where p1 - [0 - 1] amplitude, p2 - [ppm] FWHM,p3 - [ppm] freq offset with respect to water at 0
% % Amide
% poolNamesCellArr{1}='Amide pool';
% par1=      [0.025    0.5   3.5]; % initial guess
% lb1=       [0        0.4    3]; % lower limits
% ub1=       [0.2      3  4]; % upper limits
% 
% % NOE
% poolNamesCellArr{2}='NOE pool';
% par2=       [0.02  3       -3.5];
% lb2=        [0     1       -4.5];
% ub2=        [1     5       -2];

% ub2=        [1   5          2];
% 
% % Water
% poolNamesCellArr{3}='Water pool';
% par3=       [0.9     1.4     0];
% lb3=        [0.02  0.3      -1];
% ub3=        [1     10        1];
% 
% % MT
% poolNamesCellArr{4}='MT pool';
% par4=       [0.1       25      0];
% lb4=        [0        10       -4];
% ub4=        [1       100        4];
% 
% % Amine
% poolNamesCellArr{5}='Amine pool';
% par5=       [0.05      1     2];
% lb5=        [0       0.3   1.7];
% ub5=        [0.3     2.5   2.2];
% 
% % par5=       [0.01    1.5   2];
% % lb5=        [0       0.5   1];
% % ub5=        [0.2     5     3];
% % putting all parameters together
% par0=[par1 par2 par3 par4 par5];
% lb=[lb1 lb2 lb3 lb4 lb5];
% ub=[ub1 ub2 ub3 ub4 ub5];
%% fitting parameters of amplitude and linewidth is 1/10 to 10 times,the frequency offsets is +- 20%
% Amide
ampFactor = 10;
linewidthFactor = 2;
offsetFactor = 0.1;


poolNamesCellArr{1}='Amide pool';
% par1=      [0.025    0.5    3.5]; % initial guess
par1=      [0.0621    1.49    3.63]; % initial guess of in NMR

lb1=       [par1(1)/ampFactor        par1(2)/linewidthFactor    par1(3)-par1(2)*offsetFactor]; % lower limits
ub1=       [par1(1)*ampFactor        par1(2)*linewidthFactor    par1(3)+par1(2)*offsetFactor]; % upper limits
% lb1=       [0        0.4    3];       % lower limits
% ub1=       [0.2      3      4];       % upper limits
% NOE
poolNamesCellArr{2}='NOE pool';
% par2=      [0.02  3       -3.5];
par2=      [0.1505  4.28       -3.25];     % initial guess of in NMR
lb2=       [par2(1)/ampFactor        par2(2)/linewidthFactor    par2(3)-par2(2)*offsetFactor]; % lower limits
ub2=       [par2(1)*ampFactor        par2(2)*linewidthFactor    par2(3)+par2(2)*offsetFactor]; % upper limits

% lb2=        [0     1       -4.5];
% ub2=        [1     5       -2];

% Water
poolNamesCellArr{3}='Water pool';
% par3=      [0.9     1.4     0];
par3=      [0.71     1.35     0.02];   % initial guess of in NMR
lb3=       [par3(1)/ampFactor        par3(2)/linewidthFactor    par3(3)-par3(2)*offsetFactor]; % lower limits
ub3=       [par3(1)*ampFactor        par3(2)*linewidthFactor    par3(3)+par3(2)*offsetFactor]; % upper limits

% lb3=        [0.02  0.3      -1];
% ub3=        [1     10        1];


% MT
poolNamesCellArr{4}='MT pool';
% par4=       [0.1       25       0];
par4=      [0.203     27.43     -1.48];   % initial guess of in NMR
lb4=       [par4(1)/ampFactor        par4(2)/linewidthFactor    par4(3)-par4(2)*offsetFactor]; % lower limits
ub4=       [par4(1)*ampFactor        par4(2)*linewidthFactor    par4(3)+par4(2)*offsetFactor]; % upper limits


% lb4=        [0        10       -4];
% ub4=        [1       100        4];
% Amine
poolNamesCellArr{5}='Amine pool';
% par5=      [0.05      1     2];
par5=      [0.0951     2.23     2.04];   % initial guess of in NMR
lb5=       [par5(1)/ampFactor        par5(2)/linewidthFactor    par5(3)-par5(2)*offsetFactor]; % lower limits
ub5=       [par5(1)*ampFactor        par5(2)*linewidthFactor    par5(3)+par5(2)*offsetFactor]; % upper limits
% lb5=        [0       0.3   1.7];
% ub5=        [0.3     2.5   2.2];

% putting all parameters together
par0=[par1 par2 par3 par4 par5];
lb=[lb1 lb2 lb3 lb4 lb5];
ub=[ub1 ub2 ub3 ub4 ub5];

%%
Npools=numel(poolNamesCellArr);
ampMaps    =zeros(cestNormData_row,cestNormData_column,Npools);
areaMaps   =zeros(cestNormData_row,cestNormData_column,Npools);
fwhmMaps   =zeros(cestNormData_row,cestNormData_column,Npools);
offsetMaps =zeros(cestNormData_row,cestNormData_column,Npools);

% [r, c]= find(reshape(cestNormData(:,:,end),[cestNormData_row,cestNormData_column]));

% for indx=drange(1: length(r))
for i = 1 : cestNormData_row
    for j = 1 : cestNormData_column
    spectrum=1-reshape(cestNormData(i,j,:),[cestNormData_numoffset, 1]);
    options=optimset('MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-10,'TolX',1e-10,  'Display',  'off' );
    options.Algorithm = 'levenberg-marquardt';
    fitPar=lsqcurvefit(@NpoolLorFitting,par0, freq_offsets(:),spectrum(:),lb,ub,options);
    spectrumLorFit(i,j,:) = 1-NpoolLorFitting(fitPar,freq_offsets);
    temp=1;
    for poolNo=1:Npools
       ampMaps(i,j,poolNo)=fitPar(temp);
       areaMaps(i,j,poolNo)=fitPar(temp).*fitPar(temp+1).*pi;
       fwhmMaps(i,j,poolNo)=fitPar(temp+1);
       offsetMaps(i,j,poolNo)=fitPar(temp+2);
       temp=temp+3;
    end
    end
end

% for i = 1:length(poolNamesCellArr)
%     figure
%     imagesc(ampMaps(:,:,i), [0,0.07])
%     colorbar;axis off
%     colormap(gca, jet);
%     str = poolNamesCellArr{i};
%     title(str,'FontWeight','bold','FontSize',18);
% end
end

function fittedSpectrum = NpoolLorFitting(fitPar, offsets)
% fittedSpectrum = NpoolLorFitting(par, offsets)
%
% Implementation of an N-pool Lorentzian-lineshape model

count=1;
Npools=numel(fitPar)/3;
tempVar=zeros(numel(offsets),Npools);
for ii=1:Npools
    tempVar(:,ii)=fitPar(count)./(1+4 * (((offsets-fitPar(count+2))./fitPar(count+1)).^2+1e-06));
    count=count+3;
    tempVar(isnan(tempVar)) = 0;
end
fittedSpectrum=sum(tempVar,2);
end
