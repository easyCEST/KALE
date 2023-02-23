function [Lorentzian_fit,Lor_fit,Lorentzian_diff,w_offset_1,par_fit,AREX_Lorentzian_difference] = ...
                        LF_Poly(w_offset,mSigE,Voxel_Index,l1ppm,l2ppm,l3ppm,h1ppm, h2ppm, h3ppm,par0)

Nw = length(w_offset);        
w_offset_1=[min(w_offset):0.1:max(w_offset)]'; 
[mm, ~]=size(w_offset_1);
%%
% adjust l1ppm and h1ppm by x_max , x_min , x_i and par0(3).
[~,index_0] = min(abs(w_offset - 0));
mSigE_1 = mSigE([1:index_0-1,index_0+1:end],:);

[~,r10]=min(abs(w_offset_1-l3ppm)); 
[~,r625]=min(abs(w_offset_1-l2ppm));  
[~,l10]=min(abs(w_offset_1-h3ppm));  
[~,l625]=min(abs(w_offset_1-h2ppm));  
[~,r05]=min(abs(w_offset_1-l1ppm));  
[~,l05]=min(abs(w_offset_1-h1ppm)); 
%%
nROI = length(Voxel_Index);
Lorentzian_fit=zeros(mm,nROI);
Z=zeros(Nw,1) ;
ZL=zeros(mm,nROI) ;
mSigE(isnan(mSigE))=0;
for k = 1 : nROI %%%%%%%%%%需要画的ROI的个数
    Z(:,k)=mSigE(:,k);
    ZL(:,k)=spline(w_offset,Z(:,k),w_offset_1);
    options = optimset('MaxFunEvals',1000000,'TolFun',1e-10,'TolX',1e-10, 'Display',  'off' );
    options.Algorithm = 'levenberg-marquardt';
            


    % % 2.the value of data for each ppm in 1.
    toFit = ZL([r10:r625 r05:l05 l625:l10],k);
    % 3.ppm to fit
    xdata = w_offset_1([r10:r625 r05:l05 l625:l10]);
    ydata=toFit;
    xall_da= w_offset_1;
    lb = [];   
    ub = [];
    par = lsqcurvefit(@lorentz_N,par0, xdata, ydata, lb, ub, options);
    Lorentzian_fit(:,k) = lorentz_N(par, xall_da);
    par_fit(k,:) = par;
    Lorentzian_diff(:,k) =  Lorentzian_fit(:,k) - ZL(:,k);
    AREX_Lorentzian_difference(:,k) = 1./(ZL(:, k)+1e-5) - 1./(Lorentzian_fit(:, k)+1e-5);
end
for i = 1:length(w_offset)
   [~,index] = min(abs(w_offset_1 -w_offset(i)));
   index1(i,1) = index;
end
Lor_fit = Lorentzian_fit(index1,:);
end
        

function y_fit = lorentz_N(par, x)
    denum = 1+4*((x-par(1))./par(3)).^2;
    y_fit=par(4)+par(2)./denum;
end
