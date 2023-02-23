function [V_exp_Pyramid] = cestDataDowndapling(V_exp,w_offset)
% The data is subsampled to the corresponding resolution

%% down sampling layers is 4. it is used at method of Kmeans_Pyramid
V_exp_Pyramid = {};
num = length(w_offset);
w = [0.25,0.25;0.25,0.25];
for i = 1 : num
    Single_data = squeeze(V_exp(:,:,i));        % 256x256     %   64
    dsample2 = Gaussian_dsample(Single_data,2); % 128X128     %   32
    dsample4 = Gaussian_dsample(dsample2,2);    % 64x64       %   16  
    dsample8 = Gaussian_dsample(dsample4,2);    % 32X32       %   8
    dsample16 = Gaussian_dsample(dsample8,2);   % 16X16
%     dsample32 = Gaussian_dsample(dsample16,2);  % 8x8
    
    
    
    Pyramid2(:,:,i) = dsample2;
    Pyramid4(:,:,i) = dsample4;
    Pyramid8(:,:,i) = dsample8;
    Pyramid16(:,:,i) = dsample16;
%     Pyramid32(:,:,i) = dsample32;
    
    
end
V_exp_Pyramid = {Pyramid16,Pyramid8,Pyramid4,Pyramid2,V_exp}; %hunman
% V_exp_Pyramid = {Pyramid4,Pyramid2,V_exp};      % rat
end
%% 
function Idown = Gaussian_dsample(I,N) 

W = fspecial('gaussian',[3,3],1); 
G = imfilter(I, W, 'replicate');
[row,col] = size(G); 
drow = round(row/N); 
dcol = round(col/N); 
Idown = zeros(drow,dcol); 
p =1; 
q =1; 
for i = 1:N:row 
    for j = 1:N:col 
         Idown(p,q) = I(i,j); 
         q = q+1; 
    end 
    q =1; 
    p = p+1; 
end
end
%%

function img_conv = myConv(I,step,w)
    [wx,wy] = size(w);
    [Ix,Iy] = size(I);
    img_conv = [];
    for i = 1 : step :Ix-wx+1
        line = [];
        for j = 1 : step : Iy-wy+1
            line = I(i:i+wx-1,j:j+wy-1);
            line_conv = sum(sum(line.*w))/(wx*wy);
            img_conv = [img_conv,line_conv];
        end
    end
    img_conv = (reshape(img_conv,[sqrt(length(img_conv)),sqrt(length(img_conv))]))';
end








