function [mSigE,w_offset,mSigEStd] = ROIprocess_Zspecforkmeans(ROI_total,w_offset,V_exp)

%%
% insert_zero = 0;
% if insert_zero
%     if min(abs(w1)) == 0
%         w_offsetOri = [w_offset];
%     else       
%         w_offsetOri = [w_offset 0];
%     end 
%     V_exp_insertZero = cat(3, V_exp, zeros(size(V_exp,1),size(V_exp,2)));
%     Nw = length(w_offsetOri); 
%     for in = 1 : Nw
%         Vo_tmp(:,:,in) = medfilt2(V_exp_insertZero(:,:,in),[2 2]);
%     end
%     [w_offset, I] = sort(w_offsetOri);
%     V_exp = Vo_tmp(:,:,I);
% end

%% 
[x,y,z] = size(V_exp);

if length(w_offset) ~= z
    w_offset = [min(w_offset):0.1:max(w_offset)];  
end
Nw = length(w_offset);
nROI = length(ROI_total);
str=[num2str(nROI),'ROI'];  
mSigEStd=zeros(z,nROI);
for k=1:nROI
    for n_offset=1:Nw   
        temp = squeeze(V_exp(:,:,n_offset));
        mSigE(n_offset,k) = mean2(temp(ROI_total{k}));
        mSigEStd(n_offset,k) = std(temp(ROI_total{k}));
    end
end
% if strfind(computer, 'WIN')
%     xlswrite( fullfile(ROIpath, ['Zspc_',str,'.xls']),w_offset',N_echo,['B1']);
%     xlswrite( fullfile(ROIpath, ['Zspc_',str,'.xls']),mSigE,N_echo,['C1']); 
%     xlswrite( fullfile(ROIpath, ['Zspc_std',str,'.xls']),w_offset',N_echo,['B1']); 
%     xlswrite( fullfile(ROIpath, ['Zspc_std',str,'.xls']),eSigE,N_echo,['C1']);
%     save(ROIpath,'mSigE');
% elseif strfind(computer, 'MAC')    
%     csvwrite( fullfile(ROIpath, ['Zspc_',str,'.csv']),mSigE,0,0);
% end
end
