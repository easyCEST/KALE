function [V_exp_mask]=prepare(S0_use,V_exp,brainMask) 

Nw = size(V_exp,3);
Row=size(V_exp,1);        
Column=size(V_exp,2); 
V_exp_mask=zeros(Row,Column,Nw);
for n_offset=1:Nw
          m1=V_exp(:,:,n_offset)./(S0_use+1e-5);
%           m1=V_exp(:,:,n_offset);
          m1_ROI= m1.*brainMask;
%             m1_ROI= m1;

          V_exp_mask(:,:,n_offset)=m1_ROI;
end
V_exp_mask(isnan(V_exp_mask)) = 0;
% % %   save(fullfile(ROIpath,'V_exp_ROI'),'V_exp_ROI')    
% % %   save(fullfile(ROIpath,'w_offset'),'w_offset')      