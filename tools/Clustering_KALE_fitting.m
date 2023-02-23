function [cls,loli_vec_zip,index_class,index1] = Clustering_KALE_fitting(loli_vec,K,method,w_offset,offset_choice,oriindex)

if isa(loli_vec, 'cell') 
    loli_vec = cell2mat(loli_vec);
end
if isa(oriindex, 'cell') 
    oriindex = cell2mat(oriindex);
end
channel = size(loli_vec,2);
for i = 1:length(offset_choice)
   [~,index] = min(abs(w_offset -offset_choice(i)));
   index1(i,1) = index;
end
loli_vec_choice = loli_vec(:,index1);
centers = zeros(K,channel);
if (strcmp(method,'kmeans'))
    [~,cls] = k_means_2D(double(loli_vec_choice),[],K,[] );
elseif (strcmp(method,'kmeans++') || strcmp(method,'kmeanspp')) 
    [centroid,cls] = Kmeanspp(double(loli_vec_choice),K,10);
end
% diff = zeros(K,channel);
for k = 1:K
    centers(k,:)=mean(loli_vec(find(cls==k),:));
%     diff(k,:) = loli_vec(find(cls==k),:)-centers(k,:);
    index_class{k} = oriindex(find(cls == k));
end
centers(isnan(centers)) = 0;
% diff(isnan(diff)) = 0;
loli_vec_zip = centers(cls,:);

