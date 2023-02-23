function altered_colormap()
cmap = colormap;
nzeros = 2;
cmap(1:nzeros,:) = zeros(nzeros,3);
colormap(cmap); 