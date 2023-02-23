function display_image(name,brainMask,MTRname,displayrange,S0B0path)
% load(fullfile(Savepath,'brainMask.mat'));
% index_ROI = find(brainMask < 1);
% name(index_ROI) = -0.5;
[px,py]=find(brainMask);
temp = 6;
figure
%             imagesc(name,displayrange) 
            imagesc(name(min(px)-temp:max(px)+temp,min(py)-temp:max(py)+temp),displayrange);
            colormap(jet(256))
            colorbar
            altered_colormap()
            axis off
            set(gca, 'FontWeight','bold','FontSize',20)
            title(MTRname,'FontWeight','bold','FontSize',14)
%             savefig(fullfile(S0B0path, [MTRname, '.fig']))
%     function [str_name] = getvarname(var)
%         str_name = sprintf('%s',inputname(1))
%     end

end
