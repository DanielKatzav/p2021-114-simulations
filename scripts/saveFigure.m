function [fig] = saveFigure(data,fig_title,name_of_plane,graphs)
if graphs
    figureToSave = figure;
    imagesc(data);
    title(fig_title)
    colorbar();
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'),name_of_plane, ".jpg"));
    saveas(figureToSave, figFileName)
end
end