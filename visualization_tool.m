% load('tiff_results.mat')
% load('trace_seeds.mat')
function visualization_tool(image, points)
    for ii = 1:1:size(image, 3)
        % grab coords 
        pl = find(round(points(:,3))==ii);
        coord = points(pl,1:2);
        figure(1)
        imshow(image(:,:,ii))
        hold on
        plot(coord(:,2),coord(:,1),'ro')
        title(['slice = ' num2str(ii)])
        hold off
        pause
    end
end