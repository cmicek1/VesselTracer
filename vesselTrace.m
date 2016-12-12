function vesselTrace(FinalImage, seedCoord)

startEllipse = ellipsoidInit(0.25, 3);
x = startEllipse.XData;
y = startEllipse.YData;
z = startEllipse.ZData;

startSeg = Segment();
startSeg.x = x;
startSeg.y = y;
startSeg.z = z;

[numrow, numcol, numz] = size(FinalImage);
outsideWin = 30;

newSeeds = zeros(size(seedCoord, 1), 3);

for i = 1:size(seedCoord, 1)
    point = seedCoord(i, :)';
    
    % Generate temporary segment
    temp = startSeg;
    temp = temp.translate(point);
    temp = temp.intensity_est(FinalImage, outsideWin);
    
    % Next estimate new mu
    % First calculate likelihood gradient
    delL = temp.Lgrad(FinalImage) * 0.03;
    muN = delL(1:3);
    temp = temp.translate(muN);
    temp.mu(1) = temp.mu(1) * (temp.mu(1) >= 1 && temp.mu(1) <= numrow) + ...
        (temp.mu(1) < 1) + (temp.mu(1) > numrow) * numrow;
    temp.mu(2) = temp.mu(2) * (temp.mu(2) >= 1 && temp.mu(2) <= numcol) + ...
        (temp.mu(2) < 1) + (temp.mu(2) > numcol) * numcol;
    temp.mu(3) = temp.mu(3) * (temp.mu(3) >= 1 && temp.mu(3) <= numz) + ...
        (temp.mu(3) < 1) + (temp.mu(3) > numz) * numz;
    
    newSeeds(i, :) = temp.mu';
    
end

save('newSeeds.mat', 'newSeeds')

end