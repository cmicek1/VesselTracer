function vesselTrace(FinalImage, seedCoord)

dt = 0.03;

startEllipse = ellipsoidInit(0.25, 3);
x = startEllipse.XData;
y = startEllipse.YData;
z = startEllipse.ZData;

startSeg = Segment(x, y, z);

[numrow, numcol, numz] = size(FinalImage);
outsideWin = 30;

% Use subset of all points for demo
testSeeds = seedCoord(seedCoord(:, 3) >= 20 & seedCoord(:, 3) <= 30);

newSeeds = NaN(size(testSeeds, 1), 3);
L_list = NaN(size(testSeeds, 1), 1);
seg_list = startSeg;
seg_list = repmat(seg_list, length(L_list), 1);

for i = 1:length(seg_list)
    point = seedCoord(i, :)';
    
    % Generate temporary segment
    temp = seg_list(i);
    temp = temp.translate(point, size(FinalImage));
    temp = temp.intensity_est(FinalImage, outsideWin, true);
    if temp.IB == temp.IF
        continue
    end
    seg_list(i) = temp;
    newSeeds(i, :) = temp.mu';
end

cond = isnan(newSeeds);

seg_list(cond(:, 1)) = [];
L_list(cond(:, 1)) = [];
newSeeds = NaN(size(L_list, 1), 3);
epsilon = 1e-4;

prevL = zeros(10, 1);

% Update position
for i = 1:length(seg_list)
    temp = seg_list(i);
    for iter = 1:temp.initIterMax
        % Next estimate new mu
        % First calculate likelihood gradient
        [delL, L] = temp.Lgrad(FinalImage, true);
        % Check for convergence
        if all(abs(delL - prevL) < epsilon)
            break
        end
        delL = delL * dt;
        temp.L = L;
        muN = delL(1:3);
        temp = temp.translate(muN, [numrow, numcol, numz]);
        newSeeds(i, :) = temp.mu';
        L_list(i) = L;
        seg_list(i) = temp;
    end
end
% Save so if we ever get here, we don't have to do it again
save('newSeeds.mat', 'newSeeds')

% Now reestimate other parameters
maxIt = 20;

for i = 1:length(seg_list)
    temp = seg_list(i);
    for iter = 1:maxIt
        if mod(i, 5) == 0
            temp = temp.intensity_est(FinalImage, outsideWin, true);
        end
        [delL, L] = temp.Lgrad(FinalImage, true);
        % Check for convergence
        if all(abs(delL - prevL) < epsilon)
            break
        end
        delL = delL * dt;
        temp.L = L;
        muN = delL(1:3);
        temp = temp.translate(muN, [numrow, numcol, numz]);
        temp = temp.scale(delL(4:6));
        temp = temp.rotate(delL(7:9));
        temp = temp.shape(delL(10));
        seg_list(i) = temp;
    end
end
end