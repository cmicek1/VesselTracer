function vesselTrace(FinalImage, seedCoord)

dt = 0.003;

startEllipse = ellipsoidInit(0.25, 3);
x = startEllipse.XData;
y = startEllipse.YData;
z = startEllipse.ZData;

startSeg = Segment(x, y, z);

[numrow, numcol, numz] = size(FinalImage);
outsideWin = 30;

% % Use subset of all points for demo
% testSeeds = seedCoord; % (seedCoord(:, 3) >= 20 & seedCoord(:, 3) <= 30);
% newSeeds = NaN(size(testSeeds, 1), 3);
% L_list = NaN(size(testSeeds, 1), 1);
% seg_list = startSeg;
% seg_list = repmat(seg_list, length(L_list), 1);
%
% for i = 1:length(seg_list)
%     point = seedCoord(i, :)';
%
%     % Generate temporary segment
%     temp = seg_list(i);
%     temp = temp.translate(point, size(FinalImage));
%     temp = temp.intensity_est(FinalImage, outsideWin, true);
%     if temp.IB == temp.IF
%         continue
%     end
%     seg_list(i) = temp;
%     newSeeds(i, :) = temp.mu';
% end
%
% cond = isnan(newSeeds);
%
% seg_list(cond(:, 1)) = [];
% L_list(cond(:, 1)) = [];
% newSeeds = NaN(size(L_list, 1), 3);
epsilon = 1e-4;
% dt = 0.003;
% % Update position
% for i = 1:length(seg_list)
%     prevL = zeros(11, 1);
%     temp = seg_list(i);
%     for iter = 1:20
%         % Next estimate new mu
%         % First calculate likelihood gradient
%         [delL, L] = temp.Lgrad(FinalImage, true);
%         % Check for convergence
%         if all(abs(delL - prevL) < epsilon)
%             break
%         end
%         delL = delL * dt;
%         prevL = delL;
%         temp.L = L;
%         muN = delL(1:3);
%         temp = temp.translate(muN, [numrow, numcol, numz]);
%         newSeeds(i, :) = temp.mu';
%         L_list(i) = L;
%         seg_list(i) = temp;
%         if any(isinf(delL) | isnan(delL))
%             break
%         end
%     end
% end
%
% newSeeds(:, 3) = floor(newSeeds(:, 3));


% % Save so if we ever get here, we don't have to do it again
% save('newSeeds5.mat', 'newSeeds', 'seg_list', 'L_list')

% load('newSeeds5.mat')
% Now reestimate other parameters
% dt = 0.01;
% maxIt = 20;
% tic;
% for i = 1:length(seg_list)
%     if i == 265
%         continue
%     end
%     prevL = zeros(11, 1);
%     temp = seg_list(i);
%     display(i)
%     toc;
%     for iter = 1:maxIt
%         if mod(i, 5) == 0
%             temp = temp.intensity_est(FinalImage, outsideWin, false);
%         end
%         [delL, L] = temp.Lgrad(FinalImage, false);
%         Check for convergence
%         if all(abs(delL - prevL) < epsilon)
%             break
%         end
%         delL = delL * dt;
%         if any(delL >= 1e6)
%             break
%         end
%         if ~isnan(L)
%             temp.L = L;
%         else
%             L = L_list(i);
%         end
%         temp.sigma = temp.sigma + delL(4:6);
%         temp.q = temp.q + delL(7:10)';
%         temp.epsilon = temp.epsilon + delL(11);
%         Shape method applies all other changes to starting ellipsoid
%         temp = temp.shape(size(FinalImage));
%         L_list(i) = L;
%         seg_list(i) = temp;
%     end
% end
%
% save('newSeeds6.mat', 'L_list', 'seg_list')

%%
load('newSeeds6.mat')
for i = 1:length(seg_list)
    L_list(i) = seg_list(i).IB - seg_list(i).IF;
end
[~, inds] = sort(L_list);
candidates = seg_list(inds(end:-1:1));
L_list = L_list(inds(end:-1:1));
inds = [];
for i = 1:length(candidates)
    if isempty(candidates(i).x) || any(candidates(i).mu([1, 2]) < 4) ||...
            candidates(i).mu(1) > numrow - 3 ...
            || candidates(i).mu(2) > numcol - 3 % ||candidates(i).mu(3) == numz
        inds = horzcat(inds, i);
    end
end
candidates(inds) = [];
[XX, YY, ZZ] = meshgrid(1:numrow, 1:numcol, 1:numz);
visited = zeros(size(FinalImage));


numToVisit = 1;
vessels = cell(1, numToVisit);
collisions = [];
collPoints = [];
done = false;

maxIt = 30;
updateIt = 10;
for v = 5
    cur = candidates(v);
    [visited, collisions, collPoints, done] = cur.markVisited(visited,...
        {XX, YY, ZZ}, v, collisions, collPoints);
    % Make segment container
    segs = cell(1, 2 * maxIt + 1);
    segs{maxIt + 1} = cur;
    evecs = pca([cur.x', cur.y', cur.z']);
    % Should be largest eigenvector
    p_axis = evecs(:, 1);
    % Make unit length
    p_axis = p_axis / norm(p_axis);
    
    choices = [-p_axis, p_axis];
    for c = 1:2
        prevDir = choices(:, c);
        prevSeg = cur;
        for i = 1:maxIt
            newSeg = prevSeg;
            newSeg = newSeg.translate(prevDir, size(FinalImage));
            evecs = pca([newSeg.x', newSeg.y', newSeg.z']);
            newDir = evecs(:, 1) / norm(evecs(:, 1));
            angle = atan2(norm(cross(prevDir,newDir)),dot(prevDir,newDir));
            if angle > pi/2
                newDir = -newDir;
            end
            prevDir = newDir;
            prevL = zeros(11, 1);
            for it = 1:updateIt
                [delL, ~] = newSeg.Lgrad(FinalImage, false);
                % Check for convergence
                if all(abs(delL - prevL) < epsilon)
                    break
                end
                delL = delL * dt;
                if any(delL >= 1e6)
                    break
                end
                newSeg.sigma = newSeg.sigma + delL(4:6);
                newSeg.q = newSeg.q + delL(7:10)';
                newSeg.epsilon = newSeg.epsilon + delL(11);
                % Shape method applies all other changes to starting ellipsoid
                newSeg = newSeg.shape(size(FinalImage));
                prevL = delL;
                newSeg = newSeg.intensity_est(FinalImage, outsideWin, false);
            end
            
            
            [visited, collisions, collPoints, done] =...
                newSeg.markVisited(visited, {XX, YY, ZZ}, v, collisions, collPoints);
            
            if done || newSeg.tstat(FinalImage, 0.2)
                break
            end
            
            if c == 1
                segs{numToVisit + 1 + i} = newSeg;
            else
                segs{numToVisit + 1 - i} = newSeg;
            end
            
            prevSeg = newSeg;
        end
    end
    vessels{v} = segs;
    mip = min(FinalImage, [], 3);
    points = [];
    for i = 1:length(segs)
        if ~isempty(segs{i})
            points = vertcat(points, segs{i}.mu');
        end
    end
    figure(1)
    plot(points(:, 1), points(:, 2))
    imshow(mip)
    hold on
    plot(points(:, 2), points(:, 1))
    hold off
end
end