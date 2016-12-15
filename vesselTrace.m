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
% 
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
% 
% % Update position
% for i = 1:length(seg_list)
%     prevL = zeros(11, 1);
%     temp = seg_list(i);
%     for iter = 1:temp.initIterMax
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
% 
% % Save so if we ever get here, we don't have to do it again
% save('newSeeds.mat', 'newSeeds', 'seg_list', 'L_list')

load('newSeeds.mat')
% Now reestimate other parameters
maxIt = 20;
tic;
for i = 1:length(seg_list)
    prevL = zeros(11, 1);
    temp = seg_list(i);
    display(i)
    toc;
    for iter = 1:maxIt
        if mod(i, 5) == 0
            temp = temp.intensity_est(FinalImage, outsideWin, false);
        end
        [delL, L] = temp.Lgrad(FinalImage, false);
        % Check for convergence
        if all(abs(delL - prevL) < epsilon)
            break
        end
        delL = delL * dt;
        if any(delL >= 1e10)
            break
        end
        if ~isnan(L)
            temp.L = L;
        else
            L = L_list(i);
        end
        temp.sigma = temp.sigma + delL(4:6);
        temp.q = temp.q + delL(7:10)';
        temp.epsilon = temp.epsilon + delL(11);
        % Shape method applies all other changes to starting ellipsoid
        temp = temp.shape(size(FinalImage));
        L_list(i) = L;
        seg_list(i) = temp;
    end
end

save('newSeeds2.mat', 'L_list', 'seg_list')

end