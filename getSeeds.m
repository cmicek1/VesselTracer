function getSeeds(image)
load(image)

% Trim first few slices
offset = 3;
FinalImage = FinalImage(:, :, (1+offset):end);

% Estimate of average vessel width, in pixels
width = 6;
seeds = {};
[rmax, cmax, zmax] = size(FinalImage);
% Condense image planes into line grid, search for potential seed points
k_range = 1:width:(zmax - width);
j_range = 1:width:(cmax - width);
i_range = 1:width:(rmax - width);
for k = k_range
    endz = min([zmax, k + width - 1]);
    for j = j_range
        endc = min([cmax, j + width - 1]);
        for i = i_range
            endr = min([rmax, i + width - 1]);
            
            block = FinalImage(i:endr, j:endc, k:endz);
            [~, idx] = min(block(:));
            
            [seedr, seedc, seedz] = ind2sub(size(block), idx);
            seedr = i + seedr - 1;
            seedc = j + seedc - 1;
            seedz = k + seedz - 1 + offset;
            
            seeds{end + 1} = [seedr, seedc, seedz];
        end
    end
end
    
    numSeeds = numel(seeds);
    seedCoord = zeros(numSeeds,3);
    for ii = 1:1:numSeeds
        seedCoord(ii,:) = seeds{ii};
    end
    save('trace_seeds.mat', 'seeds','seedCoord')
end