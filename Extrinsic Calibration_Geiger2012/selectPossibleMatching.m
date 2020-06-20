function inlier_idx = selectPossibleMatching(matching_idx, matching_Params, images, Patterns, cameraparams)

debug = true;
global debug

for i = 1 : 3
    corners{i} = [];
    for j = 1 : numel(Patterns{i})
        corners{i} = [corners{i}; Patterns{i}{j}];
    end
end


for i = 2 : numel(images)
    
    prev_score = inf;
    for j = 1 : numel(matching_Params{i})
        
        new_score = getScore(corners{1}, corners{i}, matching_Params{i}{j});
        
        if prev_score > new_score
            fine_idx = j;
            prev_score = new_score;
        end
    end
    
    if debug 
        p3 = matching_Params{3}{1} * [corners{3} ones(size(corners{3}, 1), 1)]';
        p3 = p3(1:2, :) ./ p3(3, :);
        p3 = p3';
        figure, imshow(images{1})
        hold on 
        scatter(p3(:, 1), p3(:, 2))
        
    end
    
    [inlier_idx{i, 1}, inlier_idx{i, 2}] = find_inlier(corners{1}, corners{i});
    
    
    
    
end


[R,T] = relativeCameraPose(projective2d(Parameters{final_idx}), ...
    cameraParams{1}, cameraParams{i}, ...
    centers1(idx1, :), centers2(idx2, :));
camMatrix{i} = cameraMatrix(cameraParams{i}, R, T);

%matched 점들을 구해야함, Ax+b 로 투영하고 가장 가까운 점을 선택하는 형식

worldPoints = triangulate(matchedPoints1, matchedPoints2, ...
    camMatrix{1}, camMatrix{i})

end

function S = getScore(p1, p2, P)
[n_points1, ~] = size(p1);
[n_points2, ~] = size(p2);
p2 = [p2 ones(n_points2, 1)];
p2 = P * p2';
p2 = p2(1:2, :) ./ p2(3, :);
p2 = p2';

sum_dist = 0;
for i = 1 : n_points1
    prev_dist = inf;
    for j = 1 : n_points2
        new_dist = sum(sqrt((p1(i, :) - p2(j, :)).^2));
        if new_dist < prev_dist
            candi_dist = new_dist;
            prev_dist = new_dist;
        end
    end
    sum_dist = sum_dist + candi_dist;
end

S = sum_dist;
end

function indexPairs = find_inlier(p1, p2, im1, im2)



global debug
if debug
    matchedPoints1 = p1(indexPairs(:,1), :);
    matchedPoints2 = p2(indexPairs(:,2), :);
    figure; showMatchedFeatures(im1, im2, matchedPoints1, matchedPoints2);
end
end