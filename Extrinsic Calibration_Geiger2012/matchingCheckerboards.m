function [matching_idx, matching_Params] = matchingCheckerboards(Patterns, cameraParams)


%%
for i = 1 : 3
    nPattern(i) = numel(Patterns{i});
    for j = 1 : nPattern(i)
        centers{i}(j, :) = mean(Patterns{i}{j});
    end
end

camMatrix{1} = cameraMatrix(cameraParams{1}, diag([1 1 1]), [0 0 0]);

for i = 2 : 3
    if nPattern(1) >= nPattern(i)
        idx1 = nchoosek([1:nPattern(1)], nPattern(i));
        idx_pair1 = [];
        idx_pair2 = [];
        for j = 1 : size(idx1, 1)
            candidate = perms(idx1(j, :));
            n_candi = size(candidate, 1);
            idx_pair1 = [idx_pair1; candidate];
            idx_pair2 = [idx_pair2; repmat([1:nPattern(i)], [n_candi 1])];
        end
    else %nPattern(1) < nPattern(i)
        idxi = nchoosek([1:nPattern(i)], nPattern(1));
        idx_pair1 = [];
        idx_pair2 = [];
        for j = 1 : size(idxi, 1)
            candidate = perms(idxi(j, :));
            n_candi = size(candidate, 1);
            idx_pair1 = [idx_pair1; repmat([1:nPattern(1)], [n_candi 1])];
            idx_pair2 = [idx_pair2; candidate];
        end
    end
    
    centers1 = centers{1};
    centers2 = centers{i};
    
    global points1 points2;
    
    func = @(p)Similarity_transform(p);
    
    options = optimoptions('lsqnonlin','Display', 'off');
    tic
    for j = 1 : size(idx_pair1, 1);
        P = diag([1 1 1]);
        
        points1 = centers1(idx_pair1(j, :), :);
        points2 = centers2(idx_pair2(j, :), :);
        
        Parameters{j} = lsqnonlin(func, P, [], [], options);
        if mod(j, 100) == 0
            j / size(idx_pair1, 1)
        end
    end
    toc
    
    tic
    for j = 1 : numel(Parameters)
        p1 = centers1(idx_pair1(j, :), :);
        p2 = centers2(idx_pair2(j, :), :);
        p2 = [p2 ones(size(p2, 1) , 1)];
        s(j) = sum(getScore(p1, p2, Parameters{j}));
    end
    toc
    
    final_idx = find(s == min(s));
    for j = 1 : numel(final_idx)
        matching_idx{i, 1}{j} = idx_pair1(final_idx(j), :);
        matching_idx{i, 2}{j} = idx_pair2(final_idx(j), :);
        matching_Params{i, 1}{j} = Parameters{final_idx(j)};
    end
end

end


function E = Similarity_transform(p)
global points1 points2;

nC = size(points1, 1);

p1 = points1;
p2 = [points2 ones(nC, 1)];

E = getScore(p1, p2, p);

end

function E = getScore(p1, p2, Pmatrix)

ksi_p = Pmatrix * p2';
ksi_p = ksi_p(1:2, :) ./ ksi_p(3, :);

dist = (p1 - ksi_p').^2;
E = sqrt(sum(dist, 2));
end



















