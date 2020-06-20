function Lines3D = Line_Detection_3D_MATLAB(pc_origin, dataset_num)
rng(1);
debug = false;

%% hyperparameter
if dataset_num == 6
    roi = [-2.2 3 3 7 -1 1.5];
    maxDistance = 0.05;
    minDistance = 0.11;
    PlaneNumber = 4;
end

%% spliting point cloud
indices = findPointsInROI(pc_origin,roi);
pointcloud{1} = select(pc_origin,indices);

pointcloud{1} = pcdenoise(pointcloud{1});

%% Plane detection
for i = 2 : PlaneNumber + 1
    [model{i}, inlierIndices{i}, outlierIndices{i}] = pcfitplane(pointcloud{i-1},...
        maxDistance);
    plane{i-1} = select(pointcloud{i-1}, inlierIndices{i});
    
    %     if debug
    %         figure, pcshow(plane{i-1}), title('Plane');
    %     end
    
    remainPtCloud{i} = select(pointcloud{i-1}, outlierIndices{i});
    remainPtCloud{i} = pcdenoise(remainPtCloud{i});
    pointcloud{i} = remainPtCloud{i};
end

cnt = 0;

for i = 1 : 2
    cnt = cnt + 1;
    Plane3D{cnt} = plane{i};
end

%% pointcloud segmentation
for i = 3 : numel(plane)
    plane{i} = pcdenoise(plane{i});
    [labels,numClusters] = pcsegdist(plane{i}, minDistance);
    
    %     if debug
    %     figure, pcshow(plane{i}.Location,labels)
    %     colormap(hsv(numClusters))
    %     title('Point Cloud Clusters')
    %     end
    
    for j = 1 : numClusters
        if sum(labels == j) > 500
            cnt = cnt + 1;
            idx = find(labels == j);
            Plane3D{cnt} = pointCloud(plane{i}.Location(idx, :));
        end
    end
end

%% find bounding box
cnt = 0;
for i = 1 : numel(Plane3D)
    [coeff,score,latent] = pca(Plane3D{i}.Location);
    %     [U, S, V] = svd(Plane3D{i}.Location);
    Plane3D_proj{i} = pointCloud((coeff * Plane3D{i}.Location')');
%     if debug
%         figure, pcshow(Plane3D_proj{i}), xlabel('x'), ylabel('y'), zlabel('z');
%     end
    
    xmin = min(Plane3D_proj{i}.Location(:, 1));
    xmax = max(Plane3D_proj{i}.Location(:, 1));
    ymin = min(Plane3D_proj{i}.Location(:, 2));
    ymax = max(Plane3D_proj{i}.Location(:, 2));
    
    p1 = find(Plane3D_proj{i}.Location(:, 1) == xmin);
    p2 = find(Plane3D_proj{i}.Location(:, 1) == xmax);
    p3 = find(Plane3D_proj{i}.Location(:, 2) == ymin);
    p4 = find(Plane3D_proj{i}.Location(:, 2) == ymax);
    
    Data = [Plane3D{i}.Location(p1, :); Plane3D{i}.Location(p2, :); Plane3D{i}.Location(p3, :); Plane3D{i}.Location(p4, :)];
    K = convhull(Data(:, 1), Data(:, 2));
    
    for j = 1 : numel(K) - 1
        cnt = cnt + 1;
        Lines3D(cnt).point1 = Data(K(j), :);
        Lines3D(cnt).point2 = Data(K(j+1), :);
    end
end

%% show result
if debug
    figure, hold on, xlabel('x'), ylabel('y'), zlabel('z')
    pcshow(pc_origin);
    for i = 1 : cnt
        plot3([Lines3D(i).point1(1), Lines3D(i).point2(1)], [Lines3D(i).point1(2), Lines3D(i).point2(2)], [Lines3D(i).point1(3), Lines3D(i).point2(3)], '-o');
    end
end

end