function Plane3D = Plane_segmentation(pc_origin)

debug = false;
test = true;

if test 
    load('LiDAR_Refine.mat');
    pc_origin = pointCloud(LiDAR_Refine);
end

%% hyperparameter
maxDistance = 0.05;
minDistance = 0.11;
PlaneNumber = 5;

pointcloud{1} = pcdenoise(pc_origin);

%% Plane detection
for i = 2 : PlaneNumber + 1
    [model{i}, inlierIndices{i}, outlierIndices{i}] = pcfitplane(pointcloud{i-1},...
        maxDistance);
    plane{i-1} = select(pointcloud{i-1}, inlierIndices{i});
    
        if debug
            figure, pcshow(plane{i-1}), title('Plane');
            figure, pcshow(pointcloud{i-1}), title('Remain PC');
        end
    
    remainPtCloud{i} = select(pointcloud{i-1}, outlierIndices{i});
    remainPtCloud{i} = pcdenoise(remainPtCloud{i});
    pointcloud{i} = remainPtCloud{i};
end

cnt = 0;

for i = 1 : 2
    cnt = cnt + 1;
    Plane3D{cnt} = plane{i};
end


end