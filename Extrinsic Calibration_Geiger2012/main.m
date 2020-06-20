% --------------------------------------
% Paper: Automatic Camera and Range Sensor Calibration using a single shot(2012) - Implementation
% made by Byung-Hyun Yoon
% date : 2020-05-04
% --------------------------------------

load('dataset.mat');

%images : three image
%LiDAR : point cloud
%cameraParams : three cameraParam cell 

for i = 1 : 3
    images{i} = undistortImage(images{i}, cameraParams{i});
end

%% Camera to Camera Calibration1 (MATLAB)
nPattern = [7 7 7];
for i = 1 : 3
    Patterns{i} = cornerDetection_MATLAB(images{i}, cameraParams{i}, nPattern(i));
end

%% Camera to Camera Calibration1
% for i = 1 : 3
%     corner = cornerDetection(image);
%     subpixelRefinement
%     Patterns = structureRecovery(corner);
% end

%% Camera to Camera Calibration2
[matching_idx, matching_Params] = matchingCheckerboards(Patterns, cameraParams);
%%
inlier_idx = selectPossibleMatching(matching_idx, matching_Params, images, Patterns, cameraParams);
% cameraOpimization

Plane_C = Plane_generation(images, Patterns, cameraParams);

%% Test


%% Camera to Range Calibration

Planes_L = Plane_segmentation(LiDAR);
parameter = Global_Registration(Plane_L, Plane_C);
parameter = Fine_Registration(Plane_L, Plane_C, parameter);
