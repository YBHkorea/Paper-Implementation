%% abstract
%paper : Line-based Extrinsic Calibration of Range and Image Sensors(2013)
% this code made by Byung-Hyun Yoon

clear; clc; restoredefaultpath;
rng(1);
dataset_num = 6;

% 6 : C304_4_pattern

[images, ~, pc_origin, ~, cameraParams, ~] = datasetload(dataset_num);

for i = 1 : numel(images)
    images{i} = undistortImage(images{i}, cameraParams{i});
end
%% 3D Line Detection
Lines3D = Line_Detection_3D_MATLAB(pc_origin, dataset_num); 
%% 2D Line Detection
Lines2D = Line_Detection_2D(images);
%% 2D Line Selection(Manually)
[Lines2D, Lines3D] = human_intervation(Lines2D, Lines3D, dataset_num);
%% Line-based Registration
Parameters = Registration(Lines2D, Lines3D, cameraParams);
%% Result Visualization
visualization_result(images, pc_origin, cameraParams, Parameters);