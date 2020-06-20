function Plane_C = Plane_generation(images, Patterns, cameraParams)
addpath('funtion');
ni = numel(images);

%% preprocessing
for i = 1 : ni
    feature{i} = [];
    for j = 1 : numel(Patterns{i})
        feature{i} = [feature{i}; Patterns{i}{j}];
    end
end

%% Gnerate World 3D points
I = rgb2gray(images{1});

pointsPrev = cornerPoints(feature{1});
[featuresPrev,pointsPrev] = extractFeatures(I,pointsPrev);

vSet = viewSet;
vSet = addView(vSet, 1,'Points',pointsPrev,'Orientation',...
    orientations(:,:,1),'Location',locations(1,:));

for i = 2:images.Count
  I = rgb2gray(read(images, i));
  points = cornerPoints(feature{i});
  [features, points] = extractFeatures(I, points);
  vSet = addView(vSet,i,'Points',points,'Orientation',...
      orientations(:,:,i),'Location',locations(i,:));
  pairsIdx = matchFeatures(featuresPrev,features,'MatchThreshold',5);
  vSet = addConnection(vSet,i-1,i,'Matches',pairsIdx);
  featuresPrev = features;
end

tracks = findTracks(vSet);

cameraPoses = poses(vSet);

[xyzPoints,errors] = triangulateMultiview(tracks,cameraPoses,cameraParams);
z = xyzPoints(:,3);
idx = errors < 5 & z > 0 & z < 20;
pcshow(xyzPoints(idx, :),'VerticalAxis','y','VerticalAxisDir','down','MarkerSize',30);
hold on
plotCamera(cameraPoses, 'Size', 0.1);
hold off

end