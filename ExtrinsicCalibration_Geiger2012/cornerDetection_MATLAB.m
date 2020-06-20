function corner = cornerDetection_MATLAB(image, cameraParams, nPattern)

debug = true;

global imPoints worldPoints K optical_axis;
%% Hyper Parameter
squareSizeInMeter = 0.05; % millimeters, pattern square sizes

%% pre-setting Parameter
imageSize = [size(image, 2),size(image, 1)];
K = cameraParams.IntrinsicMatrix';
optical_axis = [K(1, 1) K(1, 3) K(2, 3)];
%% Pattern Detection
pimage = image;
for i = 1 : nPattern
    % Detect calibration pattern.
    [imagePoints{i}, boardSize{i}] = detectCheckerboardPoints(pimage);
    
    [n_points, ~] = size(imagePoints{i});
    
    xmax = max(imagePoints{i}(:, 1));
    xmin = min(imagePoints{i}(:, 1));
    ymax = max(imagePoints{i}(:, 2));
    ymin = min(imagePoints{i}(:, 2));
    
    pimage(ymin-10:ymax+10, xmin-10:xmax+10, :) = 0;
    
    % PointCloud of Camera generation
    
    worldPoints = generateCheckerboardPoints(boardSize{i}, squareSizeInMeter);
    worldPoints(:, 3) = 1;
    
    imPoints = imagePoints{i};
    
    %extrinsic calculation
    options = optimoptions('lsqnonlin', 'MaxFunctionEvaluations', 2000);
    P = [0 0 0 0 0 0 0];
    func1 = @(p)location_positioning(p);
    func2 = @(p)reprojection_error(p);
    
    RT = lsqnonlin(func1, P, [], [], options);
%     RT = lsqnonlin(func2, RT, [], [], options);
    
    
    
    %generation pointCloud
    Points3D{i} = projMatrix(RT) * [worldPoints'; ones(1, size(worldPoints, 1))];
    pc{i} = pointCloud(Points3D{i}(1:3, :)');
    
end

if debug
    figure, imshow(image);
    hold on,
    for i = 1 : nPattern
        scatter(imagePoints{i}(:, 1), imagePoints{i}(:, 2), 'r');
        points2D = projection(pc{i}.Location, K, [0 0 0 0 0 0]);
        scatter(points2D(1, :), points2D(2, :), 'g');
    end
    
    figure, hold on, xlabel('x'), ylabel('y'), zlabel('z')
    for i = 1 : nPattern
        pcshow(pc{i});
    end
    axis([-3 3 -2 2 0 12]);
end

corner = imagePoints;

end

function E = location_positioning(P)
global imPoints worldPoints K;

m = projection(worldPoints, K, P)';

xmax1 = max(imPoints(:, 1));
ymax1 = max(imPoints(:, 2));
xmin1 = min(imPoints(:, 1));
ymin1 = min(imPoints(:, 2));
center1 = mean(imPoints, 1);
var1 = var(imPoints);

xmax2 = max(m(:, 1));
ymax2 = max(m(:, 2));
xmin2 = min(m(:, 1));
ymin2 = min(m(:, 2));
center2 = mean(m, 1);
var2 = var(m);

E1 = [abs(xmax1 - xmax2) abs(ymax1 - ymax2) abs(xmin1 - xmin2) abs(ymin1 - ymin2) ...
    abs(center1(1) - center2(1)) abs(center1(2) - center2(2)) ...
    abs(var1(1) - var2(1)) abs(var1(2) - var2(2))]';

E1(1:4) = E1(1:4) * 1;
E1(5:6) = E1(5:6) * 1;
E1(7:8) = E1(7:8) * 1;

R = [cosd(45) -sind(45); sind(45) cosd(45)];
R = R(1:2, 1:2);

Ri = R * imPoints';
Ri = Ri';

Rm = R * m';
Rm = Rm';

xmax1 = max(Ri(:, 1));
ymax1 = max(Ri(:, 2));
xmin1 = min(Ri(:, 1));
ymin1 = min(Ri(:, 2));
center1 = mean(Ri, 1);
var1 = var(Ri);

xmax2 = max(Rm(:, 1));
ymax2 = max(Rm(:, 2));
xmin2 = min(Rm(:, 1));
ymin2 = min(Rm(:, 2));
center2 = mean(Rm, 1);
var2 = var(Rm);

E2 = [abs(xmax1 - xmax2) abs(ymax1 - ymax2) abs(xmin1 - xmin2) abs(ymin1 - ymin2) ...
    abs(center1(1) - center2(1)) abs(center1(2) - center2(2)) ...
    abs(var1(1) - var2(1)) abs(var1(2) - var2(2))]';



E = [E1 E2];

end

function E = reprojection_error(P)
global imPoints worldPoints K optical_axis;

m = projection(worldPoints, K, P)';

[n_points, ~] = size(imPoints);
usedIdx = ones(n_points, 1);

E1 = [];
E2 = [];
f = optical_axis(1);
cx = optical_axis(2);
cy = optical_axis(3);

center = mean(imPoints, 1);
idx_power = repmat(center, [n_points, 1]) - m;
idx_power = sum(idx_power.^2, 2);
[~, I_idx] = sort(idx_power, 'descend');

for i = I_idx' % 1 : n_points
    p1 = repmat(imPoints(i, :), [n_points 1]);
    dist1 = sum((p1 - m).^2, 2);
    dist1 = dist1 .* usedIdx;
    idx = find(dist1 == min(dist1));
    usedIdx(idx) = inf;
    E1 = [E1; dist1(idx)];

    ith_Line = [imPoints(i, :) - [cx cy] f] ./ f;
    wP = projMatrix(P) * [worldPoints'; ones(1, n_points)];
    wP = wP(1:3, :)';
    
    for j = 1 : n_points
        costheta = dot(ith_Line, wP(j, :)) ./ (norm(ith_Line) * norm(wP(j, :)));
        theta = acosd(costheta);
        dist2(j) = norm(wP(j, :) - sind(theta) * wP(j, :));
    end

    E2 = [E2; dist2(idx)];
end
E = [E1; E2];
end

function m = projection(M, K, P)
[n, ~] = size(M);
m = projMatrix(P) * [M'; ones(1, n)];
m = K * m(1:3, :);
m = m(1:2, :) ./ m(3, :);
end

function RT = projMatrix(P)
R = rod2dcm([P(1) P(2) P(3)]);
T = [P(4); P(5); P(6)];
RT = [R T];
RT = [RT; 0 0 0 1];
end