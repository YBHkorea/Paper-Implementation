function visualization_result(images, pc_origin, cameraParams, Parameters)

[IH, IW, ~] = size(images{1});

for i = 1 : 3
    K = cameraParams{i}.IntrinsicMatrix';
    E = Parameters{i};
    
    p = projection(pc_origin.Location, K, E);
    
    % projection result
    figure, imshow(images{i}), hold on, title('pointcloud projection');
    scatter(p(1, :), p(2, :));
    
    idx = zeros([4, pc_origin.Count]);
    
    xidx1 = find(p(1, :) > 1);
    idx(1, xidx1) = 1;
    
    xidx2 = find(p(1, :) < IW - 1);
    idx(2, xidx2) = 1;
    
    yidx1 = find(p(2, :) > 1);
    idx(3, yidx1) = 1;
    
    yidx2 = find(p(2, :) < IH - 1);
    idx(4, yidx2) = 1;
    
    inner_idx = sum(idx) == 4;
    outter_idx = ~inner_idx;
    
    inner_pc = pointCloud(pc_origin.Location(inner_idx, :));
    outter_pc = pointCloud(pc_origin.Location(outter_idx, :));
    
    %inner color
    n_pc = inner_pc.Count;
    p = projection(inner_pc.Location, K, E);
    color = impixel(images{i}, p(1, :), p(2, :));
    inner_pc.Color = uint8(color);
    
    %ouuter color
    n_pc = outter_pc.Count;
    color = ones([n_pc, 3]) * 255;
    outter_pc.Color = uint8(color);
    
    %combine inner and outter
    l1 = inner_pc.Location;
    l2 = outter_pc.Location;
    c1 = inner_pc.Color;
    c2 = outter_pc.Color;
    
    Location = [l1; l2];
    Color = [c1; c2];
    
    pc = pointCloud(Location, 'Color', Color);
    
    figure, pcshow(pc), title('pointcloud visualization'), xlabel('x'), ylabel('y'), zlabel('z');
end

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
