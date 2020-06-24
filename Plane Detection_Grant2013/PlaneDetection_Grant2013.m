function [Planes_result] = PlaneDetection_Grant2013(pc, n_channel)
% --------------------------------------
% Paper: Finding Planes in LiDAR Point Clouds for Real-Time Registration(2013) - Implementation
% https://github.com/YBHkorea/LiDAR-processing/tree/master/LiDARregistration-2013-Finding%20Planes%20in%20LiDAR%20Point%20Clouds%20for%20Real-Time%20Registration
% made by Byung-Hyun Yoon
% data : 2020-06-04
% --------------------------------------

debug = false;

%% (3 page)
% Clustering Lines

% 1D signal
signal_1D = sqrt(sum(pc.Location.^2, 2));
Location = pc.Location ./ signal_1D;

% splitting Lines using Clustering ( VLP-16 -> k is 16 )
idx = kmeans(Location(:, 3), n_channel);

%sorting idx
[sort_label, sort_idx] = sort(idx);
for j = 1 : n_channel
    idx = sort_idx(sort_label == j);
    Lines{j} = Location(idx, :);
    pc_Line{j} = pc.Location(idx, :);
end

%% (4 page)
%sorting idx in Lines using x-y plane distance
r_vector = [0 1];
for j = 1 : n_channel
    xy_Lines{j} = Lines{j}(:, 1:2);
    dist_ = sqrt(sum(xy_Lines{j}.^2 ,2));
    xy_Lines{j} = xy_Lines{j}./ dist_;
    
    %n : negative, p : positive
    nx = xy_Lines{j}(:, 1) <= 0;
    px = ~nx;
    
    %sort -y to y(clockwise)
    nx_points = xy_Lines{j}(nx, 2);
    [~, nx_idx] = sort(nx_points, 'ascend');
    nx_Line = pc_Line{j}(nx, :);
    nx_Line = nx_Line(nx_idx, :);
    
    %sort y to -y(clockwise)
    px_points = xy_Lines{j}(px, 2);
    [~, px_idx] = sort(px_points, 'descend');
    px_Line = pc_Line{j}(px, :);
    px_Line = px_Line(px_idx, :);
    
    pc_Line{1, j} = [nx_Line; px_Line];
end

for j = 1 : n_channel
    pc_Line{2, j} = sqrt(sum(pc_Line{1, j}.^2, 2));
end

%% (5 page, 6 page)
% gaussian and gradient filtering for break points
% pc_Line{1, :} : Lines
% pc_Line{2, :} : 1D-signal of Lines


sigma = 0.02;
Gaussian = gaussdesign(0.3, 8, 4);

for j = 1 : n_channel
    pc_Line{2, j} = conv(pc_Line{2, j}, Gaussian, 'same');
    gradient_Line{1, j} = conv(pc_Line{2, j}, [-1 1], 'same');
    gradient_Line{2, j} = conv(gradient_Line{1, j}, [-1 1], 'same');
    gradient_Line{3, j} = conv(gradient_Line{2, j}, [-1 1], 'same');
end

threshold_g2 = 0.01;
for j = 1 : n_channel
    g_idx = find(abs(gradient_Line{2, j}) > threshold_g2);
    cnt = 1;
    for k = 1 : numel(g_idx) -2
        if gradient_Line{3, j}(g_idx(k)) * gradient_Line{3, j}(g_idx(k) + 2) < 0
            breakPoints_idx{j}(cnt) = g_idx(k + 1);
            cnt = cnt + 1;
        end
    end
end

if debug
    figure, hold on;
    for j = 1 : n_channel
        intensity = ones(size(pc_Line{2, j}));
        intensity(breakPoints_idx{j}) = 2;
        pcshow(pointCloud(pc_Line{1, j}, 'Intensity', intensity), 'MarkerSize', 200);
    end
end

%% (7 page)
% Grouping close point sets
% #p : 15

% A Colloection of points
for j = 1 : n_channel
    cnt = 1;
    for k = 1 : numel(breakPoints_idx{j}) - 1
        n_points = breakPoints_idx{j}(k + 1) - breakPoints_idx{j}(k);
        if n_points > 15 % #p
            valid_Lines{j}(cnt, 1) = breakPoints_idx{j}(k);
            valid_Lines{j}(cnt, 2) = breakPoints_idx{j}(k + 1);
            cnt = cnt + 1;
        end
    end
end

%Grouping
% convert Group points to segments structure
cnt = 1;
segments = [];
for j = 1 : n_channel
    for k = 1 : size(valid_Lines{j}, 1)
        segments(cnt).points = pc_Line{1, j}(valid_Lines{j}(k, 1) : valid_Lines{j}(k, 2), :);
        segments(cnt).centroid = mean(segments(cnt).points, 1);
        
        [coeff, ~, latent] = pca(segments(cnt).points, 'Algorithm', 'eig'); %principal components analysis
        segments(cnt).v3 = coeff(1:3, 1);
        segments(cnt).v2 = coeff(1:3, 2);
        segments(cnt).v1 = coeff(1:3, 3);
        segments(cnt).lambda3 = latent(1);
        segments(cnt).lambda2 = latent(2);
        segments(cnt).lambda1 = latent(3);
        segments(cnt).row = j;
        
        segments(cnt).curvature = (segments(cnt).lambda1 + segments(cnt).lambda2) / ...
            (segments(cnt).lambda1 + segments(cnt).lambda2 + segments(cnt).lambda3);
        cnt = cnt + 1;
    end
end



% this process is ignored. in our dataset, this situation doesn't
% required.
%     t_g = 0.03; % 3 cm
%     t_d = 0.8;
%
%     converge = false;
%     while ~converge
%         converge = true;
%         n_seg = size(segments, 2);
%         for j = 1 : n_seg - 1 %search like bubble sort
%             for k = j : n_seg
%                 if j == k
%                     continue;
%                 end
%
%                 if segments(j).v3' * segments(k).v3 < t_d
%                     continue;
%                 end
%
%                 dist = pdist2(segments(j).points, segments(k).points);
%
%                 if min(min(dist)) <= t_g
%                     converge = false;
%                 end
%             end
%         end
%     end

%% (8 page)
velodyne_range = 100; % 100 meter

eta = 2; %degree

% Accumulator array
N_theta = 180;
N_pi = 90;
N_rho = 600;
Accumulator_Array = zeros(N_theta, N_pi, N_rho);
Accumulator_Counter = zeros(N_theta, N_pi, N_rho);
Accumulator_Group = cell(N_theta, N_pi, N_rho);

% Ball Accumulator



%% (9 page)
rho = velodyne_range/N_rho : velodyne_range/N_rho : velodyne_range;
for j = 1 : size(segments, 2)
    rho_cnt = 0;
    for rho_ =  rho
        rho_cnt = rho_cnt + 1;
        
        %eq. (2, 3)
        A = [segments(j).centroid; segments(j).v3'];
        Ap = pinv(A);
        
        b = [-rho_; 0];
        
        j_vector = null(A);
        k_vector = Ap * b;
        
        %quadratic equation solution (5, 6)
        % 0 = X alpha^2 + Y alpha + Z
        X = sum(j_vector.^2);
        Y = 2 * (j_vector' * k_vector);
        Z = sum(k_vector.^2) - 1;
        
        if (Y^2 - 4 * X * Z >= 0)
            alpha1 = (-Y + sqrt(Y^2 - 4 * X * Z)) / (2 * X);
            alpha2 = (-Y - sqrt(Y^2 - 4 * X * Z)) / (2 * X);
            
            normal1 = alpha1 * j_vector + k_vector;
            normal2 = alpha2 * j_vector + k_vector;
            for normal = [normal1 normal2]
                
                Phi = acos(normal(3));
                Theta = asin( normal(2)/sin(Phi) ) / pi * 180;
                Phi = Phi / pi * 180; %rad to degree
                
                if Phi < 0
                    Phi = Phi + 360;
                end
                if Phi > 360
                    Phi = Phi - 360;
                end
                if Theta < 0
                    Theta = Theta + 360;
                end
                if Theta > 360
                    Theta = Theta - 360;
                end
                
                if((Phi/4) < 0.5)
                    Phi = 360;
                end
                if((Theta/2) < 0.5)
                    Theta = 360;
                end
                
                %eq. (7)
                voting_weight = segments(j).curvature * (1 - (1 - (segments(j).v1' * normal)^2)^3) + ...
                    (1 - segments(j).curvature) * ( segments(j).curvature * ( segments(j).v1' * normal ) + 0.85 - segments(j).curvature);
                Accumulator_Array(round(Theta/2), round(Phi/4), rho_cnt) = Accumulator_Array(round(Theta/2), round(Phi/4), rho_cnt) + voting_weight;
                Accumulator_Counter(round(Theta/2), round(Phi/4), rho_cnt) = Accumulator_Counter(round(Theta/2), round(Phi/4), rho_cnt) + 1;
                Accumulator_Group{round(Theta/2), round(Phi/4), rho_cnt} = [Accumulator_Group{round(Theta/2), round(Phi/4), rho_cnt}; j];
            end
        end
    end
    
    %% compensate perturbation error (in the paper)
    rho_cnt2 = 0;
    starting_flag = false;
    for rho2_ = rho
        rho_cnt2 = rho_cnt2 + 1;
        
        cos_omega = rho2_ / norm(segments(j).centroid);
        
        c_vector = - segments(j).centroid ./ norm(segments(j).centroid);
        v_vector = segments(j).v3 ./ norm(segments(j).v3);
        
        
        alpha = c_vector(2) * v_vector(1) - c_vector(1) * v_vector(2);
        beta = c_vector(1) * v_vector(1) + c_vector(1) * v_vector(1);
        gamma = c_vector(3) * v_vector(3) - cos_omega;
        
        %solution for cos(eta)
        %quadratic equation solution (8, 9)
        A = alpha^2 + beta^2;
        B = -2 * beta * gamma;
        C = gamma^2 - alpha^2;
        
        if (B^2 - 4 * A * C < 0)
            continue;
        end
        
        cos_eta_1 = (-B + sqrt(B^2 - 4 * A * C)) / (2 * A);
        cos_eta_2 = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A);
        
        if ( abs(cos_eta_1) > 1 ) || ( abs(cos_eta_2) > 1 )
            if (starting_flag == true)
                break;
            end
            continue;
        end
        
        eta1 = acos(cos_eta_1) / pi * 180;
        eta2 = 360 - eta1;
        eta3 = acos(cos_eta_2) / pi * 180;
        eta4 = 360 - eta3;
        
        eta = [eta1 eta2 eta3 eta4];
        for e = eta
            if abs(e) > 2
                continue;
            end
            starting_flag = true;
            
            %eq. (2, 3)
            v3R = [cos(e) -sin(e) 0; sin(e) cos(e) 0; 0 0 1] * segments(j).v3;
            A = [segments(j).centroid; v3R'];
            Ap = pinv(A);
            
            b = [-rho2_; 0];
            
            j_vector = null(A);
            k_vector = Ap * b;
            
            %quadratic equation solution (5, 6)
            % 0 = X alpha^2 + Y alpha + Z
            X = sum(j_vector.^2);
            Y = 2 * (j_vector' * k_vector);
            Z = sum(k_vector.^2) - 1;
            
            if (Y^2 - 4 * X * Z >= 0)
                alpha1 = (-Y + sqrt(Y^2 - 4 * X * Z)) / (2 * X);
                alpha2 = (-Y - sqrt(Y^2 - 4 * X * Z)) / (2 * X);
                
                normal1 = alpha1 * j_vector + k_vector;
                normal2 = alpha2 * j_vector + k_vector;
                for normal = [normal1 normal2]
                    
                    Phi = acos(normal(3));
                    Theta = asin( normal(2)/sin(Phi) ) / pi * 180;
                    Phi = Phi / pi * 180; %rad to degree
                    
                    if Phi < 0
                        Phi = Phi + 360;
                    end
                    if Phi > 360
                        Phi = Phi - 360;
                    end
                    if Theta < 0
                        Theta = Theta + 360;
                    end
                    if Theta > 360
                        Theta = Theta - 360;
                    end
                    
                    if((Phi/4) < 0.5)
                        Phi = 360;
                    end
                    if((Theta/2) < 0.5)
                        Theta = 360;
                    end
                    
                    %eq. (7)
                    voting_weight = segments(j).curvature * (1 - (1 - (segments(j).v1' * normal)^2)^3) + ...
                        (1 - segments(j).curvature) * ( segments(j).curvature * ( segments(j).v1' * normal ) + 0.85 - segments(j).curvature);
                    Accumulator_Array(round(Theta/2), round(Phi/4), rho_cnt2) = Accumulator_Array(round(Theta/2), round(Phi/4), rho_cnt2) + voting_weight;
                    Accumulator_Counter(round(Theta/2), round(Phi/4), rho_cnt2) = Accumulator_Counter(round(Theta/2), round(Phi/4), rho_cnt2) + 1;
                    Accumulator_Group{round(Theta/2), round(Phi/4), rho_cnt2} = [Accumulator_Group{round(Theta/2), round(Phi/4), rho_cnt2}; j];
                end
            end
        end
    end
end

%% Filtering and Refitting ( in the paper )
%parameters
t_a =  1.0;
t_split = 0.9;
N_split = 2;
t_merge = 0.9;

% (page 10)
% by keeping all bins which pass some value, t_a, and discarding the
% others
Accumulator_Array(Accumulator_Array < t_a) = 0;

%Linearity Filtering
% All candidate planes that were voted for by less than two rows are
% removed
filtering_idx = find(Accumulator_Counter <= 2);
Accumulator_Array(filtering_idx) = 0;



%% Splitting (page 11)
% Splitting Facing Planes to other group
% (eg. left and right wall has same normals or coplanar patches)
Plane_idx = find(Accumulator_Array);
Planes_cnt = 0;
for j = Plane_idx'
    pidx = unique(Accumulator_Group{j});
    idx_used = zeros(size(pidx'));
    p_label = 0;
    for k = pidx'
        if idx_used(k == pidx') ~= 0
            continue;
        end
        
        p_label = p_label + 1;
        idx_used(k == pidx') = p_label;
        for l = pidx'
            if idx_used(l == pidx') ~= 0
                continue;
            end
            
            if segments(k).v3' * segments(l).v3 < 0
                continue;
            end
            
            dist = pdist2(segments(k).points, segments(l).points);
            
            if min(min(dist)) > t_split
                continue;
            end
            
            idx_used(l == pidx') = p_label;
        end
    end
    
    for k = 1 : p_label
        if numel(find(idx_used == k)) >= 2
            Planes_cnt = Planes_cnt + 1;
            Planes(Planes_cnt).points = [];
            Planes(Planes_cnt).group = [];
            seg_idx = pidx(idx_used == k);
            
            for l = seg_idx'
                Planes(Planes_cnt).points = [Planes(Planes_cnt).points; segments(l).points];
                Planes(Planes_cnt).group = [Planes(Planes_cnt).group; l];
            end
        end
    end
end

% Once clusters are formed, planes are re-fit to their corresponding
% points
for j = 1 : numel(Planes)
    [coeff, ~, ~] = pca(Planes(j).points);
    Planes(j).normal = coeff(:, 3);
    Planes(j).valid = 1;
end

%% Merging ( page 12 )
% normal filtering
for j = 1 : numel(segments)
    normal = [];
    for k = 1 : numel(Planes)
        if numel(find(Planes(k).group == j)) > 0
            normal = [normal Planes(k).normal];
        end
    end
    
    normal = mean(normal, 2);
    normal = normal./ norm(normal);
    if size(normal, 1)
        for k = 1 : numel(Planes)
            if numel(find(Planes(k).group == j)) > 0
                if normal' * Planes(k).normal < t_merge
                    Planes(k).valid = 0;
                end
            end
        end
    end
end

cnt = 0;
for j = 1 : numel(Planes)
    if Planes(j).valid
        cnt = cnt + 1;
        Planes_remain(cnt) = Planes(j);
    end
end

% graph based merging
graph = zeros(numel(Planes_remain), numel(segments));

for j = 1 : numel(segments)
    for k = 1 : numel(Planes_remain)
        if numel(find(Planes_remain(k).group == j)) > 0
            graph(k, j) = 1;
        end
    end
end


for j = 1 : numel(segments)
    if sum(graph(:, j)) == 0
        continue;
    end
    
    idx = find(graph(:, j));
    
    for k = idx'
        if k == idx(end)
            continue;
        end
        
        graph(idx(end), :) = graph(idx(end), :) | graph(k, :);
        graph(k, :) = 0;
    end
end

% meged result
pcnt = 0;
for j = 1 : numel(Planes_remain)
    if sum(graph(j, :)) == 0
        continue;
    end
    
    pcnt = pcnt + 1;
    
    Planes_result(pcnt).points = [];
    
    idx = find(graph(j, :));
    for k = idx
        Planes_result(pcnt).points = [Planes_result(pcnt).points; segments(k).points];
    end
    
    centroid = mean(Planes_result(pcnt).points);
    
    [coeff, ~, ~] = pca(Planes_result(pcnt).points);
    
    Planes_result(pcnt).normal = coeff(1:3, 3);
    Planes_result(pcnt).distance = centroid * coeff(1:3, 3);
    
    if Planes_result(pcnt).distance < 0
        Planes_result(pcnt).normal = - coeff(1:3, 3);
        Planes_result(pcnt).distance = centroid * (-coeff(1:3, 3));
    end
end

if debug
    figure, hold on; xlabel('x'), ylabel('y'), zlabel('z');
    for j = 1 : numel(Planes_result)
        color = uint8([randi(255) randi(255) randi(255)]);
        pcshow(pointCloud(Planes_result(j).points, 'Color', repmat(color, [size(Planes_result(j).points, 1), 1])), 'MarkerSize', 200);
    end
end

end

