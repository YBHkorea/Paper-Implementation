function [SLAM_result, tform] = LiDARSLAM_Grant2018(SLAM_result, pc, Planes_prev, Planes_new, tform)
% --------------------------------------
% Paper: Efficient Velodyne SLAM with point and plane features
% github :
% made by Byung-Hyun Yoon
% data : 2020-06-17
% --------------------------------------

%% 3.2 Plane detection and correspondence
P.t_dot = 0.5;
P.t_dist = 0.02; % 2.0 cm

correspondence = making_correspondence(Planes_prev, Planes_new, P);

%% 3.5 Computing Transforms

[T] = Translation_estimation(Planes_prev, Planes_new, correspondence);
[R] = Rotation_estimation(Planes_prev, Planes_new, correspondence);

end

function correspondence = making_correspondence(Planes_prev, Planes_new, P)
correspondence = [];
for i = 1 : numel(Planes_prev)
    dist_prev = inf;
    t_dot = 0;
    t_dist = inf;
    for j = 1 : numel(Planes_new)
        dist_new = norm(Planes_prev(i).distance * Planes_prev(i).normal - Planes_new(j).distance * Planes_new(j).normal);
        
        if dist_new < dist_prev
            dist_prev = dist_new;
            t_dot = dot(Planes_prev(i).normal, Planes_new(j).normal);
            t_dist = abs(Planes_prev(i).distance - Planes_new(j).distance) / (Planes_prev(i).distance + Planes_new(j).distance);
            
            i_idx = i;
            j_idx = j;
        end
    end
    
    if (t_dot > P.t_dot) && (t_dist < P.t_dist)
        correspondence = [correspondence; i_idx j_idx];
    end
end
end

function [T] = Translation_estimation(Planes_prev, Planes_new, correspondence)


C_d = zeros(size(correspondence, 1), size(correspondence, 1));
M = zeros(size(correspondence, 1), 3);
d = zeros(size(correspondence, 1), 1);

%% equation (9)
if size(correspondence, 1) >= 3
    for i = 1 : size(correspondence, 1)
        frame_a = Planes_prev(correspondence(i, 1));
        frame_b = Planes_new(correspondence(i,2));
        
        M(i, :) = frame_a.normal';
        d(i, :) = frame_b.distance - frame_a.distance;
        
        C_a = -pinv(frame_a.Hessian);
        C_b = -pinv(frame_b.Hessian);
        
        C_d(i, i) = abs(trace(C_a)) + abs(trace(C_b));
    end
    W = sqrt(inv(C_d));
    
    M_hat = W * M;
    d_hat = W * d;
    
    [U, S, V] = svd(M_hat);
    
    c = 400;
    tol = 10e-7;
    
    rank = 0;
    if S(1,1)^2 > tol
        rank = 1;
        rank = rank + (S(2,2)^2 > (S(1,1)^2 / c));
        rank = rank + (S(3,3)^2 > (S(1,1)^2 / c));
    end
    t = pinv(M_hat) * d_hat;
    T = [0 0 0]';
    for i = 1 : rank
        u_i_hat = U(i, :);
        u_i_hat = u_i_hat ./ norm(u_i_hat);
        v_i_hat = V(:, i);
        v_i_hat = v_i_hat ./ norm(v_i_hat);
        
        T = T + (1.0 / S(i, i)) * dot(u_i_hat, d_hat) * v_i_hat;
    end
else
    T = [0 0 0]';
end
end

function [R] = Rotation_estimation(Planes_prev, Planes_new, correspondence)
frame_a = Planes_prev;
frame_b = Planes_new;

B = 0;
Z = 0;
for i = 1 : size(correspondence, 1)
    [~, ~, latent_a] = pca(frame_a(correspondence(i, 1)).points); % sigma 값은 바꿔야함
    [~, ~, latent_b] = pca(frame_b(correspondence(i, 2)).points);
    sigma1 = latent_a(3);
    sigma2 = latent_b(3);
    
    B = B + (frame_a(correspondence(i, 1)).normal * frame_b(correspondence(i, 2)).normal') ./ (sigma1 + sigma2);
    Z = Z + cross(frame_a(correspondence(i, 1)).normal, frame_b(correspondence(i, 2)).normal) ./ (sigma1 + sigma2);
end

K = [ (B + B' - diag([1 1 1]) * trace(B)) Z; Z' trace(B)];

[V, D] = eig(K);
eigenvalue = diag(D);
idx = find(eigenvalue == max(eigenvalue));
rot_quaternion = V(:, idx);
R = quat2rotm(rot_quaternion');
end








