function [Plane] = PlaneAnalysis_Pathak2010(Plane)
% --------------------------------------
% Paper : Uncertainty analysis for optimum plane extraction from noisy 3D range-sensor point-clouds
% Paper : Fast Registration Based on Noisy Planes with Unknown Correspondences for 3D Mapping - Appendix
% reference code : Efficient Velodyne SLAM with point and plane features (RefinePlane.C code)
% 
% made by Byung-Hyun Yoon
% data : 2020-06-18
% input Params
%     @Param Plane.points
%     @Param Plane.normal
%     @Param Plane.distance
% output Params
%     @Param Plane.points
%     @Param Plane.normal
%     @Param Plane.distance
%     @Param Plane.covariance
%     @Param Plane.hessian
%     @Param Plane.eigenValues
%     @Param Plane.rhoCovariance
%     @Param Plane.normalCovariance
% --------------------------------------

n_points = size(Plane.points, 1);
centroid = zeros(3, 1);
sumOfWeights = 0;

for i = 1 : n_points
    w_i = 1.0 / (0.02 + norm(Plane.points(i, :)) );%* ( 1.0 / 1000.0 ) 
    centroid = centroid + w_i * Plane.points(i, :)';
    sumOfWeights = sumOfWeights + w_i;
    weights(i) = w_i;
end

centroid = centroid / sumOfWeights;

M = zeros(3, 3);
for i = 1 : n_points
    p_m_c = Plane.points(i, :)' - centroid;
    M = M + (weights(i) * (p_m_c * p_m_c'));
end

[V, D] = eig(M);
n = V(: , 1);

rho = dot(n, centroid);
Plane.eigenValues = [D(1,1) D(2,2) D(3,3)];

H_dd = -sumOfWeights;
H_nd = -H_dd * centroid;
H_nn = -M + H_dd * centroid * centroid' + (n' * M * n) * diag([1 1 1]);

Plane.Hessian(1:3, 1:3) = H_nn;
Plane.Hessian(1:3, 4) = H_nd;
Plane.Hessian(4, 1:3) = H_nd';
Plane.Hessian(4, 4) = H_dd;

Plane.Covariance = - pinv(Plane.Hessian);

H_nn_inv = inv(H_nn);
den = n' * H_nn_inv * H_nd;
Plane.rhoCovariance = - ( n' * H_nn_inv * n) / (den * den);

H_nn_prime = H_nn - (1.0 / H_dd) * (H_nd * H_nd');
[V, D] = eig(H_nn_prime);

test = zeros(3,3);
test = test + V(:, 3) * V(:, 3)' / D(3, 3);
test = test + V(:, 2) * V(:, 2)' / D(2, 2);
Plane.normalCovariance = test;

rho = -rho;
if ( dot(-centroid, n) < 0 ) 
    n = -n;
    rho = dot(-n, centroid);
end

Plane.distance = rho;
Plane.normal = n;
Plane.centroid = centroid;
Plane.pointCovariance = M ./ sumOfWeights;

debug = false;
if debug
    figure, pcshow(Plane.points), xlabel('x'), ylabel('y'), zlabel('z'); hold on;
    X = [centroid(1), centroid(1) + n(1)];
    Y = [centroid(2), centroid(2) + n(2)];
    Z = [centroid(3), centroid(3) + n(3)];
    line(X, Y, Z);
end
end







































