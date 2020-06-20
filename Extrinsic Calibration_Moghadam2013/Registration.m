function Parameters = Registration(Lines2D, Lines3D, cameraParams)

global M  a  b  K L;

for i = 1 : 3
    K = cameraParams{i}.IntrinsicMatrix';
    [cnt, ~] = size(Lines2D{i});
    P = [0 0 0 0 0 0];
    
    for j = 1 : cnt
        if Lines2D{i}(j, 1) == Lines2D{i}(j, 3)
            a(j) = inf;
            b(j) = nan;
        elseif Lines2D{i}(j, 2) == Lines2D{i}(j, 4)
            a(j) = 0;
            b(j) = Lines2D{i}(j, 2);
        else
            a(j) = (Lines2D{i}(j, 4) - Lines2D{i}(j, 2)) / (Lines2D{i}(j, 3) - Lines2D{i}(j, 1));
            b(j) = Lines2D{i}(j, 2) - a(j) * Lines2D{i}(j, 1);
        end
        %[a(j) -1 b(j)] %Line equation
        L(j, 1) = a(j);
        L(j, 2) = -1;
        L(j, 3) = b(j);
    end
    
    M = Lines3D{i};
    
    func = @(p)ErrorVector(p);
    
    Parameters{i} = lsqnonlin(func, P);
end

end

function E = ErrorVector(P)
global M  a  b  K L;
E = [];

m1 = projection(M(:, 1:3), K, P);
m2 = projection(M(:, 4:6), K, P);
E1 = diag(L * m1);
E2 = diag(L * m2);
E = [E; E1; E2];

end

function m = projection(M, K, P)
[n, ~] = size(M);
m = projMatrix(P) * [M'; ones(1, n)];
m = K * m(1:3, :);
m = m(1:3, :) ./ m(3, :);
end

function RT = projMatrix(P)
R = rod2dcm([P(1) P(2) P(3)]);
T = [P(4); P(5); P(6)];
RT = [R T];
RT = [RT; 0 0 0 1];
end
