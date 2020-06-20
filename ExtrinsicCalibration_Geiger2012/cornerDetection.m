function corner = cornerDetection(image)

debug = false;

addpath('function');

image = double(rgb2gray(image)) / 255;
[IH, IW] = size(image);
%% Hyper Parameter
tau_corner = 0.02;
n_nms = 5; %non-maxima-suppression
tau_nms = 0.02;


%% corner prototype
patch_size = 11; %even value

mask = zeros([patch_size patch_size]);
center = 6;
for i = 1 : patch_size
    for j = 1 : patch_size
        mask(i, j) = abs(i - center) + abs(j - center);
    end
end
% mask = 1 - mask ./ max(max(mask));

for i = 1 : 2
    for j = 1 : 4
        kernel{i, j} = zeros([patch_size patch_size]);
    end
end

kernel{1, 1}(1:5, 7:11) = mask(1:5, 7:11);
kernel{1, 2}(7:11, 1:5) = mask(7:11, 1:5);
kernel{1, 3}(1:5, 1:5) = mask(1:5, 1:5);
kernel{1, 4}(7:11, 7:11) = mask(7:11, 7:11);

for i = 1 : patch_size
    for j = 1 : patch_size
        if (i > center) && (i > j) && (i+j-1 > patch_size) %rotate B
            kernel{2, 2}(i, j) = mask(i, j);
        elseif (j < center) && (i > j) && (i+j-1 < patch_size) %rotate D
            kernel{2, 4}(i, j) = mask(i, j);
        elseif (j > center) && (i < j) && (i+j-1 > patch_size) %rotate C
            kernel{2, 3}(i, j) = mask(i, j);
        elseif (i < center) && (i < j) && (i+j-1 < patch_size)%rotate A
            kernel{2, 1}(i, j) = mask(i, j);
        end
    end
end

for i = 1 : 2
    for j = 1 : 4
        kernel{i, j} = kernel{i, j} ./ max(max(kernel{i, j}));
    end
end

if debug
    figure,
    subplot(2, 4, 1), imshow(kernel{1, 1}), title('Corner Prototype 1 - A');
    subplot(2, 4, 2), imshow(kernel{1, 2}), title('Corner Prototype 1 - B');
    subplot(2, 4, 3), imshow(kernel{1, 3}), title('Corner Prototype 1 - C');
    subplot(2, 4, 4), imshow(kernel{1, 4}), title('Corner Prototype 1 - D');
    subplot(2, 4, 5), imshow(kernel{2, 1}), title('Corner Prototype 2 - A');
    subplot(2, 4, 6), imshow(kernel{2, 2}), title('Corner Prototype 2 - B');
    subplot(2, 4, 7), imshow(kernel{2, 3}), title('Corner Prototype 2 - C');
    subplot(2, 4, 8), imshow(kernel{2, 4}), title('Corner Prototype 2 - D');
end

%% Coner likelihood c
for i = 1 : 2
    for j = 1 : 4
        f{i, j} = filter2(kernel{i, j}, image);
    end
end

u1 = 0.25 * (f{1, 1} + f{1, 2} + f{1, 3} + f{1, 4});
u2 = 0.25 * (f{2, 1} + f{2, 2} + f{2, 3} + f{2, 4});


s{1, 1} = min(min(f{1, 1}, f{1, 2}) - u1, u1 - min(f{1, 3}, f{1, 4}));
s{1, 2} = min(u1 - min(f{1, 1}, f{1, 2}), min(f{1, 3}, f{1, 4}) - u1);
s{2, 1} = min(min(f{2, 1}, f{2, 2}) - u2, u2 - min(f{2, 3}, f{2, 4}));
s{2, 2} = min(u2 - min(f{2, 1}, f{2, 2}), min(f{2, 3}, f{2, 4}) - u2);

c = max(max(s{1, 1}, s{1, 2}), max(s{2, 1}, s{2, 2}));

c(c<2.3) = 0;

if debug
    figure, imshow(c), title('Corner likelihood map');
end

%% non maxima suppression
[row, col] = nonmaxsuppts(c, 'radius', n_nms, 'thresh', tau_nms);

corner_candidate = [col row];

if debug
    figure, imshow(image);
    hold on
    scatter(col, row);
end

corner = corner_candidate;

% %% weighted orientation histogram
% % [hog, ~, ptVis] = extractHOGFeatures(image, corner_candidate, 'CellSize', [4 4], 'BlockSize', [1 1], 'NumBins', 32);
% [hog, ~, ptVis] = extractHOGFeatures(image, corner_candidate, 'CellSize', [8 8], 'BlockSize', [1 1], 'NumBins', 32);
% % [hog, ~, ptVis] = extractHOGFeatures(image, corner_candidate, 'CellSize', [12 12], 'BlockSize', [1 1], 'NumBins', 32);
% 
% if debug
%     figure, imshow(image);
%     hold on;
%     plot(ptVis,'Color','green');
% end
% 
% %% Mean-shift
% [alpha1, alpha2] = meanShift(hog);
% 
% %% normalized cross-correlation(template and gradient magnitude)
% [Gmag, ~] = imgradient(image, 'prewitt');
% if debug
%     figure, imshow(Gmag), title('gradient magnitude');
% end
% 
% [n_points, ~] = size(hog);
% 
% for i = 1 : n_points
%     if row(i) < 5 && row(i) > IH-5 && col(i) < 5 && col(i) > IW-5
%         cidx(i) = 0;
%         continue;
%     end
%     template = generateTemplate(alpha1(i), alpha2(i));
%     Similarity = normxcorr2(template, Gmag(row(i)-5:row(i)+5, col(i)-5:col(i)+5));
%     score(i) = Similarity(11, 11);
%     score(i) = score(i) * c(row(i), col(i));
%     
%     if score(i) > tau_corner
%         cidx(i) = 1;
%     else
%         cidx(i) = 0;
%     end
% end
% 
% corner = corner_candidate(find(cidx), :);
% 
% if debug 
%         figure, imshow(image);
%     hold on
%     scatter(corner_final(:, 1), corner_final(:, 2));
% end

end


function [alpha1, alpha2] = meanShift(hog)
%our goal is detecting alpha
[n_points, ~] = size(hog);

[alpha1_val, alpha1] = max(hog, [], 2, 'linear');
alpha1 = floor(alpha1 / n_points +1);
alpha1(n_points) = alpha1(n_points) - 1;

for i = 1 : n_points
    if alpha1(i) == 1
        hog(i, 1:3) = 0;
        hog(i, 31:32) = 0;
    elseif alpha1(i) == 2
        hog(i, 1:4) = 0;
        hog(i, 32) = 0;
    elseif alpha1(i) == 31
        hog(i, 1) = 0;
        hog(i, 29:32) = 0;
    elseif alpha1(i) == 32
        hog(i, 1:2) = 0;
        hog(i, 30:32) = 0;
    else
        hog(i, alpha1(i)-2 : alpha1(i)+2) = 0;
    end
end

[alpha2_val, alpha2] = max(hog, [], 2, 'linear');
alpha2 = floor(alpha2 / n_points +1);
alpha2(n_points) = alpha2(n_points) - 1;
end

function template = generateTemplate(alpha1, alpha2)
template1 = zeros([11 11]);
template2 = zeros([11 11]);

%line : y-6 = m(x-6) -> mx - y + (6 - 6m)
theta1 = (alpha1 * 180/32) - 5;
theta2 = (alpha2 + 180/32) - 5;
m1 = cos(theta1) / sin(theta1);
m2 = cos(theta2) / sin(theta2);

max_dist = sqrt(25 + 25);

for y = 1 : 11
    for x = 1 : 11
        %distance
        template1(y, x) = abs(m1 * x - y + (6 - 6 * m1))/sqrt(m1^2 + 1);
        template2(y, x) = abs(m2 * x - y + (6 - 6 * m2))/sqrt(m2^2 + 1);
    end
end

template = (template1 < 1.5) | (template2 < 1.5);

end
