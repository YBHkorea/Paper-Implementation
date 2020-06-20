function Lines2D = Line_Detection_2D(images)
rng(1);
debug = false;
%hyper parameter
threshold = 30;

for i = 1 : numel(images)
    images{i} = rgb2gray(images{i});
    
    BW{i} = edge(images{i}, 'Canny');
    
    [H{i}, theta{i}, rho{i}] = hough(BW{i});
    
    P{i}  = houghpeaks(H{i}, 30, 'threshold', ceil(0.3*max(H{i}(:))));
    
    lines{i} = houghlines(BW{i}, theta{i}, rho{i}, P{i}, 'FillGap', 10, 'MinLength', threshold);
    
    if debug
        figure, imshow(BW{i}), hold on
        max_len = 0;
        for k = 1:length(lines{i})
            xy = [lines{i}(k).point1; lines{i}(k).point2];
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
            
            % Plot beginnings and ends of lines
            plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
            plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
            
            % Determine the endpoints of the longest line segment
            len = norm(lines{i}(k).point1 - lines{i}(k).point2);
            if ( len > max_len)
                max_len = len;
                xy_long = xy;
            end
        end
    end
end

Lines2D = lines;

end