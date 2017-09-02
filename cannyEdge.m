 function E = cannyEdge(I)
%% Canny edge detector
% Author: Harshal Godhia
%  Input: A color image I = uint8(X, Y, 3), where X, Y are two dimensions of the image
%  Output: An edge map E = logical(X, Y)
%%  To DO: Write three functions findDerivatives, nonMaxSup, and edgeLink to fulfill the Canny edge detector.

%% Convert the color image to gray scale image
%  Output I = uint8(X, Y)
try
    I = rgb2gray(I);
catch
end

%% Construct 2D Gaussian filter
Gx = normpdf([-5:1:5], 0, 1);
Gy = normpdf([-5:1:5], 0, 1)';
%% Compute magnitutde and orientation of derivatives
%  J = double(X, Y), the magnitude of derivatives
%  theta = double(X, Y), the orientation of derivatives
%  Jx = double(X, Y), the magnitude of derivatives along x-axis
%  Jy = double(X, Y), the magnitude of derivatives along y-axis
[J, theta, Jx, Jy] = findDerivatives(I, Gx, Gy);
visDerivatives(I, J, Jx, Jy);

%% Detect local maximum
%  M = logical(X, Y), the edge map after non-maximal suppression
M = nonMaxSup(J, theta);
figure; imagesc(M); colormap(gray);

%% Link edges
%  E = logical(X, Y), the final edge map
E = edgeLink(M, J, theta);
figure; imagesc(E); colormap(gray);
end

 function [J, theta, Jx, Jy] = findDerivatives(I, Gx, Gy)
    I = double(I);
    Gfilter = conv2(Gx, Gy, 'full');
    dx = [1 -1];
    dy = [1 -1]';
    %Ismooth = conv2(I, Gfilter, 'same');
    Gfilterx = conv2(Gfilter, dx, 'same');
    Gfiltery = conv2(Gfilter, dy, 'same');
    
    Jx = conv2(I, Gfilterx, 'same');
    Jy = conv2(I, Gfiltery, 'same');
    
    %Jx = conv2(Jx, Gy, 'same');
    %Jx = conv2(Jx, Gx, 'same');
    
    %Jy = conv2(Jy, Gx, 'same');
    %Jy = conv2(Jy, dy, 'same');
    
    J = sqrt(Jx.*Jx + Jy.*Jy);
    theta = atan2(Jy, Jx);
 end
 
 function M = nonMaxSup(J, theta)
    [nr,nc] = size(J);
    M = false(nr,nc);
    [Ymesh1, Xmesh1] = meshgrid(1:nc, 1:nr);
    [Ymesh2, Xmesh2] = meshgrid(1:nc, 1:nr);
    
    Xco1 = Xmesh1 - sin(theta);
    Yco1 = Ymesh1 + cos(theta);
    interp_p1s = interp2(J', Xco1, Yco1);
    
    Xco2 = Xmesh2 + sin(theta);
    Yco2 = Ymesh2 - cos(theta);
    interp_p2s = interp2(J', Xco2, Yco2);
    
    for i = 1:nr
        for j = 1:nc
            pa = interp_p1s(i, j);
            pb = interp_p2s(i, j);
            if J(i,j) > pa && J(i,j) > pb
                M(i,j) = 1;
            end
        end
    end
 end
 
 function E = edgeLink(M, J, theta)
 mean_val = mean(J(:));
 J = J.*M;  %use the edgemap and magnitude to modify edge
 %fprintf('max val: %.2f\n', mean_val);
 low_th_const = 0.55*mean_val;
 high_th_const = 1.55*mean_val;
 %fprintf('lowh:%.2f\n', low_th_const);
 %fprintf('high:%.2f\n', high_th_const);
 E = false(size(J));
 low_ind = J < low_th_const;
 high_ind = J > high_th_const;
 %fprintf('no of high %i ', nnz(high_ind));
 %fprintf('no of low %i ', nnz(low_ind));
 E(low_ind) = 0;
 E(high_ind) = 1;
 
 J(low_ind) = 0;
 J(high_ind) = 1;
 %fprintf('inter elem %.2f\n', nnz(J > 0 & J ~= 1 ));
 %fprintf('beforelo:%.2f\n', nnz(E == 1));
 
 for x = 1:size(J,1)
     for y = 1:size(J,2)
         if J(x,y) == 1 %if this is a strong pixel, extend it to neighbors who are weak
             visited = zeros(size(J));
             [E, ~] = extendEdge(E, J, theta, x, y, visited);
         end
     end
 end
 %fprintf('after:%.2f\n', nnz(E == 1));
 end
 
 function p = getDiscretePoints(i, j, ang)
    if ang < pi/4 && ang >= 0
        p = [i j+1;i-1 j+1];
    elseif ang >= pi/4 && ang < pi/2
        p = [i-1 j;i-1 j+1];
    elseif ang >= pi/2 && ang < 3*pi/4
        p = [i-1 j;i-1 j-1];
    elseif ang >= 3*pi/4 && ang < pi
        p = [i j-1;i-1 j-1];
    elseif ang >= -pi/4 && ang < 0
        p = [i j+1;i+1 j+1];
    elseif ang >= -pi/2 && ang < -pi/4
        p = [i+1 j;i+1 j+1];
    elseif ang >= -3*pi/4 && ang < -pi/2
        p = [i+1 j;i+1 j-1];
    elseif ang >= -pi && ang < -3*pi/4
        p = [i+1 j-1;i j-1];
    else
        p = [];
    end
 end
 %here x, y are locations of a strong pixel whose edge we extend
 function [E, visited] = extendEdge(E, J, theta, i, j, visited)
    ang1 = theta(i,j) + pi/2;
    ang2 = theta(i,j) - pi/2;
    p1 = getDiscretePoints(i, j, ang1);
    p2 = getDiscretePoints(i, j, ang2);
    points = cat(1, p1, p2);
    for i=1:size(points, 1)
        point = points(i,:);
        if point(1) <= size(E,1) && point(1) >= 1 && point(2) <= size(E, 2) && point(2) >= 1
            %point
            cur_val = J(point(1), point(2));
            if cur_val == 1 || cur_val == 0
                return
            elseif cur_val > 0 && visited(point(1), point(2)) ~= 1
                visited(point(1), point(2)) = 1;
                E(point(1), point(2)) = 1;
                [E, visited] = extendEdge(E, J, theta, point(1), point(2), visited);
            end
        end
    end
 end