%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: calculates the distances between the sensors. 
% NOT USED IN ANY OF THE SIMULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('distances.mat');

%get distance vectors for each axis
x = distance(:,1);
y = distance(:,2);
z = distance(:,3);

%calculate the distance between each pair of sensors for every axis
%resulting matrices are symmetrical up to a minus sign around the diagonal
dstx = triu(x.'-x);
dsty = triu(y.'-y);
dstz = triu(z.'-z);

%calculate the distance sqrt(dx^2+dy^2+dz^2)
dst = sqrt(dstx.^2 + dsty.^2 + dstz.^2);
max_dist = max(dst, [], 'all');

%fill 0 with inf so we can properly get the minimum
dst_for_min = dst;
dst_for_min(dst_for_min == 0) = inf;
min_dist = min(dst_for_min, [], 'all');
