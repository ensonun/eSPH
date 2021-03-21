function nw = wallNormal(wall,kernel_type,h)
%WALLNORMAL Summary of this function goes here
%   Detailed explanation goes here

% Range search
j_id = rangesearch(wall(:,1:2),wall(:,1:2),h);

nw = zeros(size(wall(:,1:2)));

for i = 1:size(wall,1)
    for j = j_id{i}(2:end)
        rij = wall(i,1:2)-wall(j,1:2);
        nw(i,:) = nw(i,:) - wDer(norm(rij),h,kernel_type)*rij/norm(rij);        
    end
end
nw = nw./norm(nw);

if any( isnan(nw) )
    error('Nan in wall normal vector!')
end
end