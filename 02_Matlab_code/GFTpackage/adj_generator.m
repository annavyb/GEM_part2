function [ A ] = adj_generator( Coord, R )
%ADJ_GENERATOR generates the adjacency matrix for the electrode cap from the x y z coordinates of
%the electrodes.
%  INPUT:
%the coordinate matrix, where rows represent the electrodes and
%   columns represent the coordinates x y z; 
% R is the treshold distance from one electrode to another. If the distance
% between these two electrodes is smaller than a threshold -> we have 1 in
% adj matrix, otherwise we have zero.
%OUTPUT: the adjacency matrix

% calculation of a 64x64 matrix, where Distanceij is the distance from electrode i
% to electrode j. This matrix is supposed to be symmetric.
n = size(Coord,1);
Distance = zeros(n,n);
for i = 1:n
    for j = 1:n
        if (i ~= j)
            Distance(i,j) = sqrt((Coord(i,1)-Coord(j,1))^2+(Coord(i,2)-Coord(j,2))^2+(Coord(i,3)-Coord(j,3))^2);
        end
    end
    
end


% computing Adjacency matrix: Aij contains the information if electrode i
% and electrode j are sufficiently close to each other. If the distance
% between electrodes is smaller than a certain threshold -> we put 1 into
% the matrix, otherwise it is 0. If i equals j (information about the distance of electrode with itself)
%then we put a 0

A = zeros(n,n);
r = 45; % the threshold
for m = 1:n

    for l = 1:n
      if (m ~= l)
        if(Distance(m,l)<R)
         A(m,l) = 1;
        else
         A(m,l) = 0;
        end
      end
    end
end

% if(option == 'threshold')
% 
%     figure(1)
%     plot3(Coord(:,1), Coord(:,2), Coord(:,3), 'bo');
%     grid on
%     hold on 
%     [x,y,z]= sphere();
%     h = surfl(x*R+Coord(10,1), y*R+Coord(10,2), z*R+Coord(10,3));
%     set(h, 'FaceAlpha', 0.5)
%     shading interp
%     hold on
%     plot3(Coord(10,1),Coord(10,2),Coord(10,3), 'k*')
%     hold on
% end 
% end

