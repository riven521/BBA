function [P] = pgRectangle(x,y,w,h ) 
%由四要素获取四个顶点
% pgon = polyshape(p) plot(pgon)
P = zeros(4,2);
P(1,:) = [x,y];
P(2,:) = [x+w,y];
P(3,:) = [x+w,y+h];
P(4,:) = [x,y+h];
end

% leftA = A(:,1);
% bottomA = A(:,2);
% rightA = leftA + A(:,3);
% topA = bottomA + A(:,4);

%  [x y w h] 
%  rectangle('Position',[1 2 5 6])    
% rectangle
% pgRectangle(1,2,5,6)