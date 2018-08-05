% 依据LWH, 原点未基准, 获取pgon
function pgon = getPolyshape(LWH)
    m = 5;
    ox = 0; oy = 0;
    P = [NaN, NaN];
    n = size(LWH,2);
    for i=1:n
        [x,y] = getXY(ox,oy,LWH(1,i),LWH(2,i));
        P = [P; x,y; NaN, NaN];
        ox = ox + LWH(1,i) + m;
    end

    pgon = polyshape(P);    
end

function [x,y] = getXY(ox,oy,w,l)
    x=zeros(4,1);
    y=zeros(4,1);
    % point1 left-bottom 
    x(1) = ox;
    y(1) = oy;
    % point2 left-up
    x(2) = ox;
    y(2) = oy+l;
    % point3 right-up
    x(3) = ox+w;
    y(3) = oy+l;
    % point4 right-bottom
    x(4) = ox+w;
    y(4) = oy;
end