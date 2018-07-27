function h = DrawRectangleWithXY(param,style,color) %����color��ɫ hsv
% ���纯��:������ͼ ��������

if (nargin <1),
    error('Please see help for INPUT DATA.');
elseif (nargin==1)
    style='-';
    color=[1,0,0];
elseif (nargin==2)
    color=[1,0,0];
end;
[m,n] = size(param);
if(m ~= 1 || n ~= 5)
    error('param should be an 1x5 array.');
end
if(param(3) <=0 || param(4) <=0)
    error('width and height must be positive values.');
end
a = param(1);
b = param(2);
w = param(3);
h = param(4);
theta = param(5);
X = [-w/2 w/2 w/2 -w/2 -w/2];
Y = [h/2 h/2 -h/2 -h/2 h/2];
P = [X;Y];
ct = cos(theta);
st = sin(theta);
R = [ct -st;st ct];
P = R * P;
% h=plot(P(1,:)+a,P(2,:)+b,style);
h=plot(P(1,:)+a,P(2,:)+b,style,'Color',color); %������ɫ
axis equal;
end