% Item进行三维作图
function plot3DStrip(LU,Item,Veh,Type)
nIDType = unique(LU.ID);
nColors = hsv(length(nIDType)); %不同类型Item的LU赋予不同颜色

if strcmp(Type,'Item')
    yxz = Item.LWH;
    coord = Item.CoordItemStrip;
    coord(3,:) = 0;
    
    for i=1:numel(Item.Weight)
        j=Item.LID(i);
        LUColor = 0.8*nColors(nIDType==j{1}, : );
        plotcube(yxz(:,i)',coord(:,i)',0.7,LUColor);
        % Set the lable and the font size
        axis equal;         grid on;        view(60,40);
        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        xlim([0 Veh.LWH(1,1)]);         zlim([0 Veh.LWH(3,1)]);
    end
elseif strcmp(Type,'LU')
    yxz = LU.LWH;
    coord = LU.CoordLUStrip;
    
    for i=1:numel(LU.Weight)
        j=LU.ID(i);
        LUColor = 0.8*nColors(nIDType==j(1), : );
        plotcube(yxz(:,i)',coord(:,i)',0.7,LUColor);
        % Set the lable and the font size
        axis equal;         grid on;        view(120,10);  %view(60,40);
        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        xlim([0 Veh.LWH(1,1)]);         zlim([0 Veh.LWH(3,1)]);
    end
end



end