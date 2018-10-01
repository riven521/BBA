function [] = plot3DBPP(d,par)
% 作图函数:三维BPP

fields = fieldnames(par);
aField = [];
for idx = 1:length(fields), aField = [aField par.(fields{idx})];   end


nIDType = unique(d.LU.ID);
nColors = hsv(length(nIDType)); %不同类型LU赋予不同颜色
nBin = max(d.LU.LU_Bin(1,:)); %bin的个数;

for ibin=1:nBin
    figure('name',num2str([ibin, aField]));
    f = d.LU.LU_Bin(1,:)==ibin;
    
    yxz = d.LU.LWH(:,f);
    coord = d.LU.CoordLUBin(:,f);
    lid  = d.LU.ID(:,f);
    seq = d.LU.LU_Bin(2,f);
    
    for i=1:numel(seq)
        idx =find(seq==i); % LU进入Bin内的画图顺序
        LUColor = 0.8*nColors(nIDType==lid(idx), : );
        plotcube(yxz(:,idx)',coord(:,idx)',0.7,LUColor);
        
        % Set the lable and the font size
        axis equal;         grid on;        view(30,30);
        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        xlim([0 d.Veh.LWH(1,1)]);         ylim([0 d.Veh.LWH(2,1)]);       zlim([0 d.Veh.LWH(3,1)]);        
    end
end

end
