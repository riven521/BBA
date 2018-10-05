function [] = plotSolutionBBA(c,lwh,m,Veh)
% 作图函数:三维BPP

% 1: 按照LID区分颜色
 nIDType = unique(m(5,:));

nColors = hsv(length(nIDType)); %不同类型LU赋予不同颜色
nBin = max(m(2,:)); %bin的个数;

figure('name',num2str([nBin]));
for ibin=1:nBin

    subplot(2,ceil((nBin)/2),ibin);
    
    f = m(2,:)==ibin;
    yxz = lwh(:,f);
    coord = c(:,f);
    
   lid  = m(5,f);

    % LU在Bin内的顺序
    seq = m(3,f);
    
    for i=1:numel(seq)
        idx =find(seq==i); % LU进入Bin内的画图顺序
        LUColor = 0.8*nColors(nIDType==lid(idx), : );
        plotcube(yxz(:,idx)',coord(:,idx)',0.7,LUColor);
        
        % Set the lable and the font size
        axis equal;         grid on;        view(111,33); %view(60,40);
        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        
        % 车辆的长宽高调整到合适的车型
        xlim([0 Veh.LWH(1,unique(m(1,f)))]);         ylim([0 Veh.LWH(2,unique(m(1,f)))]);       zlim([0 Veh.LWH(3,unique(m(1,f)))]);        
    end
end

end
