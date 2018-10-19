function [] = plotSolutionBBA(c,lwh,m,d)
% plotSolutionBBA(output_CoordLUBin,output_LU_LWH,output_LU_Seq,daBest(bestOne).Veh);
% 作图函数:三维BPP
Veh = d.Veh;
LU = d.LU;
Item = d.Item;
Type = 'LU';
% plot3DStrip(LU,Item,Veh,Type);
% % 展示用参数 m: 1-9: [LU_Bin(1,:); LU_Bin(2,:); SID; LID; PID; LU_Item(1,:); LU_Item(2,:); hLU; LU_VehType];
% m : tmpShow =[9,1,2,3,4,6,5];  + 展示顺序  
% m1 LU_VehType  m2 LU_Bin(1,:) m3 LU_Bin(2,:) m4  SID; m5 LID;  m6
% LU_Item(1,:); m7 PID;  m8 展示甩尾 m9 展示顺序

global ISplotPause ISplotShowType
if ISplotShowType == 1
    % 1: 按照LID区分颜色
    nIDType = unique(m(5,:));
elseif ISplotShowType == 2
    % 2: 按照甩尾与否区分颜色
    nIDType = unique(m(8,:));
end

nColors = hsv(length(nIDType)); %不同类型LU赋予不同颜色
nBin = max(m(2,:)); %bin的个数;

figure('name',num2str([nBin]));
j = 1;
tmp = m(8,1);
for ibin=1:nBin

    subplot(2,ceil((nBin)/2),ibin);
    
    f = m(2,:)==ibin;
    yxz = lwh(:,f);
    coord = c(:,f);
    
    if ISplotShowType == 1
   % 1: 按照LID区分颜色 
   lid  = m(5,f);
   elseif ISplotShowType == 2
   % 2: 按照甩尾与否区分颜色
   lid  = m(8,f);
    end

    % LU在Bin内的顺序  LU_Bin(2,:)
    seq = m(3,f);

    for i=1:numel(seq)
        idx =find(seq==i); % LU进入Bin内的画图顺序
        LUColor = 0.8*nColors(nIDType==lid(idx), : );

        % 当托盘不相连出现时, 中断
        if j~=1 && tmp ~=m(8,j) 
                 if ISplotPause>=0 ,      pause(ISplotPause);   end
        end        
        tmp = m(8,j);
        j=j+1;
        
        % plot
        plotcube(yxz(:,idx)',coord(:,idx)',0.7,LUColor);
        
        % Set the lable and the font size
        axis equal;         grid on;        view(111,33); %view(60,40);
        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        
        % 车辆的长宽高调整到合适的车型
        xlim([0 Veh.LWH(1,unique(m(1,f)))]);         ylim([0 Veh.LWH(2,unique(m(1,f)))]);       zlim([0 Veh.LWH(3,unique(m(1,f)))]);        
    end
end

end
