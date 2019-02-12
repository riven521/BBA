function [LU,Veh] = cpuLULWH(LU,Veh)
% cpuLULWH 重要函数:
% 更新 LU 的LWH
        
    %% 3 更新LU：含margin的LU的LWH；更新LU的Rotaed标记
    % GET LU's LWH with margin and Rotaed or not
    % 2 Input增加间隙BUFF后的feasible的LU和BIN的长宽高转换
    LU.LWH(1,:) =  LU.LWH(1,:) +  LU.margin(1,: ) + LU.margin(2,: ); %宽度（左右）
    LU.LWH(2,:) =  LU.LWH(2,:) +  LU.margin(3,: ) + LU.margin(4,: ); %长度（上下）
    
    [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,Veh.LWH(1,1));  %第二个参数：  3按VEH车辆左右摆放的缝隙最小排序
    LU.LWH = getRotaedLWH(LU.LWH, LU.Rotaed, LU.margin);
    
    % 3 默认将LU全部采用Horizontal方向旋转（前提：该LU允许旋转）
    % NOTE: 此处将获得1: Horizontal方向的LWH和是否Rotaed标记
    % NOTE : 直接替换了原始ORIGINAL 的 LWH
    % LU.LWH

%     if pwhichSortItemOrder ==1
%         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,1); %第二个参数：1: Hori; 0: Vert；其它: 原封不动        
%     elseif pwhichSortItemOrder ==2
%         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,0); %第二个参数：1: Hori; 0: Vert；其它: 原封不动
%     elseif pwhichSortItemOrder ==3 %默认此选项
%         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,Veh.LWH(1,1)); %第二个参数：  3按VEH车辆左右摆放的缝隙最小排序
%     end
end
    