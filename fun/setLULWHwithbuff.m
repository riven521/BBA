function [LU] = setLULWHwithbuff(LU,Veh)
% setLULWHwithbuff ==> 修订LU的LWH数据 （1：增加margin；2：确定托盘旋转Rotaed）
        
    % 1 UPDATE : 更新LU的Rotaed标记( LU横放或竖放 却决于距离边界的宽度 哪个小选哪个
    if any(LU.Rotaed~=0)
        error('run此函数前，其旋转状态不应有，而由如下函数获取'); % 如果有，问题应该也不大，因为下面用不到该参数
    end
    [LU.Rotaed]= getLuRotaedBBA(LU.LWH,LU.isRota,LU.margin,Veh.LWH(1,1));  %第二个参数：  3按VEH车辆左右摆放的缝隙最小排序
    
    % 2 UPDATE : 含margin的且含Rotaed的LU的LWH   
    [LU.LWH,LU.OLWH] = LWHbuffer(LU.LWH, LU.margin, LU.Rotaed);   % V1:    LU.LWH = getRotaedLWH(LU.LWH, LU.Rotaed, LU.margin);

    % LU.OLWH = LWHunbuffer(LU.LWH, LU.margin, LU.Rotaed);


    % v1 备份：采用其它方式确认LU是否旋转： 默认将LU全部采用Horizontal方向旋转（前提：该LU允许旋转）
    % NOTE: 此处将获得1: Horizontal方向的LWH和是否Rotaed标记
    % NOTE : 直接替换了原始ORIGINAL 的 LWH

%     if pwhichSortItemOrder ==1
%         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,1); %第二个参数：1: Hori; 0: Vert；其它: 原封不动        
%     elseif pwhichSortItemOrder ==2
%         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,0); %第二个参数：1: Hori; 0: Vert；其它: 原封不动
%     elseif pwhichSortItemOrder ==3 %默认此选项
%         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,Veh.LWH(1,1)); %第二个参数：  3按VEH车辆左右摆放的缝隙最小排序
%     end
end



