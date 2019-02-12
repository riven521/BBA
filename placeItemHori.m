function [flag]  =placeItemHori(LWH,isRota,margin,HoriOrVert)
% placeItemHori ==> 更新LU的Rotaed标记( LU横放或竖放 却决于距离边界的宽度 哪个小选哪个 
%   在每次进入RunAlgorithm前进行 -> 获取LU的当前状态

        % 获取包含margin的LWH
        [LWHbuff] = LWHbuffer(LWH, margin);
        
        % 获取是否旋转的标记flag
        flag1 = LWHbuff(1,:) < LWHbuff(2,:);  %宽<高/长,       1：非Horizontal 0：是Horizontal
        flag2 = isRota == 1;                                    % 1: 允许旋转 0 不允许
        
        if HoriOrVert == 1 %按Horizontal 摆放
            error('111111111111');
            flag = flag1 & flag2;  %允许旋转且目前非Horizontal方式摆放
        elseif HoriOrVert == 0 %按Vertical 摆放
            error('111111111111');
            flag = ~flag1 & flag2;  %允许旋转且目前是Horizontal方式摆放
        elseif HoriOrVert ==2
            error('111111111111');
            flag = false(size(flag1)); %其它参数原封不动
        else
            
            wVeh = HoriOrVert;
            x=mod(wVeh,LWHbuff(1,:));   %按车辆宽度对Hori排放取余数
            y=mod(wVeh,LWHbuff(2,:));   %按车辆宽度对Vert排放取余数
            
            flag3 = x <= y; %是否Hori排放比Vert排放余数更小 1 希望Hori摆放 0 否则
            fxor = ~xor(flag1,flag3); %XOR 异或对比           
            flag =  fxor & flag2;
            
        end
end
    

% % % 将Item进行水平或竖直放置: LU或Item均可以使用
% % function [L, flag]  =placeItemHori(Item,HoriOrVert)
% %         % 获取是否旋转的标记flag
% %         flag1 = Item.LWH(1,:) < Item.LWH(2,:);  %宽<高/长, 1：非Horizontal 0：是Horizontal
% %         flag2 = Item.isRota == 1;                       % 1: 允许旋转 0 不允许
% %         if HoriOrVert == 1 %按Horizontal 摆放
% %             flag = flag1 & flag2;  %允许旋转且目前非Horizontal方式摆放
% %         elseif HoriOrVert == 0 %按Vertical 摆放
% %             flag = ~flag1 & flag2;  %允许旋转且目前是Horizontal方式摆放
% %         else
% %             flag = false(size(flag1)); %其它参数原封不动
% %         end
% % 
% %         %先赋值，后对需要的Item进行修正            
% %         L            = Item.LWH;
% %         L(1,flag) = Item.LWH(2,flag);
% %         L(2,flag) = Item.LWH(1,flag);
% % end
% %    

% % function [ItemHori,flag] = horiOrient(Item,isRota)
% % % 普通函数; 将Item转换为Horizontal Orientation摆放（长边平行与底部）
% % % NOTE : 以前是只改变需要改变的，目前是全部改变 ItemHori = zeros(size(Item));    ItemHori = Item;    
% %  
% %     ItemHori = zeros(size(Item));
% %     ItemHori(1,:) = max(Item);
% %     ItemHori(2,:) = min(Item);
% %     
% %     flag = Item(1,:) < Item(2,:);
% %     flag = double(flag);             %logical转换为double
% %  
% % %     tmpminv = min(Item);
% % %     tmpmaxv = max(Item);
% % %     ItemHori(1,flag) = tmpmaxv(1,flag);   
% % %     ItemHori(2,flag) = tmpminv(1,flag);    
% % end