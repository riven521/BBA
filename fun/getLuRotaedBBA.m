function [flag]  = getLuRotaedBBA(LWH,isRota,margin,HoriOrVert)
% getLuRotaedBBA ==> 更新LU的Rotaed标记( LU横放或竖放 却决于距离边界的宽度 哪个小选哪个 
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
            
            % 按道理是宽度W剩余小的情况作为判定托盘是否旋转依据
            flag3 = x <= y; %是否Hori排放比Vert排放余数更小 1 希望Hori摆放 0 否则
            %             flag3 = x > y; %基本无法通过, 按道理是宽度W剩余大的情况作为判定托盘是否旋转依据
            fxor = ~xor(flag1,flag3); %XOR 异或对比           
            flag =  fxor & flag2;
            
            % 默认目前均为非旋转（已在运行此函数前确认）
            % 新计算：托盘是否需要旋转
            Rotaed = false(1,length(flag));
            for i = 1: length(flag)
                if flag3(i)==true               % 如果当前排放，宽度剩余gap更小，则无需旋转
                    Rotaed(i) = false;
                else                                   % 如果当前排放，宽度剩余gap更大，则旋转（在能旋转前提下）
                    if flag2(i) == true
                        Rotaed(i) = true;
                    else
                        Rotaed(i) = false;
                    end
                end
            end

            % 两种计算Rotaed对比
            if any(Rotaed ~= flag)
                Rotaed;
                flag;
                warning('两种计算方案结果不一致');
            end
            
            % 返回的还是flage
            flag = (Rotaed);
            
        end
end
    

%% V1: 
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