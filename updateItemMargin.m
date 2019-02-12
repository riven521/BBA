function [LU,Item] = updateItemMargin(LU,Item)
% updateItemMargin 更新 堆垛 Item的 边界 Margin

nBin = max(Item.Item_Bin(1,:));
for i=1:nBin
% 找出本bin中对应LU和Item的逻辑标记值
    flagLU = LU.LU_Bin(1,:) == i;
    flagItem = Item.Item_Bin(1,:) == i; %sum后为本bin内Itme的个数
    
% 1 Update LU's LWH in Current Bin TODO : 考虑旋转？
% LU.LWH(1,flagLU) 
LU.LWH(1,flagLU) = LU.LWH(1,flagLU) - ( LU.margin(1,flagLU) + LU.margin(2,flagLU)); %左右
LU.LWH(2,flagLU) = LU.LWH(2,flagLU) - ( LU.margin(3,flagLU) + LU.margin(4,flagLU)); %上下
% LU.LWH(1,flagLU) 

% 2 Update Item's LWH in Current Bin
tmpLULWH = [LU.LWH(1:2, flagLU ); LU.LU_Item(1, flagLU ); LU.ID(1,flagLU)];
% unique: 确保thisBin内的ID, 同一个Item, 同样长宽 的只要唯一值 
[tmpU] = unique(tmpLULWH','rows','stable')' ; %unique后的顺序要务必注意 获取LU在Item内的顺序
[tmpU] = sortrows(tmpU', [3], {'ascend'})'; % 获取Item从1逐渐递增的顺序值

if sum(flagItem)~= size(tmpU,2)
    error('不应该出现的错误');
%     sum(flagItem)
%     size(tmpU,2)
%     tmpU(1:2, :);  Item.LWH(1:2, flagItem);
%     printstruct(LU);
%     printstruct(Item);
end

Item.LWH(1:2, flagItem) = tmpU(1:2, :); %flagItem = 9(Item数量=1); tmpU = 10 

% 3 Update LU's Coord in Current Bin
LU.CoordLUBin(1,flagLU)=LU.CoordLUBin(1,flagLU) + LU.margin(1,flagLU);
LU.CoordLUBin(2,flagLU)=LU.CoordLUBin(2,flagLU) + LU.margin(4,flagLU);

% 4 Update Item's Coord in Current Bin
tmpLUCoord = [LU.CoordLUBin(1:2, flagLU ); LU.LU_Item(1, flagLU ); LU.ID(1,flagLU)];
 [tmpU] = unique(tmpLUCoord','rows','stable')'; %unique后的顺序要务必注意
 [tmpU] = sortrows(tmpU', [3], {'ascend'})'; % 获取Item从1逐渐递增的顺序值
Item.CoordItemBin(1:2, flagItem ) = tmpU(1:2,:);

end

end


% % function [LU,Item] = updateItemMargin(LU,Item)
% % 
% % % Update LU's LWH
% % LU.LU_Bin
% % Item.Item_Bin
% % nBin = max(Item.Item_Bin(1,:));
% % for i=1:nBin
% %     flag = LU.LU_Bin(1,:) == i;
% % LU.LWH(1,) = LU.LWH(1,:) - ( LU.margin(1,:) + LU.margin(2,:)); %左右
% % LU.LWH(2,:) = LU.LWH(2,:) - ( LU.margin(3,:) + LU.margin(4,:)); %上下
% % % Update Item's LWH
% % tmpM = [LU.LWH(1:2, : ); LU.LU_Bin; LU.LU_Item(1, : ); LU.ID];
% % tmpU = unique(tmpM','rows')';
% % Item.LWH(1:2, : ) = tmpU(1:2,:);
% % 
% % % Update LU's Coord
% % LU.CoordLUBin(1,:)=LU.CoordLUBin(1,:) + LU.margin(1,:);
% % LU.CoordLUBin(2,:)=LU.CoordLUBin(2,:) + LU.margin(4,:);
% % % Update Item's Coord
% % tmpM = [LU.CoordLUBin(1:2, : ); LU.LU_Bin; LU.LU_Item(1, : ); LU.ID];
% % tmpU = unique(tmpM','rows')'; %顺序不能变啊....
% % Item.CoordItemBin(1:2, : ) = tmpU(1:2,:);
% % 
% % end
% % 
% % end