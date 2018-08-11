function [LU,Item] = updateItemMargin(LU,Item)


nBin = max(Item.Item_Bin(1,:));
for i=1:nBin
% �ҳ���bin�ж�ӦLU��Item���߼����ֵ
    flagLU = LU.LU_Bin(1,:) == i;
    flagItem = Item.Item_Bin(1,:) == i;
% Update LU's LWH in Current Bin
LU.LWH(1,flagLU) = LU.LWH(1,flagLU) - ( LU.margin(1,flagLU) + LU.margin(2,flagLU)); %����
LU.LWH(2,flagLU) = LU.LWH(2,flagLU) - ( LU.margin(3,flagLU) + LU.margin(4,flagLU)); %����
% Update Item's LWH in Current Bin
tmpM = [LU.LWH(1:2, flagLU ); LU.LU_Item(1, flagLU ); LU.ID(1,flagLU)];
[tmpU] = unique(tmpM','rows','stable')' ; %unique���˳��Ҫ���ע�� ��ȡLU��Item�ڵ�˳��
[tmpU] = sortrows(tmpU', [3], {'ascend'})'; % ��ȡItem��1�𽥵�����˳��ֵ
Item.LWH(1:2, flagItem) = tmpU(1:2, :); 


% Update LU's Coord in Current Bin
LU.CoordLUBin(1,flagLU)=LU.CoordLUBin(1,flagLU) + LU.margin(1,flagLU);
LU.CoordLUBin(2,flagLU)=LU.CoordLUBin(2,flagLU) + LU.margin(4,flagLU);
% Update Item's Coord in Current Bin
tmpM = [LU.CoordLUBin(1:2, flagLU ); LU.LU_Item(1, flagLU ); LU.ID(1,flagLU)];
 [tmpU] = unique(tmpM','rows','stable')'; %unique���˳��Ҫ���ע��
 [tmpU] = sortrows(tmpU', [3], {'ascend'})'; % ��ȡItem��1�𽥵�����˳��ֵ
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
% % LU.LWH(1,) = LU.LWH(1,:) - ( LU.margin(1,:) + LU.margin(2,:)); %����
% % LU.LWH(2,:) = LU.LWH(2,:) - ( LU.margin(3,:) + LU.margin(4,:)); %����
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
% % tmpU = unique(tmpM','rows')'; %˳���ܱ䰡....
% % Item.CoordItemBin(1:2, : ) = tmpU(1:2,:);
% % 
% % end
% % 
% % end