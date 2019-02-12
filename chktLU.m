%% 功能：核验是否数据是上轻下重；是否Z高度与ItemSEQ一致；是否ITEMID与XY坐标一致;
function chktLU(t) % t必定是托盘类型的table或struct
%% 0 结构体转换为table 预处理
if isstruct(t)
    t = struct2table(structfun(@(x) x', t,'UniformOutput',false));  end

% t.Properties.VariableNames
                                                                            % tLU1 = t(:,{'ID','LWH','LID','SID','PID','Weight',,'LU_Bin','LU_Item','isShuaiWei','Rotaed'});
%% 1 CHECK Weight 上轻下重 通过相同ITEMID 其X Y 坐标判断 应是金标准
% 1.1 初步核查LU_Item 如不存在获取子ITEMID 
if any(strcmp('LU_Item', t.Properties.VariableNames))
        t.ITEMID = t.LU_Item(:,1);
        t.ITEMSEQ = t.LU_Item(:,2); end

                                                                                                                            % if any(strcmp('ITEMID', t.Properties.VariableNames))
                                                                                                                            %     if any(t.LU_Item(:,1) ~= t.ITEMID) || any(t.LU_Item(:,2) ~= t.ITEMSEQ) %是否存在及检验 可能是老版本, 建议以LU_Item为标准;
                                                                                                                            %             error('LU_Item与ITEMID或ITEMSEQ不同'); end
if any(strcmp('CoordLUBin', t.Properties.VariableNames))
    t.X = t.CoordLUBin(:,1);    t.Y = t.CoordLUBin(:,2);   t.Z = t.CoordLUBin(:,3);   end
    
%% 2 相同ITEMID下的CHEK
if ~any(strcmp('LU_Item', t.Properties.VariableNames))
    return;
else
uniItemID = unique(t.ITEMID(:));
% 2.1 如果坐标CoordLUBin存在,用坐标系判断与ITEMID的差异(1:XY值(坐标是否相同); 2:重量和Z值)
for iItem = 1:length(uniItemID) %对ITEM进行循环
    flagIdx = t.ITEMID==uniItemID(iItem); %对单一ItemId去逻辑值;
    if any(strcmp('Y', t.Properties.VariableNames))
        vX = t{flagIdx,'X'};
        vY = t{flagIdx,'Y'};
        if any(vX ~= vX(1)) || any(vY ~= vY(1))
            error('相同ITEM,但X或Y坐标错位'); end
        
        v = t(flagIdx,{'Z','Weight','ITEMSEQ'});
    else
        v = t(flagIdx,{'Weight','ITEMSEQ'});
    end
    
    v = sortrows(v,'ITEMSEQ');
    
    if any(strcmp('Y', t.Properties.VariableNames))
        % 如果重量不是递减或Z高度不是递增
        if ~issorted(v.Z,'ascend') || ~issorted(v.Weight,'descend')
            issorted(v.Z,'ascend');  issorted(v.Weight,'descend');
            error('相同ITEM,但重量不是递减或Z高度不是递增'); end
    else
        if ~issorted(v.Weight,'descend')
            uniItemID(iItem)
            v;
            error('相同ITEM,但重量不是递减'); end
        end
end
end



%% V1 先判断后循环;
% % %% 2 相同ITEMID下的CHEK
% % uniItemID = unique(t.ITEMID(:));
% % % 2.1 如果坐标CoordLUBin存在,用坐标系判断与ITEMID的差异(1:XY值(坐标是否相同); 2:重量和Z值)
% % if any(strcmp('CoordLUBin', t.Properties.VariableNames))
% %     if ~any(strcmp('Y', t.Properties.VariableNames))
% %         t.X = t.CoordLUBin(:,1);    t.Y = t.CoordLUBin(:,2);   t.Z = t.CoordLUBin(:,3);   end
% % if any(strcmp('Y', t.Properties.VariableNames))
% %     for iItem = 1:length(uniItemID) %对ITEM进行循环
% %         % a 对坐标判断
% %         flagIdx = t.ITEMID==uniItemID(iItem); %对单一ItemId去逻辑值;
% %         vX = t{flagIdx,'X'};
% %         vY = t{flagIdx,'Y'};
% %         if any(vX ~= vX(1)) || any(vY ~= vY(1))
% %             error('相同ITEM,但X或Y坐标错位'); end
% %         % b 对重量按坐标判断
% %         v = t(flagIdx,{'Z','Weight','ITEMSEQ'});   v = sortrows(v,'ITEMSEQ');
% %         % 如果重量不是递减或Z高度不是递增
% %         if ~issorted(v.Z,'ascend') || ~issorted(v.Weight,'descend')
% %             issorted(v.Z,'ascend')
% %             issorted(v.Weight,'descend')
% %             error('相同ITEM,但重量不是递减或Z高度不是递增'); end
% %     end
% % end
% % % 2.2 如果坐标CoordLUBin不存在,用序号判断ITEMID的差异(2:重量)
% % else 
% %     for iItem = 1:length(uniItemID)
% %         flagIdx = t.ITEMID==uniItemID(iItem); %对单一ItemId去逻辑值;
% %         v = t(flagIdx,{'Weight','ITEMSEQ'});   v = sortrows(v,'ITEMSEQ');
% %         % 如果重量不是递减或Z高度不是递增
% %         if ~issorted(v.Weight,'descend')
% %             issorted(v.Weight,'descend')
% %             error('相同ITEM,但重量不是递减'); end
% %     end
% % end