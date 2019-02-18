%% 功能：核验是否数据是上轻下重；是否Z高度与ItemSEQ一致；是否ITEMID与XY坐标一致;
function chktLU(T) % t必定是托盘类型的table或struct
%% 0 结构体转换为table 预处理
if isstruct(T), T = getTableLU(T); end

%% 1 CHECK LU (TODO 或许后期把initCheck放到此处,仅chk托盘)

%% 2 CHECK ITEM (1:Weight 上轻下重 
if ~any(strcmp('LU_Item', T.Properties.VariableNames)) %如没有,表明还未堆垛,可直接返回
    return;
end

% 2.1 CHECK 堆垛与托盘 关系是否正确 LU_Item
chkA_B(T.LU_Item')

uniItem = unique(T.ITEMID);

% 2.1 CHECK 是否满足LW和ID号一致
if any(strcmp('L', T.Properties.VariableNames)) && any(strcmp('ID', T.Properties.VariableNames))
    
        for iItem = 1:length(uniItem)
        flagLUIdx = T.ITEMID==uniItem(iItem);
        
        subT = T(flagLUIdx,{'L','W','ID'});        
        
        if ~isscalar(unique(subT.L)) || ~isscalar(unique(subT.W)) || ~isscalar(unique(subT.ID))
            error('相同ITEM,但L或W或ID不一致'); end        
        
        end
    
end

% 2.2 CHECK 堆垛重量 满足上轻下重
for iItem = 1:length(uniItem) %对ITEM进行循环
    flagLUIdx = T.ITEMID==uniItem(iItem); %对单一ItemId去逻辑值;
    
    subT = T(flagLUIdx,{'Weight','ITEMSEQ'});
    
    if  isUpDownWeight(subT) %表明不满足
        error('相同ITEM,但重量不是递减');  % T.isWeightUpDown(flagLUIdx) = 1;
    end
    
    %     subT = sortrows(subT,'ITEMSEQ');    
    %     if ~issorted(subT.Weight,'descend')            
    %             error('相同ITEM,但重量不是递减');  % T.isWeightUpDown(flagLUIdx) = 1;
    %     end
end

% 2.3 CHECK 堆垛坐标 是否满足 X Y Z 正确 
% 2.1 如果坐标CoordLUBin存在,用坐标系判断与ITEMID的差异(1:XYZ值)
if any(strcmp('Y', T.Properties.VariableNames))
    
    for iItem = 1:length(uniItem)
        flagLUIdx = T.ITEMID==uniItem(iItem);
        
        subT = T(flagLUIdx,{'X','Y','Z','ITEMSEQ'});
        
        subT = sortrows(subT,'ITEMSEQ');
        
        if ~isscalar(unique(subT.X)) || ~isscalar(unique(subT.Y))
            error('相同ITEM,但X或Y坐标错位'); end
        
        if ~issorted(subT.Z,'ascend')
            error('相同ITEM,但高度不是递增');
        end
        
    end
end



%% 3 LWH1(含margin) LWH2（不含margin） 和 OLWH（含Rotaed，但不含mragin） 核验
    % v1 : 老版本,后期可被下面的V2取代
    OLWH = T.OLWH';    
    Rotaed = T.Rotaed';
    margin= T.margin';
    
    OLWHtmp = OLWH;
    OLWH(1,Rotaed) = OLWHtmp(2,Rotaed);
    OLWH(2,Rotaed) = OLWHtmp(1,Rotaed);
    
    LWH1 = T.LWH';
    LWH2 = LWHunbuffer(LWH1, margin);
    
    if ~isequal(OLWH,LWH1) && ~isequal(OLWH,LWH2)
        error('t的LWH在变化前后有错误，需要核验');
    end
    
    % v2: 也是核验LWH是否正确的
    LWHOLD = T.OLWH';
    LWHNEW = LWHunbuffer(T.LWH', T.margin', T.Rotaed');
    if ~isequal(LWHOLD,LWHNEW) 
        error('t的LWH在变化前后有错误，需要核验');
    end
end
  


%% 局部函数
function chkA_B(array)
if size(array,1) ~= 2, error('输入错误'); end

a = array(1,:);
b = array(2,:);
uniA = unique(a);

if ~issorted(uniA,'strictascend') || uniA(1)~= 1
    error('第一行非从1开始严格递增');
end

for i=1:length(uniA)
    uniB = sort(b(a == uniA(i) ));
    if isempty(uniB) || ~issorted(uniB,'strictascend') || uniB(1)~= 1
        error('第而行非从1开始严格递增');
    end
end

end



%% v2:基于repair的check
% % function chktLU(t) % t必定是托盘类型的table或struct
% % %% 0 结构体转换为table 预处理
% % if isstruct(t)
% %     t = struct2table(structfun(@(x) x', t,'UniformOutput',false));  end
% % 
% % % t.Properties.VariableNames
% %                                                                             % tLU1 = t(:,{'ID','LWH','LID','SID','PID','Weight',,'LU_Bin','LU_Item','isShuaiWei','Rotaed'});
% % %% 1 CHECK Weight 上轻下重 通过相同ITEMID 其X Y 坐标判断 应是金标准
% % % 1.1 初步核查LU_Item 如不存在获取子ITEMID 
% % if any(strcmp('LU_Item', t.Properties.VariableNames))
% %         t.ITEMID = t.LU_Item(:,1);
% %         t.ITEMSEQ = t.LU_Item(:,2); end
% % 
% %                                                                                                                             % if any(strcmp('ITEMID', t.Properties.VariableNames))
% %                                                                                                                             %     if any(t.LU_Item(:,1) ~= t.ITEMID) || any(t.LU_Item(:,2) ~= t.ITEMSEQ) %是否存在及检验 可能是老版本, 建议以LU_Item为标准;
% %                                                                                                                             %             error('LU_Item与ITEMID或ITEMSEQ不同'); end
% % if any(strcmp('CoordLUBin', t.Properties.VariableNames))
% %     t.X = t.CoordLUBin(:,1);    t.Y = t.CoordLUBin(:,2);   t.Z = t.CoordLUBin(:,3);   end
% %     
% % %% 2 相同ITEMID下的CHEK
% % if ~any(strcmp('LU_Item', t.Properties.VariableNames))
% %     return;
% % else
% % uniItemID = unique(t.ITEMID(:));
% % % 2.1 如果坐标CoordLUBin存在,用坐标系判断与ITEMID的差异(1:XY值(坐标是否相同); 2:重量和Z值)
% % for iItem = 1:length(uniItemID) %对ITEM进行循环
% %     flagIdx = t.ITEMID==uniItemID(iItem); %对单一ItemId去逻辑值;
% %     if any(strcmp('Y', t.Properties.VariableNames))
% %         vX = t{flagIdx,'X'};
% %         vY = t{flagIdx,'Y'};
% %         if any(vX ~= vX(1)) || any(vY ~= vY(1))
% %             error('相同ITEM,但X或Y坐标错位'); end
% %         
% %         v = t(flagIdx,{'Z','Weight','ITEMSEQ'});
% %     else
% %         v = t(flagIdx,{'Weight','ITEMSEQ'});
% %     end
% %     
% %     v = sortrows(v,'ITEMSEQ');
% %     
% %     if any(strcmp('Y', t.Properties.VariableNames))
% %         % 如果重量不是递减或Z高度不是递增
% %         if ~issorted(v.Z,'ascend') || ~issorted(v.Weight,'descend')
% %             issorted(v.Z,'ascend');  issorted(v.Weight,'descend');
% %             error('相同ITEM,但重量不是递减或Z高度不是递增'); end
% %     else
% %         if ~issorted(v.Weight,'descend')
% %             uniItemID(iItem)
% %             v;
% %             error('相同ITEM,但重量不是递减'); end
% %         end
% % end
% % end
% % 
% % %% 3 LWH1(含margin) LWH2（不含margin） 和 OLWH（含Rotaed，但不含mragin） 核验
% %     OLWH = t.OLWH';
% %     Rotaed = t.Rotaed';
% %     margin= t.margin';
% %     
% %     
% %     OLWHtmp = OLWH;
% %     OLWH(1,Rotaed) = OLWHtmp(2,Rotaed);
% %     OLWH(2,Rotaed) = OLWHtmp(1,Rotaed);
% %     
% %     LWH1 = t.LWH';
% %     LWH2 = LWHunbuffer(LWH1, margin);
% %     
% %     if ~isequal(OLWH,LWH1) && ~isequal(OLWH,LWH2)
% %         error('t的LWH在变化前后有错误，需要核验');
% %     end
% %     
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