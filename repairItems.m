%% ����1 : isWeightFine: �ж�LU�Ƿ��������ع���
% v3 repairItems
function ITEMSEQ = repairItems(LU)  % t�ض����������͵�table��struct
    
T = getTableLU(LU);

uniItemID = unique(T.ITEMID(:));
for iItem = 1:length(uniItemID)         %��ITEM����ѭ��
    
    flagLUIdx = T.ITEMID==uniItemID(iItem); %�Ե�һItemIdȥ�߼�ֵ;
    
    subT = T(flagLUIdx,{'Weight','ITEMSEQ'});
    
    if isUpDownWeight(subT) 
        T.ITEMSEQ(flagLUIdx) = repairItemWeight(T,flagLUIdx); % ����t��LU_Item
    end
    
end

ITEMSEQ = T.ITEMSEQ';

% LU = getSturctT(T); %���ò����ؽṹ��LU
end

%% ����2 : repairItemWeight: LU�粻����������,�����޸�
% V3: ����struct������table�޸� �����޸�LU.LU_Item�ĵڶ��е�ֵ LU.LU_Item(2,flagLU)
function b = repairItemWeight(T,flagLUIdx)
    tmpWeight=T.Weight(flagLUIdx);
    [~,b] = sort(tmpWeight,'descend');
    [~,b] = sort(b);
    T.ITEMSEQ(flagLUIdx) = b;
end











%% v2 repairItems
% % function LU = repairItems(t)  % t�ض����������͵�table��struct
% %     chktLU(t)
% % % % 0 �ṹ��ת��Ϊtable Ԥ����
% % T = getTableLU(t)
% % T.Properties.VariableNames
% %    t=T
% % 
% % %% 2 ��ͬITEMID�µ�CHEK
% % 
% % uniItemID = unique(t.ITEMID(:));
% % % 2.1 �������CoordLUBin����,������ϵ�ж���ITEMID�Ĳ���(1:XYֵ(�����Ƿ���ͬ); 2:������Zֵ)
% % for iItem = 1:length(uniItemID) %��ITEM����ѭ��
% %     flagIdx = t.ITEMID==uniItemID(iItem); %�Ե�һItemIdȥ�߼�ֵ;
% %     if any(strcmp('Y', t.Properties.VariableNames))
% %         vX = t{flagIdx,'X'};
% %         vY = t{flagIdx,'Y'};
% %         if any(vX ~= vX(1)) || any(vY ~= vY(1))
% %             error('��ͬITEM,��X��Y�����λ'); end
% %         
% %         v = t(flagIdx,{'Z','Weight','ITEMSEQ'});
% %     else
% %         v = t(flagIdx,{'Weight','ITEMSEQ'});
% %     end
% %     
% %     v = sortrows(v,'ITEMSEQ');
% %     
% %     if any(strcmp('Y', t.Properties.VariableNames))
% %         % ����������ǵݼ���Z�߶Ȳ��ǵ���
% %         if ~issorted(v.Z,'ascend') || ~issorted(v.Weight,'descend')
% %             issorted(v.Z,'ascend')  
% %             issorted(v.Weight,'descend')
% %             error('��ͬITEM,���������ǵݼ���Z�߶Ȳ��ǵ���'); end
% %         else
% %         if ~issorted(v.Weight,'descend')
% %                                                             %             uniItemID(iItem)
% %                                                             %             v
% %             t = repairItemWeight(t,uniItemID(iItem)); % ����t��LU_Item
% %                                                             %             v = t(flagIdx,{'Weight','LU_Item'});
% %                                                             %             v = sortrows(v,'LU_Item');
% %         end
% %     end
% % end
% % LU = t;
% % 
% % % ɾ����ʼ����ȡ����
% % if any(strcmp('ITEMID', LU.Properties.VariableNames))
% %         LU.ITEMID = []; end
% % if any(strcmp('ITEMSEQ', LU.Properties.VariableNames))
% %         LU.ITEMSEQ = []; end
% %     if any(strcmp('X', LU.Properties.VariableNames))
% %         LU.X = []; end
% %     if any(strcmp('Y', LU.Properties.VariableNames))
% %         LU.Y = []; end
% %     if any(strcmp('Z', LU.Properties.VariableNames))
% %         LU.Z = []; end
% %     
% % if istable(LU)
% %     LU = table2struct(LU,'ToScalar',true);
% %     LU = (structfun(@(x) x',LU,'UniformOutput',false));
% % end
% % 
% % end






%%
% % %% ����2 : repairItemWeight: LU�粻����������,�����޸�
% % % V2: �����޸�LU.LU_Item�ĵڶ��е�ֵ LU.LU_Item(2,flagLU)
% % function LU = repairItemWeight(oLU,itemIdx)
% %     if istable(oLU)
% %         LU = table2struct(oLU,'ToScalar',true);
% %         LU = (structfun(@(x) x',LU,'UniformOutput',false));
% %     else
% %         LU = oLU;
% %     end
% % 
% %     flagLU = (LU.LU_Item(1,:)==itemIdx);   %�ҳ���item��Ӧ��lu��index
% %     tmpWeight=LU.Weight(:,flagLU);
% %     [~,b] = sort(tmpWeight,'descend');
% %     [~,b] = sort(b);
% % 
% %     LU.LU_Item(2,flagLU) = b;
% % 
% %     if istable(oLU) %���ص�ҲҪ��table
% %         LU = struct2table(structfun(@(x) x',LU,'UniformOutput',false));
% %     end
% % end


%% V1�� isWeightUpDown �޷���ȫ�����������ص�case
% % function Item = isWeightUpDown(Item,LU)
% % for iItem = 1:max(LU.LU_Item(1,:)) %��ITEM����ѭ��
% %     [~,idx] = find(LU.LU_Item(1,:)==iItem);
% %     if iItem==21
% %         1
% %     end
% %     if isempty(idx), 
% %         1
% %     end
% %     nbLUinItem = length(idx);    
% %     % ��ITME�ں�2������LU�Ľ����ж�
% %     if nbLUinItem > 1 %Item������ֻһ��Item,��Ҫ�ж��Ƿ������صı仯
% %         currLUWeight = zeros(1,nbLUinItem);
% %         for iIdx = 1:nbLUinItem
% %             currIdx = idx(LU.LU_Item(2,idx) == iIdx);
% %             currLUWeight(iIdx) = LU.Weight(:,currIdx);          %                 currLUHight(iIdx) = LU.LWH(3,currIdx);
% %         end
% %         diff(currLUWeight)
% %         all(diff(currLUWeight) > 0) 
% %         if all(diff(currLUWeight) > 0) % ������������
% %             % �޸�LU.LU_Item��ֵ             1 5: idx  Ϊ 1 2: LU.LU_Item(2,idx) == iIdx��Ϊ 2 1
% %             Item.isWeightFine(1,iItem) = 0;
% %         else
% %             Item.isWeightFine(1,iItem) = 1;
% %         end
% %     else  %ITEM��ֻ��1��LU, �ض���������
% %         Item.isWeightFine(1,iItem) = 1; 
% %     end
% % end
% % end

%% V1�� repairItemWeight ����insert2Item�и����ʼֵ, �˴��ٽ����޸�
% ****************** �Խ��߿����޸� ************ �ر�
% % % V1: ����insert2Item�и����ʼֵ, �˴��ٽ����޸�
% % % repairItemFull: ������ڷ������case, ����΢��:Ϊ0�ĸ�Ϊ1,��ʣ��߶�С��ITEM�Խ��ߵĸ߶�, ��ΪFULL ����
% % if ~all(Item.isHeightFull)
% %     [~,b] = find(Item.isHeightFull == 0);
% %     for i=1:length(b)
% %            Item = repairItemFull(Item,hVeh,b(i)); %DONE 
% %     end
% % end