%% GET STRIP �������
% 5 Strip.nbItem % ����������ֵ, ����ITEM�ĶѶ���� ��ͷ�ڷ����� -1�����strip

%% ����
function   [Strip] = cpuStripnbItem(Strip,Item,Veh)
Strip.nbItem = ones(size(Strip.Weight))*-1;   %��STRIP�ڲ�ITEM���͸���, �����Ĭ��Ϊ-1

%% 5,6,7
%Strip.nbItem: ���:-1; ����: ��ӦStrip�ڲ���Item��nbLID���͸���,��ֵԽ��,����LU����Խ��
LIDinItemsArray = cellfun(@(x) x(1), Item.LID); % arrayAllLID: ����ITEM��Ӧ��LIDֵ ������ʽ

uniItem = unique(Item.Item_Strip(1,:));
for i=1:length(Strip.nbItem)
    if ~Strip.isMixed(1,i) %���ǵ�����        
        cellLID = Item.LID(Item.Item_Strip(1,:) == uniItem(i)); % cellLID: ��Strip�ڵ�ITEM��Ӧ��LIDֵ
        LIDinThisItemArray = cellfun(@(x)x(1), cellLID);
        if isscalar(unique(LIDinThisItemArray))
            
            Strip.nbItem(1,i) = sum(LIDinItemsArray == unique(LIDinThisItemArray));
            
        else
             error('������STRIP�ڵ�ITEM�����Ͳ�ͬ'); %arrayLID            
        end
    end
end
end
