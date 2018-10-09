%% GET STRIP �������
% 5 Strip.nbItem % ����������ֵ, ����ITEM�ĶѶ���� ��ͷ�ڷ����� -1�����strip
% TODO ����Strip.nbLU

%% ����
function   [Strip] = cpuStripnbItem(Strip,Item,LU)
Strip.nbItem = ones(size(Strip.Weight))*-1;    %��STRIP�ڲ�ITEM���͸���, �����Ĭ��Ϊ-1
Strip.nbLU = ones(size(Strip.Weight))*-1;       %��STRIP�ڲ�LU���͸���, �����Ĭ��Ϊ-1

%% 7 Strip.nbItem Strip.nbLU
% Strip.nbItem: ���:-1; ����: ��ӦStrip�ڲ���Item��nbLID���͸���,��ֵԽ��,����LU����Խ��
% tmpItemLID = cell2mat(Item.LID);
% uniStripIdx = unique(Item.Item_Strip(1,:));

% �˴�Stripһ��43��,����Item_Strip�е�Strip����,��Ϊ��ʱStrip�ǲ���ĳ��bin�ڵ�Strip
for iStrip=1:length(Strip.Weight)
    
    if ~Strip.isMixed(1,iStrip)  %���ǵ�����
        
        if isfield(Strip,'f') && Strip.f(iStrip)~=1 % �����re������,��Ҫ�ж���
               continue;
        end
        
        tmpLID = cell2mat(Strip.LID(iStrip));  % tmpLID������LID��; ��������˳��; LID��ֻ��1��,����˳����ܶ��
        if ~isscalar(tmpLID), error('������STRIP�ڵ�LID�����Ͳ�ͬ'); end
        
        % 1 GET Strip.nbLU
        flagLID = LU.ID == tmpLID;
        Strip.nbLU(1,iStrip) = unique(LU.nbLID(flagLID)); % ��ͬLID�� ��Ӧ�� nbLID һ����ͬ
        
        % 2 GET Strip.nbItem
        flagItem = ismember(Item.Item_Strip(1,:),iStrip);
        if ~any(flagItem),
            error('ismember(Item.Item_Strip(1,:),iStrip)����,Item_Strip��������� '); 
        end
        Strip.nbItem(1,iStrip) = unique(Item.nbItem(flagItem));
         
        end
end

% LIDinItemsArray = cell2mat(Item.LID); %����ITEM��Ӧ��LIDֵ ������ʽ
% 
% for iStrip=1:length(Strip.Weight)
%     if ~Strip.isMixed(1,iStrip) %���ǵ�����        
%         LU.ID
%         cellLID = Item.LID(Item.Item_Strip(1,:) == uniStripIdx(iStrip)); % cellLID: ��Strip�ڵ�ITEM��Ӧ��LIDֵ
%         LIDinThisItemArray = cellfun(@(x)x(1), cellLID);
%         if isscalar(unique(LIDinThisItemArray))
%             Strip.nbLU(1,iStrip) = LU.nbLID(LU.ID(unique(LIDinThisItemArray)))
%             Strip.nbItem(1,iStrip) = sum(LIDinItemsArray == unique(LIDinThisItemArray));         
%             if Strip.nbItem(1,iStrip) > Strip.nbLU(1,iStrip)
%                 error('strip.nbItem(1,iItem)< Strip.nbLU((1,iItem)'); 
%             end
%             
%         else
%              error('������STRIP�ڵ�ITEM�����Ͳ�ͬ'); %arrayLID            
%         end
%     end
% end
end
