%% GET STRIP �������
% 5 Strip.nbItem % ����������ֵ, ����ITEM�ĶѶ���� ��ͷ�ڷ����� -1�����strip
% TODO ����Strip.nbLU

%% ����
% V2: �ƺ���Ϊ�˻���Strip����ԭ��������, ���Ӷ�ʣ��Strip���ж�, ���fֵ 
    % �˴�Strip����,����Item_Strip�е�Strip����,��Ϊ��ʱStrip�ǲ���ĳ��bin�ڵ�Strip
    function   [nbItem, nbLU, nbLULID] = cpuStripnbItem(Strip,Item,LU)
    nbItem = ones(size(Strip.Weight))*-1;    %��STRIP�ڲ�ITEM���͸���, �����Ĭ��Ϊ-1
    nbLU = ones(size(Strip.Weight))*-1;       %��STRIP�ڲ�LU���͸���, �����Ĭ��Ϊ-1
    nbLULID = ones(size(Strip.Weight))*-1;       %��STRIP�ڲ�LULID���͸���, �����Ĭ��Ϊ-1   

    if isfield(Strip,'f')
        % UPDATE LU.nbID AND  Item.nbItem �ر���Ҫ,����re���ȵ���˵
        tmpItemLID = cell2mat(Item.LID);
            LU.nbLID = ones(size(LU.nbLID))*-1;   %�粻����,���ܵ��±���ԭ����
            LU.nbID = ones(size(LU.nbID))*-1;   %�粻����,���ܵ��±���ԭ����
            Item.nbItem= ones(size(Item.nbItem))*-1;  %�粻����,���ܵ��±���ԭ����
        Item.nbItem(Item.f) = sum(tmpItemLID(Item.f) == tmpItemLID(Item.f)');
        LU.nbID(LU.f) = sum(LU.ID(LU.f) == LU.ID(LU.f)');
        LU.nbLID(LU.f) = sum(LU.LID(LU.f) == LU.LID(LU.f)');
    end
    
    for iStrip=1:length(Strip.Weight)
        if ~Strip.isMixed(1,iStrip)  %���ǵ�����
            if isfield(Strip,'f') && Strip.f(iStrip)~=1 % �����re������,������f�ֶ�, ��Ҫ�жϸ�Strip�Ƿ�δ1
                continue; 
            end
            
            % 1 GET Strip.nbLU
            %  TODO : flagLID��ؿ���ɾ��
                    iStripLID = cell2mat(Strip.LID(iStrip));  if ~isscalar(iStripLID), error('������STRIP�ڵ�LID�����Ͳ�ͬ'); end% tmpLID������LID��; ��������˳��; LID��ֻ��1��,����˳����ܶ��
                    if isfield(Strip,'f'),        flagLID(LU.f) = LU.ID(LU.f) == iStripLID;
                    else  flagLID = LU.ID == iStripLID; end         
            flagLID2 = ismember(LU.LU_Strip(1,:),iStrip);
                    if ~any(flagLID2),    error('ismember(LU.LU_Strip(1,:),iStrip)����,LU_Strip��������� ');     end
                  if unique(LU.nbID(flagLID))~= unique(LU.nbID(flagLID2)),    error('����,LU_Strip��������� ');     end
            nbLU(1,iStrip) = unique(LU.nbID(flagLID2)); % ��ͬLID�� ��Ӧ�� nbLID һ����ͬ NOTE: Ҫ����LU.nbLID
            
            % 1.1 GET Strip.nbLLU
            flagLULID = ismember(LU.LU_Strip(1,:),iStrip);
                 if ~any(flagLULID),    error('ismember(LU.LU_Strip(1,:),iStrip)����,LU_Strip��������� ');     end
                 if length(unique(LU.nbLID(flagLULID)))>1
                     nbLULID(1,iStrip) = -1;  % �����1��,������LID��ϵ�ID����ϵ�Strip,
                 else
                     nbLULID(1,iStrip) = unique(LU.nbLID(flagLULID)); % ��ͬLID�� ��Ӧ�� nbLID һ����ͬ NOTE: Ҫ����LU.nbLID
                 end
            
            
            % 2 GET Strip.nbItem
            flagItem = ismember(Item.Item_Strip(1,:),iStrip);
                  if ~any(flagItem),    error('ismember(Item.Item_Strip(1,:),iStrip)����,Item_Strip��������� ');     end
            nbItem(1,iStrip) = unique(Item.nbItem(flagItem)); % NOTE: Ҫ����Item.nbItem
        end
    end
            
    

% % for iStrip=1:length(Strip.Weight)    
% %     if ~Strip.isMixed(1,iStrip)  %���ǵ�����
% %         
% %         if isfield(Strip,'f') && Strip.f(iStrip)~=1 % �����re������,������f�ֶ�, ��Ҫ�жϸ�Strip�Ƿ�δ1
% %                continue;
% %         end
% %         
% %         iStripLID = cell2mat(Strip.LID(iStrip));  % tmpLID������LID��; ��������˳��; LID��ֻ��1��,����˳����ܶ��
% %                 if ~isscalar(iStripLID), error('������STRIP�ڵ�LID�����Ͳ�ͬ'); end
% %         
% %         % 1 GET Strip.nbLU
% %         LU.nbLID = sum(LU.ID==LU.ID');
% %         flagLID = LU.ID == iStripLID;
% %         flagLID2 = ismember(LU.LU_Strip(1,:),iStrip);
% %             if ~any(flagLID2),    error('ismember(LU.LU_Strip(1,:),iStrip)����,LU_Strip��������� ');     end
% %         Strip.nbLU(1,iStrip) = unique(LU.nbLID(flagLID)); % ��ͬLID�� ��Ӧ�� nbLID һ����ͬ NOTE: Ҫ����LU.nbLID
% %         
% %         % 2 GET Strip.nbItem
% %         flagItem = ismember(Item.Item_Strip(1,:),iStrip);
% %             if ~any(flagItem),    error('ismember(Item.Item_Strip(1,:),iStrip)����,Item_Strip��������� ');     end
% %         Strip.nbItem(1,iStrip) = unique(Item.nbItem(flagItem)); % NOTE: Ҫ����Item.nbItem
% %         
% % 
% %          %Strip.nbLU(1,iStrip) = unique(LU.nbLID(  LU.ID(:) == LIDinThisItemArray))
% %          %Strip.nbItem(1,iStrip) = sum(LIDinItemsArray == LIDinThisItemArray);
% %         end
% % end

% 
% function   [nbItem, nbLU] = cpuStripnbItem(Strip,Item,LU)
% Strip.nbItem = ones(size(Strip.Weight))*-1;    %��STRIP�ڲ�ITEM���͸���, �����Ĭ��Ϊ-1
% Strip.nbLU = ones(size(Strip.Weight))*-1;       %��STRIP�ڲ�LU���͸���, �����Ĭ��Ϊ-1

% % %% 7 Strip.nbItem Strip.nbLU
% % % Strip.nbItem: ���:-1; ����: ��ӦStrip�ڲ���Item��nbLID���͸���,��ֵԽ��,����LU����Խ��
% % % tmpItemLID = cell2mat(Item.LID);
% % 
% % % V2: �ƺ���Ϊ�˻���Strip����ԭ��������, ���Ӷ�ʣ��Strip���ж�, ���fֵ 
% % % �˴�Strip����,����Item_Strip�е�Strip����,��Ϊ��ʱStrip�ǲ���ĳ��bin�ڵ�Strip
% % for iStrip=1:length(Strip.Weight)
% %     
% %     if ~Strip.isMixed(1,iStrip)  %���ǵ�����
% %         
% %         if isfield(Strip,'f') && Strip.f(iStrip)~=1 % �����re������,������f�ֶ�, ��Ҫ�жϸ�Strip�Ƿ�δ1
% %                continue;
% %         end
% %         
% %         tmpLID = cell2mat(Strip.LID(iStrip));  % tmpLID������LID��; ��������˳��; LID��ֻ��1��,����˳����ܶ��
% %                 if ~isscalar(tmpLID), error('������STRIP�ڵ�LID�����Ͳ�ͬ'); end
% %         
% %         % 1 GET Strip.nbLU
% %         flagLID = LU.ID == tmpLID;
% %         Strip.nbLU(1,iStrip) = unique(LU.nbLID(flagLID)); % ��ͬLID�� ��Ӧ�� nbLID һ����ͬ
% %         
% %         % 2 GET Strip.nbItem
% %         flagItem = ismember(Item.Item_Strip(1,:),iStrip);
% %         if ~any(flagItem),
% %             error('ismember(Item.Item_Strip(1,:),iStrip)����,Item_Strip��������� '); 
% %         end
% %         Strip.nbItem(1,iStrip) = unique(Item.nbItem(flagItem));
% %          
% %         end
% % end


% V1: ԭʼ���, �����ǲ���Strip, ����ƺ���ȷ.
% % uniStripIdx = unique(Item.Item_Strip(1,:));
% % 
% % LIDinItemsArray = cell2mat(Item.LID); %����ITEM��Ӧ��LIDֵ ������ʽ
% % 
% % unique(LIDinItemsArray)
% % for i=1:length(unique(LIDinItemsArray))
% %     fprintf('LID %1.0f �� %1.0f �� \n', i,sum(LIDinItemsArray(:) == i));    
% % end
% % 
% % for iStrip=1:length(Strip.Weight)
% %     if ~Strip.isMixed(1,iStrip) %���ǵ�����        
% %         LU.ID;
% %         cellLID = Item.LID(Item.Item_Strip(1,:) == uniStripIdx(iStrip)); % cellLID: ��Strip�ڵ�ITEM��Ӧ��LIDֵ
% %         LIDinThisItemArray = cellfun(@(x)x(1), cellLID);
% %         LIDinThisItemArray = unique(LIDinThisItemArray)
% %         if isscalar(LIDinThisItemArray)
% %            
% %             %             Strip.nbLU(1,iStrip) = LU.nbLID(LU.ID(LIDinThisItemArray))
% %             a = unique(LU.nbLID(  LU.ID(:) == LIDinThisItemArray))
% %             Strip.nbLU(1,iStrip) = a
% %             Strip.nbItem(1,iStrip) = sum(LIDinItemsArray == LIDinThisItemArray);
% %             if Strip.nbItem(1,iStrip) > Strip.nbLU(1,iStrip)
% %                 error('strip.nbItem(1,iItem)< Strip.nbLU((1,iItem)'); 
% %             end
% %             
% %         else
% %              error('������STRIP�ڵ�ITEM�����Ͳ�ͬ'); %arrayLID            
% %         end
% %     end
% % end

end
