%% V5 : ���ӿɶ��� V4 ��EID�����ݽṹ������д
% �Ȱ�SID����, �ٰ�EID����, �ٸ���priorityΪSID�ڲ�˳��
% �ص㣺63-65 �� tmpM��˳��

function StripPriorityArray = getOrderofLID(SIDord, EIDord, Strip)

%% Ŀ�� -> ��ȡallPriority, �ɶ��perPriority��Ϲ���
% �ֱ����ͬSID/EID�µ�STRIP����priorityֵ, ÿ�������priority���ظ�,����1��ʼ
% ��ͬSIDord��ȡ��ͬ��˳��, ��SID/EID�ĵ�һ����ʼ

% ��ʼ��
StripPriorityArray = zeros(1,length(SIDord));        

uniSIDord = unique(SIDord);     if ~issorted(uniSIDord,'strictascend'), error('���ϸ����'); end

for i=1:length(uniSIDord)
    
    fsid = SIDord==uniSIDord(i);
    
    uniEIDord = unique(EIDord(fsid));  if ~issorted(uniEIDord,'strictascend'), error('���ϸ����'); end %��SID�µ�EID, Ҳ�������ϸ������
    
    for j=1:length(uniEIDord)
        
        feid = EIDord==uniEIDord(j);
        
        f = fsid & feid;   if ~any(f), error('��SID��EID�����κ�STRIP'); end %����������Ը�Ϊwarning,��Ϊ���ܲ�ͬSID�µ�EID��ͬ
        
        priority = 1;  % 555 priority: ����ֵ, ��Զ����,��Сֵ=1,���ֵΪLU����
        
        perPriority = zeros(1,sum(f));     % perPriority: ����, ÿ��SID/EID�µ�f������priorityֵ,��Ҫ����ֵ        
        
        %% TO DEL �����õ�ĳЩֵ ��Ϊ����STRIP���� �������ƺ��ò������������ο���
        % UNUSED ��ȡ��SID��EID�ڵ�STRIP��Ӧ��ֵtLW;tLL;tL;tID;celltID.
        %         tLW = SLW(:,f);   tLL = SLoadingRateLimit(:,f);  tPured = Spured(:,f);   tSingle = Ssingle(:,f);   %����δ֪        
        % FINALLY USED ���ȳ�ͷ�ڷŵ�˳��
        %         tnbLU = SnbLU(:,f);     tnbLULID = SnbLULID(:,f); tnbItem = SnbItem(:,f);
        %         tMixed = SisMixed(:,f);     tFull = SHeightfullfull(:,f);   tLR = SLoadingrate(:,f);        % max(tL)  %  min(tL)  % mean(tL)

        %% ��ʼ��
        tID = Strip.LID(:,f);    
        celltID = tID;                   %celltID �� ͬһSID/EID�µ���������  %Strip.LIDʵ��ΪLU��ID,ƴ���ж�
        
        if length(tID) == 1    % ���ֻ��1��  == �ȼ���  if sum(f)==1
            mID = tID{1};
            isID = all(mID,2);
        else
            [mID, isID] = padcat(tID{:});  
        end
        
         if iscolumn(mID), mID = mID'; end    % mID��ÿ�� f (����) ��Ӧ��ID���ͺ�; isID: ÿ�� f (����) ��Ӧ��ID���ͺ��Ƿ����(�߼��ж�)
         if iscolumn(isID), isID = isID'; end
        
        % nbIDperStrip = sum(isID,1);
        
        %% 1 ��ȡperPriority: ��ͬSID/EID��, ����STRIP������˳��; �̶������ڰڷŽǶȳ���, ϸ��ÿ���ڲ�STRIP��priority
        
        if sum(f)==1 % 1.1 �����STRIP��ֻ��1��STRIP����Ϊ������STRIP�����ͣ�����ֵȨ��Ϊ1.
            
            perPriority = 1;     if ~isscalar(isID),    warning('�ǻ����strip.  Ҳ��ֵΪ1, Ӧ�����ⲻ��,�Ƿ�����������Ҫ�������;');   end
            
        else  % 1.2 �����STRIP�ĺ��ж��STRIP, �����STRIP���������ȡorder
            
            tmpM = [Strip.nbLU(:,f);  Strip.nbLULID(:,f);  Strip.nbItem(:,f); Strip.isMixed(:,f);  Strip.isHeightFull(:,f);  Strip.loadingrate(:,f); ];
            % ����ļ��������⣬������ todo �����޸�. 
%              tmpM = [Strip.nLUID(:,f);  Strip.nLULID(:,f);  Strip.nbItem(:,f); Strip.isMixed(:,f);  Strip.isHeightFull(:,f);  Strip.loadingrate(:,f); ];
          % �ǻ�����ȣ������������ݼ�������LID�ݼ����Ѷ����ݼ�����ͬ����������������ȣ��� 
          [~,order] = sortrows(tmpM',[4,1,2,3,5,6],{'ascend','descend', 'descend','descend','descend','descend'});          if ~isrow(order), order=order'; end
%             [~,order] = sortrows(tmpM',[1,2,3,4,5,6],{'descend', 'descend', 'descend','ascend','descend','descend'});          if ~isrow(order), order=order'; end
              
%             [~,order] = sortrows(tmpM',[4,1,2,5,6],{'ascend','descend', 'descend','descend','descend'});          if ~isrow(order), order=order'; end          

            %% 2 ����1�İڷ�˳��, ��������STRIP��priority��torder
            while any(perPriority==0) %����SID�����priorityδ��ȫʱ, ��������ѭ��(δȫ��ԭ�����Ҳ������ڵ�strip��,�ٻ���order����һ����ѡstrip)
                
                % SIDorder�з�0��order�������Ǹ���Ϊ��ѡSTRIP
                
                % 2.1 �ҳ�������( ��ϻ�ǻ�� )ǰ����,��һ��o����Ӧ��order(o)
                [~,o]=find(perPriority(order) == 0,1,'first');
                
                perPriority(order(o)) = priority;
                
                priority=priority+1;
                
                % 2.2 �ҳ���order(o)λ��tLID��Ӧ������Strip.
                tnbItem = celltID{:,order(o)};
                
                % NOTE �״�tLID��Ӧ�ó��ֻ����STRIP, ���簴�߶�����, �ǿ��ܳ��ֵ�; Ŀǰ��Ϊ��������򣬷ǻ���ǰ��Ӧ�ò�����ֵ��ˡ�
                if Strip.isMixed(order(o))
                    warning('�״β�ϣ�����ֻ����STRIP,�����ֲ������ҽ�����������STRIP');   % ֮ǰΪwarning
                end
                
                [perPriority,priority] = getAdjPriority(priority,order,perPriority,mID,tnbItem);
                
            end
        end
        
        StripPriorityArray(f) = perPriority;
        
    end % EOF EID 
    
end % EOF SID

% �Ŵ�CODE
 if ~all(StripPriorityArray), error('allPriority ����δ��ֵ��'); end

end

%% �ֲ����� %%
%% V2 : 
%����1: Ѱ��tLID��Ӧ������STRIP,����һ���߼�
function [SIDorder,priority] = getAdjPriority(priority,order,SIDorder,tID,tLID)
     
      % ��tLID�ж����ѭ��Ѱ��
     for p=1:length(tLID)
        atPID = tLID(p);
        getorderofatLID(atPID);
     end
    
    % ��Ե���tPID, Ѱ������STRIP
    function getorderofatLID(tLID)
        if isscalar(tLID)  %���ֻ��1��
            %�ҳ��Ƿ��ж�Ӧ���ID
            % Ѱ��tID�е���: ���н���1��LID; ���ж�Ӧ��SIDorderΪ0�����а��������ڣ�tLID
            [~,nbcol1] = find(tID(:,:)==tLID & sum(~isnan( tID(:,:) ),1) ==1 &  SIDorder==0);
            if isrow(nbcol1),  nbcol1=nbcol1'; end
            % ���ж��nbcol1, ��order��˳���SIDorder��priority��ֵ
            if ~isempty(nbcol1)
                [~,b]= ismember(nbcol1,order);
                tmp= [nbcol1,b];
                [a,~] = sortrows(tmp, 2 ,{'ascend'});
                for i1=1:length(nbcol1)
                    SIDorder(a(i1,1)) = priority;
                    priority=priority+1;
                end
            end
            
            % Ѱ��tID�е���: ������2������LID; ���ж�Ӧ��SIDorderΪ0�����а����������ڣ�tLID
            [~,nbcol2] = find(tID(:,:)==tLID & sum(~isnan( tID(:,:) ),1) ~=1 &  SIDorder==0);
            if ~isempty(nbcol2)
                if ~isscalar(nbcol2)
%                     tID(:,:)
                    % tID(:,:)=tLID
                    warning('nbcol2������ֻһ����: ���ڶ����ͬ��ϵ�STRIP; -> ͬһ����Strip, �����ϵ�����Strip, ��������2��������strip�Ļ��strip,�޴�ѡ���ĸ���Ϊ����strip��ì�� ');   
%                     error('nbcol2������ֻһ����: ���ڶ����ͬ��ϵ�STRIP; -> ͬһ����Strip, �����ϵ�����Strip, ��������2��������strip�Ļ��strip,�޴�ѡ���ĸ���Ϊ����strip��ì�� ');   
                end
                SIDorder(nbcol2) = priority;
                priority=priority+1;
            end
            % Ѱ��nbcol2�г�tLID�������,���ݹ���������STRIP,ֱ��������STRIP
            hyLID = tID(:,nbcol2);
            hyLID(hyLID==tLID) =[];
            if ~isempty(hyLID)
                [SIDorder,priority] = getAdjPriority(priority,order,SIDorder,tID,hyLID);
            end
        else
            error('tLID�������ֵ');
        end
    end
end


%% V1 ����1: Ѱ��tLID��Ӧ������STRIP,����һ���߼�
% % function [SIDorder,priority] = getAdjPriority(priority,order,SIDorder,tID,tLID)
% %      % ��tLID�ж����ѭ��Ѱ��
% %      for p=1:length(tLID)
% %         atPID = tLID(p);
% %         getorderofatLID(atPID);
% %      end
% %     
% %     % ��Ե���tPID, Ѱ������STRIP
% %     function getorderofatLID(tLID)
% %         if isscalar(tLID)  %���ֻ��1��
% %             %�ҳ��Ƿ��ж�Ӧ���ID
% %             % Ѱ��tID�е���: ���н���1��LID; ���ж�Ӧ��SIDorderΪ0�����а��������ڣ�tLID
% %             [~,nbcol1] = find(tID(:,:)==tLID & sum(~isnan( tID(:,:) ),1) ==1 &  SIDorder==0);
% %             if isrow(nbcol1),  nbcol1=nbcol1'; end
% %             % ���ж��nbcol1, ��order��˳���SIDorder��priority��ֵ
% %             if ~isempty(nbcol1)
% %                 [~,b]= ismember(nbcol1,order);
% %                 tmp= [nbcol1,b];
% %                 [a,~] = sortrows(tmp,[2],{'ascend'});
% %                 for i1=1:length(nbcol1)
% %                     SIDorder(a(i1,1)) = priority;
% %                     priority=priority+1;
% %                 end
% %             end
% %             
% %             % Ѱ��tID�е���: ������2������LID; ���ж�Ӧ��SIDorderΪ0�����а����������ڣ�tLID
% %             [~,nbcol2] = find(tID(:,:)==tLID & sum(~isnan( tID(:,:) ),1) ~=1 &  SIDorder==0);
% %             if ~isempty(nbcol2)
% %                 if ~isscalar(nbcol2),  
% %                     tID(:,:)
% %                     tID(:,:)=tLID
% %                     error('nbcol2������ֻһ����: ���ڶ����ͬ��ϵ�STRIP');   
% %                 end
% %                 SIDorder(nbcol2) = priority;
% %                 priority=priority+1;
% %             end
% %             % Ѱ��nbcol2�г�tLID�������,���ݹ���������STRIP,ֱ��������STRIP
% %             hyLID = tID(:,nbcol2);
% %             hyLID(hyLID==tLID) =[];
% %             if ~isempty(hyLID)
% %                 [SIDorder,priority] = getAdjPriority(priority,order,SIDorder,tID,hyLID);
% %             end
% %         else
% %             error('tLID�������ֵ');
% %         end
% %     end
% % end


%% %% V3 : ����cell��ʽ �ҷֿ�ITEM��STRIP���� �ص��ų�order���ظ���ֵ�����
% % % �Ȱ�SID����, �ٸ���priorityΪSID�ڲ�˳��
% % % �ص㣺27-29��tmpM��˳��
% % 
% % % IDorder = getOrderofLID(SIDorder, Strip.isSingleItem, Strip.isAllPured, Strip.nbItem, Strip.isHeightFull,...
% % %                                         Strip.isMixed, Strip.LID, Strip.LW(1:2,:), Strip.loadingrateLimit, Strip.loadingrate);
% %                                     
% % function allPriority = getOrderofLID(SIDord,Ssingle, Spured,SnbItem,SnbLU,SnbLULID, ...
% %     SHeightfullfull, SisMixed, SLID,SLW,SLoadingRateLimit,SLoadingrate)
% % % ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';
% % if ~iscell(SLID),error('�����cell');end
% % 
% % allPriority = zeros(1,size(SLID,2));
% % uniOrd = unique(SIDord);
% % 
% % % ��ͬSID��ȡ��ͬ��˳��, ��SID�����1��ʼ
% % for i=1:length(uniOrd)    
% %     idxSID = SIDord==uniOrd(i);
% %     
% %     % UNUSED ��ȡ��SID�ڵ�STRIP��Ӧ��ֵtLW;tLL;tL;tID;celltID.
% %     tLW = SLW(:,idxSID);   tLL = SLoadingRateLimit(:,idxSID);
% %     tID = SLID(:,idxSID);    celltID = tID;
% %     tPured = Spured(:,idxSID);
% %     tSingle = Ssingle(:,idxSID);   %����δ֪
% %     
% %     % FINALLY USED ���ȳ�ͷ�ڷŵ�˳��
% %     tnbItem = SnbItem(:,idxSID);
% %     tnbLU = SnbLU(:,idxSID);
% %     tnbLULID = SnbLULID(:,idxSID);
% %     
% %     
% %     tMixed = SisMixed(:,idxSID);
% %     tFull = SHeightfullfull(:,idxSID);
% %     tLR = SLoadingrate(:,idxSID);        % max(tL)  %  min(tL)  % mean(tL)
% %     
% %     % 555 priority: ����ֵ, ��Զ����,��Сֵ=1,���ֵΪLU����
% %     priority = 1;
% %     SIDpriority = zeros(1,size(tID,2)); % SIDpriority: ����, ÿ��SID�µ�priorityֵ
% %     szRow = cellfun(@(x)size(x,1), tID);
% %     if isscalar(szRow) %�����STRIP��ֻ��1STRIP����Ϊ������STRIP����ֵΪ1
% %         SIDpriority = 1;
% %     else
% %         [tID, ~] = padcat(tID{:});   if iscolumn(tID), tID = tID'; end
% %         
% %         % 1 �������������ڵ�Strip�ڷ�˳�� ->             
% %         % ���ú�Item2Strip���Ƶ��Ը߶����ȵ�����ʽ, ���ΪLoadingRate
% %         %         tmpM = [tLL;tL;tLW;];
% %         % %         [~,order] = sortrows(tmpM',[1,2,4],{'descend','descend','descend'}); %order is index vector ����ԭ����ֵ ���һ��Ϊ3����ʾԭarray�е�3��Ŀǰ�ǵ�1��
% %         %          tmpM = [tLL;tL;tLW;];[~,order] = sortrows(tmpM',[4,1],{'descend','descend'}); %order is index vector ����ԭ����ֵ ���һ��Ϊ3����ʾԭarray�е�3��Ŀǰ�ǵ�1��
% % 
% %         % 1 ���� > ���� > ���� > LoadingRate
% %         %         tmpM = [tLID; tMixed; tFull; tLL;tL;tLW;];  [~,order] = sortrows(tmpM',[2,3,1,5],{'ascend','descend','descend','descend'}); 
% %         % 2 ���� > ���� > ���� > LoadingRate        
% %         tmpM = [tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
% %         [~,order] = sortrows(tmpM',[1,2,3,4],{'descend','ascend','descend','descend'}); 
% %         tmpM = [tnbLU; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
% %         [~,order] = sortrows(tmpM',[1,2,3,4,5],{'descend', 'descend','ascend','descend','descend'}); 
% %         tmpM = [tnbLU; tnbLULID; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
% %         [~,order] = sortrows(tmpM',[1,2,3,4,5,6],{'descend', 'descend', 'descend','ascend','descend','descend'});         
% %         
% % %         tmpM = [tPured;tnbLU; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
% % %         [~,order] = sortrows(tmpM',[1,2,3,4,5,6],{'descend', 'descend','descend','ascend','descend','descend'}); 
% %         
% %         % 3 ȫ�� > ���� > ���� > ���� > LoadingRate
% %         %         tmpM = [tPured; tLID; tMixed; tFull; tLL;tL;tLW;];  [~,order] = sortrows(tmpM',[1,2,3,4,6],{'descend','descend','ascend','descend','descend'}); 
% %         %         tmpM = [tLID; tMixed; tFull; tLL;tL;tLW;];        [~,order] = sortrows(tmpM',[2,3],{'ascend','descend'});
% %         
% %         if ~isrow(order), order=order'; end
% % 
% %         % 2 ����1�İڷ�˳��, ��������STRIP˳��torder
% %         while any(SIDpriority==0) %����SID�����priorityδ��ȫʱ, ��������ѭ��
% %             % SIDorder�з�0��order�������Ǹ���Ϊ��ѡSTRIP
% %             % 2.1 �ҳ�������ǰ����,��һ��o����Ӧ��order(o)
% %             [~,o]=find(SIDpriority(order) == 0,1,'first');
% %             SIDpriority(order(o)) = priority;  
% %             priority=priority+1;
% %             
% %             % 2.2 �ҳ���order(o)λ��tLID��Ӧ������Strip.
% %             tnbItem = celltID{:,order(o)};
% %             % NOTE �״�tLID��Ӧ�ó��ֻ����STRIP, ���簴�߶�����, �ǿ��ܳ��ֵ�
% %             [SIDpriority,priority] = getAdjPriority(priority,order,SIDpriority,tID,tnbItem);
% %         end
% %     end
% %     allPriority(idxSID) = SIDpriority;
% % 
% % end
% % 
% % % ������䣺
% % if any(allPriority==0), error('allPriority ����δ��ֵ��'); end
% % 
% % end



%% V4 : ��EID�����ݽṹ������д
% % % �Ȱ�SID����, �ٰ�EID����, �ٸ���priorityΪSID�ڲ�˳��
% % % �ص㣺27-29��tmpM��˳��                                    
% % % function allPriority = getOrderofLID(SIDord, EIDord, Ssingle, Spured,SnbItem,SnbLU,SnbLULID, ...
% % %     SHeightfullfull, SisMixed, SLID,SLW,SLoadingRateLimit,SLoadingrate)
% % function allPriority = getOrderofLID(SIDord, EIDord, Strip)
% % 
% % %% Ŀ�� -> ��ȡallPriority, �ɶ��perPriority��Ϲ���
% % % �ֱ����ͬSID/EID�µ�STRIP����priorityֵ, ÿ�������priority���ظ�,����1��ʼ
% % % ��ͬSIDord��ȡ��ͬ��˳��, ��SID/EID�ĵ�һ����ʼ
% % allPriority = zeros(1,length(SIDord));        
% % uniSIDord = unique(SIDord); if ~issorted(uniSIDord,'strictascend'), error('���ϸ����'); end
% % for i=1:length(uniSIDord)
% %     fsid = SIDord==uniSIDord(i);
% %     uniEIDord = unique(EIDord(fsid));  if ~issorted(uniEIDord,'strictascend'), error('���ϸ����'); end %��SID�µ�EID, Ҳ�������ϸ������
% %     for j=1:length(uniEIDord)
% %         feid = EIDord==uniEIDord(j);
% %         f = fsid & feid;   if ~any(f), error('��SID��EID�����κ�STRIP'); end %����������Ը�Ϊwarning,��Ϊ���ܲ�ͬSID�µ�EID��ͬ
% %         
% %         priority = 1;  % 555 priority: ����ֵ, ��Զ����,��Сֵ=1,���ֵΪLU����
% %         perPriority = zeros(1,sum(f));     % perPriority: ����, ÿ��SID/EID�µ�f������priorityֵ,��Ҫ����ֵ
% %         
% %         
% %         %% TO DEL �����õ�ĳЩֵ ��Ϊ����STRIP����
% %         % UNUSED ��ȡ��SID��EID�ڵ�STRIP��Ӧ��ֵtLW;tLL;tL;tID;celltID.
% % %         tLW = SLW(:,f);   tLL = SLoadingRateLimit(:,f);  tPured = Spured(:,f);   tSingle = Ssingle(:,f);   %����δ֪        
% %         % FINALLY USED ���ȳ�ͷ�ڷŵ�˳��
% % %         tnbLU = SnbLU(:,f);     tnbLULID = SnbLULID(:,f); tnbItem = SnbItem(:,f);
% % %         tMixed = SisMixed(:,f);     tFull = SHeightfullfull(:,f);   tLR = SLoadingrate(:,f);        % max(tL)  %  min(tL)  % mean(tL)
% % 
% %         %% 0 ��ʼ��     tID = SLID(:,f);    celltID = tID;
% %         tID = Strip.LID(:,f);    celltID = tID;  %Strip.LIDʵ��ΪLU��ID,ƴ���ж�
% %         if length(tID) == 1
% %             mID = tID{1};
% %             isID = all(mID,2);
% %         else
% %             [mID, isID] = padcat(tID{:});  
% %         end
% %          if iscolumn(mID), mID = mID'; end;  if iscolumn(isID), isID = isID'; end
% %         
% %         % nbIDperStrip = sum(isID,1);
% %         
% %         %% 1 ��ȡperPriority: ��ͬSID/EID��, ����STRIP������˳��; �̶������ڰڷŽǶȳ���, ϸ��ÿ���ڲ�STRIP��priority
% %         % 1.1 �����STRIP��ֻ��1��STRIP����Ϊ������STRIP�����ͣ�����ֵȨ��Ϊ1.
% %         if sum(f)==1
% %             perPriority = 1;     if ~isscalar(isID), warning('�ǻ����stripҲ��ֵΪ1, Ӧ�����ⲻ��,�Ƿ�����������Ҫ�������;'); end
% %             % 1.2 �����STRIP�ĺ��ж��STRIP, �����STRIP���������ȡorder
% %         else
% %             %             tmpM = [tnbLU; tnbLULID; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
% % 
% %             tmpM = [Strip.nbLU(:,f);  Strip.nbLULID(:,f);  Strip.nbItem(:,f); Strip.isMixed(:,f);  Strip.isHeightFull(:,f);  Strip.loadingrate(:,f); ];
% %             [~,order] = sortrows(tmpM',[4,1,2,3,5,6],{'ascend','descend', 'descend', 'descend','descend','descend'});          if ~isrow(order), order=order'; end
% % %             [~,order] = sortrows(tmpM',[1,2,3,4,5,6],{'descend', 'descend', 'descend','ascend','descend','descend'});          if ~isrow(order), order=order'; end
% % %             [~,order] = sortrows(tmpM',[4],{'descend'});          if ~isrow(order), order=order'; end
% %             
% %             %% 2 ����1�İڷ�˳��, ��������STRIP��priority��torder
% %             while any(perPriority==0) %����SID�����priorityδ��ȫʱ, ��������ѭ��(δȫ��ԭ�����Ҳ������ڵ�strip��,�ٻ���order����һ����ѡstrip)
% %                 % SIDorder�з�0��order�������Ǹ���Ϊ��ѡSTRIP
% %                 % 2.1 �ҳ�������ǰ����,��һ��o����Ӧ��order(o)
% %                 [~,o]=find(perPriority(order) == 0,1,'first');
% %                 perPriority(order(o)) = priority;
% %                 priority=priority+1;
% %                 
% %                 % 2.2 �ҳ���order(o)λ��tLID��Ӧ������Strip.
% %                 tnbItem = celltID{:,order(o)};
% %                 % NOTE �״�tLID��Ӧ�ó��ֻ����STRIP, ���簴�߶�����, �ǿ��ܳ��ֵ�
% %                 if Strip.isMixed(order(o))
% %                     warning('�״β�ϣ�����ֻ����STRIP,�����ֲ������ҽ�����������STRIP'); 
% %                 end
% %                 [perPriority,priority] = getAdjPriority(priority,order,perPriority,mID,tnbItem);
% %             end
% %         end
% %         allPriority(f) = perPriority;
% %     end
% % end
% % 
% % % �Ŵ�CODE
% %  if ~all(allPriority), error('allPriority ����δ��ֵ��'); end
% % 
% % end
