%% V4 : SID is cell��ʽ �ҷֿ�ITEM��STRIP���� �ص��ų�order���ظ���ֵ�����
function SIDallPriority = getOrderofSID(SID)

    if ~iscell(SID),error('�����cell');end
    
    [mSID, isSID] = padcat(SID{:});     
    if iscolumn(mSID), mSID = mSID'; end 
    if iscolumn(isSID), isSID = isSID'; end % ����padcat: ��cell�ڵľ���ת��Ϊ�������
    nbSIDperStrip = sum(isSID,1);     %nbSIDperStrip: ����,ÿ��strip��SID�ĸ���     nbSIDperStrip = cellfun(@(x)size(x,1), SID);

    % CASE 1: ����STRIP��Ϊ��һ�� �����1��Strip�ĵ�һ��(isSID��һάscalar=1)������(isSID��һά������)
    if isvector(isSID) || isscalar(isSID)
        % mSID = cell2mat(SID);
        SIDallPriority = mSID;  return;

    % CASE 2: ���ڰ������SID��STRIP
    else %��STRIPΪ�����-���������������ϵĻ��
        if max(nbSIDperStrip)>=3,  error('ͬһSTRIP��ϴ�������������SID');    end
        if max(nbSIDperStrip)==1,  error('STRIP���ǻ����,Ӧ���Ѿ����ų���');    end
        warning('ͬһSTRIP����ʹ�������SID');
        
        uniSID = unique(mSID(~isnan(mSID)));          if ~issorted(uniSID,'strictascend'), error('��BIN��STRIP��SID��δ��ָ��SID(�ϸ�)��������.'); end
        priority=1;  
        SIDallPriority = zeros(1,length(nbSIDperStrip));
        
        for idxSID=1:length(uniSID)
         
            fsid = any(ismember(mSID,uniSID(idxSID)),1);     % ������sid��trip�߼��ж�
            fprio = SIDallPriority==0;                                       % δ��Ȩ��ֵ�߼��ж�
            fmix = nbSIDperStrip>1;                                         % ������߼��ж�
            
            f1 = fsid & fprio & ~fmix;       % �ǻ��
            f2 = fsid & fprio & fmix;          % ���
            
            % ��һ��strip��idx�ǿ� �� ��priorityδȫ�������� �Ĵ���
            if any(f1)
                SIDallPriority(f1) = priority;
                priority=priority+1;
            end
            
            % �����strip��idx�ǿ� �� ��priorityδȫ�������� �Ĵ���
            if any(f2)
                if sum(f2)>1, error('ͬһSID,������ڶ�������strip! '); end
                if ~issorted(mSID(:,f2),'strictascend'), error('ͬһSID,����strip��SID�Ƿǵ�������ֵ(��SID����)'); end % TODO , �ų�mSID(:,f2)����NaN�����
                SIDallPriority(f2) = priority;
                priority=priority+1;
                
                % �Ŵ�CODE
                nextSID = mSID(:,f2);  nextSID(nextSID(:)==uniSID(idxSID))=[];   if ~isscalar(nextSID) || isempty(nextSID),   error('ͬһSTRIP����ʹ�������������SID��Ϊ�գ�');    end
                if uniSID(idxSID+1) ~= nextSID, error('������ڶ��SID,��һ��SID�����Ǻͱ�SID��ϵĹ�Ӧ��,�����Ƿ�����������'); end
            end
            
            % �Ŵ�CODE
            if ~any(f1) && ~any(f2), warning('�����ڰ�����SID��strip��,�����Ǹ�SID��������,�����ڻ��case�а���'); end
        end
        
        % �Ŵ�CODE
        if ~all(SIDallPriority),  error('����δȨ�ص�SIDֵ'); end
        
    end
end

            

%         while ~done
%             if all(SIDallPriority),  done=false;   continue;  end   %���orderȫ������0,��ʾ�Ѿ�ȫ������
   
            % V2: idx1(��stripΪ��һ��)�����ж�� idx2(��stripΪ�����) ֻ����һ��
%             [~,nbcol1] = find(mID(:,:)==uniSID(SID) & sum(~isnan( mID(:,:) ),1) ==1 & allPriority==0 ); %�ҳ�mID�а���SID����,�Ҹ��н���1��
%             [~,nbcol2] = find( mID(:,:)==uniSID(SID) & sum(~isnan( mID(:,:) ),1) ~=1 & allPriority==0); %�ҳ�mID�а���SID����,�Ҹ��в�ֻ1��
%             
            % V1: idx1(stripΪ��һ��)�����ж�� idx2(stripΪ�����) ֻ����һ��
%             [~,nbcol1] = find(mSID(:,:)==idxSID & sum(~isnan( mSID(:,:) ),1) ==1 & SIDallPriority==0 ); %�ҳ�mID�а���SID����,�Ҹ��н���1��
%             [~,nbcol2] = find( mSID(:,:)==idxSID & sum(~isnan( mSID(:,:) ),1) ~=1 & SIDallPriority==0); %�ҳ�mID�а���SID����,�Ҹ��в�ֻ1��
             
%             % ��һ��strip��idx1�ǿ� �� ��orderδȫ��������
%             if ~isempty(nbcol1)
%                 SIDallPriority(nbcol1) = priority;
%                 priority=priority+1;                
%                 % ������һ��SID: �ǻ���Ͱ��������� 
%                 nextSID = idxSID + 1; %���SID�ж�,SID�м䲻�ܷ�����
%             end
            
%             % �����strip��idx2�ǿ� �� ��orderδȫ��������
%             if ~isempty(nbcol2)
%                 if any(diff(mSID(:,nbcol2))~=1)  %������ڻ����SID���һ�ϵķ�����ֵ, ֻ���Ǵ���
%                     error('�����strip�Ĵ��ڻ�Ϸ�����ֵ��');
%                 end                               
%                 if numel(nbcol2)>1
%                     error('�����strip��δ��ֵ�Ĵ��ڶ����');
%                 end
%                 SIDallPriority(nbcol2) = priority;
%                 priority=priority+1;
            
%                 % ������һ��SID: �ҳ�idx2���ų���SID�����һ��SID ����Ϊ���
%                 hyLID = mSID(:, nbcol2);
%                 hyLID(hyLID(:)==idxSID)=[];
%                 if ~isscalar(hyLID),   error('ͬһSTRIP����ʹ�������������SID��');    end
%                 nextSID = hyLID;

            
%             % �����ʩ - ��Ӧ����
%             if isempty(nbcol1) && isempty(nbcol2)
%                 warning('�����ڰ�����SID��strip��,�����Ǹ�SID��������,���ڻ��case�а���');
%                 if ~exist('nextSID'),
%                     error('������nextSID');
%                     nextSID=idxSID; end;
%                 nextSID=nextSID+1;
%             end
            
%             idxSID = nextSID;
            



%% 
% ��ʼ��
% % mID = sort(mID,1); %�ȶ�ID��������->��С�ķ�ǰ��
% % lastTID = -1;
% % nbID = max(mID(:));
% % typeID = zeros(1,nbID);
% % tmpidx=zeros(1,numel(szRow));
% % while 1
% %     if all(typeID), break; end %���aIDȫ��=1,��ʾ�Ѿ�ȫ������
% %     idx3 = find(sum(~isnan( mID(:,:) ),1) >=3); %�ҳ���t����t����idx �ҷ�NaNֵ�ж���1������θ�Ȩ
% %     if ~isempty(idx3), warning('��������ͬ��һ��strip,��Ҫ��ע');end
% %     
% % 
% %     1
% %     if ~isempty(idx1)
% %         priority=priority+1;
% %         order(idx1) = priority;
% %         typeID(tID) = 1; %tID�Ѿ�ʹ��
% %     end
% %     if ~isempty(idx2)
% %         priority=priority+1;
% %         order(idx2) = priority;
% %         typeID(tID) = 1; %tID�Ѿ�ʹ��
% %     end    
% %     if isempty(idx1) && isempty(idx2)
% %         tmptID = find(typeID(:) == 0); %��һ��������������ŷ���
% %         if ~isempty(tmptID)
% %             for i=1:numel(tmptID)
% %               tID = find(typeID(:) == 0,i); 
% %                if tID ~= lastTID,
% %                    lastTID = tID;
% %                    continue;
% %                end
% %             end
% %         end
% %     end
% %     if ~isempty(idx2)
% %         tID =  mID(2,idx2); %ʹ�ñ�strip��Ӧ�Ķ����1��
% %     elseif ~isempty(idx1)
% %         tID = find(typeID(:) == 0,1); %��һ��������������ŷ���
% %     end
% %     
% % end
% % 
% % end


%% %% V3 : ����cell��ʽ �ҷֿ�ITEM��STRIP���� �ص��ų�order���ظ���ֵ�����
% % % SID: ����SID��С�����˳������
% % function allPriority = getOrderofSID(ID)
% % % ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';
% % 
% % if ~iscell(ID),error('�����cell');end
% % szRow = cellfun(@(x)size(x,1), ID);
% % 
% %     if (max(szRow)==min(szRow) && ~isscalar(szRow)) %����STRIP��Ϊ��һ��
% %         mID = cell2mat(ID);
% %         allPriority = mID;
% %         return;
% %     else %��STRIPΪ�����-���������������ϵĻ��
% %         if max(szRow)>=3
% %             error('ͬһSTRIP����ʹ�������������SID');
% %         end
% %         warning('ͬһSTRIP����ʹ�������SID');
% %         
% %         
% %         [mID, ~] = padcat(ID{:}); % �������padcat:��cell�ڵľ���ת��Ϊ�������        
% %         uniSID = unique(mID(~isnan(mID)));
% %         allPriority = zeros(1,size(mID,2));
% %         priority=1;        
% %         SID = 1;  %��ʵSID,��1��ʼ        
% %         while 1
% %             if all(allPriority)
% %                  %mID
% %                  %allPriority
% %                 break;
% %             end %���orderȫ������0,��ʾ�Ѿ�ȫ������
% %             
% %             % V2: idx1(��stripΪ��һ��)�����ж�� idx2(��stripΪ�����) ֻ����һ��
% % %             [~,nbcol1] = find(mID(:,:)==uniSID(SID) & sum(~isnan( mID(:,:) ),1) ==1 & allPriority==0 ); %�ҳ�mID�а���SID����,�Ҹ��н���1��
% % %             [~,nbcol2] = find( mID(:,:)==uniSID(SID) & sum(~isnan( mID(:,:) ),1) ~=1 & allPriority==0); %�ҳ�mID�а���SID����,�Ҹ��в�ֻ1��
% % %             
% %             % V1: idx1(��stripΪ��һ��)�����ж�� idx2(��stripΪ�����) ֻ����һ��
% %             [~,nbcol1] = find(mID(:,:)==SID & sum(~isnan( mID(:,:) ),1) ==1 & allPriority==0 ); %�ҳ�mID�а���SID����,�Ҹ��н���1��
% %             [~,nbcol2] = find( mID(:,:)==SID & sum(~isnan( mID(:,:) ),1) ~=1 & allPriority==0); %�ҳ�mID�а���SID����,�Ҹ��в�ֻ1��
% %             
% %             % ��һ��strip��idx1�ǿ� �� ��orderδȫ��������
% %             if ~isempty(nbcol1)
% %                 allPriority(nbcol1) = priority;
% %                 priority=priority+1;                
% %                 % ������һ��SID: �ǻ���Ͱ��������� 
% %                 nextSID = SID + 1; %���SID�ж�,SID�м䲻�ܷ�����
% %             end
% %             
% %             % �����strip��idx2�ǿ� �� ��orderδȫ��������
% %             if ~isempty(nbcol2)
% %                 if any(diff(mID(:,nbcol2))~=1)  %������ڻ����SID���һ�ϵķ�����ֵ, ֻ���Ǵ���
% %                     error('�����strip�Ĵ��ڻ�Ϸ�����ֵ��');
% %                 end                               
% %                 if numel(nbcol2)>1
% %                     error('�����strip��δ��ֵ�Ĵ��ڶ����');
% %                 end
% %                 allPriority(nbcol2) = priority;
% %                 priority=priority+1;
% %             
% %                 % ������һ��SID: �ҳ�idx2���ų���SID�����һ��SID ����Ϊ���
% %                 hyLID = mID(:, nbcol2);
% %                 hyLID(hyLID(:)==SID)=[];
% %                 if ~isscalar(hyLID),   error('ͬһSTRIP����ʹ�������������SID��');    end
% %                 nextSID = hyLID;
% %             end
% %             
% %             % �����ʩ - ��Ӧ����
% %             if isempty(nbcol1) && isempty(nbcol2)
% %                 warning('�����ڰ�����SID��strip��,�����Ǹ�SID��������,���ڻ��case�а���');
% %                 if ~exist('nextSID'),
% %                     error('������nextSID');
% %                     nextSID=SID; end;
% %                 nextSID=nextSID+1;
% %             end
% %             
% %             SID = nextSID;
% %             
% %         end
% %     end
% % end

%% 
% ��ʼ��
% % mID = sort(mID,1); %�ȶ�ID��������->��С�ķ�ǰ��
% % lastTID = -1;
% % nbID = max(mID(:));
% % typeID = zeros(1,nbID);
% % tmpidx=zeros(1,numel(szRow));
% % while 1
% %     if all(typeID), break; end %���aIDȫ��=1,��ʾ�Ѿ�ȫ������
% %     idx3 = find(sum(~isnan( mID(:,:) ),1) >=3); %�ҳ���t����t����idx �ҷ�NaNֵ�ж���1������θ�Ȩ
% %     if ~isempty(idx3), warning('��������ͬ��һ��strip,��Ҫ��ע');end
% %     
% % 
% %     1
% %     if ~isempty(idx1)
% %         priority=priority+1;
% %         order(idx1) = priority;
% %         typeID(tID) = 1; %tID�Ѿ�ʹ��
% %     end
% %     if ~isempty(idx2)
% %         priority=priority+1;
% %         order(idx2) = priority;
% %         typeID(tID) = 1; %tID�Ѿ�ʹ��
% %     end    
% %     if isempty(idx1) && isempty(idx2)
% %         tmptID = find(typeID(:) == 0); %��һ��������������ŷ���
% %         if ~isempty(tmptID)
% %             for i=1:numel(tmptID)
% %               tID = find(typeID(:) == 0,i); 
% %                if tID ~= lastTID,
% %                    lastTID = tID;
% %                    continue;
% %                end
% %             end
% %         end
% %     end
% %     if ~isempty(idx2)
% %         tID =  mID(2,idx2); %ʹ�ñ�strip��Ӧ�Ķ����1��
% %     elseif ~isempty(idx1)
% %         tID = find(typeID(:) == 0,1); %��һ��������������ŷ���
% %     end
% %     
% % end
% % 
% % end


%% V2 : ����cell��ʽ  ��ȡSTRIP�ڲ�����ID�������˳�� BUG: δ���ǵ�strip��LUID,�����Ǹ�LUID�ڵ�strip��Ϊ
% % function order = getOrderofID(ID)
% % % ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';
% % 
% % if ~iscell(ID),error('�����cell');end
% % sz = cellfun(@(x)size(x,1), ID);
% % 
% % % ��cellΪ�ȳ�,��ֱ��ת��;�������padcat
% % if numel(unique(sz)) == 1
% %     ID = cell2mat(ID);
% % else
% %     % ��cell�ڵľ���ת��Ϊ������� padcat
% %     [ID, ~] = padcat(ID{:});
% % end
% % 
% % % ��ʼ��
% % ID = sort(ID,1); %�ȶ�ID��������->��С�ķ�ǰ��
% % priority=0;
% % order = zeros(1,size(ID,2));
% % 
% % for t=1:max(ID(:)) %��PID/SID�ĵ�1��ֱ�����1��,���Ѱ�Ҳ�����˳��
% %     idx3 = find(sum(~isnan( ID(:,:) ),1) >=3); %�ҳ���t����t����idx �ҷ�NaNֵ�ж���1������θ�Ȩ
% %     if ~isempty(idx3), warning('��������ͬ��һ��strip,��Ҫ��ע');end
% %     
% %     idx1 = find(ID(1,:)==t & sum(~isnan( ID(:,:) ),1) ==1); %�ҳ���t����t����idx �ҷ�NaNֵ����1���������ȸ�Ȩ
% %     idx2 = find(ID(1,:)==t & sum(~isnan( ID(:,:) ),1) ~=1); %�ҳ���t����t����idx �ҷ�NaNֵ�ж���1������θ�Ȩ
% % 
% %     if ~isempty(idx1)
% %         priority=priority+1;
% %         order(idx1) = priority;
% %     end
% %     if ~isempty(idx2)
% %         priority=priority+1;
% %         order(idx2) = priority;
% %     end    
% % end
% % 
% % end


%% V1 : ֱ�Ӳ��þ���ʽ
% % % function z = getOrderofID(ID)
% % % t = ID;
% % %  ss = sum(t,1);  %ÿ��STRIP�ڰ�����SID����
% % % 
% % %  z = zeros(1,size(t,2)); 
% % %  k=1;
% % %  for i=1:numel(ss)
% % %     if ss(i) ==1
% % %         if i>1 && ss(i-1) ==1 && find(t(:, i)==1) ~= find(t(:, i-1)==1) %�жϵ�ǰ��ǰһ��strip�Ƿ�����ͬ��SID
% % %              k=k+1;
% % %         end     
% % %         z(i) = k;
% % % %     elseif  i>1 && ss(i) >1 %ֻҪ����STRIP����2�������ϵ�STRIPʱ������˳��                
% % % %         k=k+1;
% % % %         z(i)=k;
% % % %         k=k+1;
% % %     elseif ss(i) >1 %ֻҪ����STRIP����2�������ϵ�STRIPʱ������˳��      
% % %         k=k+1;
% % %         z(i)=k;
% % %         k=k+1;
% % %  end
% % % end


%  z = zeros(1,size(t,2)); 
%  k=1;
%  for i=1:numel(ss)
%     if ss(i) ==1
%         if i>1 && ss(i-1) ==1 && find(t(:, i)==1) ~= find(t(:, i-1)==1) %�жϵ�ǰ��ǰһ��strip�Ƿ�����ͬ��SID
%             k=k+1;
%         end
%         z(i) = k;
%     elseif ss(i) >1 %ֻҪ����STRIP����2�������ϵ�STRIPʱ������˳��        
%         k=k+1;
%         z(i)=k;
%         k=k+1;
%     end
%  end
% end