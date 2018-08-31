%% V3 : ����cell��ʽ �ҷֿ�ITEM��STRIP���� �ص��ų�order���ظ���ֵ�����
% �Ȱ�SID����, �ٸ���priorityΪSID�ڲ�˳��
% �ص㣺27-29��tmpM��˳��
function allPriority = getOrderofLID(SIDord,ID,LW,LoadingRateLimit,Loadingrate)
% ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';

if ~iscell(ID),error('�����cell');end

allPriority = zeros(1,size(ID,2));
uniOrd = unique(SIDord);
for i=1:length(uniOrd)
    idxSID = SIDord==uniOrd(i);
    
    tLW = LW(:,idxSID);   tLL = LoadingRateLimit(:,idxSID);     tL = Loadingrate(:,idxSID);
    
    tID = ID(:,idxSID);   celltID = tID;
    
    priority = 1;
    SIDorder = zeros(1,size(tID,2));
    szRow = cellfun(@(x)size(x,1), tID);
    if isscalar(szRow)
        SIDorder = 1;
    else
        [tID, ~] = padcat(tID{:});   if iscolumn(tID), tID = tID'; end
        
        % 1 �������������ڵ�Strip�ڷ�˳��
        tmpM = [tLL;tL;tLW;];
            %         [~,order] = sortrows(tmpM',[1,2,4],{'descend','descend','descend'}); %order is index vector ����ԭ����ֵ ���һ��Ϊ3����ʾԭarray�е�3��Ŀǰ�ǵ�1��
        % ���ú�Item2Strip���Ƶ��Ը߶����ȵ�����ʽ, ���δLoadingRate
        [~,order] = sortrows(tmpM',[4,1],{'descend','descend'}); %order is index vector ����ԭ����ֵ ���һ��Ϊ3����ʾԭarray�е�3��Ŀǰ�ǵ�1��
        if ~isrow(order), order=order'; end
        
        % 2 ����1�İڷ�˳��, ��������STRIP˳��torder
        while any(SIDorder==0)
            % SIDorder�з�0��order�������Ǹ���Ϊ��ѡSTRIP
            % 2.1 �ҳ�������ǰ����,��һ��o����Ӧ��order(o)
            [~,o]=find(SIDorder(order) == 0,1,'first');
            SIDorder(order(o)) = priority;  priority=priority+1;
            
            % 2.2 �ҳ���order(o)λ��tLID��Ӧ������Strip.
            tLID = celltID{:,order(o)};
            % �״�tLID��Ӧ�ó��ֻ����STRIP
            [SIDorder,priority] = getAdjPriority(priority,order,SIDorder,tID,tLID);
        end
    end
    allPriority(idxSID) = SIDorder;
end
end

% Ѱ��tLID��Ӧ������STRIP,����һ���߼�
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
                [a,~] = sortrows(tmp,[2],{'ascend'});
                for i1=1:length(nbcol1)
                    SIDorder(a(i1,1)) = priority;
                    priority=priority+1;
                end
            end
            
            % Ѱ��tID�е���: ������2������LID; ���ж�Ӧ��SIDorderΪ0�����а����������ڣ�tLID
            [~,nbcol2] = find(tID(:,:)==tLID & sum(~isnan( tID(:,:) ),1) ~=1 &  SIDorder==0);
            if ~isempty(nbcol2)
                if ~isscalar(nbcol2),  error('nbcol2������ֻһ����');   end
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
            error('eeee');
        end
    end
end


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