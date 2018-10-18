%% V3 : ����cell��ʽ �ҷֿ�ITEM��STRIP���� �ص��ų�order���ظ���ֵ�����
% �Ȱ�SID����, �ٸ���priorityΪSID�ڲ�˳��
% �ص㣺27-29��tmpM��˳��

% IDorder = getOrderofLID(SIDorder, Strip.isSingleItem, Strip.isAllPured, Strip.nbItem, Strip.isHeightFull,...
%                                         Strip.isMixed, Strip.LID, Strip.LW(1:2,:), Strip.loadingrateLimit, Strip.loadingrate);
                                    
function allPriority = getOrderofLID(SIDord,Ssingle, Spured,SnbItem,SnbLU,SnbLULID, ...
    SHeightfullfull, SisMixed, SLID,SLW,SLoadingRateLimit,SLoadingrate)
% ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';
if ~iscell(SLID),error('�����cell');end

allPriority = zeros(1,size(SLID,2));
uniOrd = unique(SIDord);

% ��ͬSID��ȡ��ͬ��˳��, ��SID�����1��ʼ
for i=1:length(uniOrd)    
    idxSID = SIDord==uniOrd(i);
    
    % UNUSED ��ȡ��SID�ڵ�STRIP��Ӧ��ֵtLW;tLL;tL;tID;celltID.
    tLW = SLW(:,idxSID);   tLL = SLoadingRateLimit(:,idxSID);
    tID = SLID(:,idxSID);    celltID = tID;
    tPured = Spured(:,idxSID);
    tSingle = Ssingle(:,idxSID);   %����δ֪
    
    % FINALLY USED ���ȳ�ͷ�ڷŵ�˳��
    tnbItem = SnbItem(:,idxSID);
    tnbLU = SnbLU(:,idxSID);
    tnbLULID = SnbLULID(:,idxSID);
    
    
    tMixed = SisMixed(:,idxSID);
    tFull = SHeightfullfull(:,idxSID);
    tLR = SLoadingrate(:,idxSID);        % max(tL)  %  min(tL)  % mean(tL)
    
    % 555 priority: ����ֵ, ��Զ����,��Сֵ=1,���ֵΪLU����
    priority = 1;
    SIDpriority = zeros(1,size(tID,2)); % SIDpriority: ����, ÿ��SID�µ�priorityֵ
    szRow = cellfun(@(x)size(x,1), tID);
    if isscalar(szRow) %�����STRIP��ֻ��1STRIP����Ϊ������STRIP����ֵΪ1
        SIDpriority = 1;
    else
        [tID, ~] = padcat(tID{:});   if iscolumn(tID), tID = tID'; end
        
        % 1 �������������ڵ�Strip�ڷ�˳�� ->             
        % ���ú�Item2Strip���Ƶ��Ը߶����ȵ�����ʽ, ���ΪLoadingRate
        %         tmpM = [tLL;tL;tLW;];
        % %         [~,order] = sortrows(tmpM',[1,2,4],{'descend','descend','descend'}); %order is index vector ����ԭ����ֵ ���һ��Ϊ3����ʾԭarray�е�3��Ŀǰ�ǵ�1��
        %          tmpM = [tLL;tL;tLW;];[~,order] = sortrows(tmpM',[4,1],{'descend','descend'}); %order is index vector ����ԭ����ֵ ���һ��Ϊ3����ʾԭarray�е�3��Ŀǰ�ǵ�1��

        % 1 ���� > ���� > ���� > LoadingRate
        %         tmpM = [tLID; tMixed; tFull; tLL;tL;tLW;];  [~,order] = sortrows(tmpM',[2,3,1,5],{'ascend','descend','descend','descend'}); 
        % 2 ���� > ���� > ���� > LoadingRate        
        tmpM = [tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
        [~,order] = sortrows(tmpM',[1,2,3,4],{'descend','ascend','descend','descend'}); 
        tmpM = [tnbLU; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
        [~,order] = sortrows(tmpM',[1,2,3,4,5],{'descend', 'descend','ascend','descend','descend'}); 
        tmpM = [tnbLU; tnbLULID; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
        [~,order] = sortrows(tmpM',[1,2,3,4,5,6],{'descend', 'descend', 'descend','ascend','descend','descend'});         
        
%         tmpM = [tPured;tnbLU; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
%         [~,order] = sortrows(tmpM',[1,2,3,4,5,6],{'descend', 'descend','descend','ascend','descend','descend'}); 
        
        % 3 ȫ�� > ���� > ���� > ���� > LoadingRate
        %         tmpM = [tPured; tLID; tMixed; tFull; tLL;tL;tLW;];  [~,order] = sortrows(tmpM',[1,2,3,4,6],{'descend','descend','ascend','descend','descend'}); 
        %         tmpM = [tLID; tMixed; tFull; tLL;tL;tLW;];        [~,order] = sortrows(tmpM',[2,3],{'ascend','descend'});
        
        if ~isrow(order), order=order'; end

        % 2 ����1�İڷ�˳��, ��������STRIP˳��torder
        while any(SIDpriority==0) %����SID�����priorityδ��ȫʱ, ��������ѭ��
            % SIDorder�з�0��order�������Ǹ���Ϊ��ѡSTRIP
            % 2.1 �ҳ�������ǰ����,��һ��o����Ӧ��order(o)
            [~,o]=find(SIDpriority(order) == 0,1,'first');
            SIDpriority(order(o)) = priority;  
            priority=priority+1;
            
            % 2.2 �ҳ���order(o)λ��tLID��Ӧ������Strip.
            tnbItem = celltID{:,order(o)};
            % NOTE �״�tLID��Ӧ�ó��ֻ����STRIP, ���簴�߶�����, �ǿ��ܳ��ֵ�
            [SIDpriority,priority] = getAdjPriority(priority,order,SIDpriority,tID,tnbItem);
        end
    end
    allPriority(idxSID) = SIDpriority;

end

% ������䣺
if any(allPriority==0), error('allPriority ����δ��ֵ��'); end

end

%% �ֲ����� %%
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
                [a,~] = sortrows(tmp,[2],{'ascend'});
                for i1=1:length(nbcol1)
                    SIDorder(a(i1,1)) = priority;
                    priority=priority+1;
                end
            end
            
            % Ѱ��tID�е���: ������2������LID; ���ж�Ӧ��SIDorderΪ0�����а����������ڣ�tLID
            [~,nbcol2] = find(tID(:,:)==tLID & sum(~isnan( tID(:,:) ),1) ~=1 &  SIDorder==0);
            if ~isempty(nbcol2)
                if ~isscalar(nbcol2),  error('nbcol2������ֻһ����: ���ڶ����ͬ��ϵ�STRIP');   end
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