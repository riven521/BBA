%% ����: Strip˦β v3 ����˦β��Χ
function   [tempStrip_Bin, StripisShuaiWei,LUisShuaiWei,TF] = HStripSW(Strip,LU)
% 1 ��Щ˦β: ��Ȳ���isWidthFull��߶Ȳ���isHeightFull�� ��Ϊȫ��˦β��������Χ
% 2 ��Щ�ǿ�Ȳ���: 
% 3 ��Щ�Ǹ߶Ȳ���:     

    tempStrip_Bin = Strip.Strip_Bin(2,:);
    TF = false;    
    %     TF = true;    
    LUisShuaiWei = zeros(size(LU.Weight));      % �ж�LU�Ƿ�˦β������,������ͼʱ��������
    StripisShuaiWei = zeros(size(Strip.Weight));       % StripseqShuaiWei = zeros(size(Strip.Weight)); % seqShuaiWeiԽ��,����Խ��˦β; ԽС, Խ��˦β, ��Խ�����ڳ�ͷ.  
       
% % %% 1: ********************** ˦β ********************************** 
% % % 1 ��Щ˦β: ��Ȳ���isWidthFull��߶Ȳ���isHeightFull��
% % if any(~Strip.isWidthFull | ~Strip.isHeightFull) % | ~Strip.isHeightBalance)
% %     
% %     fprintf(1,'       Exsiting ˦β���� in HStripSW (Strip��߲���)...\n');
% %     TF = true;
% %     
% %     % Get b : strip index to be move to end of Vehicle
% %     [~,bNOTheightfull] = find(Strip.isHeightFull == 0);
% %     [~,bNOTwidthfull] = find(Strip.isWidthFull == 0);
% %     [~,bNOTheightbalance] = find(Strip.isHeightBalance == 0);
% %     
% %     b = unique([bNOTheightfull, bNOTwidthfull],'stable');    % ���ڷų�β��˳����ȫ��order
% % % b = unique([bNOTheightfull, bNOTwidthfull, bNOTheightbalance],'stable');   
% % 
% % % 2 ���������˦β��Strip, Ҫ�������? ��order
% %    if ~isempty(b)       
% %        %��1 Sort b by LoadingRate,Height etc.
% %             %tmpM = [Strip.loadingrate(b); Strip.maxHeight(b);Strip.loadingrateLimit(b);];
% %             %[~,order] = sortrows(tmpM',[1,2,3],{'descend','descend','descend'});
% %        %��2 Sort b by 1 ʣ���ȴ�ĺ��; 2 ʣ��߶ȶ�ĺ����
% %             %tmpM = [Strip.LW(1,b); Strip.maxHeight(b)];
% %             %[~,order] = sortrows(tmpM',[1,2],{'ascend','ascend'});
% % 
% %         %��3 Sort b 1 �ǵ����ķ�����; �ǻ��isMixed�ķ�����; �߶Ⱦ���isHeightBalance�ķ�����; ��͸߶�lowestHeight�ݼ��ķ�����
% %              %tmpM = [Strip.isSingleItem(b);Strip.isMixed(b);Strip.isHeightBalance(b); Strip.lowestHeight(b); Strip.maxHeight(b);Strip.meanHeight(b)];
% %             %[~,order] = sortrows(tmpM',[1,2,3,4,5],{'ascend','ascend','descend','descend','descend'});      
% %        
% %        %��4 Sort b 1 �ǻ��isMixed�ķ�����;  ƽ���߶�t�ݼ��ķ�����
% %        tmpM = [Strip.isMixed(b);Strip.meanHeight(b)];
% %        [~,order] = sortrows(tmpM',[1,2],{'ascend','descend'});   
% %        
% %        b = b(order);
% %        
% %        StripisShuaiWei(b) = 1;        % StripseqShuaiWei(b) = order;
% % 
% %        %%% *********** LU˦β���ҳ���,��ͼ�� ***************
% %        LUisShuaiWei( ismember(LU.LU_Strip(1,:), b) ) = 1;  % ��strip�ڵ�ֵ
% %    end   
% %   
% %        
% %     %%%  ***************** �Ƿ�˦β�Ŀ��� *************
% %     for i=1:length(b)
% %         %     Strip = repairStripPlace(Strip,b(i));    %V1 ������SID/EID��˦β Strip.Strip_Bin
% %         %     Strip = repairStripPlace2(Strip,b(i));    %V2 ����SID/EDI��˦β Strip.Strip_Bin
% %         tempStrip_Bin = repairStripPlace2(Strip,b(i));    %V3 ����Strip.Strip_Bin
% %         Strip.Strip_Bin(2,:) = tempStrip_Bin;  % ������ Strip.Strip_Bin(2,:) ��Ҫѭ����ȥ
% %     end
% %        
% % end

% % %% 1: ********************** ˦β ********************************** 
% % % 1 ��Щ˦β: ���ж�˦β��������ͷ�󣬶Գ��ڽ�������ڷš�
       %��5 Sort b 1 �ǻ��isMixed�ķ�����;  ƽ���߶�t�ݼ��ķ�����
       abnStrip = zeros(1,length(Strip.Weight));
       idx = ~Strip.isWidthFull | ~Strip.isHeightFull  | Strip.isMixed;  %todo isMixed�Ƿ�˦β�أ�
       abnStrip(idx) = 1;
       
       %        [Strip.nbItem, Strip.nbLU, Strip.nbLULID]
       
       nLU = Strip.nLUIDBin;
       nLULID = Strip.nLULIDBin;
       
       nLU(idx) = 0;
       nLULID(idx) = 0;              

       tmpM = [abnStrip; Strip.isMixed; nLU; nLULID; Strip.meanHeight];
       [~,order] = sortrows(tmpM',[1,2,3,4,5],{'ascend','ascend','descend','descend','descend'});   
    
       b = 1:length(Strip.Weight);
       b = b(order);
       
       StripisShuaiWei(b) = 1;        % StripseqShuaiWei(b) = order;

       %%% *********** LU˦β���ҳ���,��ͼ�� ***************
       LUisShuaiWei( ismember(LU.LU_Strip(1,:), b(idx) ) ) = 1;  % ��strip�ڵ�ֵ    
        %%%  ***************** �Ƿ�˦β�Ŀ��� *************
    for i=1:length(b)
        %     Strip = repairStripPlace(Strip,b(i));    %V1 ������SID/EID��˦β Strip.Strip_Bin
        %     Strip = repairStripPlace2(Strip,b(i));    %V2 ����SID/EDI��˦β Strip.Strip_Bin
        tempStrip_Bin = repairStripPlace2(Strip,b(i));    %V3 ����Strip.Strip_Bin
        Strip.Strip_Bin(2,:) = tempStrip_Bin;  % ������ Strip.Strip_Bin(2,:) ��Ҫѭ����ȥ
    end
    
end

%% �ֲ�����
%% V3 ���ٷ��ز��� 555 �ص����Ѻ���
function Strip_Bin_Seq = repairStripPlace2(Strip,stripidx)
    Strip_Bin_Seq = Strip.Strip_Bin(2,:);
    
    % 1 �ж����stripidx��SID�Ļ�ϻ�EID�Ļ��,����˦β,ֱ��return
    if ~isscalar(Strip.SID{stripidx}) || ~isscalar(Strip.EID{stripidx})
        return; end    
    
    % 2 �ж����stripidx����SID/EID��ϵ�,
    % ����Ƿǻ�ϵ�SID/EID,����˦β,������˦����β,��ʵ˦����stripidx��ͬ��SID/EIDβ��(��˦SID��,����˦EID��)
    sid = Strip.SID{stripidx};
    eid = Strip.EID{stripidx};
    bid = Strip.Strip_Bin(1,stripidx); 
    
    [tSID, ~] = padcat(Strip.SID{:});    if iscolumn(tSID), tSID = tSID'; end
    [tEID, ~] = padcat(Strip.EID{:});    if iscolumn(tEID), tEID = tEID'; end
    
    % 3 �ҵ�stripidx��Ӧ��BIN�µ�����Strip�������߼�ֵ
    fsid = Strip.isMixedSID == 0  &  any(ismember(tSID, sid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
    feid = Strip.isMixedEID == 0  &  any(ismember(tEID, eid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
    fbid = Strip.Strip_Bin(1,:) == bid;
    fbig =  Strip.Strip_Bin(2,:) > Strip.Strip_Bin(2,stripidx);       % ���аڷ�˳������stripidx���߼��ж�    
    f = fbid&fsid&feid&fbig;
    
    % 4 ����˳��Strip.Strip_Bin
    maxstripidx = max(Strip.Strip_Bin(2,f));        % ������(���ͻ���ǰ)��strip����˳������ֵ,���豾stripidx
    if ~isempty(maxstripidx)
        % 2.1 �������ڱ�Bin�� & �Ұڷ�˳������stripidx ��˳���1, ����ǰ�ڷ�
        Strip_Bin_Seq(f)  = Strip.Strip_Bin(2,f)  - 1;
        % 2.2 ��ǰstripidx�ڷŵ�(��ͬSID/EIDβ��), ��˳�����õ����
        Strip_Bin_Seq(stripidx) = maxstripidx;
    end
    
end

%% V2 ����SID��EID��˦β
% % function Strip = repairStripPlace2(Strip,stripidx)
% %     % 1 �ж����stripidx��SID�Ļ�ϻ�EID�Ļ��,����˦β,ֱ��return
% %     if ~isscalar(Strip.SID{stripidx}) || ~isscalar(Strip.EID{stripidx})
% %         return; end    
% %     
% %     % 2 �ж����stripidx����SID/EID��ϵ�,
% %     % ����Ƿǻ�ϵ�SID/EID,����˦β,������˦����β,��ʵ˦����stripidx��ͬ��SID/EIDβ��(��˦SID��,����˦EID��)
% %     sid = Strip.SID{stripidx};
% %     eid = Strip.EID{stripidx};
% %     bid = Strip.Strip_Bin(1,stripidx); 
% %     
% %     [tSID, ~] = padcat(Strip.SID{:});    if iscolumn(tSID), tSID = tSID'; end
% %     [tEID, ~] = padcat(Strip.EID{:});    if iscolumn(tEID), tEID = tEID'; end
% %     
% %     % 3 �ҵ�stripidx��Ӧ��BIN�µ�����Strip�������߼�ֵ
% %     fsid = Strip.isMixedSID == 0  &  any(ismember(tSID, sid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
% %     feid = Strip.isMixedEID == 0  &  any(ismember(tEID, eid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
% %     fbid = Strip.Strip_Bin(1,:) == bid;
% %     fbig =  Strip.Strip_Bin(2,:) > Strip.Strip_Bin(2,stripidx);       % ���аڷ�˳������stripidx���߼��ж�    
% %     f = fbid&fsid&feid&fbig;
% %     
% %     % 4 ����˳��Strip.Strip_Bin
% %     maxstripidx = max(Strip.Strip_Bin(2,f));        % ������(���ͻ���ǰ)��strip����˳������ֵ,���豾stripidx
% %     if ~isempty(maxstripidx)
% %         % 2.1 �������ڱ�Bin�� & �Ұڷ�˳������stripidx ��˳���1, ����ǰ�ڷ�
% %         Strip.Strip_Bin(2,f)  = Strip.Strip_Bin(2,f)  - 1;
% %         % 2.2 ��ǰstripidx�ڷŵ�(��ͬSID/EIDβ��), ��˳�����õ����
% %         Strip.Strip_Bin(2,stripidx) = maxstripidx;
% %     end
% %     
% % end

%% V1 repairStripPlace
% % function Strip = repairStripPlace(Strip,stripidx)
% %     % 1 �ҵ�stripidx��Ӧ��BIN�µ�����Strip�������߼�ֵ
% %     flagStrip   =  Strip.Strip_Bin(1,:) == Strip.Strip_Bin(1,stripidx);    % ����ͬ����stripidx��Bin�ڵ�strip�߼��ж�
% %     flagBigIdx = Strip.Strip_Bin(2,:) > Strip.Strip_Bin(2,stripidx);       % ���аڷ�˳������stripidx���߼��ж�
% % 
% %     % 2 �������ڱ�Bin�� & �Ұڷ�˳������stripidx ��˳���1, ����ǰ�ڷ�
% %     Strip.Strip_Bin(2,flagBigIdx & flagStrip)  = Strip.Strip_Bin(2,flagBigIdx & flagStrip)  - 1;
% %     Strip.Strip_Bin(2,stripidx) = sum(flagStrip); % ��ǰstripidx�ڷŵ���β, ��˳�����õ����
% % 
% %     % [~,maxSeq]=max(Strip.Strip_Bin(2,Strip.Strip_Bin(1,:) == Strip.Strip_Bin(1,stripidx) ));
% % end



% % %% ����: Strip˦β v2
% % function   [tempStrip_Bin, StripisShuaiWei,LUisShuaiWei,TF] = HStripSW(Strip,LU)
% % % 1 ��Щ˦β: ��Ȳ���isWidthFull��߶Ȳ���isHeightFull��
% % % 2 ��Щ�ǿ�Ȳ���: 
% % % 3 ��Щ�Ǹ߶Ȳ���: 
% %     
% % 
% %     tempStrip_Bin = Strip.Strip_Bin(2,:);
% %     TF = false;    
% %     LUisShuaiWei = zeros(size(LU.Weight));      % �ж�LU�Ƿ�˦β������,������ͼʱ��������
% %     StripisShuaiWei = zeros(size(Strip.Weight));       % StripseqShuaiWei = zeros(size(Strip.Weight)); % seqShuaiWeiԽ��,����Խ��˦β; ԽС, Խ��˦β, ��Խ�����ڳ�ͷ.  
% %     
% % %% 1: ********************** ˦β ********************************** 
% % % 1 ��Щ˦β: ��Ȳ���isWidthFull��߶Ȳ���isHeightFull��
% % if any(~Strip.isWidthFull | ~Strip.isHeightFull) % | ~Strip.isHeightBalance)
% %     
% %     fprintf(1,'       Exsiting ˦β���� in HStripSW (Strip��߲���)...\n');
% %     TF = true;
% %     
% %     % Get b : strip index to be move to end of Vehicle
% %     [~,bNOTheightfull] = find(Strip.isHeightFull == 0);
% %     [~,bNOTwidthfull] = find(Strip.isWidthFull == 0);
% %     [~,bNOTheightbalance] = find(Strip.isHeightBalance == 0);
% %     
% %     b = unique([bNOTheightfull, bNOTwidthfull],'stable');    % ���ڷų�β��˳����ȫ��order
% % % b = unique([bNOTheightfull, bNOTwidthfull, bNOTheightbalance],'stable');   
% % 
% % % 2 ���������˦β��Strip, Ҫ�������? ��order
% %    if ~isempty(b)       
% %        %��1 Sort b by LoadingRate,Height etc.
% %             %tmpM = [Strip.loadingrate(b); Strip.maxHeight(b);Strip.loadingrateLimit(b);];
% %             %[~,order] = sortrows(tmpM',[1,2,3],{'descend','descend','descend'});
% %        %��2 Sort b by 1 ʣ���ȴ�ĺ��; 2 ʣ��߶ȶ�ĺ����
% %             %tmpM = [Strip.LW(1,b); Strip.maxHeight(b)];
% %             %[~,order] = sortrows(tmpM',[1,2],{'ascend','ascend'});
% % 
% %         %��3 Sort b 1 �ǵ����ķ�����; �ǻ��isMixed�ķ�����; �߶Ⱦ���isHeightBalance�ķ�����; ��͸߶�lowestHeight�ݼ��ķ�����
% %              %tmpM = [Strip.isSingleItem(b);Strip.isMixed(b);Strip.isHeightBalance(b); Strip.lowestHeight(b); Strip.maxHeight(b);Strip.meanHeight(b)];
% %             %[~,order] = sortrows(tmpM',[1,2,3,4,5],{'ascend','ascend','descend','descend','descend'});      
% %        
% %        %��4 Sort b 1 �ǻ��isMixed�ķ�����;  ƽ���߶�t�ݼ��ķ�����
% %        tmpM = [Strip.isMixed(b);Strip.meanHeight(b)];
% %        [~,order] = sortrows(tmpM',[1,2],{'ascend','descend'});   
% %        
% %        b = b(order);
% %        
% %        StripisShuaiWei(b) = 1;        % StripseqShuaiWei(b) = order;
% % 
% %        %%% *********** LU˦β���ҳ���,��ͼ�� ***************
% %        LUisShuaiWei( ismember(LU.LU_Strip(1,:), b) ) = 1;  % ��strip�ڵ�ֵ
% %    end   
% %    
% % %%%  ***************** �Ƿ�˦β�Ŀ��� *************
% % 
% %     for i=1:length(b)
% %         %     Strip = repairStripPlace(Strip,b(i));    %V1 ������SID/EID��˦β Strip.Strip_Bin
% %         %     Strip = repairStripPlace2(Strip,b(i));    %V2 ����SID/EDI��˦β Strip.Strip_Bin
% %         tempStrip_Bin = repairStripPlace2(Strip,b(i));    %V3 ����Strip.Strip_Bin
% %         Strip.Strip_Bin(2,:) = tempStrip_Bin;  % ������ Strip.Strip_Bin(2,:) ��Ҫѭ����ȥ
% %     end
% % 
% %        
% % end
% % 
% % %% �ֲ�����
% % %% V3 ���ٷ��ز��� 555 �ص����Ѻ���
% % function Strip_Bin_Seq = repairStripPlace2(Strip,stripidx)
% %     Strip_Bin_Seq = Strip.Strip_Bin(2,:);
% %     
% %     % 1 �ж����stripidx��SID�Ļ�ϻ�EID�Ļ��,����˦β,ֱ��return
% %     if ~isscalar(Strip.SID{stripidx}) || ~isscalar(Strip.EID{stripidx})
% %         return; end    
% %     
% %     % 2 �ж����stripidx����SID/EID��ϵ�,
% %     % ����Ƿǻ�ϵ�SID/EID,����˦β,������˦����β,��ʵ˦����stripidx��ͬ��SID/EIDβ��(��˦SID��,����˦EID��)
% %     sid = Strip.SID{stripidx};
% %     eid = Strip.EID{stripidx};
% %     bid = Strip.Strip_Bin(1,stripidx); 
% %     
% %     [tSID, ~] = padcat(Strip.SID{:});    if iscolumn(tSID), tSID = tSID'; end
% %     [tEID, ~] = padcat(Strip.EID{:});    if iscolumn(tEID), tEID = tEID'; end
% %     
% %     % 3 �ҵ�stripidx��Ӧ��BIN�µ�����Strip�������߼�ֵ
% %     fsid = Strip.isMixedSID == 0  &  any(ismember(tSID, sid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
% %     feid = Strip.isMixedEID == 0  &  any(ismember(tEID, eid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
% %     fbid = Strip.Strip_Bin(1,:) == bid;
% %     fbig =  Strip.Strip_Bin(2,:) > Strip.Strip_Bin(2,stripidx);       % ���аڷ�˳������stripidx���߼��ж�    
% %     f = fbid&fsid&feid&fbig;
% %     
% %     % 4 ����˳��Strip.Strip_Bin
% %     maxstripidx = max(Strip.Strip_Bin(2,f));        % ������(���ͻ���ǰ)��strip����˳������ֵ,���豾stripidx
% %     if ~isempty(maxstripidx)
% %         % 2.1 �������ڱ�Bin�� & �Ұڷ�˳������stripidx ��˳���1, ����ǰ�ڷ�
% %         Strip_Bin_Seq(f)  = Strip.Strip_Bin(2,f)  - 1;
% %         % 2.2 ��ǰstripidx�ڷŵ�(��ͬSID/EIDβ��), ��˳�����õ����
% %         Strip_Bin_Seq(stripidx) = maxstripidx;
% %     end
% %     
% % end
% % 
% % %% V2 ����SID��EID��˦β
% % % % function Strip = repairStripPlace2(Strip,stripidx)
% % % %     % 1 �ж����stripidx��SID�Ļ�ϻ�EID�Ļ��,����˦β,ֱ��return
% % % %     if ~isscalar(Strip.SID{stripidx}) || ~isscalar(Strip.EID{stripidx})
% % % %         return; end    
% % % %     
% % % %     % 2 �ж����stripidx����SID/EID��ϵ�,
% % % %     % ����Ƿǻ�ϵ�SID/EID,����˦β,������˦����β,��ʵ˦����stripidx��ͬ��SID/EIDβ��(��˦SID��,����˦EID��)
% % % %     sid = Strip.SID{stripidx};
% % % %     eid = Strip.EID{stripidx};
% % % %     bid = Strip.Strip_Bin(1,stripidx); 
% % % %     
% % % %     [tSID, ~] = padcat(Strip.SID{:});    if iscolumn(tSID), tSID = tSID'; end
% % % %     [tEID, ~] = padcat(Strip.EID{:});    if iscolumn(tEID), tEID = tEID'; end
% % % %     
% % % %     % 3 �ҵ�stripidx��Ӧ��BIN�µ�����Strip�������߼�ֵ
% % % %     fsid = Strip.isMixedSID == 0  &  any(ismember(tSID, sid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
% % % %     feid = Strip.isMixedEID == 0  &  any(ismember(tEID, eid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
% % % %     fbid = Strip.Strip_Bin(1,:) == bid;
% % % %     fbig =  Strip.Strip_Bin(2,:) > Strip.Strip_Bin(2,stripidx);       % ���аڷ�˳������stripidx���߼��ж�    
% % % %     f = fbid&fsid&feid&fbig;
% % % %     
% % % %     % 4 ����˳��Strip.Strip_Bin
% % % %     maxstripidx = max(Strip.Strip_Bin(2,f));        % ������(���ͻ���ǰ)��strip����˳������ֵ,���豾stripidx
% % % %     if ~isempty(maxstripidx)
% % % %         % 2.1 �������ڱ�Bin�� & �Ұڷ�˳������stripidx ��˳���1, ����ǰ�ڷ�
% % % %         Strip.Strip_Bin(2,f)  = Strip.Strip_Bin(2,f)  - 1;
% % % %         % 2.2 ��ǰstripidx�ڷŵ�(��ͬSID/EIDβ��), ��˳�����õ����
% % % %         Strip.Strip_Bin(2,stripidx) = maxstripidx;
% % % %     end
% % % %     
% % % % end
% % 
% % %% V1 repairStripPlace
% % % % function Strip = repairStripPlace(Strip,stripidx)
% % % %     % 1 �ҵ�stripidx��Ӧ��BIN�µ�����Strip�������߼�ֵ
% % % %     flagStrip   =  Strip.Strip_Bin(1,:) == Strip.Strip_Bin(1,stripidx);    % ����ͬ����stripidx��Bin�ڵ�strip�߼��ж�
% % % %     flagBigIdx = Strip.Strip_Bin(2,:) > Strip.Strip_Bin(2,stripidx);       % ���аڷ�˳������stripidx���߼��ж�
% % % % 
% % % %     % 2 �������ڱ�Bin�� & �Ұڷ�˳������stripidx ��˳���1, ����ǰ�ڷ�
% % % %     Strip.Strip_Bin(2,flagBigIdx & flagStrip)  = Strip.Strip_Bin(2,flagBigIdx & flagStrip)  - 1;
% % % %     Strip.Strip_Bin(2,stripidx) = sum(flagStrip); % ��ǰstripidx�ڷŵ���β, ��˳�����õ����
% % % % 
% % % %     % [~,maxSeq]=max(Strip.Strip_Bin(2,Strip.Strip_Bin(1,:) == Strip.Strip_Bin(1,stripidx) ));
% % % % end
% % 
% % end

%% ����: Strip˦β V1
% function   [Strip,LUisShuaiWei,TF] = HStripSW(Strip,LU)
% % 1 ��Щ˦β: ��Ȳ���isWidthFull��߶Ȳ���isHeightFull��
% % 2 ��Щ�ǿ�Ȳ���: 
% % 3 ��Щ�Ǹ߶Ȳ���: 
%     
%     TF = false;
%     LUisShuaiWei = zeros(size(LU.Weight));      % �ж�LU�Ƿ�˦β������,������ͼʱ��������
%     Strip.isShuaiWei = zeros(size(Strip.Weight));
%     
%     % StripseqShuaiWei = zeros(size(Strip.Weight)); % seqShuaiWeiԽ��,����Խ��˦β; ԽС, Խ��˦β, ��Խ�����ڳ�ͷ.  
%     
% %% 1: ********************** ˦β ********************************** 
% % 1 ��Щ˦β: ��Ȳ���isWidthFull��߶Ȳ���isHeightFull��
% % Strip.isHeightFull
% % Strip.isWidthFull
% if any(~Strip.isWidthFull | ~Strip.isHeightFull)
%     fprintf(1,'       Exsiting ˦β���� in HStripSW (Strip��߲���)...\n');
%             Strip.isWidthFull 
%             Strip.isHeightFull
%             Strip.isHeightBalance
%             ~Strip.isWidthFull | ~Strip.isHeightFull
%     
%     TF = true;
%     
%     % Get b : strip index to be move to end of Vehicle
%     [~,bNOTheightfull] = find(Strip.isHeightFull == 0);
%     [~,bNOTwidthfull] = find(Strip.isWidthFull == 0);
%     b = unique([bNOTheightfull, bNOTwidthfull],'stable');    % ���ڷų�β��˳����ȫ��order
% 
% % 2 ���������˦β��Strip, Ҫ�������? ��order
%    if ~isempty(b)       
%        %��1 Sort b by LoadingRate,Height etc.
%             %tmpM = [Strip.loadingrate(b); Strip.maxHeight(b);Strip.loadingrateLimit(b);];
%             %[~,order] = sortrows(tmpM',[1,2,3],{'descend','descend','descend'});
%        %��2 Sort b by 1 ʣ���ȴ�ĺ��; 2 ʣ��߶ȶ�ĺ����
%             %tmpM = [Strip.LW(1,b); Strip.maxHeight(b)];
%             %[~,order] = sortrows(tmpM',[1,2],{'ascend','ascend'});
% 
%         %��3 Sort b 1 �ǵ����ķ�����; �ǻ��isMixed�ķ�����; �߶Ⱦ���isHeightBalance�ķ�����; ��͸߶�lowestHeight�ݼ��ķ�����
%              %tmpM = [Strip.isSingleItem(b);Strip.isMixed(b);Strip.isHeightBalance(b); Strip.lowestHeight(b); Strip.maxHeight(b);Strip.meanHeight(b)];
%             %[~,order] = sortrows(tmpM',[1,2,3,4,5],{'ascend','ascend','descend','descend','descend'});      
%        
%        %��4 Sort b 1 �ǻ��isMixed�ķ�����;  ƽ���߶�t�ݼ��ķ�����
%        tmpM = [Strip.isMixed(b);Strip.meanHeight(b)];
%        [~,order] = sortrows(tmpM',[1,2],{'ascend','descend'});   
%        
%        b = b(order);
%        Strip.isShuaiWei(b) = 1;        % StripseqShuaiWei(b) = order;
%        
%        tmpLIDcellarray = Strip.LID(b);  tmpLIDmatarray=vertcat(tmpLIDcellarray{:});
%        LUisShuaiWei(tmpLIDmatarray)=1; 
%        
%    end   
%    
% %%%  ***************** �Ƿ�˦β�Ŀ��� *************
% for i=1:length(b)
% %     Strip = repairStripPlace(Strip,b(i));    %V1 ������SID/EID��˦β Strip.Strip_Bin
% %     Strip = repairStripPlace2(Strip,b(i));    %V2 ����SID/EDI��˦β Strip.Strip_Bin
% % Strip.Strip_Bin(2,:) 
%         tempStrip_Bin = repairStripPlace2(Strip,b(i));    %V3 ����Strip.Strip_Bin
%         Strip.Strip_Bin(2,:) = tempStrip_Bin;
% end
% 
% end
% 
% %% �ֲ�����
% %% V3 ���ٷ��ز���
% function Strip_Bin_Seq = repairStripPlace2(Strip,stripidx)
%     Strip_Bin_Seq = Strip.Strip_Bin(2,:);
%     
%     % 1 �ж����stripidx��SID�Ļ�ϻ�EID�Ļ��,����˦β,ֱ��return
%     if ~isscalar(Strip.SID{stripidx}) || ~isscalar(Strip.EID{stripidx})
%         return; end    
%     
%     % 2 �ж����stripidx����SID/EID��ϵ�,
%     % ����Ƿǻ�ϵ�SID/EID,����˦β,������˦����β,��ʵ˦����stripidx��ͬ��SID/EIDβ��(��˦SID��,����˦EID��)
%     sid = Strip.SID{stripidx};
%     eid = Strip.EID{stripidx};
%     bid = Strip.Strip_Bin(1,stripidx); 
%     
%     [tSID, ~] = padcat(Strip.SID{:});    if iscolumn(tSID), tSID = tSID'; end
%     [tEID, ~] = padcat(Strip.EID{:});    if iscolumn(tEID), tEID = tEID'; end
%     
%     % 3 �ҵ�stripidx��Ӧ��BIN�µ�����Strip�������߼�ֵ
%     fsid = Strip.isMixedSID == 0  &  any(ismember(tSID, sid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
%     feid = Strip.isMixedEID == 0  &  any(ismember(tEID, eid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
%     fbid = Strip.Strip_Bin(1,:) == bid;
%     fbig =  Strip.Strip_Bin(2,:) > Strip.Strip_Bin(2,stripidx);       % ���аڷ�˳������stripidx���߼��ж�    
%     f = fbid&fsid&feid&fbig;
%     
%     % 4 ����˳��Strip.Strip_Bin
%     maxstripidx = max(Strip.Strip_Bin(2,f));        % ������(���ͻ���ǰ)��strip����˳������ֵ,���豾stripidx
%     if ~isempty(maxstripidx)
%         % 2.1 �������ڱ�Bin�� & �Ұڷ�˳������stripidx ��˳���1, ����ǰ�ڷ�
%         Strip_Bin_Seq(f)  = Strip.Strip_Bin(2,f)  - 1;
%         % 2.2 ��ǰstripidx�ڷŵ�(��ͬSID/EIDβ��), ��˳�����õ����
%         Strip_Bin_Seq(stripidx) = maxstripidx;
%     end
%     
% end
% 
% %% V2 ����SID��EID��˦β
% % % function Strip = repairStripPlace2(Strip,stripidx)
% % %     % 1 �ж����stripidx��SID�Ļ�ϻ�EID�Ļ��,����˦β,ֱ��return
% % %     if ~isscalar(Strip.SID{stripidx}) || ~isscalar(Strip.EID{stripidx})
% % %         return; end    
% % %     
% % %     % 2 �ж����stripidx����SID/EID��ϵ�,
% % %     % ����Ƿǻ�ϵ�SID/EID,����˦β,������˦����β,��ʵ˦����stripidx��ͬ��SID/EIDβ��(��˦SID��,����˦EID��)
% % %     sid = Strip.SID{stripidx};
% % %     eid = Strip.EID{stripidx};
% % %     bid = Strip.Strip_Bin(1,stripidx); 
% % %     
% % %     [tSID, ~] = padcat(Strip.SID{:});    if iscolumn(tSID), tSID = tSID'; end
% % %     [tEID, ~] = padcat(Strip.EID{:});    if iscolumn(tEID), tEID = tEID'; end
% % %     
% % %     % 3 �ҵ�stripidx��Ӧ��BIN�µ�����Strip�������߼�ֵ
% % %     fsid = Strip.isMixedSID == 0  &  any(ismember(tSID, sid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
% % %     feid = Strip.isMixedEID == 0  &  any(ismember(tEID, eid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
% % %     fbid = Strip.Strip_Bin(1,:) == bid;
% % %     fbig =  Strip.Strip_Bin(2,:) > Strip.Strip_Bin(2,stripidx);       % ���аڷ�˳������stripidx���߼��ж�    
% % %     f = fbid&fsid&feid&fbig;
% % %     
% % %     % 4 ����˳��Strip.Strip_Bin
% % %     maxstripidx = max(Strip.Strip_Bin(2,f));        % ������(���ͻ���ǰ)��strip����˳������ֵ,���豾stripidx
% % %     if ~isempty(maxstripidx)
% % %         % 2.1 �������ڱ�Bin�� & �Ұڷ�˳������stripidx ��˳���1, ����ǰ�ڷ�
% % %         Strip.Strip_Bin(2,f)  = Strip.Strip_Bin(2,f)  - 1;
% % %         % 2.2 ��ǰstripidx�ڷŵ�(��ͬSID/EIDβ��), ��˳�����õ����
% % %         Strip.Strip_Bin(2,stripidx) = maxstripidx;
% % %     end
% % %     
% % % end
% 
% %% V1 repairStripPlace
% % % function Strip = repairStripPlace(Strip,stripidx)
% % %     % 1 �ҵ�stripidx��Ӧ��BIN�µ�����Strip�������߼�ֵ
% % %     flagStrip   =  Strip.Strip_Bin(1,:) == Strip.Strip_Bin(1,stripidx);    % ����ͬ����stripidx��Bin�ڵ�strip�߼��ж�
% % %     flagBigIdx = Strip.Strip_Bin(2,:) > Strip.Strip_Bin(2,stripidx);       % ���аڷ�˳������stripidx���߼��ж�
% % % 
% % %     % 2 �������ڱ�Bin�� & �Ұڷ�˳������stripidx ��˳���1, ����ǰ�ڷ�
% % %     Strip.Strip_Bin(2,flagBigIdx & flagStrip)  = Strip.Strip_Bin(2,flagBigIdx & flagStrip)  - 1;
% % %     Strip.Strip_Bin(2,stripidx) = sum(flagStrip); % ��ǰstripidx�ڷŵ���β, ��˳�����õ����
% % % 
% % %     % [~,maxSeq]=max(Strip.Strip_Bin(2,Strip.Strip_Bin(1,:) == Strip.Strip_Bin(1,stripidx) ));
% % % end
% 
% end


%% V1 HStripSW v2: �򻯲���Ҫ�����
% % %% ����: Strip˦β
% % function   [Strip,LUisShuaiWei,TF] = HStripSW(Strip,LU)
% % % 1 ��Щ˦β: ��Ȳ���isWidthFull��߶Ȳ���isHeightFull��
% % % 2 ��Щ�ǿ�Ȳ���: 
% % % 3 ��Щ�Ǹ߶Ȳ���: 
% % 
% %     Strip.isShuaiWei = zeros(size(Strip.Weight));
% %     Strip.seqShuaiWei = zeros(size(Strip.Weight)); % seqShuaiWeiԽ��,����Խ��˦β; ԽС, Խ��˦β, ��Խ�����ڳ�ͷ.
% %     LUisShuaiWei = zeros(size(LU.Weight));      % �ж�LU�Ƿ�˦β������,������ͼʱ��������
% %     TF = false;
% % %% 1: ********************** ˦β ********************************** 
% % % 1 ��Щ˦β: ��Ȳ���isWidthFull��߶Ȳ���isHeightFull��
% % % Strip.isHeightFull
% % % Strip.isWidthFull
% % if any(~Strip.isWidthFull | ~Strip.isHeightFull)
% %     fprintf(1,'       Exsiting ˦β���� in HStripSW (Strip��߲���)...\n');
% %             Strip.isWidthFull 
% %             Strip.isHeightFull
% %             Strip.isHeightBalance
% %             ~Strip.isWidthFull | ~Strip.isHeightFull
% %     
% %     TF = true;
% %     % Get b : strip index to be move to end of Vehicle
% %     [~,bNOTheightfull] = find(Strip.isHeightFull == 0);
% %     [~,bNOTwidthfull] = find(Strip.isWidthFull == 0);
% %     b = unique([bNOTheightfull, bNOTwidthfull],'stable');    % ���ڷų�β��˳����ȫ��order
% % 
% % % 2 ���������˦β��Strip, Ҫ�������? ��order
% %    if ~isempty(b)       
% %        %��1 Sort b by LoadingRate,Height etc.
% %             %tmpM = [Strip.loadingrate(b); Strip.maxHeight(b);Strip.loadingrateLimit(b);];
% %             %[~,order] = sortrows(tmpM',[1,2,3],{'descend','descend','descend'});
% %        %��2 Sort b by 1 ʣ���ȴ�ĺ��; 2 ʣ��߶ȶ�ĺ����
% %             %tmpM = [Strip.LW(1,b); Strip.maxHeight(b)];
% %             %[~,order] = sortrows(tmpM',[1,2],{'ascend','ascend'});
% % 
% %         %��3 Sort b 1 �ǵ����ķ�����; �ǻ��isMixed�ķ�����; �߶Ⱦ���isHeightBalance�ķ�����; ��͸߶�lowestHeight�ݼ��ķ�����
% %              %tmpM = [Strip.isSingleItem(b);Strip.isMixed(b);Strip.isHeightBalance(b); Strip.lowestHeight(b); Strip.maxHeight(b);Strip.meanHeight(b)];
% %             %[~,order] = sortrows(tmpM',[1,2,3,4,5],{'ascend','ascend','descend','descend','descend'});      
% %        
% %        %��4 Sort b 1 �ǻ��isMixed�ķ�����;  ƽ���߶�t�ݼ��ķ�����
% %        tmpM = [Strip.isMixed(b);Strip.meanHeight(b)];
% %        [~,order] = sortrows(tmpM',[1,2],{'ascend','descend'});   
% %        
% %        b = b(order);
% %        Strip.isShuaiWei(b) = 1;
% %        Strip.seqShuaiWei(b) = order;
% %        
% %        tmpLIDcellarray = Strip.LID(b);  tmpLIDmatarray=vertcat(tmpLIDcellarray{:});
% %        LUisShuaiWei(tmpLIDmatarray)=1;      
% %        
% %                     %    Strip.seqSW                                                                     (b) = 1:length(b);
% %                     %     Strip.loadingrateLimit(b)
% %    end   
% %    
% % %%%  ***************** �Ƿ�˦β�Ŀ��� *************
% % for i=1:length(b)
% % %     Strip = repairStripPlace(Strip,b(i));    %V1 ������SID/EID��˦β Strip.Strip_Bin
% %     Strip = repairStripPlace2(Strip,b(i));    %V2 ����SID/EDI��˦β Strip.Strip_Bin
% % end
% % end
% % 
% % %% �ֲ�����
% % %% V2 ����SID��EID��˦β
% % function Strip = repairStripPlace2(Strip,stripidx)
% %     % 1 �ж����stripidx��SID�Ļ�ϻ�EID�Ļ��,����˦β,ֱ��return
% %     if ~isscalar(Strip.SID{stripidx}) || ~isscalar(Strip.EID{stripidx})
% %         return; end    
% %     
% %     % 2 �ж����stripidx����SID/EID��ϵ�,
% %     % ����Ƿǻ�ϵ�SID/EID,����˦β,������˦����β,��ʵ˦����stripidx��ͬ��SID/EIDβ��(��˦SID��,����˦EID��)
% %     sid = Strip.SID{stripidx};
% %     eid = Strip.EID{stripidx};
% %     bid = Strip.Strip_Bin(1,stripidx); 
% %     
% %     [tSID, ~] = padcat(Strip.SID{:});    if iscolumn(tSID), tSID = tSID'; end
% %     [tEID, ~] = padcat(Strip.EID{:});    if iscolumn(tEID), tEID = tEID'; end
% %     
% %     % 3 �ҵ�stripidx��Ӧ��BIN�µ�����Strip�������߼�ֵ
% %     fsid = Strip.isMixedSID == 0  &  any(ismember(tSID, sid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
% %     feid = Strip.isMixedEID == 0  &  any(ismember(tEID, eid),1);  % ��һstrip��sidֻҪ������stripidx��sid,��=1
% %     fbid = Strip.Strip_Bin(1,:) == bid;
% %     fbig =  Strip.Strip_Bin(2,:) > Strip.Strip_Bin(2,stripidx);       % ���аڷ�˳������stripidx���߼��ж�    
% %     f = fbid&fsid&feid&fbig;
% %     
% %     % 4 ����˳��Strip.Strip_Bin
% %     maxstripidx = max(Strip.Strip_Bin(2,f));        % ������(���ͻ���ǰ)��strip����˳������ֵ,���豾stripidx
% %     if ~isempty(maxstripidx)
% %         % 2.1 �������ڱ�Bin�� & �Ұڷ�˳������stripidx ��˳���1, ����ǰ�ڷ�
% %         Strip.Strip_Bin(2,f)  = Strip.Strip_Bin(2,f)  - 1;
% %         % 2.2 ��ǰstripidx�ڷŵ�(��ͬSID/EIDβ��), ��˳�����õ����
% %         Strip.Strip_Bin(2,stripidx) = maxstripidx;
% %     end
% %     
% % end
% % 
% % %% V1 repairStripPlace
% % % % function Strip = repairStripPlace(Strip,stripidx)
% % % %     % 1 �ҵ�stripidx��Ӧ��BIN�µ�����Strip�������߼�ֵ
% % % %     flagStrip   =  Strip.Strip_Bin(1,:) == Strip.Strip_Bin(1,stripidx);    % ����ͬ����stripidx��Bin�ڵ�strip�߼��ж�
% % % %     flagBigIdx = Strip.Strip_Bin(2,:) > Strip.Strip_Bin(2,stripidx);       % ���аڷ�˳������stripidx���߼��ж�
% % % % 
% % % %     % 2 �������ڱ�Bin�� & �Ұڷ�˳������stripidx ��˳���1, ����ǰ�ڷ�
% % % %     Strip.Strip_Bin(2,flagBigIdx & flagStrip)  = Strip.Strip_Bin(2,flagBigIdx & flagStrip)  - 1;
% % % %     Strip.Strip_Bin(2,stripidx) = sum(flagStrip); % ��ǰstripidx�ڷŵ���β, ��˳�����õ����
% % % % 
% % % %     % [~,maxSeq]=max(Strip.Strip_Bin(2,Strip.Strip_Bin(1,:) == Strip.Strip_Bin(1,stripidx) ));
% % % % end
% % 
% % end
