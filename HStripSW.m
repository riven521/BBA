%% ����: Strip˦β
function   [Strip,LUisShuaiWei] = HStripSW(Strip,LU)
%% ��ʼ��

    Strip.isShuaiWei = zeros(size(Strip.Weight));
    Strip.seqShuaiWei = zeros(size(Strip.Weight)); % seqShuaiWeiԽ��,����Խ��˦β; ԽС, Խ��˦β, ��Խ�����ڳ�ͷ.
    LUisShuaiWei = zeros(size(LU.Weight));
    
%% 1: ********************** ˦β ********************************** 
% 1 ��Щ˦β: ��Ȳ���isWidthFull��߶Ȳ���isHeightFull��
Strip.isHeightFull
Strip.isWidthFull
if any(~Strip.isWidthFull | ~Strip.isHeightFull)
    % Get b : strip index to be move to end of Vehicle
    [~,bNOTheightfull] = find(Strip.isHeightFull == 0);
    [~,bNOTwidthfull] = find(Strip.isWidthFull == 0);
    b = unique([bNOTheightfull, bNOTwidthfull],'stable');    % ���ڷų�β��˳����ȫ��order

% 2 ���������˦β��Strip, Ҫ�������? ��order
   if ~isempty(b)       
       %��1 Sort b by LoadingRate,Height etc.
        tmpM = [Strip.loadingrate(b); Strip.maxHeight(b);Strip.loadingrateLimit(b);];
        [~,order] = sortrows(tmpM',[1,2,3],{'descend','descend','descend'});
       %��2 Sort b by 1 ʣ���ȴ�ĺ��; 2 ʣ��߶ȶ�ĺ����
       tmpM = [Strip.LW(1,b); Strip.maxHeight(b)];
       [~,order] = sortrows(tmpM',[1,2],{'ascend','ascend'});

        %��3 Sort b 1 �ǵ����ķ�����; �ǻ��isMixed�ķ�����; �߶Ⱦ���isHeightBalance�ķ�����; ��͸߶�lowestHeight�ݼ��ķ�����
       tmpM = [Strip.isSingleItem(b);Strip.isMixed(b);Strip.isHeightBalance(b); Strip.lowestHeight(b); Strip.maxHeight(b);Strip.meanHeight(b)];
       [~,order] = sortrows(tmpM',[1,2,3,4,5],{'ascend','ascend','descend','descend','descend'});      
       
       %��4 Sort b 1 �ǻ��isMixed�ķ�����;  ƽ���߶�t�ݼ��ķ�����
       tmpM = [Strip.isMixed(b);Strip.meanHeight(b)];
       [~,order] = sortrows(tmpM',[1,2],{'ascend','descend'});   
       
       b = b(order);
       Strip.isShuaiWei(b) = 1;
       Strip.seqShuaiWei(b) = order;
       
       tmpLIDcellarray = Strip.LID(b);  tmpLIDmatarray=vertcat(tmpLIDcellarray{:});
       LUisShuaiWei(tmpLIDmatarray)=1;      
       
                    %    Strip.seqSW(b) = 1:length(b);
                    %     Strip.loadingrateLimit(b)
   end   
   
%%%  ***************** �Ƿ�˦β�Ŀ��� *************
for i=1:length(b)
    Strip = repairStripPlace(Strip,b(i));    % Strip.Strip_Bin
end
end

%% �ֲ�����
function Strip = repairStripPlace(Strip,stripidx)
    % 1 �ҵ�stripidx��Ӧ��BIN�µ�����Strip�������߼�ֵ
    flagStrip   =  Strip.Strip_Bin(1,:) == Strip.Strip_Bin(1,stripidx);    % ����ͬ����stripidx��Bin�ڵ�strip�߼��ж�
    flagBigIdx = Strip.Strip_Bin(2,:) > Strip.Strip_Bin(2,stripidx);       % ���аڷ�˳������stripidx���߼��ж�

    % 2 �������ڱ�Bin�� & �Ұڷ�˳������stripidx ��˳���1, ����ǰ�ڷ�
    Strip.Strip_Bin(2,flagBigIdx & flagStrip)  = Strip.Strip_Bin(2,flagBigIdx & flagStrip)  - 1;
    Strip.Strip_Bin(2,stripidx) = sum(flagStrip); % ��ǰstripidx�ڷŵ���β, ��˳�����õ����

    % [~,maxSeq]=max(Strip.Strip_Bin(2,Strip.Strip_Bin(1,:) == Strip.Strip_Bin(1,stripidx) ));
end

end
