function [LU,Veh] = cpuLUVeh(LU,Veh)
% cpuLUVeh ��Ҫ����:
% ���� Veh ��Volume/order
% ���� LU ��ID/LU��LWH/LU��maxL/LU��maxHLayer/LU��nbID/LU��isNonMixed/isMixedTile
        
     % global  ISisGpreprocLU1

    %% 1 ���³���Veh˳�� ����order
    % Veh�������->С   Ĭ��˳��
    Veh.Volume = prod(Veh.LWH);
    [~,order] = sortrows(Veh.Volume', [1],{'descend'});    
    Veh = structfun(@(x) x(:,order),Veh,'UniformOutput',false);
    if ~isrow(order)
        order = order';
    end
    Veh.order = order;
    
     %% 2 ����LU����ID�ţ�����Old��
     % 1 �����ͬID���£���ӦSID��Ҫ���벻ͬ������ͬ�ı�ID�ţ�ֱ����������ͬ��ID�ڲ�ͬSID��;
     % TODO : ���Ӳ�ͬSID/��ͬEID�µ�ID�Ų����ظ�
     if isrepeated(LU.ID,LU.SID)
        warning('��������ID���ڲ�ͬ SID���ظ�, ��Ҫ����'); 
        LU.ID = reviseID(LU.ID,LU.SID);
     end

     if isrepeated(LU.EID,LU.SID)
         warning('��������EP LOCATION ID���ڲ�ͬSID���ظ�, ��Ҫ����');
         LU.EID = reviseID(LU.EID,LU.SID);
     end
     
     if isrepeated(LU.PID,LU.SID)
         warning('��������PID ID���ڲ�ͬSID���ظ�, ��Ҫ����');
         LU.PID = reviseID(LU.PID,LU.SID);
     end
     
     % �Ƿ���Ҫ��fixme 
      if isrepeated(LU.LID,LU.SID)
         warning('��������LID ID���ڲ�ͬSID���ظ�, ��Ҫ����');
         LU.LID = reviseID(LU.LID,LU.SID);
      end    
     
     if isrepeated(LU.ID,LU.EID)
         error('��������ID���ڲ�ͬ EID���ظ�, ��Ҫ����');
     end
     
    % 1 ID ת��Ϊ��1��ʼ������� ������������ID����Ϣ    TODO ����LID��OLID��
    if isfield(LU, 'ID'),    [LU.ID,LU.OID]      = idExchange(LU.ID);    end
    if isfield(LU, 'PID'),  [LU.PID,LU.OPID]  = idExchange(LU.PID); end
    if isfield(LU, 'EID'),  [LU.EID,LU.OEID]  = idExchange(LU.EID);  end
    if isfield(LU, 'SID'),  [LU.SID,LU.OSID]  = idExchange(LU.SID);  end
    
    % 1 �����ͬID��PID/EID/LID�ȣ����£���ӦSID��Ҫ���벻ͬ������ͬ�ı�ID�ţ�ֱ����������ͬ��ID�ڲ�ͬSID��;
    if isrepeated(LU.ID,LU.SID),  error('����ID���ظ�, Ӧ����checkʱ�ѵ���');   end
    
    %% 3 ����LU����margin��LU��LWH������LU��Rotaed���
% %     % GET LU's LWH with margin and Rotaed or not
% %     % 2 Input���Ӽ�϶BUFF���feasible��LU��BIN�ĳ����ת��
% %     LU.OLWH = LU.LWH;
% %     LU.LWH(1,:) =  LU.LWH(1,:) +  LU.margin(1,: ) + LU.margin(2,: ); %��ȣ����ң�
% %     LU.LWH(2,:) =  LU.LWH(2,:) +  LU.margin(3,: ) + LU.margin(4,: ); %���ȣ����£�
% %     
% %     [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,Veh.LWH(1,1));  %�ڶ���������  3��VEH�������Ұڷŵķ�϶��С����
% %     LU.LWH = getRotaedLWH(LU.LWH, LU.Rotaed, LU.margin);
% %     
% %     % 3 Ĭ�Ͻ�LUȫ������Horizontal������ת��ǰ�᣺��LU������ת��
% %     % NOTE: �˴������1: Horizontal�����LWH���Ƿ�Rotaed���
% %     % NOTE : ֱ���滻��ԭʼORIGINAL �� LWH
% %     % LU.LWH
% % 
% % %     if pwhichSortItemOrder ==1
% % %         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,1); %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��        
% % %     elseif pwhichSortItemOrder ==2
% % %         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,0); %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
% % %     elseif pwhichSortItemOrder ==3 %Ĭ�ϴ�ѡ��
% % %         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,Veh.LWH(1,1)); %�ڶ���������  3��VEH�������Ұڷŵķ�϶��С����
% % %     end
        
    %% 4 ����LU���ڵ�ǰ�����µ���󳤿�߲��� ����margin
    if ~isfield(LU, 'maxL') 
        LU.maxL(1,:) =  floor(Veh.LWH(1,1)./LU.LWH(1,:));
        LU.maxL(2,:) =  floor(Veh.LWH(2,1)./LU.LWH(2,:));
        LU.maxL(3,:) =  floor(Veh.LWH(3,1)./LU.LWH(3,:));   %����ÿ������LU�ĸ߶ȵ�������
    end

    if ~isfield(LU, 'maxHLayer') % �粻����,ֱ����maxLָʾ; �����,����Сֵȡ��
        LU.maxHLayer = LU.maxL(3,:); 
    else
        LU.maxHLayer = min( [LU.maxHLayer',LU.maxL(3,:)'], [], 2 )';  % ͬ�����¶�ӦLU��maxHLayer    
    end
    
    if ~isfield(LU, 'Rotaed')
        LU.Rotaed = zeros(size(LU.ID));
    end

    %% 5 ����LUͬ��ID/�ɶѶ�ID�µĸ���
    % V1"
    %     for i=1:length(LU.Weight)
    %         LU.nbID(i) = sum(LU.ID == LU.ID(i));
    %     end
    % V2: ֱ�Ӽ���
    LU.nbID = sum(LU.ID==LU.ID');
    LU.nbLID = sum(LU.LID==LU.LID');
    
    %% 6��555 ITEMBALANCE���� 7 V2: ����LU�����isNonMixed/isMixedTile�Ƿ�Ϊ����Ҫ  ���� ��ƴ/��ƴ���� 
    % 555 ����ΪLU��ͬ�߶ȣ����ݲ����жϣ� �� LU���ָ߶ȣ����ݸ߶��ж� TODO ���Ӳ����жϣ���������
    [LU] = cpuLU(LU,Veh);
    
  end

function [exID,ID] = idExchange(ID)
        uniID = unique(ID);
        exID=ID; %�м����
        for i=1:length(uniID)
            exID(ID(:)==uniID(i)) = i;
        end
end









%%
% % %      LU.isNonMixed = ones(1,length(LU.Weight))*-1;     %1 LU�����γ����� 0 �ض��з��������� LU��������֮һ(�����γɷ������LU��ǰ��ȡ�Ѷ�,Ŀ�ľ������ٻ�϶Ѷ������)
% % %      LU.isMixedTile = zeros(1,length(LU.Weight));          %1 ��isNonMixed=1ʱ,��������,isMixedTile��Ϊ0; ��isNonMixed=0ʱ, �з�����,�����µĽ���Ϊ�������LU��isMixedTile��ֵΪ1
% % %     % GET LU.isNonMixed: ����ÿ��LU�Ƿ�Ϊ����Ҫ��ƴ�Ŀ���
% % %     % ѭ��: ID����
% % %     tLU = getTableLU(LU);
% % %     %     [tLUsorted,torder] = sortrows(tLU,{'SID','ID','LID','H'},{'ascend','ascend','ascend','descend'});
% % %     [tLUsorted,torder] = sortrows(tLU,{'SID','EID','ID','LID','H'},{'ascend','ascend','ascend','ascend','descend'});  %V2: ����EID˳��
% % %                         %     CateOrder = tLUsorted(:,{'SID','ID','LID','H'});  %����maxHLayer?
% % %                         %     [GID,g] = findgroups(CateOrder)
% % %     VehHeight = Veh.LWH(3,1);  % �����߶�
% % %     VehWidth = Veh.LWH(1,1);  % �������
% % %     
% % %     %% �ԶѶ���ID����LU��isNonMixed��isMixedTile
% % %     uniID = unique(tLUsorted.ID);
% % %     for iLu=1:length(uniID)
% % %         % Item i ���ڵ�LU flag���
% % %         flagLU = tLUsorted.ID == uniID(iLu);
% % %         idxLU = find(flagLU);
% % %         unimaxHLayerperLUID = unique(tLUsorted.maxHLayer(flagLU));
% % %         
% % %         % 1.1 �����ӦID��LU�ĸ߶�H��һ��,���ò����ж�
% % %         if isscalar(unimaxHLayerperLUID) %�ƺ���ʹ�����÷ǵ���Ҳ����,�����ڻ�ϸ߶ȵ�Ҳ�ǿ����õ�; �������Բ���Ϊ����,�ر��ǿ���ƽ��
% % %             nb = sum(flagLU);
% % %             nbmod = mod(nb,unimaxHLayerperLUID);       if nb ==0 || nbmod>nb, error('Gpreproc�м���isNonMixed����'); end    
% % %             
% % %             % �ж�isNonMixed
% % %             if nbmod == 0 %modΪ0������ID������LU���������ķŵ��߶�����,�������и߶Ȳ�����Item����. ����Ҫ��� ����ϵ���ǰ��order����ǰ
% % %                 tLUsorted.isNonMixed(flagLU) = 1;
% % %             else
% % %                 tLUsorted.isNonMixed(flagLU)= 0;
% % %                 % ����tLUsorted.isMixedTile
% % %                 idxMod = idxLU(end-nbmod+1:end);        % ������������е�isMixedTile��ֵΪ1.%��Ϊ�Ѿ�����,���Կ��԰�����滻        
% % %                 tLUsorted.isMixedTile(idxMod)=1;
% % %             end
% % %             
% % %         % 1.2 �����ӦID��LU�ĸ߶�H��һ��,���÷ǲ���(�߶�)�ж�    (һ����ISisGpreprocLU1�Ѷ�߶Ⱦ���ʹ��)
% % %         else            
% % %             nbLU = length(idxLU);
% % %             idxVeh = zeros(nbLU,1);
% % %             tmpsum = 0;
% % %             tmpidx = 1;
% % %             
% % %             for idx=1:nbLU
% % %                 tmpsum = tmpsum + tLUsorted.H(idxLU(idx));
% % %                 if tmpsum > VehHeight
% % %                     tmpsum = tLUsorted.H(idxLU(idx));
% % %                     tmpidx=tmpidx+1;
% % %                 end
% % %                 idxVeh(idx) =tmpidx;
% % %             end
% % %             
% % %             % �ж�isNonMixed
% % %             if ( VehHeight - tmpsum ) > min(tLUsorted.H)  %��߶ȼ�϶ > ��СLU�߶�
% % %                 % ???? �����LU�Ŀ�ȿ��Է�2������, ��Ŀǰ�����ɵ�Ҳ��2������,����󼸸����и߶Ⱦ�������,
% % %                 maxWLayerLUID = unique(tLUsorted.maxL(flagLU,1));   if ~isscalar(maxWLayerLUID), error('Gpreproc�м���maxWLayerLUID����'); end
% % %                 if ISisGpreprocLU1==1 && maxWLayerLUID>1 && max(idxVeh) >1 && mod(max(idxVeh),maxWLayerLUID)~=1 %������1��Item������1��
% % %                     tLUsorted.isNonMixed(flagLU)= 1;
% % %                 else
% % %                     tLUsorted.isNonMixed(flagLU)= 0;
% % %                     idxMod = idxLU(idxVeh==max(idxVeh));
% % %                     tLUsorted.isMixedTile(idxMod)=1;
% % %                 end
% % %             else
% % %                 tLUsorted.isNonMixed(flagLU)= 1;
% % %             end
% % %             
% % %         end
% % %         
% % %     end
% % %     
% % %     % ��˳�򷵻�
% % %     [~,sorder] = sort(torder);
% % %     LU.isNonMixed = tLUsorted.isNonMixed(sorder)';
% % %     LU.isMixedTile = tLUsorted.isMixedTile(sorder)';
    
    
    %% 7 V1: BUG MAYBE ����LU�����isNonMixed/isMixedTile�Ƿ�Ϊ����Ҫ  ���� ��ƴ/��ƴ����
% % % %      LU.isNonMixed = ones(1,length(LU.Weight))*-1;    %Item�Ƿ����Ҫ����ж�,��ż������Item��ǰ����Strip����
% % % %      LU.isMixedTile = zeros(1,length(LU.Weight));          %Item���,�ҳ����������Item��β�и�ֵΪ1
% % % %     % GET LU.isNonMixed: ����ÿ��LU�Ƿ�Ϊ����Ҫ��ƴ�Ŀ���
% % % %     % ѭ��: LID����
% % % %     uniLID = unique(LU.LID);
% % % % %     uniLID = unique(LU.ID);
% % % %     for iLu=1:length(uniLID)
% % % %         % Item i ���ڵ�LU flag���
% % % %         VehHeight = Veh.LWH(3,1);  % �����߶�
% % % %         
% % % %         flagLU = LU.LID(:) == uniLID(iLu);
% % % % %         flagLU = LU.ID(:) == uniLID(iLu);
% % % %         LUHeight = unique(LU.LWH(3,flagLU));     % ͬ��LULIDʱ�� LU�߶�
% % % %             if length(unique(LUHeight)) > 2, warning('ͬ��LULIDʱ��LU�߶�ֵ>2,��Ԥ�ڴ���'); end %������������ɵĿ������������
% % % %         if length(unique(LUHeight)) >= 2, LUHeight = max(LUHeight); end % ��LU���ֵ��Ϊ�Ѷ��ж�ֵ
% % % % 
% % % %         LU.maxL(3,flagLU)
% % % %         LU.maxHLayer(flagLU)
% % % %         maxHeightLayer= floor(VehHeight/LUHeight); %LU�߶Ȳ���
% % % %         
% % % %         nb = sum(flagLU);
% % % %         nbmod = mod(nb,maxHeightLayer);
% % % %             if nb ==0 || nbmod>nb, error('Gpreproc�м���isNonMixed����'); end
% % % %             
% % % %         if nbmod == 0 %modΪ0���� ����Ҫ��� ����ϵ���ǰ��order����ǰ
% % % %             LU.isNonMixed(flagLU) = 1;
% % % %         else
% % % %             LU.isNonMixed(flagLU)= 0;            
% % % %             % ����LU��isMixedTile
% % % %             tmpSort=[LU.SID; LU.LWH(3,:); LU.PID; ]; %LU.Weight LU.maxL;            
% % % %             [~, order]=sortrows(tmpSort(:,flagLU)', [1,2,3],{'ascend','ascend','ascend'});
% % % %             flagLUIdx = find(flagLU);
% % % %             flagmodIdx = flagLUIdx(order(1:nbmod));
% % % %             LU.isMixedTile(flagmodIdx)=1;
% % % %         end
% % % %     end
    
    
%     LUID = getLUIDArray(LU); %% ���㣺LU����������� ��ʱ����
%      printInput();
%% Ƕ�׺���:
% function printInput()
%     fprintf('������ֻ��һ������ ��=%1.0f ��=%1.0f ��=%1.0f \n', unique(Veh.LWH','rows')');
%     fprintf('������ֻ��һ������ ���϶=%1.0f ����϶=%1.0f �߼�϶=%1.0f \n', unique(Veh.buff','rows')');
%     fprintf('�������� %d ����Ʒ,����߷ֱ�Ϊ \n',numel(LU.ID(:)));
%     fprintf('%1.1f %1.1f %1.1f \n',LU.LWH);
% end





% % %%  ��ȡLUID�����������(ͬ����ID��������������Ƿ����ת)
% % function LUID = getLUIDArray(LU)
% % LUID.ID = unique(LU.ID);
% % for iID = 1:length(LUID.ID)
% %     LUID.Weight(iID) = sum(LU.Weight .* (LU.ID == LUID.ID(iID)) );
% %     LUID.Volume(iID) = sum(prod(LU.LWH) .* (LU.ID == LUID.ID(iID)) );
% %     LUID.isRota(iID) =  unique(LU.isRota(LU.ID == LUID.ID(iID)));
% %     LUID.maxL(iID) =  unique(LU.maxL(LU.ID == LUID.ID(iID)));
% %     LUID.yID(iID) =  unique(LU.yID(LU.ID == LUID.ID(iID)));
% %     if ~isscalar(LUID.isRota(iID))||~isscalar(LUID.maxL(iID))||~isscalar(LUID.yID(iID)), error('��������'); end
% % end
% % end
        

    %% 
% %     function inputExchange()
% % % %         % 1 ���������ж�ת��+�����ж�ת��
% % % %         % d.LU.ID,d.LU.Weight: ������row������        
% % % %         % d.Veh.LWH: ������column������
% % % %         % d.Veh.BUFF: ������column������
% % % %         if ~isrow(d.LU.ID),              d.LU.ID=d.LU.ID';              end
% % % %         if ~isrow(d.LU.Weight),              d.LU.Weight=d.LU.Weight';              end
% % % %         if ~iscolumn(d.Veh.LWH),    d.Veh.LWH=d.Veh.LWH';    end
% % % %         if ~iscolumn(d.Veh.BUFF),   d.Veh.BUFF=d.Veh.BUFF';   end        
% % % %         % d.LU.BUFF: ������column������,�̶�ת��Ϊmatrix����,����Ϊ��������,����Ϊ3;
% % % %         if ~iscolumn(d.LU.BUFF),    d.LU.BUFF=d.LU.BUFF';    end
% % % %         % d.LU.LWH: ������matrix����,����Ϊ��������,����Ϊ3; δ����3��3������
% % % %         if size(d.LU.LWH,1)~=3,     d.LU.LWH=d.LU.LWH';     end
% %         
% %         % ������û�и���ֵ
% % %         if ~isfield(d.LU,'isRota')
% % %             d.LU.isRota = ones(size(d.LU.ID));
% % %         end
% %             
% %         % ��whichRotation��0��1ʱ���ж�, ��Ϊ2ʱ���򲻱�
% % %         if ParaArray.whichRotation == 1 % ȫ��������ת
% % %             d.LU.isRota = ones(size(d.LU.isRota));
% % %         elseif ParaArray.whichRotation == 0    % ȫ����ֹ��ת
% % %             d.LU.isRota = zeros(size(d.LU.isRota));
% % %         elseif ParaArray.whichRotation == 2    % ����������ת
% % %             1;
% % %         end
% % 
% %         % 2 Input���Ӽ�϶BUFF���feasible��LU��BIN�ĳ����ת��
% %         d.Veh.LWH = d.Veh.LWH - d.Veh.BUFF;
% %         d.LU.BUFF = [d.LU.BUFF; 0]; %�û�������BuffΪÿ���������ӵĳߴ�(�ܼ�϶)
% %         d.LU.BUFF = repmat(d.LU.BUFF,1,numel(d.LU.ID));
% %         d.LU.LWH = d.LU.LWH + d.LU.BUFF;
% %         %TODO: �˴����Ӽ�϶ΪȨ��֮�ʣ�δ����rotation��ı仯�����ڿ������㷨�����Ӽ�϶
% %         
% %         % 3 Ĭ�Ͻ�LUȫ������Horizontal������ת��ǰ�᣺��LU������ת��
% %         % NOTE: �˴������1: Horizontal�����LWH���Ƿ�Rotaed���
% %         % NOTE : ֱ���滻��ԭʼORIGINAL �� LWH
% %         [d.LU.Rotaed]= placeItemHori(d.LU.LWH,d.LU.isRota,1); %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
% %         d.LU.LWH = getRotaedLWH(d.LU.LWH, d.LU.Rotaed, d.LU.BUFF); 
% %         
% % end

    