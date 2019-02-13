function   [LU] = cpuLU(LU,Veh)
% cpuLU ==> ��ȡLU��isNonMixed��LU��isMixedTile

%% 1 ��ȡ���µ�LU.isNonMixed��LU.isMixedTile

    %1.1 �����µ�������������LU������Ѷ����
    % ����Ԥ�����ص㣺��ȡLU.Rotaed,�����Ƿ�����
    LU.isNonMixed = ones(1,length(LU.Weight))*-1;     %1 LU�����γ����� 0 �ض��з��������� LU��������֮һ(�����γɷ������LU��ǰ��ȡ�Ѷ�,Ŀ�ľ������ٻ�϶Ѷ������)
    LU.isMixedTile = zeros(1,length(LU.Weight));            %1 ��isNonMixed=1ʱ,��������,isMixedTile��Ϊ0; ��isNonMixed=0ʱ, �з�����,�����µĽ���Ϊ�������LU��isMixedTile��ֵΪ1
    % GET LU.isNonMixed: ����ÿ��LU�Ƿ�Ϊ����Ҫ��ƴ�Ŀ���
    % ѭ��: ID����
    tLU = getTableLU(LU);

             %  V1:   [tLUsorted,torder] = sortrows(tLU,{'SID','ID','LID','H'},{'ascend','ascend','ascend','descend'});
    [tLUsorted,torder] = sortrows(tLU,{'SID','EID','ID','LID','H'},{'ascend','ascend','ascend','ascend','descend'});  %V2: ����EID˳��
            %     CateOrder = tLUsorted(:,{'SID','ID','LID','H'});  %����maxHLayer?
            %     [GID,g] = findgroups(CateOrder)
    VehHeight = Veh.LWH(3,1);  % �����߶�
    VehWidth = Veh.LWH(1,1);  % �������

    %% 1.2  �ԶѶ���ID����LU��isNonMixed��isMixedTile
    uniID = unique(tLUsorted.ID);
    for iLu=1:length(uniID)
        % Item i ���ڵ�LU flag���
        flagLU = tLUsorted.ID == uniID(iLu);
        idxLU = find(flagLU);
        unimaxHLayerperLUID = unique(tLUsorted.maxHLayer(flagLU));

        % 1.1 �����ӦID��LU�ĸ߶�H��һ��,���ò����ж�
        if isscalar(unimaxHLayerperLUID) %�ƺ���ʹ�����÷ǵ���Ҳ����,�����ڻ�ϸ߶ȵ�Ҳ�ǿ����õ�; �������Բ���Ϊ����,�ر��ǿ���ƽ��
            nb = sum(flagLU);
            nbmod = mod(nb,unimaxHLayerperLUID);       if nb ==0 || nbmod>nb, error('Gpreproc�м���isNonMixed����'); end

            % �ж�isNonMixed
            if nbmod == 0 %modΪ0������ID������LU���������ķŵ��߶�����,�������и߶Ȳ�����Item����. ����Ҫ��� ����ϵ���ǰ��order����ǰ
                tLUsorted.isNonMixed(flagLU) = 1;
            else
                tLUsorted.isNonMixed(flagLU)= 0;
                % ����tLUsorted.isMixedTile
                idxMod = idxLU(end-nbmod+1:end);        % ������������е�isMixedTile��ֵΪ1.%��Ϊ�Ѿ�����,���Կ��԰�����滻
                tLUsorted.isMixedTile(idxMod)=1;
            end

            % 1.2 �����ӦID��LU�ĸ߶�H��һ��,���÷ǲ���(�߶�)�ж�    (һ����ISisGpreprocLU1�Ѷ�߶Ⱦ���ʹ��)
        else
            nbLU = length(idxLU);
            idxVeh = zeros(nbLU,1);
            tmpsum = 0;
            tmpidx = 1;

            for idx=1:nbLU
                tmpsum = tmpsum + tLUsorted.H(idxLU(idx));
                if tmpsum > VehHeight
                    tmpsum = tLUsorted.H(idxLU(idx));
                    tmpidx=tmpidx+1;
                end
                idxVeh(idx) =tmpidx;
            end

            % �ж�isNonMixed
            if ( VehHeight - tmpsum ) > min(tLUsorted.H)  %��߶ȼ�϶ > ��СLU�߶�
                % ???? �����LU�Ŀ�ȿ��Է�2������, ��Ŀǰ�����ɵ�Ҳ��2������,����󼸸����и߶Ⱦ�������,
                maxWLayerLUID = unique(tLUsorted.maxL(flagLU,1));   if ~isscalar(maxWLayerLUID), error('Gpreproc�м���maxWLayerLUID����'); end
                if maxWLayerLUID>1 && max(idxVeh) >1 && mod(max(idxVeh),maxWLayerLUID)~=1 %������1��Item������1��
                    tLUsorted.isNonMixed(flagLU)= 1;
                else
                    tLUsorted.isNonMixed(flagLU)= 0;
                    idxMod = idxLU(idxVeh==max(idxVeh));
                    tLUsorted.isMixedTile(idxMod)=1;
                end
            else
                tLUsorted.isNonMixed(flagLU)= 1;
            end

        end

    end

    % ��˳�򷵻�
    [~,sorder] = sort(torder);
    LU.isNonMixed = tLUsorted.isNonMixed(sorder)';
    LU.isMixedTile = tLUsorted.isMixedTile(sorder)';
    
    %% 2 ����LU������ָ
    
end