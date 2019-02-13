function   [LU] = cpuLU(LU,Veh)
% cpuLU ==> 获取LU的isNonMixed和LU的isMixedTile

%% 1 获取最新的LU.isNonMixed和LU.isMixedTile

    %1.1 依据新的最大层数再重新LU的排序堆垛过程
    % 数据预处理：重点：获取LU.Rotaed,托盘是否排序
    LU.isNonMixed = ones(1,length(LU.Weight))*-1;     %1 LU可以形成满垛 0 必定有非满垛生成 LU排序依据之一(不会形成非满垛的LU提前获取堆垛,目的尽量减少混合堆垛的情形)
    LU.isMixedTile = zeros(1,length(LU.Weight));            %1 当isNonMixed=1时,都是满垛,isMixedTile都为0; 当isNonMixed=0时, 有非满垛,将余下的将称为非满垛的LU的isMixedTile赋值为1
    % GET LU.isNonMixed: 计算每个LU是否为不需要混拼的可能
    % 循环: ID个数
    tLU = getTableLU(LU);

             %  V1:   [tLUsorted,torder] = sortrows(tLU,{'SID','ID','LID','H'},{'ascend','ascend','ascend','descend'});
    [tLUsorted,torder] = sortrows(tLU,{'SID','EID','ID','LID','H'},{'ascend','ascend','ascend','ascend','descend'});  %V2: 增加EID顺序
            %     CateOrder = tLUsorted(:,{'SID','ID','LID','H'});  %增加maxHLayer?
            %     [GID,g] = findgroups(CateOrder)
    VehHeight = Veh.LWH(3,1);  % 车辆高度
    VehWidth = Veh.LWH(1,1);  % 车辆宽度

    %% 1.2  以堆垛标记ID区分LU的isNonMixed和isMixedTile
    uniID = unique(tLUsorted.ID);
    for iLu=1:length(uniID)
        % Item i 对于的LU flag标记
        flagLU = tLUsorted.ID == uniID(iLu);
        idxLU = find(flagLU);
        unimaxHLayerperLUID = unique(tLUsorted.maxHLayer(flagLU));

        % 1.1 如果对应ID的LU的高度H都一样,采用层数判断
        if isscalar(unimaxHLayerperLUID) %似乎即使单个用非单的也可以,但对于混合高度的也是可以用的; 但还是以层数为重心,特别是考虑平铺
            nb = sum(flagLU);
            nbmod = mod(nb,unimaxHLayerperLUID);       if nb ==0 || nbmod>nb, error('Gpreproc中计算isNonMixed错误'); end

            % 判断isNonMixed
            if nbmod == 0 %mod为0表明此ID的托盘LU可以完美的放到高度满垛,不会有有高度不满的Item出现. 不需要混合 不混合的提前在order中提前
                tLUsorted.isNonMixed(flagLU) = 1;
            else
                tLUsorted.isNonMixed(flagLU)= 0;
                % 更新tLUsorted.isMixedTile
                idxMod = idxLU(end-nbmod+1:end);        % 将非满垛的余托的isMixedTile赋值为1.%因为已经排序,所以可以按这个替换
                tLUsorted.isMixedTile(idxMod)=1;
            end

            % 1.2 如果对应ID的LU的高度H不一样,采用非层数(高度)判断    (一般结合ISisGpreprocLU1堆垛高度均衡使用)
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

            % 判断isNonMixed
            if ( VehHeight - tmpsum ) > min(tLUsorted.H)  %如高度间隙 > 最小LU高度
                % ???? 如果该LU的宽度可以放2层以上, 且目前能生成的也是2个以上,把最后几个进行高度均衡设置,
                maxWLayerLUID = unique(tLUsorted.maxL(flagLU,1));   if ~isscalar(maxWLayerLUID), error('Gpreproc中计算maxWLayerLUID错误'); end
                if maxWLayerLUID>1 && max(idxVeh) >1 && mod(max(idxVeh),maxWLayerLUID)~=1 %不能余1个Item单独放1层
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

    % 按顺序返回
    [~,sorder] = sort(torder);
    LU.isNonMixed = tLUsorted.isNonMixed(sorder)';
    LU.isMixedTile = tLUsorted.isMixedTile(sorder)';
    
    %% 2 计算LU的其它指
    
end