function [LU,Veh] = cpuLUVeh(LU,Veh)
% cpuLUVeh 重要函数:
% 更新 Veh 的Volume/order
% 更新 LU 的ID/LU的LWH/LU的maxL/LU的maxHLayer/LU的nbID/LU的isNonMixed/isMixedTile
        
     % global  ISisGpreprocLU1

    %% 1 更新车辆Veh顺序 增加order
    % Veh从体积大->小   默认顺序
    Veh.Volume = prod(Veh.LWH);
    [~,order] = sortrows(Veh.Volume', [1],{'descend'});    
    Veh = structfun(@(x) x(:,order),Veh,'UniformOutput',false);
    if ~isrow(order)
        order = order';
    end
    Veh.order = order;
    
     %% 2 更新LU：的ID号（保留Old）
     % 1 如果相同ID号下，对应SID号要必须不同；如相同改变ID号，直到不存在相同的ID在不同SID内;
     % TODO : 增加不同SID/不同EID下的ID号不可重复
     if isrepeated(LU.ID,LU.SID)
        warning('存在托盘ID号在不同 SID下重复, 需要更正'); 
        LU.ID = reviseID(LU.ID,LU.SID);
     end

     if isrepeated(LU.EID,LU.SID)
         warning('存在托盘EP LOCATION ID号在不同SID下重复, 需要更正');
         LU.EID = reviseID(LU.EID,LU.SID);
     end
     
     if isrepeated(LU.PID,LU.SID)
         warning('存在托盘PID ID号在不同SID下重复, 需要更正');
         LU.PID = reviseID(LU.PID,LU.SID);
     end
     
     % 是否需要？fixme 
      if isrepeated(LU.LID,LU.SID)
         warning('存在托盘LID ID号在不同SID下重复, 需要更正');
         LU.LID = reviseID(LU.LID,LU.SID);
      end    
     
     if isrepeated(LU.ID,LU.EID)
         error('存在托盘ID号在不同 EID下重复, 需要更正');
     end
     
    % 1 ID 转换为从1开始的类序号 方便刘工输入ID类信息    TODO 增加LID的OLID？
    if isfield(LU, 'ID'),    [LU.ID,LU.OID]      = idExchange(LU.ID);    end
    if isfield(LU, 'PID'),  [LU.PID,LU.OPID]  = idExchange(LU.PID); end
    if isfield(LU, 'EID'),  [LU.EID,LU.OEID]  = idExchange(LU.EID);  end
    if isfield(LU, 'SID'),  [LU.SID,LU.OSID]  = idExchange(LU.SID);  end
    
    % 1 如果相同ID（PID/EID/LID等）号下，对应SID号要必须不同；如相同改变ID号，直到不存在相同的ID在不同SID内;
    if isrepeated(LU.ID,LU.SID),  error('存在ID号重复, 应该在check时已调整');   end
    
    %% 3 更新LU：含margin的LU的LWH；更新LU的Rotaed标记
% %     % GET LU's LWH with margin and Rotaed or not
% %     % 2 Input增加间隙BUFF后的feasible的LU和BIN的长宽高转换
% %     LU.OLWH = LU.LWH;
% %     LU.LWH(1,:) =  LU.LWH(1,:) +  LU.margin(1,: ) + LU.margin(2,: ); %宽度（左右）
% %     LU.LWH(2,:) =  LU.LWH(2,:) +  LU.margin(3,: ) + LU.margin(4,: ); %长度（上下）
% %     
% %     [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,Veh.LWH(1,1));  %第二个参数：  3按VEH车辆左右摆放的缝隙最小排序
% %     LU.LWH = getRotaedLWH(LU.LWH, LU.Rotaed, LU.margin);
% %     
% %     % 3 默认将LU全部采用Horizontal方向旋转（前提：该LU允许旋转）
% %     % NOTE: 此处将获得1: Horizontal方向的LWH和是否Rotaed标记
% %     % NOTE : 直接替换了原始ORIGINAL 的 LWH
% %     % LU.LWH
% % 
% % %     if pwhichSortItemOrder ==1
% % %         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,1); %第二个参数：1: Hori; 0: Vert；其它: 原封不动        
% % %     elseif pwhichSortItemOrder ==2
% % %         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,0); %第二个参数：1: Hori; 0: Vert；其它: 原封不动
% % %     elseif pwhichSortItemOrder ==3 %默认此选项
% % %         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,Veh.LWH(1,1)); %第二个参数：  3按VEH车辆左右摆放的缝隙最小排序
% % %     end
        
    %% 4 计算LU：在当前车型下的最大长宽高层数 包含margin
    if ~isfield(LU, 'maxL') 
        LU.maxL(1,:) =  floor(Veh.LWH(1,1)./LU.LWH(1,:));
        LU.maxL(2,:) =  floor(Veh.LWH(2,1)./LU.LWH(2,:));
        LU.maxL(3,:) =  floor(Veh.LWH(3,1)./LU.LWH(3,:));   %具体每个托盘LU的高度的最大层数
    end

    if ~isfield(LU, 'maxHLayer') % 如不存在,直接以maxL指示; 如存在,以最小值取代
        LU.maxHLayer = LU.maxL(3,:); 
    else
        LU.maxHLayer = min( [LU.maxHLayer',LU.maxL(3,:)'], [], 2 )';  % 同步更新对应LU的maxHLayer    
    end
    
    if ~isfield(LU, 'Rotaed')
        LU.Rotaed = zeros(size(LU.ID));
    end

    %% 5 计算LU同样ID/可堆垛ID下的个数
    % V1"
    %     for i=1:length(LU.Weight)
    %         LU.nbID(i) = sum(LU.ID == LU.ID(i));
    %     end
    % V2: 直接计算
    LU.nbID = sum(LU.ID==LU.ID');
    LU.nbLID = sum(LU.LID==LU.LID');
    
    %% 6：555 ITEMBALANCE必须 7 V2: 计算LU下面的isNonMixed/isMixedTile是否为不需要  层面 混拼/混拼排序 
    % 555 区分为LU相同高度（依据层数判断） 或 LU多种高度（依据高度判断 TODO 增加层数判断）两种类型
    [LU] = cpuLU(LU,Veh);
    
  end

function [exID,ID] = idExchange(ID)
        uniID = unique(ID);
        exID=ID; %中间变量
        for i=1:length(uniID)
            exID(ID(:)==uniID(i)) = i;
        end
end









%%
% % %      LU.isNonMixed = ones(1,length(LU.Weight))*-1;     %1 LU可以形成满垛 0 必定有非满垛生成 LU排序依据之一(不会形成非满垛的LU提前获取堆垛,目的尽量减少混合堆垛的情形)
% % %      LU.isMixedTile = zeros(1,length(LU.Weight));          %1 当isNonMixed=1时,都是满垛,isMixedTile都为0; 当isNonMixed=0时, 有非满垛,将余下的将称为非满垛的LU的isMixedTile赋值为1
% % %     % GET LU.isNonMixed: 计算每个LU是否为不需要混拼的可能
% % %     % 循环: ID个数
% % %     tLU = getTableLU(LU);
% % %     %     [tLUsorted,torder] = sortrows(tLU,{'SID','ID','LID','H'},{'ascend','ascend','ascend','descend'});
% % %     [tLUsorted,torder] = sortrows(tLU,{'SID','EID','ID','LID','H'},{'ascend','ascend','ascend','ascend','descend'});  %V2: 增加EID顺序
% % %                         %     CateOrder = tLUsorted(:,{'SID','ID','LID','H'});  %增加maxHLayer?
% % %                         %     [GID,g] = findgroups(CateOrder)
% % %     VehHeight = Veh.LWH(3,1);  % 车辆高度
% % %     VehWidth = Veh.LWH(1,1);  % 车辆宽度
% % %     
% % %     %% 以堆垛标记ID区分LU的isNonMixed和isMixedTile
% % %     uniID = unique(tLUsorted.ID);
% % %     for iLu=1:length(uniID)
% % %         % Item i 对于的LU flag标记
% % %         flagLU = tLUsorted.ID == uniID(iLu);
% % %         idxLU = find(flagLU);
% % %         unimaxHLayerperLUID = unique(tLUsorted.maxHLayer(flagLU));
% % %         
% % %         % 1.1 如果对应ID的LU的高度H都一样,采用层数判断
% % %         if isscalar(unimaxHLayerperLUID) %似乎即使单个用非单的也可以,但对于混合高度的也是可以用的; 但还是以层数为重心,特别是考虑平铺
% % %             nb = sum(flagLU);
% % %             nbmod = mod(nb,unimaxHLayerperLUID);       if nb ==0 || nbmod>nb, error('Gpreproc中计算isNonMixed错误'); end    
% % %             
% % %             % 判断isNonMixed
% % %             if nbmod == 0 %mod为0表明此ID的托盘LU可以完美的放到高度满垛,不会有有高度不满的Item出现. 不需要混合 不混合的提前在order中提前
% % %                 tLUsorted.isNonMixed(flagLU) = 1;
% % %             else
% % %                 tLUsorted.isNonMixed(flagLU)= 0;
% % %                 % 更新tLUsorted.isMixedTile
% % %                 idxMod = idxLU(end-nbmod+1:end);        % 将非满垛的余托的isMixedTile赋值为1.%因为已经排序,所以可以按这个替换        
% % %                 tLUsorted.isMixedTile(idxMod)=1;
% % %             end
% % %             
% % %         % 1.2 如果对应ID的LU的高度H不一样,采用非层数(高度)判断    (一般结合ISisGpreprocLU1堆垛高度均衡使用)
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
% % %             % 判断isNonMixed
% % %             if ( VehHeight - tmpsum ) > min(tLUsorted.H)  %如高度间隙 > 最小LU高度
% % %                 % ???? 如果该LU的宽度可以放2层以上, 且目前能生成的也是2个以上,把最后几个进行高度均衡设置,
% % %                 maxWLayerLUID = unique(tLUsorted.maxL(flagLU,1));   if ~isscalar(maxWLayerLUID), error('Gpreproc中计算maxWLayerLUID错误'); end
% % %                 if ISisGpreprocLU1==1 && maxWLayerLUID>1 && max(idxVeh) >1 && mod(max(idxVeh),maxWLayerLUID)~=1 %不能余1个Item单独放1层
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
% % %     % 按顺序返回
% % %     [~,sorder] = sort(torder);
% % %     LU.isNonMixed = tLUsorted.isNonMixed(sorder)';
% % %     LU.isMixedTile = tLUsorted.isMixedTile(sorder)';
    
    
    %% 7 V1: BUG MAYBE 计算LU下面的isNonMixed/isMixedTile是否为不需要  层面 混拼/混拼排序
% % % %      LU.isNonMixed = ones(1,length(LU.Weight))*-1;    %Item是否非需要混合判定,将偶数个的Item提前进行Strip生成
% % % %      LU.isMixedTile = zeros(1,length(LU.Weight));          %Item混合,找出奇数个混合Item的尾托赋值为1
% % % %     % GET LU.isNonMixed: 计算每个LU是否为不需要混拼的可能
% % % %     % 循环: LID个数
% % % %     uniLID = unique(LU.LID);
% % % % %     uniLID = unique(LU.ID);
% % % %     for iLu=1:length(uniLID)
% % % %         % Item i 对于的LU flag标记
% % % %         VehHeight = Veh.LWH(3,1);  % 车辆高度
% % % %         
% % % %         flagLU = LU.LID(:) == uniLID(iLu);
% % % % %         flagLU = LU.ID(:) == uniLID(iLu);
% % % %         LUHeight = unique(LU.LWH(3,flagLU));     % 同样LULID时的 LU高度
% % % %             if length(unique(LUHeight)) > 2, warning('同样LULID时的LU高度值>2,非预期错误'); end %但对于随机生成的可能有这个错误
% % % %         if length(unique(LUHeight)) >= 2, LUHeight = max(LUHeight); end % 用LU最大值作为堆垛判断值
% % % % 
% % % %         LU.maxL(3,flagLU)
% % % %         LU.maxHLayer(flagLU)
% % % %         maxHeightLayer= floor(VehHeight/LUHeight); %LU高度层数
% % % %         
% % % %         nb = sum(flagLU);
% % % %         nbmod = mod(nb,maxHeightLayer);
% % % %             if nb ==0 || nbmod>nb, error('Gpreproc中计算isNonMixed错误'); end
% % % %             
% % % %         if nbmod == 0 %mod为0表明 不需要混合 不混合的提前在order中提前
% % % %             LU.isNonMixed(flagLU) = 1;
% % % %         else
% % % %             LU.isNonMixed(flagLU)= 0;            
% % % %             % 计算LU的isMixedTile
% % % %             tmpSort=[LU.SID; LU.LWH(3,:); LU.PID; ]; %LU.Weight LU.maxL;            
% % % %             [~, order]=sortrows(tmpSort(:,flagLU)', [1,2,3],{'ascend','ascend','ascend'});
% % % %             flagLUIdx = find(flagLU);
% % % %             flagmodIdx = flagLUIdx(order(1:nbmod));
% % % %             LU.isMixedTile(flagmodIdx)=1;
% % % %         end
% % % %     end
    
    
%     LUID = getLUIDArray(LU); %% 计算：LU类型相关数据 暂时无用
%      printInput();
%% 嵌套函数:
% function printInput()
%     fprintf('本算例只有一个箱型 宽=%1.0f 长=%1.0f 高=%1.0f \n', unique(Veh.LWH','rows')');
%     fprintf('本算例只有一个箱型 宽间隙=%1.0f 长间隙=%1.0f 高间隙=%1.0f \n', unique(Veh.buff','rows')');
%     fprintf('本算例有 %d 个物品,其宽长高分别为 \n',numel(LU.ID(:)));
%     fprintf('%1.1f %1.1f %1.1f \n',LU.LWH);
% end





% % %%  获取LUID类型相关数据(同类型ID的体积，重量，是否可旋转)
% % function LUID = getLUIDArray(LU)
% % LUID.ID = unique(LU.ID);
% % for iID = 1:length(LUID.ID)
% %     LUID.Weight(iID) = sum(LU.Weight .* (LU.ID == LUID.ID(iID)) );
% %     LUID.Volume(iID) = sum(prod(LU.LWH) .* (LU.ID == LUID.ID(iID)) );
% %     LUID.isRota(iID) =  unique(LU.isRota(LU.ID == LUID.ID(iID)));
% %     LUID.maxL(iID) =  unique(LU.maxL(LU.ID == LUID.ID(iID)));
% %     LUID.yID(iID) =  unique(LU.yID(LU.ID == LUID.ID(iID)));
% %     if ~isscalar(LUID.isRota(iID))||~isscalar(LUID.maxL(iID))||~isscalar(LUID.yID(iID)), error('致命错误'); end
% % end
% % end
        

    %% 
% %     function inputExchange()
% % % %         % 1 行列向量判断转换+矩阵判断转换
% % % %         % d.LU.ID,d.LU.Weight: 必须是row行向量        
% % % %         % d.Veh.LWH: 必须是column列向量
% % % %         % d.Veh.BUFF: 必须是column列向量
% % % %         if ~isrow(d.LU.ID),              d.LU.ID=d.LU.ID';              end
% % % %         if ~isrow(d.LU.Weight),              d.LU.Weight=d.LU.Weight';              end
% % % %         if ~iscolumn(d.Veh.LWH),    d.Veh.LWH=d.Veh.LWH';    end
% % % %         if ~iscolumn(d.Veh.BUFF),   d.Veh.BUFF=d.Veh.BUFF';   end        
% % % %         % d.LU.BUFF: 必须是column列向量,继而转换为matrix矩阵,列数为托盘数量,行数为3;
% % % %         if ~iscolumn(d.LU.BUFF),    d.LU.BUFF=d.LU.BUFF';    end
% % % %         % d.LU.LWH: 必须是matrix矩阵,列数为托盘数量,行数为3; 未考虑3行3列情形
% % % %         if size(d.LU.LWH,1)~=3,     d.LU.LWH=d.LU.LWH';     end
% %         
% %         % 如输入没有该数值
% % %         if ~isfield(d.LU,'isRota')
% % %             d.LU.isRota = ones(size(d.LU.ID));
% % %         end
% %             
% %         % 对whichRotation在0和1时的判断, 如为2时，则不变
% % %         if ParaArray.whichRotation == 1 % 全部允许旋转
% % %             d.LU.isRota = ones(size(d.LU.isRota));
% % %         elseif ParaArray.whichRotation == 0    % 全部禁止旋转
% % %             d.LU.isRota = zeros(size(d.LU.isRota));
% % %         elseif ParaArray.whichRotation == 2    % 部分允许旋转
% % %             1;
% % %         end
% % 
% %         % 2 Input增加间隙BUFF后的feasible的LU和BIN的长宽高转换
% %         d.Veh.LWH = d.Veh.LWH - d.Veh.BUFF;
% %         d.LU.BUFF = [d.LU.BUFF; 0]; %用户给定的Buff为每个托盘增加的尺寸(总间隙)
% %         d.LU.BUFF = repmat(d.LU.BUFF,1,numel(d.LU.ID));
% %         d.LU.LWH = d.LU.LWH + d.LU.BUFF;
% %         %TODO: 此处增加间隙为权宜之际，未考虑rotation后的变化；后期考虑在算法中增加间隙
% %         
% %         % 3 默认将LU全部采用Horizontal方向旋转（前提：该LU允许旋转）
% %         % NOTE: 此处将获得1: Horizontal方向的LWH和是否Rotaed标记
% %         % NOTE : 直接替换了原始ORIGINAL 的 LWH
% %         [d.LU.Rotaed]= placeItemHori(d.LU.LWH,d.LU.isRota,1); %第二个参数：1: Hori; 0: Vert；其它: 原封不动
% %         d.LU.LWH = getRotaedLWH(d.LU.LWH, d.LU.Rotaed, d.LU.BUFF); 
% %         
% % end

    