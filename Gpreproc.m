%% GPREPROC 重要函数:输入数据da转换+输入数据da核对
%% Form
%    [d] = Gpreproc(d)
%% Descripttion
%    输入数据预处理
%
%% Inputs
%    d                      (1,1)    结构体：用户输入数据
%
%% Outputs
%    d                     (1,1)
%

function [LU,Veh] = Gpreproc(LU,Veh,pwhichSortItemOrder)
    % 1 ID 转换为从1开始的类序号 方便刘工输入ID类信息    
    if isfield(LU, 'ID'),  [LU.ID,LU.OID] = idExchange(LU.ID); end 
    if isfield(LU, 'PID'),  [LU.PID,LU.OPID] = idExchange(LU.PID); end
    if isfield(LU, 'SID'),  [LU.SID,LU.OSID] = idExchange(LU.SID); end
%     if isfield(LU, 'UID'),  LU.UID = idExchange(LU.UID); end
    
    % 2 Input增加间隙BUFF后的feasible的LU和BIN的长宽高转换
    %     Veh.LWH = Veh.LWH - Veh.buff;  %Veh的buff无用，输入直接给可用内径长宽高
                %         LU.margin([2],:) = 2 %测试增加不同margin的效果
                %         LU.margin([3],:) = 3
    LU.LWH(1,:) =  LU.LWH(1,:) +  LU.margin(1,: ) + LU.margin(2,: ); %宽度（左右）
    LU.LWH(2,:) =  LU.LWH(2,:) +  LU.margin(3,: ) + LU.margin(4,: ); %长度（上下）
    
    % V1 : buff: 托盘间的间隙
    %     LU.buff = [LU.buff; 0]; %用户给定的Buff为每个托盘增加的尺寸(总间隙)
    %     LU.buff = repmat(LU.buff,1,numel(LU.ID));
%     LU.LWH = LU.LWH + LU.buff;
    %TODO: 此处增加间隙为权宜之际，未考虑rotation后的变化；后期考虑在算法中增加间隙

    % 3 默认将LU全部采用Horizontal方向旋转（前提：该LU允许旋转）
    % NOTE: 此处将获得1: Horizontal方向的LWH和是否Rotaed标记
    % NOTE : 直接替换了原始ORIGINAL 的 LWH
%        LU.LWH

    if pwhichSortItemOrder ==1
        [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,1); %第二个参数：1: Hori; 0: Vert；其它: 原封不动        
    elseif pwhichSortItemOrder ==2
        [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,0); %第二个参数：1: Hori; 0: Vert；其它: 原封不动
    elseif pwhichSortItemOrder ==3 %默认此选项
        [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,Veh.LWH(1,1)); %第二个参数：  3按VEH车辆左右摆放的缝隙最小排序
    end
    LU.LWH = getRotaedLWH(LU.LWH, LU.Rotaed, LU.margin);
        
    % 4 Veh从体积大->小   默认顺序
    Veh.Volume = prod(Veh.LWH);
    [~,order] = sortrows(Veh.Volume', [1],{'descend'});    
    Veh = structfun(@(x) x(:,order),Veh,'UniformOutput',false);
    if ~isrow(order),
        order = order';
    end
    Veh.order = order;
    
    % 5 计算LU在当前车型下的最大长宽高层数 TODO 考虑margin
    for i=1:length(LU.Weight)
        LU.maxL(1,i) =  floor(Veh.LWH(1,1)/LU.LWH(1,i));
        LU.maxL(2,i) =  floor(Veh.LWH(2,1)/LU.LWH(2,i));
        LU.maxL(3,i) =  floor(Veh.LWH(3,1)/LU.LWH(3,i));   %具体每个托盘LU的高度的最大层数
    end

    % 6 计算LU同样ID/可堆垛ID下的个数
    % V1"
    for i=1:length(LU.Weight)
        LU.nbID(i) = sum(LU.ID == LU.ID(i));
    end
    % V2: 直接计算
    LU.nbID = sum(LU.ID==LU.ID');
    LU.nbLID = sum(LU.LID==LU.LID');

    
    % 7 计算LU下面的isNonMixed/isMixedTile是否为不需要混拼/混拼排序
     LU.isNonMixed = ones(1,length(LU.Weight))*-1;    %Item是否非需要混合判定,将偶数个的Item提前进行Strip生成
     LU.isMixedTile = zeros(1,length(LU.Weight));    %Item混合,找出奇数个混合Item的尾托赋值为1
    % GET LU.isNonMixed: 计算每个LU是否为不需要混拼的可能
    % 循环: LID个数
    uniLID = unique(LU.LID)
    for iLu=1:length(uniLID)
        % Item i 对于的LU flag标记
        VehHeight = Veh.LWH(3,1);  % 车辆宽度
        
        flagLU = LU.LID(:) == uniLID(iLu);
        LUHeight = unique(LU.LWH(3,flagLU));     % 同样LULID时的 LU高度
            if length(unique(LUHeight)) > 2, warning('同样LULID时的LU高度值>2,非预期错误'); end %但对于随机生成的可能有这个错误
        if length(unique(LUHeight)) >= 2, LUHeight = max(LUHeight); end % 用LU最大值作为堆垛判断值

        maxHeightLayer= floor(VehHeight/LUHeight); %LU高度层数
        
        nb = sum(flagLU);
        nbmod = mod(nb,maxHeightLayer);
            if nb ==0 || nbmod>nb, error('Gpreproc中计算isNonMixed错误'); end
            
        if nbmod == 0 %mod为0表明 不需要混合 不混合的提前在order中提前
            LU.isNonMixed(flagLU) = 1;
        else
            LU.isNonMixed(flagLU)= 0;            
            % 计算LU的isMixedTile
            tmpSort=[LU.SID; LU.LWH(3,:); LU.PID; ]; %LU.Weight LU.maxL;            
            [~, order]=sortrows(tmpSort(:,flagLU)', [1,2,3],{'ascend','ascend','ascend'});
            flagLUIdx = find(flagLU);
            flagmodIdx = flagLUIdx(order(1:nbmod));
            LU.isMixedTile(flagmodIdx)=1;
        end
    end

    
%     LUID = getLUIDArray(LU); %% 计算：LU类型相关数据 暂时无用

%      printInput();
%% 嵌套函数:
function printInput()
    fprintf('本算例只有一个箱型 宽=%1.0f 长=%1.0f 高=%1.0f \n', unique(Veh.LWH','rows')');
    fprintf('本算例只有一个箱型 宽间隙=%1.0f 长间隙=%1.0f 高间隙=%1.0f \n', unique(Veh.buff','rows')');
    fprintf('本算例有 %d 个物品,其宽长高分别为 \n',numel(LU.ID(:)));
    fprintf('%1.1f %1.1f %1.1f \n',LU.LWH);
end

end


function [exID,ID] = idExchange(ID)
        uniID = unique(ID);
        exID=ID; %中间变量
        for i=1:length(uniID)
            exID(ID(:)==uniID(i)) = i;
        end
end



%%  获取LUID类型相关数据(同类型ID的体积，重量，是否可旋转)
function LUID = getLUIDArray(LU)
LUID.ID = unique(LU.ID);
for iID = 1:length(LUID.ID)
    LUID.Weight(iID) = sum(LU.Weight .* (LU.ID == LUID.ID(iID)) );
    LUID.Volume(iID) = sum(prod(LU.LWH) .* (LU.ID == LUID.ID(iID)) );
    LUID.isRota(iID) =  unique(LU.isRota(LU.ID == LUID.ID(iID)));
    LUID.maxL(iID) =  unique(LU.maxL(LU.ID == LUID.ID(iID)));
    LUID.yID(iID) =  unique(LU.yID(LU.ID == LUID.ID(iID)));
    if ~isscalar(LUID.isRota(iID))||~isscalar(LUID.maxL(iID))||~isscalar(LUID.yID(iID)), error('致命错误'); end
end
end
        

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

    