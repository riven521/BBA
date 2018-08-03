%% GCHECKINPUT 重要函数:输入数据da转换+输入数据da核对
%% Form
%    [d] = GcheckInput(d,ParaArray)
%% Descripttion
%    输入数据检验
%
%% Inputs
%    d                      (1,1)    结构体：用户输入数据
%    ParaArray          (1,1)    结构体：算法测试参数
%
%% Outputs
%    d                     (1,1)
%

function [d] = GcheckInput(d,ParaArray)
%  行数:长宽高(row);  列数:托盘数量(coloum);

initCheck();
inputExchange(); %333
inputLUIDExchange(); %555
deepCheck(); %333
getLUIDArray(); %% 计算：LU类型相关数据
printInput();

    %% 嵌套函数
    %% initCheck 输入数据初步检验    
    function initCheck()
        % 初步判断输入: d.Veh
        if numel(d.Veh.LWH(:))~=numel(d.Veh.BUFF(:)) || numel(d.Veh.LWH(:)) ~= 3 ||...
                ~isscalar(d.Veh.Weight)
            error('本算例da.BinArray错误,超出预期 \n');
        end
        if any(d.Veh.LWH(:)<=0) || any(d.Veh.Weight(:)<=0) || any(d.Veh.BUFF(:)<0)
            error('本算例da.BinArray存在异常输入,请核验数据 \n');
        end
        % 初步判断输入: d.LU
        if numel(d.LU.ID(:))~=numel(d.LU.Weight(:)) || numel(d.LU.BUFF(:)) ~= 2 || ...
                numel(d.LU.ID(:))*3~=numel(d.LU.LWH(:)) || ~ismatrix(d.LU.LWH)
                        %             numel(d.LU.ID(:))~=numel(d.LU.Weight(:)) 
                        %             numel(d.LU.BUFF(:)) ~= 2 
                        %             numel(d.LU.ID(:))*3~=numel(d.LU.LWH(:)) 
                        %             ~ismatrix(d.LU.LWH)           
            error('本算例da.LUArray错误,超出预期 \n');
        end
        if any(d.LU.LWH(:)<=0) || any(d.LU.Weight(:)<=0)  || any(d.LU.BUFF(:)<0) || ...
                any(d.LU.isRota(:)<0) || any(d.LU.isRota(:)>1)
            error('本算例da.BinArray存在异常输入,请核验数据 \n');
        end
    end

    %% Input输入转换
    function inputExchange()
        % 1 行列向量判断转换+矩阵判断转换
        % d.LU.ID,d.LU.Weight: 必须是row行向量        
        % d.Veh.LWH: 必须是column列向量
        % d.Veh.BUFF: 必须是column列向量
        if ~isrow(d.LU.ID),              d.LU.ID=d.LU.ID';              end
        if ~isrow(d.LU.Weight),              d.LU.Weight=d.LU.Weight';              end
        if ~iscolumn(d.Veh.LWH),    d.Veh.LWH=d.Veh.LWH';    end
        if ~iscolumn(d.Veh.BUFF),   d.Veh.BUFF=d.Veh.BUFF';   end        
        % d.LU.BUFF: 必须是column列向量,继而转换为matrix矩阵,列数为托盘数量,行数为3;
        if ~iscolumn(d.LU.BUFF),    d.LU.BUFF=d.LU.BUFF';    end
        % d.LU.LWH: 必须是matrix矩阵,列数为托盘数量,行数为3; 未考虑3行3列情形
        if size(d.LU.LWH,1)~=3,     d.LU.LWH=d.LU.LWH';     end
        
        % 如输入没有该数值
        if ~isfield(d.LU,'isRota')
            d.LU.isRota = ones(size(d.LU.ID));
        end
            
        % 对whichRotation在0和1时的判断, 如为2时，则不变
        if ParaArray.whichRotation == 1 % 全部允许旋转
            d.LU.isRota = ones(size(d.LU.isRota));
        elseif ParaArray.whichRotation == 0    % 全部禁止旋转
            d.LU.isRota = zeros(size(d.LU.isRota));
        elseif ParaArray.whichRotation == 2    % 部分允许旋转
            1;
        end

        % 2 Input增加间隙BUFF后的feasible的LU和BIN的长宽高转换
        d.Veh.LWH = d.Veh.LWH - d.Veh.BUFF;
        d.LU.BUFF = [d.LU.BUFF; 0]; %用户给定的Buff为每个托盘增加的尺寸(总间隙)
        d.LU.BUFF = repmat(d.LU.BUFF,1,numel(d.LU.ID));
        d.LU.LWH = d.LU.LWH + d.LU.BUFF;
        %TODO: 此处增加间隙为权宜之际，未考虑rotation后的变化；后期考虑在算法中增加间隙
        
        % 3 默认将LU全部采用Horizontal方向旋转（前提：该LU允许旋转）
        % NOTE: 此处将获得1: Horizontal方向的LWH和是否Rotaed标记
        % NOTE : 直接替换了原始ORIGINAL 的 LWH
        [d.LU.Rotaed]= placeItemHori(d.LU.LWH,d.LU.isRota,1); %第二个参数：1: Hori; 0: Vert；其它: 原封不动
        d.LU.LWH = getRotaedLWH(d.LU.LWH, d.LU.Rotaed, d.LU.BUFF); 
        
        
        % 3 Input增加LUArray是否Rotation的判断flag TOBE DELE
%         d.LU.RotaFlag = zeros(1,numel(d.LU.ID)); TOBE DELE
    end

    %% LUID输入转换: d.LU.ID 转换为从1开始的类序号 方便刘工输入ID信息
    function inputLUIDExchange()
        uniLUID = unique(d.LU.ID);
        nLUid = size(uniLUID,2);
        tmpLUID=d.LU.ID; %中间变量
        for i=1:nLUid
            tmpLUID(d.LU.ID(:)==uniLUID(i)) = i;
        end
        d.LU.ID=tmpLUID;
    end

    %% 深入判断
    function deepCheck()        
        % 基础输入获取
        dBin = unique(d.Veh.LWH','rows')'; %获取非相同车型数据
         dLU = d.LU.LWH;
        dLUisRota = d.LU.isRota;
        dLUid = d.LU.ID;
        nLUid = size(unique(dLUid),2);
        % 高度约束 罕见
        if any(dBin(3) < dLU(3,:)),     error('错误: 存在托盘高度 大于 本车型可用高度数据 \n');     end
        
        % 长宽约束 罕见 (托盘长宽大于)
        flag1 = dLUisRota == 1;  %允许rota的LU标记
        flag0 = dLUisRota == 0;  %不允许rota的LU标记
        % 测试不允许旋转的LU
        if any(dBin(1) < dLU(1,flag0)),  error('错误: 存在托盘宽度 大于 本车型可用宽度数据 \n');    end
        if any(dBin(2) < dLU(2,flag0)),  error('错误: 存在托盘长度 大于 本车型可用长度数据 \n');    end
        % 测试允许旋转的LU
        if ~all(~flag1) && min(dBin(1:2)) < max(min(dLU(1:2, flag1)))  %flag1不是全0 且满足后面条件
            error('错误: 存在托盘最短边的最大值 大于 本车型长或宽 \n');
        end
        if ~all(~flag1) && max(dBin(1:2)) < max(max(dLU(1:2, flag1 ))) %flag1不是全0 且满足后面条件
            error('错误: 存在托盘最长边的最大值 大于 本车型长或宽 \n');
        end

                            % 长宽约束 罕见 (托盘长宽大于) TOBE DELE
                    %         if ParaArray.whichRotation == 1
                    %             if min(dBin(1:2)) < max(min(dLU(1:2,:))) || max(dBin(1:2)) < max(max(dLU(1:2,:)))
                    %                 error('错误: 存在托盘长或宽 大于 本车型长或宽 \n');
                    %             end
                    %         else %不准rotation
                    %             if any(dBin(1) < dLU(1,:)),  error('错误: 存在托盘宽度 大于 本车型可用宽度数据 \n');    end
                    %             if any(dBin(2) < dLU(2,:)),  error('错误: 存在托盘长度 大于 本车型可用长度数据 \n');    end
                    %         end

        % 托盘ID（材质，类型）对应长宽约束  常见555
        for iLU = 1:nLUid
            tmp = dLU(1:2,dLUid==iLU)';   %获取宽和长,不要高
            if numel(unique(tmp,'rows')) > 2,  error('错误: 存在托盘ID类型相同 但其长宽不同数据 \n');    end
        end
    end
    
    %%  获取LUID类型相关数据(同类型ID的体积，重量，是否可旋转)
    function getLUIDArray()
%         printstruct(d);
        d.LUID.ID = unique(d.LU.ID);
        nLUID = numel(d.LUID.ID);        

        for iID = 1:nLUID
        d.LUID.Weight(iID) = sum(d.LU.Weight .* (d.LU.ID == d.LUID.ID(iID)) );
        d.LUID.Volume(iID) = sum(prod(d.LU.LWH) .* (d.LU.ID == d.LUID.ID(iID)) );
        d.LUID.isRota(iID) =  unique(d.LU.isRota(d.LU.ID == d.LUID.ID(iID))); if ~isscalar(d.LUID.isRota(iID)), error('致命错误'); end
        end

    end

   %% 打印Input数据
    function printInput()
        fprintf('本算例只有一个箱型 宽=%1.0f 长=%1.0f 高=%1.0f \n', unique(d.Veh.LWH','rows')');
        fprintf('本算例只有一个箱型 宽间隙=%1.0f 长间隙=%1.0f 高间隙=%1.0f \n', unique(d.Veh.BUFF','rows')');
        fprintf('本算例有 %d 个物品,其宽长高分别为 \n',numel(d.LU.ID(:)));
        fprintf('%1.1f %1.1f %1.1f \n',d.LU.LWH);
    end

end





