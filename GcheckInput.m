%% GCHECKINPUT 重要函数:输入数据da转换+输入数据da核对
%% Form
%    [da] = GcheckInput(da,ParaArray)
%% Descripttion
%    输入数据检验
%
%% Inputs
%    da                      (1,1)    结构体：用户输入数据
%    ParaArray          (1,1)    结构体：算法测试参数
%
%% Outputs
%    da                     (1,1)
%

function [da] = GcheckInput(da,ParaArray)
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
        % 初步判断输入: da.BinArray
        if numel(da.BinArray.LWH(:))~=numel(da.BinArray.BUFF(:)) || numel(da.BinArray.LWH(:)) ~= 3 ||...
                ~isscalar(da.BinArray.Weight)
            error('本算例da.BinArray错误,超出预期 \n');
        end
        if any(da.BinArray.LWH(:)<=0) || any(da.BinArray.Weight(:)<=0) || any(da.BinArray.BUFF(:)<0)
            error('本算例da.BinArray存在异常输入,请核验数据 \n');
        end
        % 初步判断输入: da.LUArray
        if numel(da.LUArray.ID(:))~=numel(da.LUArray.Weight(:)) || numel(da.LUArray.BUFF(:)) ~= 2 || ...
                numel(da.LUArray.ID(:))*3~=numel(da.LUArray.LWH(:)) || ~ismatrix(da.LUArray.LWH)
                        %             numel(da.LUArray.ID(:))~=numel(da.LUArray.Weight(:)) 
                        %             numel(da.LUArray.BUFF(:)) ~= 2 
                        %             numel(da.LUArray.ID(:))*3~=numel(da.LUArray.LWH(:)) 
                        %             ~ismatrix(da.LUArray.LWH)           
            error('本算例da.LUArray错误,超出预期 \n');
        end
        if any(da.LUArray.LWH(:)<=0) || any(da.LUArray.Weight(:)<=0)  || any(da.LUArray.BUFF(:)<0) || ...
                any(da.LUArray.isRota(:)<0) || any(da.LUArray.isRota(:)>1)
            error('本算例da.BinArray存在异常输入,请核验数据 \n');
        end
    end

    %% Input输入转换
    function inputExchange()
        % 1 行列向量判断转换+矩阵判断转换
        % da.LUArray.ID,da.LUArray.Weight: 必须是row行向量        
        % da.BinArray.LWH: 必须是column列向量
        % da.BinArray.BUFF: 必须是column列向量
        if ~isrow(da.LUArray.ID),              da.LUArray.ID=da.LUArray.ID';              end
        if ~isrow(da.LUArray.Weight),              da.LUArray.Weight=da.LUArray.Weight';              end
        if ~iscolumn(da.BinArray.LWH),    da.BinArray.LWH=da.BinArray.LWH';    end
        if ~iscolumn(da.BinArray.BUFF),   da.BinArray.BUFF=da.BinArray.BUFF';   end        
        % da.LUArray.BUFF: 必须是column列向量,继而转换为matrix矩阵,列数为托盘数量,行数为3;
        if ~iscolumn(da.LUArray.BUFF),    da.LUArray.BUFF=da.LUArray.BUFF';    end
        % da.LUArray.LWH: 必须是matrix矩阵,列数为托盘数量,行数为3; 未考虑3行3列情形
        if size(da.LUArray.LWH,1)~=3,     da.LUArray.LWH=da.LUArray.LWH';     end
        
        % 如输入没有该数值
        if ~isfield(da.LUArray,'isRota')
            da.LUArray.isRota = ones(size(da.LUArray.ID));
        end
            
        % 对whichRotation在0和1时的判断, 如为2时，则不变
        if ParaArray.whichRotation == 1 % 全部允许旋转
            da.LUArray.isRota = ones(size(da.LUArray.isRota));
        elseif ParaArray.whichRotation == 0    % 全部禁止旋转
            da.LUArray.isRota = zeros(size(da.LUArray.isRota));
        elseif ParaArray.whichRotation == 2    % 部分允许旋转
            1;
        end

        % 2 Input增加间隙BUFF后的feasible的LU和BIN的长宽高转换
        da.BinArray.LWH = da.BinArray.LWH - da.BinArray.BUFF;
        da.LUArray.BUFF = [da.LUArray.BUFF; 0]; %用户给定的Buff为每个托盘增加的尺寸(总间隙)
        da.LUArray.BUFF = repmat(da.LUArray.BUFF,1,numel(da.LUArray.ID));
        da.LUArray.LWH = da.LUArray.LWH + da.LUArray.BUFF;
        %TODO: 此处增加间隙为权宜之际，未考虑rotation后的变化；后期考虑在算法中增加间隙
        
        % 3 默认将LU全部采用Horizontal方向旋转（前提：该LU允许旋转）
        % NOTE: 此处将获得1: Horizontal方向的LWH和是否Rotaed标记
        % NOTE : 直接替换了原始ORIGINAL 的 LWH
        [da.LUArray.Rotaed]= placeItemHori(da.LUArray.LWH,da.LUArray.isRota,1); %第二个参数：1: Hori; 0: Vert；其它: 原封不动
        da.LUArray.LWH = getRotaedLWH(da.LUArray.LWH, da.LUArray.Rotaed, da.LUArray.BUFF); 
        
        
        % 3 Input增加LUArray是否Rotation的判断flag TOBE DELE
%         da.LUArray.RotaFlag = zeros(1,numel(da.LUArray.ID)); TOBE DELE
    end

    %% LUID输入转换: da.LUArray.ID 转换为从1开始的类序号 方便刘工输入ID信息
    function inputLUIDExchange()
        uniLUID = unique(da.LUArray.ID);
        nLUid = size(uniLUID,2);
        tmpLUID=da.LUArray.ID; %中间变量
        for i=1:nLUid
            tmpLUID(da.LUArray.ID(:)==uniLUID(i)) = i;
        end
        da.LUArray.ID=tmpLUID;
    end

    %% 深入判断
    function deepCheck()        
        % 基础输入获取
        dBin = unique(da.BinArray.LWH','rows')'; %获取非相同车型数据
         dLU = da.LUArray.LWH;
        dLUisRota = da.LUArray.isRota;
        dLUid = da.LUArray.ID;
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
%         printstruct(da);
        da.LUIDArray.ID = unique(da.LUArray.ID);
        nLUID = numel(da.LUIDArray.ID);        

        for iID = 1:nLUID
        da.LUIDArray.Weight(iID) = sum(da.LUArray.Weight .* (da.LUArray.ID == da.LUIDArray.ID(iID)) );
        da.LUIDArray.Volume(iID) = sum(prod(da.LUArray.LWH) .* (da.LUArray.ID == da.LUIDArray.ID(iID)) );
        da.LUIDArray.isRota(iID) =  unique(da.LUArray.isRota(da.LUArray.ID == da.LUIDArray.ID(iID))); if ~isscalar(da.LUIDArray.isRota(iID)), error('致命错误'); end
        end

    end

   %% 打印Input数据
    function printInput()
        fprintf('本算例只有一个箱型 宽=%1.0f 长=%1.0f 高=%1.0f \n', unique(da.BinArray.LWH','rows')');
        fprintf('本算例只有一个箱型 宽间隙=%1.0f 长间隙=%1.0f 高间隙=%1.0f \n', unique(da.BinArray.BUFF','rows')');
        fprintf('本算例有 %d 个物品,其宽长高分别为 \n',numel(da.LUArray.ID(:)));
        fprintf('%1.1f %1.1f %1.1f \n',da.LUArray.LWH);
    end

end





