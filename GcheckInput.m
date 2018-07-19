function [da] = GcheckInput(da,ParaArray)
% 先对输入数据进行转换,后对数据进行check
% 输入转换
    if size(da.LUArray.ID,1)~=1 && size(da.LUArray.ID,1)>size(da.LUArray.ID,2) %行向量
        da.LUArray.ID=da.LUArray.ID';
    end
    if size(da.BinArray.LWH,2)~=1 && size(da.BinArray.LWH,2)>size(da.BinArray.LWH,1) %列向量
        da.BinArray.LWH=da.BinArray.LWH';
    end
    if size(da.LUArray.LWH,2)~=3  %确保用户输入是三维列向量;
        da.LUArray.LWH=da.LUArray.LWH';        
    end
    da.LUArray.LWH = da.LUArray.LWH';  %转换三维列向量为三维行向量
% 输入间隙的转换
    if size(da.BinArray.BUFF,2)~=1  %BinArray.BUFF必须为列向量
        da.BinArray.BUFF=da.BinArray.BUFF';
    end
    if size(da.LUArray.BUFF,2)~=1   %LUArray.BUFF必须为列向量
        da.LUArray.BUFF=da.LUArray.BUFF';
    end

    
% 输入转换: da.LUArray.ID 转换为从1开始的类序号 方便刘工输入ID信息
    nLUid = size(unique(da.LUArray.ID),2);
    uniLUID = unique(da.LUArray.ID);
    tmpLUID=da.LUArray.ID; %中间变量
    for i=1:nLUid
        tmpLUID(da.LUArray.ID(:)==uniLUID(i)) = i;
    end
    da.LUArray.ID=tmpLUID;
% 输入转换2: 增加间隙后的转换
    da.LUArray.BUFF = [da.LUArray.BUFF;0];
    da.BinArray.LWH = da.BinArray.LWH - da.BinArray.BUFF;
    da.LUArray.LWH = da.LUArray.LWH + da.LUArray.BUFF;
    

% 输入Check
    nDim = size(da.LUArray.LWH,1); 
%     for i=1:nDim
%         fprintf('%1.0f', da.LUArray.LWH(i,:));
%         fprintf('\n');
%     end
    if nDim ~=3, error('本算例非三维算例,超出预期 \n'); end
    dLU = da.LUArray.LWH(1:nDim,:);
    dLUid = da.LUArray.ID;
    nLUid = size(unique(dLUid),2);
    tmpbinDataMatrix = da.BinArray.LWH(1:nDim,:);
    dBin = unique(tmpbinDataMatrix','rows')';
    nItem = size(dLU,2);    
    if size(dBin,2)==1
        fprintf('本算例只有一个箱型 宽=%1.0f 长=%1.0f 高=%1.0f \n', dBin);
        fprintf('本算例有 %d 个物品,其宽长高分别为 \n',nItem);
        fprintf('%1.0f %1.0f %1.0f \n',dLU);
    else
        error('本算例有多个箱型,超出期望 \n');
    end
    if numel(da.BinArray.BUFF)~=3, error('本算例箱型间隙参数不是3个,超出期望 \n'); end
    if numel(da.LUArray.BUFF)~=3, error('本算例箱型间隙参数不是3个,超出期望 \n'); end

    %% 高度约束
    if any(dBin(3) < dLU(3,:))
        error('错误: 存在托盘高度 大于 本车型可用高度数据 \n');
    end 
    %% 长宽约束
    if ParaArray.whichRotation == 1
        if min(dBin(1:2)) < max(min(dLU(1:2,:))) || max(dBin(1:2)) < max(max(dLU(1:2,:)))
            error('错误: 存在托盘长或宽 大于 本车型长或宽 \n');
        end
    else %不准rotation 
        if any(dBin(1) < dLU(1,:))
            error('错误: 存在托盘宽度 大于 本车型可用宽度数据 \n');
        end
        if any(dBin(2) < dLU(2,:))
            error('错误: 存在托盘长度 大于 本车型可用长度数据 \n');
        end
    end
    %% 托盘ID（材质，类型）与长宽约束  
    for iLU = 1:nLUid
%         dLULW = dLU(1:2,:); %获取宽和长 不要高
        tmp = dLU(1:2,dLUid==iLU)'; %获取
        if numel(unique(tmp,'rows')) > 2
            error('错误: 存在托盘ID相同 但其长宽不同数据 \n');
        end
    end    
end





