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

function [d] = GcheckInput(d)
    % Check if INPUT is legal with minor modification
     [d.LU, d.Veh] = initCheck(d.LU,d.Veh); %DO More time

    % Deep Check if INPUT is logical legal
    deepCheck(d.LU,d.Veh);
end

function deepCheck(LU,Veh)
    % 高度约束 罕见
    if min(Veh.LWH(3, :)) < max(LU.LWH(3, : )),     error('错误: 存在托盘高度 大于 本车型可用高度数据 \n');     end

    % 长宽约束 罕见 (托盘长宽大于)
    flag1 = LU.isRota == 1;  %允许rota的LU标记
    flag0 = LU.isRota == 0;  %不允许rota的LU标记
    % 测试不允许旋转的LU
    if min(Veh.LWH(1, :)) < max(LU.LWH(1,flag0)), error('错误: 存在托盘宽度 大于 本车型可用宽度数据 \n');    end
    if min(Veh.LWH(2, :)) < max(LU.LWH(2,flag0)), error('错误: 存在托盘长度 大于 本车型可用长度数据 \n');    end
    % 测试允许旋转的LU
    if ~all(~flag1) && min(Veh.LWH(1:2)) < max(min(LU.LWH(1:2, flag1)))  %flag1不是全0 且满足后面条件
        error('错误: 存在托盘最短边的最大值 大于 本车型长或宽 \n');
    end
    if ~all(~flag1) && max(Veh.LWH(1:2)) < max(max(LU.LWH(1:2, flag1 ))) %flag1不是全0 且满足后面条件
        error('错误: 存在托盘最长边的最大值 大于 本车型长或宽 \n');
    end

    %% V2 : 测试LU：存在托盘ID类型相同 （且不同供应商） 但其长宽不同数据
% %     uniSID = unique(LU.SID);
% %     uniID = unique(LU.ID);
% %     % 托盘ID - 对应长宽/isRota/yID/maxL  常见555
% %     for iSID = 1:length(uniSID)
% %         fSID = LU.SID==uniSID(iSID);
% %         for iLU = 1:length(uniID)
% %             fID = LU.ID==uniID(iLU);
% %             f = fSID & fID;
% %             tmpLULW = LU.LWH(1:2, f)';                       if isempty(tmpLULW),  error('不应该出现的错误');  end %获取宽和长,不要高           
% %             if numel(unique(tmpLULW,'rows')) > 2,  error('错误: 存在托盘ID类型相同 但其长宽不同数据 \n');    end
% %             
% %             tmpLULW = LU.isRota(:,f)';                           if isempty(tmpLULW),      error('不应该出现的错误');       end
% %             if numel(unique(tmpLULW)) > 2,      error('错误: 存在托盘ID类型相同 但其允许旋转类型不同数据 \n');    end
% %             
% %             %         tmp = LU.maxL(:,LU.ID==uniID(iLU));
% %             %         if isempty(tmp),      error('不应该出现的错误');       end
% %             %         if numel(unique(tmp)) > 2,  error('错误: 存在托盘ID类型相同 但其最大层数maxL类型不同数据 \n');    end
% %             %         tmp = LU.yID(:,LU.ID==uniID(iLU));
% %             %         if isempty(tmp),      error('不应该出现的错误');       end
% %             %         if numel(unique(tmp)) > 2,  error('错误: 存在托盘ID类型相同 但其yID类型不同数据 \n');    end
% %         end
% %     end
    

    %% V1 : 测试LU：存在托盘ID类型相同 但其长宽不同数据 (通过修改ID,确保不同SID下的ID是不同的)
    uniID = unique(LU.ID);
    % 托盘ID - 对应长宽/isRota/yID/maxL  常见555
    for iLU = 1:length(uniID)
        tmp = LU.LWH(1:2,LU.ID==uniID(iLU))';   %获取宽和长,不要高 
        if isempty(tmp),      error('不应该出现的错误');       end
                        %x=[LU.LWH;LU.ID]
                %         x=[LU.LWH;LU.ID;LU.SID]
                %         f = LU.ID==uniID(iLU)
                %         y=x(:,f)'
                %         unique(y,'rows')
        if numel(unique(tmp,'rows')) > 2,  error('错误: 存在托盘ID类型相同 但其长宽不同数据 \n');    end        
        tmp = LU.isRota(:,LU.ID==uniID(iLU));
        if isempty(tmp),      error('不应该出现的错误');       end
        if numel(unique(tmp)) > 2,  error('错误: 存在托盘ID类型相同 但其允许旋转类型不同数据 \n');    end
%         tmp = LU.maxL(:,LU.ID==uniID(iLU));
%         if isempty(tmp),      error('不应该出现的错误');       end
%         if numel(unique(tmp)) > 2,  error('错误: 存在托盘ID类型相同 但其最大层数maxL类型不同数据 \n');    end
%         tmp = LU.yID(:,LU.ID==uniID(iLU));
%         if isempty(tmp),      error('不应该出现的错误');       end
%         if numel(unique(tmp)) > 2,  error('错误: 存在托盘ID类型相同 但其yID类型不同数据 \n');    end
    end
    



end
%                         tmpVeh =tmpM(1:2,LU.ID==uniID(iLU))';    %获取宽和长,不要高

    
%%
                                    %         fLU = structfun(@isrow, rmfield(d.LU,'LWH'), 'UniformOutput', false);
                                    %         fVeh = structfun(@isrow, rmfield(d.Veh,{'LWH','buff'}), 'UniformOutput', false);
                                    %         idx = find(cell2mat( struct2cell(fVeh) )== 0)

                                    % 初步判断输入: d.Veh
                                        %         if numel(d.Veh.LWH(:))~=numel(d.Veh.BUFF(:)) || numel(d.Veh.LWH(:)) ~= 3 ||...
                                        %                 ~isscalar(d.Veh.Weight)
                                        %             error('本算例da.BinArray错误,超出预期 \n');
                                        %         end
                                            % 初步判断输入: d.LU
                                    %         if numel(d.LU.ID(:))~=numel(d.LU.Weight(:)) || numel(d.LU.BUFF(:)) ~= 2 || ...
                                    %                 numel(d.LU.ID(:))*3~=numel(d.LU.LWH(:)) || ~ismatrix(d.LU.LWH)    
                                    %             error('本算例da.LUArray错误,超出预期 \n');

% if any(dBin(1) < LU.LWH(1,flag0)),  error('错误: 存在托盘宽度 大于 本车型可用宽度数据 \n');    end
% if any(dBin(2) < LU.LWH(2,flag0)),  error('错误: 存在托盘长度 大于 本车型可用长度数据 \n');    end
% 长宽约束 罕见 (托盘长宽大于) TOBE DELE
%         if ParaArray.whichRotation == 1
%             if min(dBin(1:2)) < max(min(dLU(1:2,:))) || max(dBin(1:2)) < max(max(dLU(1:2,:)))
%                 error('错误: 存在托盘长或宽 大于 本车型长或宽 \n');
%             end
%         else %不准rotation
%             if any(dBin(1) < dLU(1,:)),  error('错误: 存在托盘宽度 大于 本车型可用宽度数据 \n');    end
%             if any(dBin(2) < dLU(2,:)),  error('错误: 存在托盘长度 大于 本车型可用长度数据 \n');    end
%         end


    %% 嵌套函数
%     %% 行列转换检验
%     function initCheck()
%         % 判断小于0
%         fLU = structfun(@(x) any( x(:) < 0), d.LU, 'UniformOutput', false);
%         fVeh = structfun(@(x) any( x(:) < 0), d.LU, 'UniformOutput', false);
%         if any(cell2mat(struct2cell(fLU))) || any(cell2mat(struct2cell(fVeh)))
%             error('输入数据存在小于0的数'); end
%         
%         nVeh = length(d.Veh.ID);
%         nLU = length(d.LU.ID);
%         % 判断Veh向量是否行向量; 矩阵列是否一直
%         fields  = fieldnames(d.Veh);
%         for idx = 1:length(fields)
%             aField = d.Veh.(fields{idx});
%             if ~ismatrix(aField) || ~isnumeric(aField), error('输入数据存在错误'); end
%             if isvector(aField)
%                 if ~isrow(aField)
%                     d.Veh.(fields{idx}) = aField';    end
%                 if length(aField) ~= nVeh
%                     error('Veh数据输入存在错误');   end
%             else % ismatrix
%                 if size(aField,1) == nVeh && size(aField,2) ~= nVeh
%                     d.Veh.(fields{idx}) = aField';
%                 if size(aField,2)  ~= nVeh
%                     error('Veh数据输入存在错误');   end                    
%                 end
%             end                
%         end
%         % 判断LU向量是否行向量; 矩阵列是否一直
%         fields  = fieldnames(d.LU);
%         for idx = 1:length(fields)
%             aField = d.LU.(fields{idx});
%             if ~ismatrix(aField) || ~isnumeric(aField), error('输入数据存在错误'); end
%             if isvector(aField)
%                 if ~isrow(aField)
%                     d.LU.(fields{idx}) = aField';     end
%                 if length(aField) ~= nLU
%                     error('LU数据输入存在错误');   end
%             else % ismatrix
%                 if size(aField,1) == nLU && size(aField,2) ~= nLU
%                     d.LU.(fields{idx}) = aField';    end
%                 if size(aField,2)  ~= nLU
%                     error('LU数据输入存在错误');   end
%             end
%         end
%     end




    %% LUID输入转换: d.LU.ID 转换为从1开始的类序号 方便刘工输入ID信息
%     function idExchange()
% 
%         
% %         uniLUID = unique(d.LU.ID);
% %         nLUid = length(uniLUID);
% %         tmpLUID=d.LU.ID; %中间变量
% %         for i=1:nLUid
% %             tmpLUID(d.LU.ID(:)==uniLUID(i)) = i;
% %         end
% %         d.LU.ID=tmpLUID;
%     end