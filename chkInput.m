function [d] = chkInput(d)
% GCHECKINPUT 重要函数:输入数据da转换(不增加新数据）+输入数据da核对
% 函数1：initCheck
%   1：确保LU和Veh均为数值型矩阵
%   2：对颠倒的Veh进行转置
%   3：对不同SID下有相同ID的进行ID修改（REMOVE)
%   4：对颠倒的LU进行转置
% 函数2：deepCheck
%   1：LU高度约束
%   2：LU长宽约束
%   3：LU可堆垛的一致性约束
     
     %  Check if INPUT is legal with minor modification
     [d.LU, d.Veh] = initCheck(d.LU,d.Veh); 

    % Deep Check if INPUT is logical legal
    deepCheck(d.LU,d.Veh);
end

% 局部函数
function [LU,Veh] = initCheck(LU,Veh)
% initCheck if INPUT is legal AND 转置颠倒的 INPUT 
% 1：确保LU和Veh均为数值型矩阵
% 2：对颠倒的Veh进行转置
% 3：对不同SID下有相同ID的进行ID修改
% 4：对颠倒的LU进行转置

nVeh = length(Veh.Weight);
nLU = length(LU.Weight);

%% *************** 1 Lu和Veh的综合判断：确保为数值型矩阵 ***************
% 1: 均必须为非负/二维数组，包括标量、向量、矩阵和空数组的double型数组
try
    structfun(@(x) validateattributes(x,{'double'},{'nonnegative','2d'}), LU, 'UniformOutput', false);
    structfun(@(x) validateattributes(x,{'double'},{'nonnegative','2d'}), Veh, 'UniformOutput', false);
catch
    fLU = structfun(@(x) any( x(:) < 0), LU, 'UniformOutput', false);
    fVeh = structfun(@(x) any( x(:) < 0), LU, 'UniformOutput', false);
    if any(cell2mat(struct2cell(fLU))) || any(cell2mat(struct2cell(fVeh)))
        error('输入数据存在小于0的数'); end
end

%% *************** 2 Veh向量判断 ***************
% 判断Veh的每个field数组判定：向量/矩阵 是否不一致，是否列数不为nVeh
fields  = fieldnames(Veh);
for idx = 1:length(fields)
    aField = Veh.(fields{idx});
    try
        % V2 版本判定
        % 向量判定
        if ~strcmp('LWH',fields{idx}) % && ~strcmp('buff',fields{idx}) 目前没有buff变量
            if ~isrow(aField),  Veh.(fields{idx}) = aField';     end %如果车辆数据非行向量,修改
            validateattributes(Veh.(fields{idx}),{'double'},{'nonnegative','vector','ncols', nVeh});
        % 矩阵判定
        else
            if size(aField,1) == nVeh && size(aField,2) ~= nVeh %如果矩阵的行数==车辆数 且 列数~=车辆数
                Veh.(fields{idx}) = aField'; %行列转置
            end
            
            % 特殊矩阵判定（方阵且行列为3个或其他的判定）
            if size(aField,1) == nVeh && size(aField,2) == nVeh
                
                if nVeh==3
                    Veh.(fields{idx}) = aField';
                    warning('正在处理3个车辆的情况,未予以判断,需谨慎判断以下数据是否正确');
                    fprintf('车辆的的长宽高为: %d \n ',Veh.(fields{idx}));
                    tmpVeh = Veh.(fields{idx});
                    if ~(tmpVeh(1)*1.5 <= tmpVeh(2)) %如果第一辆车辆的宽度的1.5倍大于车辆的长度,明现为长宽高放反了
                        Veh.(fields{idx}) = aField;
                    end
                    
                else
                    warning('不是3个托盘且是方阵,仅可能在非长宽高时出现');
                end
                
            end
            
            validateattributes(Veh.(fields{idx}),{'double'},{'nonnegative','2d','ncols', nVeh});
            
        end
    catch    
        error('try error in Veh向量判断');
        % V1 版本判定
        % 判断1: 必须是数值型
        if ~ismatrix(aField) || ~isnumeric(aField), error('车辆数据存在错误(非数值型矩阵)'); end
        % 判断2: 列数如不是nVeh,必须转换(区分两种情况,对于可能是矩阵类型的LWH和buff分开判断)
        if isvector(aField) && ~strcmp('LWH',fields{idx}) && ~strcmp('buff',fields{idx}) %如果车辆是向量且不是后面两个值(此处不判断车辆长宽高LWH)
            if ~isrow(aField) %如果车辆数据非行向量,修改
                Veh.(fields{idx}) = aField';         end
        else % ismatrix %此处判断车辆长宽高是否写反 LWH buff等
            if size(aField,1) == nVeh && size(aField,2) ~= nVeh %如果矩阵的行数==车辆数 且 列数~=车辆数
                Veh.(fields{idx}) = aField'; %行列转置
            end
            if size(aField,1) == nVeh && size(aField,2) == nVeh

                if nVeh==3
                    Veh.(fields{idx}) = aField';
                    warning('正在处理3个车辆的情况,未予以判断,需谨慎判断以下数据是否正确');
                    fprintf('车辆的的长宽高为: %d \n ',Veh.(fields{idx}));
                    tmpVeh = Veh.(fields{idx});
                    if ~(tmpVeh(1)*1.5 <= tmpVeh(2)) %如果第一辆车辆的宽度的1.5倍大于车辆的长度,明现为长宽高放反了
                        Veh.(fields{idx}) = aField;
                    end             

                else
                    warning('不是3个托盘且是方阵,仅可能在非长宽高时出现');  
                end

            end        
        end
        % 判断3: 转换后必须列数是nVeh
        if size(Veh.(fields{idx}),2)  ~= nVeh %5555 重要判断
            error('Veh数据存在列数不等于车辆数情况,输入错误2');
        end
    
    end
end

%% *************** 3 LU的ID号不可在不同SID下重复(Milkrun环境) - 改变ID号 放到预处理函数Gpreproc *************** 
     % 2 如果相同ID（PID/EID/LID等）号下，对应SID号要必须不同；如相同改变ID号，直到不存在相同的ID在不同SID内;
%      if isrepeated(LU.ID,LU.SID)
%         warning('存在托盘ID号重复, 需要更正'); 
%         LU.ID = reviseID(LU.ID,LU.SID);
%      end

    % V1: 可能造成更新后ID与其它ID号相同
%     done = false;
%     while ~done
%     flag = 0;
%     uniID = unique(LU.ID);                  %不同的ID号
%     for iID = 1:length(uniID)
%         fID = LU.ID==uniID(iID);
%         uniSID = unique(LU.SID(fID));
%         if length(uniSID) > 1                %如果存在，改变ID号（ID+SID号）
%             for iSID = 1:length(uniSID)
%                 fSID = LU.SID==uniSID(iSID);
%                 f = fSID & fID;
%                 LU.ID(f) = LU.ID(f) + LU.SID(f);               
%             end
%             flag = 1;                                   % 如果flag==1，表明有修改，但还要while循环，直到不循环
%         end
%     end
%     if flag==0, done=true; end
%     end
    
%% *************** 4 LU向量判断 ***************
% 判断LU向量是否行向量; 矩阵列是否一致; 不一致进行旋转
fields  = fieldnames(LU);
for idx = 1:length(fields)
    aField = LU.(fields{idx});
      try
        % V2 版本判定
        % 向量判定
        if ~strcmp('LWH',fields{idx}) && ~strcmp('maxL',fields{idx}) && ~strcmp('margin',fields{idx}) %目前没有buff变量
            if ~isrow(aField),  LU.(fields{idx}) = aField';    end %如果车辆数据非行向量,修改
            validateattributes(LU.(fields{idx}),{'double','logical'},{'nonnegative','vector','ncols', nLU});
        % 矩阵判定
        else
            if size(aField,1) == nLU && size(aField,2) ~= nLU
                LU.(fields{idx}) = aField';   end
            
            % 特殊矩阵判定（方阵且行列为3个或其他的判定）
            if size(aField,1) == nLU && size(aField,2) == nLU
                
                if nLU==3
                    LU.(fields{idx}) = aField';
                    warning('正在处理3个托盘的情况,需谨慎处理');
                    fprintf('托盘的的长宽高为: %d \n ',LU.(fields{idx}));
                    uniID = unique(LU.ID);
                    % 托盘ID - 对应长宽LW/isRota/yID/maxL  etc... 矩阵表述的转置
                    for iLU = 1:length(uniID)
                        tmpM = LU.(fields{idx});
                        tmpVeh =tmpM(1:2,LU.ID==uniID(iLU))';    %获取宽和长,不要高
                        if isempty(tmpVeh),   error('不应该出现的错误');        end
                        if numel(unique(tmpVeh,'rows')) > 2
                            LU.(fields{idx}) = aField;
                            warning('原托盘长宽高等进行了转置 \n');
                            break;
                        end
                    end
                else
                    warning('不是3个托盘且是方阵,仅可能在非长宽高时出现');
                end
                
            end
            
            validateattributes(LU.(fields{idx}),{'double'},{'nonnegative','2d','ncols', nLU});
            
        end
    catch
         error('try error in LU向量判断');
        % V1 版本判定
        % 判断1: 必须是数值型
        if ~ismatrix(aField) || ~isnumeric(aField), error('托盘数据存在错误(非数值型矩阵)'); end
        % 判断2: 列数如不是nLU,必须转换    (区分两种情况,对于可能是矩阵类型的LWH和buff和margin分开判断)
        if isvector(aField) && ~strcmp('LWH',fields{idx}) && ~strcmp('maxL',fields{idx}) && ~strcmp('margin',fields{idx})%&& ~strcmp('buff',fields{idx} ) %如果托盘是向量且不是后面两个值
            if ~isrow(aField)
                LU.(fields{idx}) = aField';     end
        else % ismatrix %此处判断托盘长宽高是否写反
            if size(aField,1) == nLU && size(aField,2) ~= nLU
                LU.(fields{idx}) = aField';   end
            
            if size(aField,1) == nLU && size(aField,2) == nLU
                
                if nLU==3
                    LU.(fields{idx}) = aField';
                    warning('正在处理3个托盘的情况,需谨慎处理');
                    fprintf('托盘的的长宽高为: %d \n ',LU.(fields{idx}));
                    uniID = unique(LU.ID);
                    % 托盘ID - 对应长宽LW/isRota/yID/maxL  etc... 矩阵表述的转置
                    for iLU = 1:length(uniID)
                        tmpM = LU.(fields{idx});
                        tmpVeh =tmpM(1:2,LU.ID==uniID(iLU))';    %获取宽和长,不要高
                        if isempty(tmpVeh),   error('不应该出现的错误');        end
                        if numel(unique(tmpVeh,'rows')) > 2
                            LU.(fields{idx}) = aField;
                            warning('原托盘长宽高等进行了转置 \n');
                            break;
                        end
                    end
                else
                    warning('不是3个托盘且是方阵,仅可能在非长宽高时出现');
                end
                
            end
        end
        % 判断3: 转换后必须列数是nLU
        if size(LU.(fields{idx}),2)  ~= nLU
            error('LU数据存在列数不等于托盘数情况,存在不可原谅的错误2');
        end
    end
end 

end

function deepCheck(LU,Veh)
    % LU高度约束 罕见
    if min(Veh.LWH(3, :)) < max(LU.LWH(3, : )),     error('错误: 存在托盘高度 大于 本车型可用高度数据 \n');     end

    % LU长宽约束 罕见 (托盘长宽大于)
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
    % (确定：不同SID供应商下的托盘不能拼载)
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
    

%% INPUT CHECK: 不允许相同LID，相同SID，相同EID，却不能拼载，即ID不同。似乎没有这个问题
% %     tlu = array2table([LU.SID;LU.EID;LU.LID]')
% %     tluid = array2table([LU.ID]')
% %     [a,b]=findgroups(tlu);
% %     x = splitapply(@(x) length(unique(x)),tluid.Var1,a)
% % %     Num = splitapply(@(w) length(unique(w)),T2.DUEDATA,GID) %取值个数计算
% %    


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