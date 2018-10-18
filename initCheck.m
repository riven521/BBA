%% 行列转换检验
function [LU,Veh] = initCheck(LU,Veh)
% 判断小于0
fLU = structfun(@(x) any( x(:) < 0), LU, 'UniformOutput', false);
fVeh = structfun(@(x) any( x(:) < 0), LU, 'UniformOutput', false);
if any(cell2mat(struct2cell(fLU))) || any(cell2mat(struct2cell(fVeh)))
    error('输入数据存在小于0的数'); end

nVeh = length(Veh.ID);
nLU = length(LU.ID);
%% *************** Veh向量判断 ***************
% 判断Veh向量是否行向量; 矩阵列是否一致; 不一致进行旋转
fields  = fieldnames(Veh);
% 对Veh矩阵中每个字段分别判断, 每个字段的列数必须是车辆数nVeh
for idx = 1:length(fields)
    aField = Veh.(fields{idx});
    % 判断1: 必须是数值型
    if ~ismatrix(aField) || ~isnumeric(aField), error('车辆数据存在错误(非数值型矩阵)'); end
    % 判断2: 列数如不是nVeh,必须转换(区分两种情况,对于可能是矩阵类型的LWH和buff分开判断)
    if isvector(aField) && ~strcmp('LWH',fields{idx}) && ~strcmp('buff',fields{idx}) %如果车辆是向量且不是后面两个值(此处不判断车辆长宽高LWH)
        if ~isrow(aField) %如果车辆数据非行向量,修改
            Veh.(fields{idx}) = aField';         end
                                    %         if length(aField) ~= nVeh %如有对vector才可以用length判断
                                    %             error('Veh数据输入存在错误1');            end
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
%% *************** LU向量判断 ***************
% 判断LU向量是否行向量; 矩阵列是否一致; 不一致进行旋转
fields  = fieldnames(LU);
for idx = 1:length(fields)
    aField = LU.(fields{idx});
    % 判断1: 必须是数值型
    if ~ismatrix(aField) || ~isnumeric(aField), error('托盘数据存在错误(非数值型矩阵)'); end
    % 判断2: 列数如不是nLU,必须转换    (区分两种情况,对于可能是矩阵类型的LWH和buff和margin分开判断)
    if isvector(aField) && ~strcmp('LWH',fields{idx}) && ~strcmp('maxL',fields{idx}) && ~strcmp('margin',fields{idx})%&& ~strcmp('buff',fields{idx} ) %如果托盘是向量且不是后面两个值
        if ~isrow(aField)
            LU.(fields{idx}) = aField';     end
                                            %         if length(aField) ~= nLU %如有对vector才可以用length判断
                                            %             error('LU数据输入存在不可原谅的错误1');   end
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
        size(LU.(fields{idx}),2)
        nLU
        error('LU数据存在列数不等于托盘数情况,存在不可原谅的错误2');   end
end 

end