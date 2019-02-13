%% chkStrip Strip结构体核对
%% Form
%    chkStrip(Strip)

function chkStrip(Strip)
     if ~isstruct(Strip), error('NOT a STRUCT!'); end
     if isstruct(Strip)
     T = struct2table(structfun(@(x) x', Strip,'UniformOutput',false));  end
    
     % 多维度判断：Strip混合与否 Strip内LID不能超2个
     % if any(~iscell(T.LID(fmix))),   error('混合Strip非cell或LID数量不为2'); end
     fmix = T.isMixed == 1;
     if any(T.nbLID(fmix) >=3) || any(T.nbLID(fmix) <= 1)
         error('混合Strip内LID数量不为2'); end     
     if any(T.nbLID(~fmix) ~= 1)
         error('非混合Strip内LID数量不为1'); end
     
     % 判断：Width Not Full 有可能时 Mixed
     fwid = T.isWidthFull == 0;
     if any(fwid&fmix), 
         fwid
         fmix
         warning('宽度不满的Strip不能为混合Strip'); end

end

% TO 修改使用
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
end