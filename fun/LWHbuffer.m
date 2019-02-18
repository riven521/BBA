function [LWHBuff,OLWH] = LWHbuffer(LWH, margin, rota)
% add LWH的buffer 且增加旋转状态
% 参数3必须有，默认为则不恢复旋转
    narginchk(2,3);

    OLWH = LWH;
    
    LWHBuff = LWH;    
    if nargin > 2
        rota = logical(rota);
        % 不旋转的不动（包括buff）；旋转的长宽替换（包括buff）
        LWHBuff(1, rota) = LWH(2, rota) + (margin(1, rota) +margin(2, rota)); %宽度（左右）
        LWHBuff(1, ~rota) = LWH(1, ~rota) + (margin(1, ~rota) + margin(2, ~rota));
        
        LWHBuff(2, rota) = LWH(1, rota) + (margin(3, rota) +margin(4, rota)); %长度（上下）
        LWHBuff(2, ~rota) = LWH(2, ~rota) + (margin(3, ~rota) + margin(4, ~rota));
    else
        LWHBuff(1, :) = LWH(1, :) + (margin(1, :) +margin(2, :));  %宽度（左右）
        LWHBuff(2, :) = LWH(2, :) + (margin(3, :) + margin(4, :));  %长度（上下）
    end

end