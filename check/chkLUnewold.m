% 每次RunAlgorithm后，判断d.LU 输入与输出的差距
function chkLUnewold(lui,luo)
% checkLUnewold: 核验 lui 和 luo 均为LU
% lui is a struct of LU
% luo is a struct of LU

    if ~isstruct(lui) || ~isstruct(luo)
        error('lui or luo is not a struct');
    end

    lui = struct2table(structfun(@(x) x',lui,'UniformOutput',false));
    luo = struct2table(structfun(@(x) x',luo,'UniformOutput',false));   
    
    luo.LWH(luo.Rotaed,1) = lui.LWH(luo.Rotaed,1);
    luo.LWH(luo.Rotaed,2) = lui.LWH(luo.Rotaed,2);
    
    % 1.1 LWH CHECK
    if any(luo.LWH ~= lui.LWH)
       % any(tout.LWH ~= tin.LWH);
        error('LWH');
    end
    
    % 1.2 Weight CHECK
    if any(luo.Weight ~= lui.Weight) || any(luo.OPID ~= lui.PID) || any(luo.OSID ~= lui.SID)  || any(luo.OEID ~= lui.EID)
       % any(tout.Weight ~= tin.Weight),
        error('other');
    end
end