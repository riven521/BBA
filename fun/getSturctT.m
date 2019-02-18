function [LU] = getSturctT(TLU)
%% 转变为struct为table
if ~istable(TLU)
    error('input is not a table');
end
   
if any(strcmp('LWH', TLU.Properties.VariableNames))
    TLU.L = []; 
    TLU.W = [];
    TLU.H = [];
end
if any(strcmp('LU_Bin', TLU.Properties.VariableNames))
    TLU.BINID = []; 
    TLU.BINSEQ = []; 
end
if any(strcmp('LU_Item', TLU.Properties.VariableNames))
    TLU.ITEMID = []; 
    TLU.ITEMSEQ = []; 
end
if any(strcmp('CoordLUBin', TLU.Properties.VariableNames))
    TLU.X = [];
    TLU.Y = [];
    TLU.Z = [];
end

if istable(TLU)
    LU = table2struct(TLU,'ToScalar',true);
    LU = (structfun(@(x) x',LU,'UniformOutput',false));
end

end