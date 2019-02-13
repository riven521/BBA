function [TLU] = getTableLU(d)
%% ת��ΪTable: TLU
if isfield(d, 'LU')
    TLU = struct2table(structfun(@(x) x',d.LU,'UniformOutput',false)); % TVEH = struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false));
else
    TLU = struct2table(structfun(@(x) x',d,'UniformOutput',false)); 
end
                                                                            % TLU = TLU(:,{'CoordLUBin','LWH','LID','SID','PID','Weight','LU_VehType','LU_Bin','LU_Item','isShuaiWei','Rotaed'});
% splitvar ��ȡ��ά�б���->������
TLU.L = TLU.LWH(:,1); 
TLU.W = TLU.LWH(:,2); 
TLU.H = TLU.LWH(:,3);  %���ӳ���ߵ���

TLU.Properties.VariableNames
if any(strcmp('LU_Bin', TLU.Properties.VariableNames)), 
    TLU.BINID = TLU.LU_Bin(:,1);  
    TLU.BINSEQ = TLU.LU_Bin(:,2); 
end
if any(strcmp('LU_Item', TLU.Properties.VariableNames)), 
    TLU.ITEMID = TLU.LU_Item(:,1);  
    TLU.ITEMSEQ = TLU.LU_Item(:,2);
end
if any(strcmp('CoordLUBin', TLU.Properties.VariableNames)), 
    TLU.X = TLU.CoordLUBin(:,1);  
    TLU.Y = TLU.CoordLUBin(:,2); 
    TLU.Z = TLU.CoordLUBin(:,3); 
end
end