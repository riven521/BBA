function [TLU] = getTableLU(d)
%% ת��ΪTable: TLU
                                                                            % TVEH = struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false));
TLU = struct2table(structfun(@(x) x',d.LU,'UniformOutput',false));
                                                                            % TLU = TLU(:,{'CoordLUBin','LWH','LID','SID','PID','Weight','LU_VehType','LU_Bin','LU_Item','isShuaiWei','Rotaed'});
% splitvar ��ȡ��ά�б���->������
TLU.L = TLU.LWH(:,1); 
TLU.W = TLU.LWH(:,2); 
TLU.H = TLU.LWH(:,3);  %���ӳ���ߵ���

TLU.BINID = TLU.LU_Bin(:,1);  
TLU.BINSEQ = TLU.LU_Bin(:,2);
TLU.ITEMID = TLU.LU_Item(:,1); 
TLU.ITEMSEQ = TLU.LU_Item(:,2);
end