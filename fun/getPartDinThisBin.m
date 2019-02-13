function partd = getPartDinThisBin(do,ibin,d)
% getPartdinThisBin ==> GET �ṹ��d�е����ڱ�bin�Ĳ�������
%   ��nagrin����3����ȫ����do��ȡ���ݣ����򣬴�d��ȡ����
%   ָ��ibin,��ȡ��bin�ڵ�LU,Veh����Ϊ��������,�ص���LU����

luIdx = do.LU.LU_Bin(1,:) == ibin;
stripidx = do.Strip.Strip_Bin(1,:) == ibin;
itemidx = do.Item.Item_Bin(1,:) == ibin;

if nargin < 3   %ȫ����do��ȡ���ݣ�
    partd.Veh = do.Veh;
    partd.Par = do.Par;
    partd.LU = structfun(@(x) x(:,luIdx),do.LU,'UniformOutput',false);
    partd.Bin = structfun(@(x) x(:,ibin),do.Bin,'UniformOutput',false);
    partd.Strip = structfun(@(x) x(:,stripidx),do.Strip,'UniformOutput',false);
    partd.Item = structfun(@(x) x(:,itemidx),do.Item,'UniformOutput',false);
else    %ȫ����d��ȡ���ݣ�
    partd.Veh = d.Veh;
    partd.Par = d.Par;
    partd.LU = structfun(@(x) x(:,luIdx),d.LU,'UniformOutput',false);
end






%% V1 : V2 ����Ϊ3��������˫����;
% function partd = getPartDinThisBin(d,luIdx)
% % getPartdinThisBin ==> GET �ṹ��d�е����ڱ�bin�Ĳ�����������
% %   ָ��ibin,��ȡ��bin�ڵ�LU,Veh����Ϊ��������,�ص���LU����
%         
%         partd.Veh = d.Veh;        
%         partd.Par = d.Par;
%         partd.LU = structfun(@(x) x(:,luIdx),d.LU,'UniformOutput',false);
%         
% end

%% ע��
        % 2 bin�ڵ�LU
%         luIdx = tmpd.LU.LU_Bin(1,:) == ibin;    %tmpd.LU.LU_Strip(1,:) == istrip
        
        % thisd.Veh = rmfield(thisd.Veh,{'Volume','order'});
%         thisd.LU.LWH([1,2], thisd.LU.Rotaed ) = flipud(thisd.LU.LWH([1,2], thisd.LU.Rotaed)); %LU.LWH ����ת,��ָ�ԭ��
%         thisd.LU.PID = thisd.LU.OPID;     thisd.LU.SID = thisd.LU.OSID;  %  thisd.LU.LID = thisd.LU.OLID;
%         thisd.LU = rmfield(thisd.LU,{'Rotaed','order','LU_Item','DOC','LU_Strip',...
%             'LU_Bin','CoordLUBin','CoordLUStrip','LU_VehType','OPID','OSID'});
        
        
        
%     % tmpd�е�Bin��������, ����С�Ŀ�ʼ��
%     tmpusedVehIdx = max(tmpd.LU.LU_Bin(1,:)); %tmpusedVehIdx: ���һ��Bin��indexֵ
%     flagusedLUIdx = tmpd.LU.LU_Bin(1,:)==tmpusedVehIdx; % flagused: �ҳ����һ��Bin��Ӧ��LUindexֵ
%     if isSameCol(tmpd.LU)
%         % ��ȡ�����һ��Bin����������
%         lastd.LU = structfun(@(x) x(:,flagusedLUIdx),tmpd.LU,'UniformOutput',false);  %��ȡ���һ�����ڵ�LU
%         lastd.LU.LWH([1,2], lastd.LU.Rotaed ) = flipud(lastd.LU.LWH([1,2], lastd.LU.Rotaed)); %LU.LWH ����ת,��ָ�ԭ��
%         lastd.LU = rmfield(lastd.LU,{'Rotaed','order','LU_Item','DOC','LU_Strip','LU_Bin','CoordLUBin','maxL','CoordLUStrip'}); 
%         lastd.Par = tmpd.Par;
%     else
%         error('����ʹ��structfun');
%     end
