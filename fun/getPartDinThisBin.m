function partd = getPartDinThisBin(do,ibin,d)
% getPartdinThisBin ==> GET 结构体d中的属于本bin的部分数据
%   若nagrin不到3个，全部从do中取数据；否则，从d中取数据
%   指定ibin,获取该bin内的LU,Veh等作为输入数据,重点是LU数据

luIdx = do.LU.LU_Bin(1,:) == ibin;
stripidx = do.Strip.Strip_Bin(1,:) == ibin;
itemidx = do.Item.Item_Bin(1,:) == ibin;

if nargin < 3   %全部从do中取数据；
    partd.Veh = do.Veh;
    partd.Par = do.Par;
    partd.LU = structfun(@(x) x(:,luIdx),do.LU,'UniformOutput',false);
    partd.Bin = structfun(@(x) x(:,ibin),do.Bin,'UniformOutput',false);
    partd.Strip = structfun(@(x) x(:,stripidx),do.Strip,'UniformOutput',false);
    partd.Item = structfun(@(x) x(:,itemidx),do.Item,'UniformOutput',false);
else    %全部从d中取数据；
    partd.Veh = d.Veh;
    partd.Par = d.Par;
    partd.LU = structfun(@(x) x(:,luIdx),d.LU,'UniformOutput',false);
end






%% V1 : V2 增加为3个参数，双重用途
% function partd = getPartDinThisBin(d,luIdx)
% % getPartdinThisBin ==> GET 结构体d中的属于本bin的部分托盘数据
% %   指定ibin,获取该bin内的LU,Veh等作为输入数据,重点是LU数据
%         
%         partd.Veh = d.Veh;        
%         partd.Par = d.Par;
%         partd.LU = structfun(@(x) x(:,luIdx),d.LU,'UniformOutput',false);
%         
% end

%% 注释
        % 2 bin内的LU
%         luIdx = tmpd.LU.LU_Bin(1,:) == ibin;    %tmpd.LU.LU_Strip(1,:) == istrip
        
        % thisd.Veh = rmfield(thisd.Veh,{'Volume','order'});
%         thisd.LU.LWH([1,2], thisd.LU.Rotaed ) = flipud(thisd.LU.LWH([1,2], thisd.LU.Rotaed)); %LU.LWH 如旋转,则恢复原形
%         thisd.LU.PID = thisd.LU.OPID;     thisd.LU.SID = thisd.LU.OSID;  %  thisd.LU.LID = thisd.LU.OLID;
%         thisd.LU = rmfield(thisd.LU,{'Rotaed','order','LU_Item','DOC','LU_Strip',...
%             'LU_Bin','CoordLUBin','CoordLUStrip','LU_VehType','OPID','OSID'});
        
        
        
%     % tmpd中的Bin是排序后的, 从最小的开始试
%     tmpusedVehIdx = max(tmpd.LU.LU_Bin(1,:)); %tmpusedVehIdx: 最后一个Bin的index值
%     flagusedLUIdx = tmpd.LU.LU_Bin(1,:)==tmpusedVehIdx; % flagused: 找出最后一个Bin对应的LUindex值
%     if isSameCol(tmpd.LU)
%         % 获取仅最后一个Bin的输入数据
%         lastd.LU = structfun(@(x) x(:,flagusedLUIdx),tmpd.LU,'UniformOutput',false);  %仅取最后一辆车内的LU
%         lastd.LU.LWH([1,2], lastd.LU.Rotaed ) = flipud(lastd.LU.LWH([1,2], lastd.LU.Rotaed)); %LU.LWH 如旋转,则恢复原形
%         lastd.LU = rmfield(lastd.LU,{'Rotaed','order','LU_Item','DOC','LU_Strip','LU_Bin','CoordLUBin','maxL','CoordLUStrip'}); 
%         lastd.Par = tmpd.Par;
%     else
%         error('不能使用structfun');
%     end
