function [nLUID, nLULID,nLUIDBin,nLULIDBin] = getStripnID(LU)
% 本期望替换cpuStripnbItem函数，但似乎不可以，主要在量大车头下
% 通过LU计算：Strip对应的LUID数量、Strip对应相应bin内的ID数量

uniStrip = unique(LU.LU_Strip(1,:));
nStrip = length(uniStrip);

[nLUID,nLULID,nLUIDBin,nLULIDBin] = deal(zeros(1,nStrip));

for i=1:nStrip
    
    istrip = uniStrip(i);
    
    idx = LU.LU_Strip(1,:) == istrip;
    
    %    unique(LU.ID(idx))
    
    if length(unique(LU.ID(idx))) == 1  %只对该strip内LUID相同的赋值        
        nLUID(istrip) = unique(LU.nID(idx));       % 赋值Strip的托盘ID号数量 = 唯一数
    else
        nLUID(istrip) = -1;
    end
        
    if length(unique(LU.LID(idx))) == 1
        nLULID(istrip) = unique(LU.nLID(idx));
    else
        nLULID(istrip) = -1;
    end
    
    if isfield(LU, 'LU_Bin')
        if length(unique(LU.ID(idx))) == 1
            nLUIDBin(istrip) = unique(LU.nIDBin(idx));
        else
            nLUIDBin(istrip) = -1;
        end
        
        if length(unique(LU.LID(idx))) == 1
            nLULIDBin(istrip) = unique(LU.nLIDBin(idx));
        else
            nLULIDBin(istrip) = -1;
        end
    end
    
end

end



% % 
% % function [nLUID, nLULID,nLUIDBin,nLULIDBin] = getStripnID(LU)
% % % 本期望替换cpuStripnbItem函数，但似乎不可以，主要在量大车头下
% % % 通过LU计算：Strip对应的LUID数量、Strip对应相应bin内的ID数量
% % 
% % uniStrip = unique(LU.LU_Strip(1,:));
% % nStrip = length(uniStrip);
% % 
% % [nLUID,nLULID,nLUIDBin,nLULIDBin] = deal(zeros(1,nStrip));
% % 
% % for i=1:nStrip
% %     
% %     istrip = uniStrip(i);
% %     
% %     idx = LU.LU_Strip(1,:) == istrip;
% %     
% %     %    unique(LU.ID(idx))
% %     
% %     if length(unique(LU.ID(idx))) == 1  %只对该strip内LUID相同的赋值        
% %         nLUID(istrip) = unique(LU.nID(idx));       % 赋值Strip的托盘ID号数量 = 唯一数
% %     else
% %         nLUID(istrip) = -1;
% %     end
% %         
% %     if length(unique(LU.LID(idx))) == 1
% %         nLULID(istrip) = unique(LU.nLID(idx));
% %     else
% %         nLULID(istrip) = -1;
% %     end
% %     
% %     if isfield(LU, 'LU_Bin')
% %         if length(unique(LU.ID(idx))) == 1
% %             nLUIDBin(istrip) = unique(LU.nIDBin(idx));
% %         else
% %             nLUIDBin(istrip) = -1;
% %         end
% %         
% %         if length(unique(LU.LID(idx))) == 1
% %             nLULIDBin(istrip) = unique(LU.nLIDBin(idx));
% %         else
% %             nLULIDBin(istrip) = -1;
% %         end
% %     end
% %     
% % end
% % 
% % end
% % 
