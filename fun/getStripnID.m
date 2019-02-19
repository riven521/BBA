function [nLUID, nLULID,nLUIDBin,nLULIDBin] = getStripnID(LU)

uniStrip = unique(LU.LU_Strip(1,:));
nStrip = length(uniStrip);

[nLUID,nLULID,nLUIDBin,nLULIDBin] = deal(zeros(1,nStrip));

for i=1:nStrip
    
    istrip = uniStrip(i);
    
    idx = LU.LU_Strip(1,:) == istrip;
    
    %    unique(LU.ID(idx))
    
    if length(unique(LU.ID(idx))) == 1
        nLUID(istrip) = unique(LU.nID(idx));
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

