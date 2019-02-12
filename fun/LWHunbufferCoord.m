function smallCoordLUBin = LWHunbufferCoord(CoordLUBin, margin, idx)
% Update LWHµÄCoordLUBin ¿¼ÂÇmargin ²»¿¼ÂÇÐý×ª

    narginchk(2,2);
    
    smallCoordLUBin = CoordLUBin;        
    
    if nargin > 2
        smallCoordLUBin(1,idx) = CoordLUBin(1,idx) + margin(1,idx);
        smallCoordLUBin(2,idx) = CoordLUBin(2,idx) + margin(4,idx);
    else
        smallCoordLUBin(1,:) = CoordLUBin(1,:) + margin(1,:);
        smallCoordLUBin(2,:) = CoordLUBin(2,:) + margin(4,:);
    end

end

