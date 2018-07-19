function [LWHItemSortHori,idxHori] = horiOrient(LWHItemSort)
% 普通函数;获取
    LWHItemSortHori = LWHItemSort;
    idxHori = LWHItemSort(1,:) < LWHItemSort(2,:);
    tmpminv = min(LWHItemSort);
    tmpmaxv = max(LWHItemSort);
    LWHItemSortHori(1,idxHori) = tmpmaxv(1,idxHori);
    LWHItemSortHori(2,idxHori) = tmpminv(1,idxHori);
end