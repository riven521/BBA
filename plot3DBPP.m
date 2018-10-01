function [] = plot3DBPP(d,par)
% 作图函数:三维BPP

fields = fieldnames(par);
aField = [];
for idx = 1:length(fields), aField = [aField par.(fields{idx})];   end

% 1: 按照LID区分颜色
% 2: 按照Item内Lu isHeightFull区分颜色
% 3: 按照Strip内Lu isWidthFull区分颜色
% 23: 按照Strip内Lu isWidthFull和isHeightFull区分颜色
colorType = 1;  

if colorType == 1
 nIDType = unique(d.LU.ID);
elseif colorType == 2
    d.LU.isHeightFull = ones(size(d.LU.ID))*1;
    nIDType = unique(d.Item.isHeightFull);
    f = find(d.Item.isHeightFull == 0);    
    for i=1:numel(f)
        f2 = d.LU.LU_Item(1,:) == f(i);
        d.LU.isHeightFull(f2) = 0;
    end
elseif colorType == 3
    d.LU.isWidthFull = ones(size(d.LU.ID))*1;
    nIDType = unique(d.Strip.isWidthFull);
    f = find(d.Strip.isWidthFull == 0);    
    for i=1:numel(f)
        f2 = d.LU.LU_Strip(1,:) == f(i);
        %         f22 = find(d.LU.LU_Strip(1,:) == f(i))
        d.LU.isWidthFull(f2) = 0;
    end
elseif colorType == 23
     d.LU.isHeightFull = ones(size(d.LU.ID))*1;
    d.LU.isWidthFull = ones(size(d.LU.ID))*1;
    nIDType = unique(d.Item.isHeightFull);
    f = find(d.Item.isHeightFull == 0);    
    for i=1:numel(f)
        f2 = d.LU.LU_Item(1,:) == f(i);
        d.LU.isHeightFull(f2) = 0;
    end
    f = find(d.Strip.isWidthFull == 0);
    for i=1:numel(f)
        f2 = d.LU.LU_Strip(1,:) == f(i);
        d.LU.isWidthFull(f2) = 0;
    end
end

% d.Strip.LW(:,15)
nColors = hsv(length(nIDType)); %不同类型LU赋予不同颜色
nBin = max(d.LU.LU_Bin(1,:)); %bin的个数;

for ibin=1:nBin
    figure('name',num2str([ibin, aField]));
    f = d.LU.LU_Bin(1,:)==ibin;
    yxz = d.LU.LWH(:,f);
    coord = d.LU.CoordLUBin(:,f);
    
    if colorType == 1
        lid  = d.LU.ID(:,f);
    elseif colorType == 2   
        lid  = d.LU.isHeightFull(:,f);
    elseif colorType == 3
        lid  = d.LU.isWidthFull(:,f);
    elseif colorType == 23
        lid1  = logical(d.LU.isHeightFull(:,f));
        lid2  = logical(d.LU.isWidthFull(:,f));
        lid = lid1 & lid2;
    end
    % LU在Bin内的顺序
    seq = d.LU.LU_Bin(2,f);
    
    for i=1:numel(seq)
        idx =find(seq==i); % LU进入Bin内的画图顺序
        LUColor = 0.8*nColors(nIDType==lid(idx), : );
        plotcube(yxz(:,idx)',coord(:,idx)',0.7,LUColor);
        
        % Set the lable and the font size
        axis equal;         grid on;        view(130,30);
        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        xlim([0 d.Veh.LWH(1,1)]);         ylim([0 d.Veh.LWH(2,1)]);       zlim([0 d.Veh.LWH(3,1)]);        
    end
end

end
