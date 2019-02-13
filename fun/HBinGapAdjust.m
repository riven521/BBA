function [LU] = HBinGapAdjust(LU,VEH)
% HBinGapAdjust ==> Gap Adjust for Each Bin
 
TLU = struct2table(structfun(@(x) x',LU,'UniformOutput',false));
TVEH = struct2table(structfun(@(x) x',VEH,'UniformOutput',false));

typeVeh = unique(TLU.LU_Bin(:,1));
numVeh = length(typeVeh);      if numVeh>1, error('多个车辆需要调整Gap'); end
% 循环每个VEH，找到对应子托盘subTLU和子车辆subVeh
for idxVeh = 1:numVeh
    subTLU = TLU(TLU.LU_Bin(:,1) == typeVeh(idxVeh), : );    
    subVeh = TVEH(unique(subTLU.LU_VehType), :); %TODO CHECK 是否LU_VehType就是车辆排序的第一个?
        if height(subVeh)>1,  error('NOT LU in the same Veh type'); end
        
    % 555  子托盘调整重要函数
    subTLUNew = HGapAdjust(subTLU,subVeh);   
    TLU(TLU.LU_Bin(:,1) == typeVeh(idxVeh), : ) = subTLUNew;    
end

% TLU to LU
LU = table2struct(TLU,'ToScalar',true); % struct2table(structfun(@(x) x',dd.LU,'UniformOutput',false));
LU = (structfun(@(x) x',LU,'UniformOutput',false));
end

% 通用函数 ： 获取托盘（集）的多边形 
%   TLU：托盘（集合）
%   isRota : 是否需要旋转，改变长宽
%   varargin：是否采用固定坐标（即LU的起始位置）
function [pgon] = getPgLU(TLU,isRota,varargin)
% polygon of LUs ( Note: only get polyshape whose hight = 0 )
TLU = sortrows(TLU,'CoordLUBin');
P = [];
for idxl=1:height(TLU)
    % 采用固定坐标, 若有可变输入参数
    if nargin>2
        x = varargin{1};
        y = varargin{2};
    else
        x=TLU.CoordLUBin(idxl,1)-TLU.margin(idxl,1);
        y=TLU.CoordLUBin(idxl,2)-TLU.margin(idxl,4);
    end
    
    w = TLU.LWH(idxl,1) + TLU.margin(idxl,1 ) + TLU.margin(idxl,2 );
     l = TLU.LWH(idxl,2) + TLU.margin(idxl,3 ) + TLU.margin(idxl,4 );
    
    % 旋转宽长, 若isRota==1
    if isRota
        P =  [P;pgRectangle(x,y,l,w);[NaN,NaN]];
    else
        P =  [P;pgRectangle(x,y,w,l);[NaN,NaN]];
    end
end
pgon = polyshape(P);
end

% 获取待判断的内点, 由Boundary和Centroid共同构成
function [P] = getBoundaryCentroid(pgon)
[x,y] = boundary(pgon);
[xcenter,ycenter] = centroid(pgon);
P = [[x;xcenter],[y;ycenter]];
end

%% 555  子托盘调整重要函数
% 均在同一bin内的LU, VEH是一个车
function [LU] = HGapAdjust(LU,VEH)
global parMulipleGap ISplotShowGapAdjust

% INITILIZE
flagGap=0; % 1: 调整 0: 未调整

bottomLU = LU(LU.CoordLUBin(:,3)==0, : );  % 底层的托盘 % plotSolutionT(LU,VEH);

% Get pgGap: VEH内的剩余空白区域  
pgLU = getPgLU(bottomLU,0);  % plot(pgLU);
pgVEH = polyshape(pgRectangle(0,0,VEH.LWH(1,1),VEH.LWH(1,2)));	    
pgGap = subtract(pgVEH,pgLU);    if pgGap.NumRegions  > 1,  warning('Exsit %d Regions in this pgon', pgGap.NumRegions);  end
                        % TODO pgGap排序? ,gBlanks = sortregions(pgBlanks,'area','ascend');

% CoordGap : pgGap 的 坐标 先X后Y 
[XBound,YBound]=boundary(pgGap);
coordGapArray = fliplr(sortrows([YBound,XBound]));  % 先Y后X排序, 继而左右改变顺序

% 循环每个Boundary顶点
for iVertex=1:size(coordGapArray,1)
    % Get Each GapBound Coord X and Y
    % 每个边界顶点的坐标值
    coordX = coordGapArray(iVertex,1);
    coordY = coordGapArray(iVertex,2);
            
    % 获取: 满足条件托盘序号 Get idxLUs i.e. Y of LU is larger than coordY of bound vertex 
    flagLeqLU = bottomLU.CoordLUBin(:,2) - bottomLU.margin(:,3) > coordY;
    idxLUs = find(flagLeqLU);         % 找出比yPoint值 高的且是bottom的LU索引号    
    
    % 如果没有该类型托盘，继续.
    if isempty(idxLUs), continue; end       
           
    % 节约时间: 如果该Boundary的顶点不能放下最小的LU边(即超出车辆的边界), 则跳出该Boundary到下一个
    minLW = min(min(LU.LWH(idxLUs,[1,2])));
    if ~isinterior(pgVEH,coordX+minLW,coordY) || ~isinterior(pgVEH,coordX,coordY+minLW) || ~isinterior(pgVEH,coordX+minLW/2, coordY+minLW/2),    continue;    end
    
    % 循环每个Boundary顶点下每个符合条件的LU    
    %   排序: 托盘排序
    [~,ff] = sortrows([bottomLU.CoordLUBin(idxLUs,[2,1])]);  idxLUs = idxLUs(ff);
    for ii=1:length(idxLUs)
        idxLU = idxLUs(ii); % bottomLU的序号
        thisLU = bottomLU(idxLU,:);
        
        % pgLU1 待合并的LU pgGapNew 合并后的新Gap
        pgLU1 = getPgLU(thisLU, 0);        
        pgGapNew = union(pgGap,pgLU1);
        
        % pgLU2 待新摆放的LU位置
        pgLU2 = getPgLU(thisLU, 0, coordX,coordY);        %         pgLU2 =  polyshape(pgRectangle(coordX,coordY,wLU,lLU));  
        flagLU = all(isinterior(pgGapNew, getBoundaryCentroid(pgLU2)));        
        
        pgIntersect = intersect(pgLU2,subtract(pgVEH,pgGapNew));
        if pgIntersect.NumRegions ~= 0, flagLU = 0; end
        
        % pgLU2Rota 待新摆放的旋转后LU位置
        flagLURota=0;
        if thisLU.isRota     %bottomLU.isRota(idxLU)
            pgLU2Rota =  getPgLU(thisLU, 1, coordX,coordY); %polyshape(pgRectangle(coordX,coordY,lLU,wLU));
            flagLURota = all(isinterior(pgGapNew, getBoundaryCentroid(pgLU2Rota))); 
            
            pgIntersectRota = intersect(pgLU2Rota,subtract(pgVEH,pgGapNew));
            if pgIntersectRota.NumRegions ~= 0, flagLURota = 0; end
        end        
        
        % 若都可放下, 二者选择一个，横放优先
        if flagLU && flagLURota
            if thisLU.LWH(1) > thisLU.LWH(2) %wLU>lLU
                flagLURota=0;
            else
                flagLU=0;
            end
        end
            
        if ISplotShowGapAdjust
        pausetime = 0.0;   
        % plot 某些vertex的尝试过程
        plot(pgVEH,'FaceColor','white','FaceAlpha',0.01);
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
        plot(pgLU,'FaceColor','green','FaceAlpha',0.2)
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
        plot(pgGap,'FaceColor','blue','FaceAlpha',0.2)
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
        plot(coordGapArray(:,1),coordGapArray(:,2),'.', 'MarkerSize', 8);
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
        plot(coordX,coordY,'.', 'MarkerSize', 20);
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
        plot(pgLU1,'FaceColor','green','FaceAlpha',0.5)
        %         plot(pgGapNew,'FaceColor','green','FaceAlpha',0.8)
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
            plot(pgLU2Rota,'FaceColor','red','FaceAlpha',0.5)
            hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
            plot(pgLU2,'FaceColor','red','FaceAlpha',0.5)
            if flagLU || flagLURota
                axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime*6);  hold off;
            else                
        axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime*1.5);  hold off;
            end
        
%              plot(pgGapNew,'FaceColor','green','FaceAlpha',0.8)
%             hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
        end
        
        if flagLU || flagLURota
            if ISplotShowGapAdjust
                % plot 仅调整的vertex的过程
                plot(pgVEH,'FaceColor','white','FaceAlpha',0.01);
                hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
                plot(pgLU,'FaceColor','green','FaceAlpha',0.2)
                hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
                plot(pgGap,'FaceColor','blue','FaceAlpha',0.2)
                hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);     pause(pausetime/100);
                plot(coordGapArray(:,1),coordGapArray(:,2),'.', 'MarkerSize', 8);
                hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
                plot(coordX,coordY,'.', 'MarkerSize', 20);
                hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
                plot(pgLU1,'FaceColor','green','FaceAlpha',0.5)
                hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
                if flagLURota
                    plot(pgLU2Rota,'FaceColor','red','FaceAlpha',0.5)
                else
                    plot(pgLU2,'FaceColor','red','FaceAlpha',0.5)
                end
                axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime*6);  hold off;
            end
            fprintf(1,'       Exsiting 混装间隙 in HBinGapAdjust (do2/do3)...\n');    
            % fLU: 需要调整的LU标记  从botttomLU计算, 但需返回到LU替换, 即非底部LU也要调整
            fLU = LU.CoordLUBin(:,1) == thisLU.CoordLUBin(1) & LU.CoordLUBin(:,2) == thisLU.CoordLUBin(2);
            
            LU.CoordLUBin(fLU,1) = coordX + LU.margin(idxLU,1);
            LU.CoordLUBin(fLU,2) = coordY + LU.margin(idxLU,4);
            if flagLURota
                tmpLWH = LU.LWH(fLU,1);
                LU.LWH(fLU,1) = LU.LWH(fLU,2);
                LU.LWH(fLU,2) = tmpLWH;
                LU.Rotaed(fLU)=~LU.Rotaed(fLU);
            end
            flagGap=1;            
        end
        if flagGap,   break;    end
    end
    if flagGap,   break;    end
end

% 如果本次调整Gap成功，则需要递归, 重新对该bin进行调整
if flagGap && parMulipleGap
     [LU] = HGapAdjust(LU,VEH);
end

% 如果调整Gap不成功, 则做些后处理
if ISplotShowGapAdjust
    pausetime = 0.0;
    % plot 某些vertex的尝试过程
    plot(pgVEH,'FaceColor','white','FaceAlpha',0.01);
    hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
    plot(pgLU,'FaceColor','green','FaceAlpha',0.2)
    hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
    plot(pgGap,'FaceColor','blue','FaceAlpha',0.2)
    hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
    plot(coordGapArray(:,1),coordGapArray(:,2),'.', 'MarkerSize', 8);
    hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime*1.5);
    hold on;
end

% 190108 优化
%  1 ：去除整层间隔间隙
gapY=sort(pgGap.Vertices(:,2));
gapY=unique(gapY);
% 当多个pgon时,会出现NaN值, 需要排除
gapY = gapY(~isnan(gapY));
for g=1:length(gapY)-1
    pgRect=polyshape([0 gapY(g);  VEH.LWH(1,1) gapY(g); VEH.LWH(1,1) gapY(g+1); 0 gapY(g+1)]); % 矩阵
    
    if ISplotShowGapAdjust
        plot(pgRect,'FaceColor','red','FaceAlpha',0.2)
        axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime*1.5);  hold on;
    end
    
    % 如果存在整层间隔且后面有堆垛，就往前移动
    if all(isinterior(pgGap, getBoundaryCentroid(pgRect)))
        warning('车内存在整层间隔');
        fYLU = LU.CoordLUBin(:,2)>gapY(g);
        if any(fYLU)
            LU.CoordLUBin(fYLU,2) =  LU.CoordLUBin(fYLU,2) - (gapY(g+1) - gapY(g));  %向前移动一定距离
        end
    end
end

if ISplotShowGapAdjust
    hold off;
end


% n1=fliplr(sortrows(n))

% LU.CoordLUBin(:,1)
% LU.CoordLUBin(:,2)
% sLU = sortrows(LU,{'LU_Bin'},{'ascend'})
% 1

end