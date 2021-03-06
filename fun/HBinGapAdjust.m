function [LU] = HBinGapAdjust(LU,VEH)
% HBinGapAdjust ==> Gap Adjust for Each Bin
global ISplotShowGapAdjust

TLU = getTableLU(LU);
TVEH = getTableLU(VEH);

typeVeh = unique(TLU.LU_Bin(:,1));
numVeh = length(typeVeh);      
if numVeh>1, error('多个车辆需要调整Gap'); end

% 循环每个VEH，找到对应子托盘subTLU和子车辆subVeh
for idxVeh = 1:numVeh
    subTLU = TLU(TLU.LU_Bin(:,1) == typeVeh(idxVeh), : );    
    subVeh = TVEH(unique(subTLU.LU_VehType), :); %TODO CHECK 是否LU_VehType就是车辆排序的第一个?
    
    if height(subVeh)>1,  error('NOT LU in the same Veh type'); end
        
    % 555  子托盘调整重要函数
   
    subTLUNew = HGapAdjust(subTLU,subVeh);
    
    
    %% V2   稳定性判定   
    
    isStablility = chkHorizontalStability(subTLUNew,subVeh);
    
    if ~isStablility,         subTLUNew=subTLU;     end
    
    if ISplotShowGapAdjust
        plotSolutionT(subTLUNew,subVeh,0,0,0,1,3,'subTLUNew');
        plotSolutionT(subTLU,subVeh,0,0,0,1,3,'subTLU');
    end
    

%% V1    稳定性判定(canbe removed)
% % % % %%  待增加：对比subTLU和subTLUNew的gap调整前后; 找一个指标表明gap调整后是可行的，不能以strip为维度了，难，暂时没有解决。
% % % %     % 设置条件：某些调整后状态不佳时，整个gap调整都不做。即gap调整单个堆垛允许,但最终结果不允许, 改回原LU
% % % %     pgVEH12 = polyshape(pgRectangle(0,0,subVeh.LWH(1,1),subVeh.LWH(1,2)));
% % % %     
% % % %     bottomLU1 = subTLU(subTLU.CoordLUBin(:,3)==0, : );  % 底层的托盘 % plotSolutionT(LU,VEH);
% % % %     pgLU1 = getPgLU(bottomLU1,0);  % plot(pgLU);    % pgLU1: 俯视 polygon
% % % %     pgGap1 = subtract(pgVEH12,pgLU1);    % pgGap1：俯视gap间隙
% % % %     
% % % %     bottomLU2 = subTLUNew(subTLUNew.CoordLUBin(:,3)==0, : );  % 底层的托盘 % plotSolutionT(LU,VEH);
% % % %     pgLU2 = getPgLU(bottomLU2,0);  % plot(pgLU);
% % % %     pgGap2 = subtract(pgVEH12,pgLU2);    
% % % % 
% % % %     % 对比pgLU1/2 or pgGap1/2
% % % %     [x1,y1]=centroid(pgLU1);
% % % %     [x2,y2]=centroid(pgLU2);
% % % %     % 2 Get CoordGap : pgGap 的 坐标 先X后Y 
% % % %     [XBound,YBound]=boundary(pgGap1);
% % % %     coordGapArray1 = fliplr(sortrows([YBound,XBound]));  % 先Y后X排序, 继而左右改变顺序
% % % %     [XBound,YBound]=boundary(pgGap2);
% % % %     coordGapArray2 = fliplr(sortrows([YBound,XBound]));  % 先Y后X排序, 继而左右改变顺序
% % % % 
% % % %     % 对比最大的Y值 小好 数量，少好 数量一致 最大Y
% % % %     % 如何不要呢？hole的多少，hole的大小
% % % %     %  依据每个LU的Coord与Veh的边界; 或依据Gap形状   
% % % %     if size(coordGapArray2,1)>size(coordGapArray1,1) % (st1: ) gap或lu的分布
% % % %         size(coordGapArray1,1)
% % % %         size(coordGapArray2,1)
% % % %         x2,x1
% % % %         
% % % %         area(pgLU1)
% % % %         area(pgLU2)
% % % %         subTLUNew = subTLU;
% % % %     end


    TLU(TLU.LU_Bin(:,1) == typeVeh(idxVeh), : ) = subTLUNew;

end

[LU] = getSturctT(TLU);

end

%% 555  子托盘调整重要函数
% 均在同一bin内的LU, VEH是一个车
function [LU] = HGapAdjust(LU,VEH)
global parMulipleGap ISplotShowGapAdjust ISplotGapCompare ISplotShowType

% 1 INITILIZE
flagGap=0; % 1: 调整 0: 未调整

% 2 Get bottomLU
bottomLU = LU(LU.CoordLUBin(:,3)==0, : );  % 底层的托盘 % plotSolutionT(LU,VEH);

% 2 Get pgGap: VEH内的剩余空白区域  
pgLU = getPgLU(bottomLU,0);  % plot(pgLU);
pgVEH = polyshape(pgRectangle(0,0,VEH.LWH(1,1),VEH.LWH(1,2)));	    
pgGap = subtract(pgVEH,pgLU);    if pgGap.NumRegions  > 1,  warning('Exsit %d Regions in this pgon', pgGap.NumRegions);  end
                        % TODO pgGap排序? ,gBlanks = sortregions(pgBlanks,'area','ascend');

% 2 Get CoordGap : pgGap 的 坐标 先X后Y 
[XBound,YBound]=boundary(pgGap);
coordGapArray = fliplr(sortrows([YBound,XBound]));  % 先Y后X排序, 继而左右改变顺序

% 3 循环每个Boundary顶点
for iVertex=1:size(coordGapArray,1)
    % Get Each GapBound Coord X and Y
    % 每个边界顶点的坐标值
    coordX = coordGapArray(iVertex,1);
    coordY = coordGapArray(iVertex,2);
            
    % 获取: 满足条件托盘序号 Get idxLUs i.e. Y of LU is larger than coordY of bound vertex 
%     flagLeqLU = bottomLU.CoordLUBin(:,2) - bottomLU.margin(:,3) > coordY;
%     flagLeqLU = bottomLU.CoordLUBin(:,2) - bottomLU.margin(:,3) >= coordY;
    
    flagLeqLU = (bottomLU.CoordLUBin(:,2) - bottomLU.margin(:,3)) >= coordY;
    flagEqLU = abs((bottomLU.CoordLUBin(:,2) - bottomLU.margin(:,3)) - coordY) < 0.1;  % NEW ADD 190420    
%     a = abs((bottomLU.CoordLUBin(:,2) - bottomLU.margin(:,3)) - coordY) < 0.1
    
    idxLUs = find(flagLeqLU);         % 找出比yPoint值 高的且是bottom的LU索引号    
    
    % 如果没有该类型托盘，继续.
    if isempty(idxLUs), continue; end       
           
    % 节约时间: 如果该Boundary的顶点不能放下最小的LU边(即超出车辆的边界), 则跳出该Boundary到下一个
    minLW = min(min(LU.LWH(idxLUs,[1,2])));
    if ~isinterior(pgVEH,coordX+minLW,coordY) || ~isinterior(pgVEH,coordX,coordY+minLW) || ~isinterior(pgVEH,coordX+minLW/2, coordY+minLW/2),    
        continue;   
    end
    
    % 循环每个Boundary顶点下每个符合条件的LU    
    %   排序: 托盘排序，
%     [~,ff] = sortrows([bottomLU.CoordLUBin(idxLUs,[2,1])],[1,2],{'ascend' 'descend'} );  % V2 先按y递增; 再按x递减
    [~,ff] = sortrows([bottomLU.CoordLUBin(idxLUs,[2,1])],[1,2],{'ascend' 'ascend'} ); % V1
    idxLUs = idxLUs(ff);
    for ii=1:length(idxLUs)
        idxLU = idxLUs(ii); % bottomLU的序号
        thisLU = bottomLU(idxLU,:);
        
        % pgLU1 待合并的LU pgGapNew 合并后的新Gap
        pgLU1 = getPgLU(thisLU, 0);        % pgLU1-待移动堆垛(不含margin)
        
        % pgGapNew pgLU1拿走，原Gap面积增加
        pgGapNew = union(pgGap,pgLU1);
        
        % pgLU2 待新摆放的LU位置           % pgLU2-待移动堆垛位置(含margin,图示不会存在空隙)
        pgLU2 = getPgLU(thisLU, 0, coordX,coordY);        %         pgLU2 =  polyshape(pgRectangle(coordX,coordY,wLU,lLU));  
        
        % 判定flagLU：是否pgLU2可以放置
        flagLU = all(isinterior(pgGapNew, getBoundaryCentroid(pgLU2)));        
        
        % 特殊情况1flagLU：是否pgLU2可以放置
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
        
        % 特殊情况2flagLU/flagLURota：是否pgLU2可以放置, 在同Y值移动情况; 2.1 仅允许沿X变小移/左移; 2.2 左移不允许旋转移
        % 190420 NEW ADD 如果相同Y：flagEqLU 且挪动增加了X值 则flagLU/flagLURota = 0 即不允许右移
                                            if thisLU.CoordLUBin(1) - thisLU.margin(1) <0, error('The margin minus is wrong'); end %放置margin减错
        if flagEqLU(idxLU)  % 如果同Y值移动
            % coordX: 新X位置; thisLU.CoordLUBin(1) - thisLU.margin(1): 原X位置
            % if coordX>thisLU.CoordLUBin(1) - thisLU.margin(1)  % 如果同Y值移动 如右移 不允许
            % if coordX>=(thisLU.CoordLUBin(1) - thisLU.margin(1))  % 如果同Y值移动 如右移或不动 均不允许
            if coordX - (thisLU.CoordLUBin(1) - thisLU.margin(1)) >= -0.1 % 如果同Y值移动 如右移或不动 均不允许
            flagLU = 0;
            flagLURota = 0;
            else  % 如果同Y值移动 如左移 不允许旋转
                flagLURota = 0;
            end
        end
        
        
        % 特殊情况3：isShuaiWei的thisLU不调整（不可以，甩尾的也有部分可以填充进来，使得内部充实，）
% %         if thisLU.isShuaiWei
% % %             flagLU = 0;         flagLURota = 0;
% %         end
        

        %         % 特殊情况4：Y值前进不足1/3不给推  % 改变不了特殊
%         if flagLU || flagLURota 
% %             if ~flagEqLU(idxLU)
%             if (thisLU.CoordLUBin(2) - thisLU.margin(3)) - coordY  <= 1/3*thisLU.LWH(2)
% %                 flagLU = 0;   flagLURota = 0;
%             end
% %             end
%         end
            
        % 特殊情况均排除后: 若都可放下, 二者选择一个，横放优先 ??? 确定更好??? todo
        if flagLU && flagLURota
            if thisLU.LWH(1) > thisLU.LWH(2) %wLU>lLU
                flagLURota=0;
            else
                flagLU=0;
            end
        end
            
        % 1此处，plot每个可能的gap顶点,即使可能转换fail,当然转换成功也在此
        if ISplotShowGapAdjust            
            plotGapPgon(VEH,pgLU,pgGap,pgVEH,coordGapArray,coordX,coordY,pgLU1,pgLU2,pgLU2Rota,flagLU,flagLURota);
        end
        
      % 允许调整 -> 修改坐标 / (旋转: 长宽和旋转状态) -> 修改BINSEQ/LU_Bin
        if flagLU || flagLURota
            
            fprintf(1,'       Exsiting 混装间隙 in HBinGapAdjust (do2/do3)...\n');    
            
            % fLU: 需要调整的LU标记  从botttomLU计算, 但需返回到LU替换, 即非底部LU也要调整
            fLU = LU.CoordLUBin(:,1) == thisLU.CoordLUBin(1) & LU.CoordLUBin(:,2) == thisLU.CoordLUBin(2);
            
            % bug fix : 基于LU的magin改为基于bottomLU的margin. 因为idxLU为bottomLU的索引号
            LU.CoordLUBin(fLU,1) = coordX + bottomLU.margin(idxLU,1);%LU.CoordLUBin(fLU,1) = coordX + LU.margin(idxLU,1);
            LU.CoordLUBin(fLU,2) = coordY + bottomLU.margin(idxLU,4);%LU.CoordLUBin(fLU,2) = coordY + LU.margin(idxLU,4);
            
            if flagLURota
                tmpLWH = LU.LWH(fLU,1);
                LU.LWH(fLU,1) = LU.LWH(fLU,2);
                LU.LWH(fLU,2) = tmpLWH;
                LU.Rotaed(fLU)=~LU.Rotaed(fLU);
            end
            
            % 190421 ADD 增加修改LU.LU_Bin -> LU为该bin内的托盘 不考虑序号
            [~,b2] = sortrows(LU.CoordLUBin,[2,1,3],{'ascend' 'ascend' 'ascend'} ); 
            [~,LU.LU_Bin(:,2)]=sort(b2);
            
            
            flagGap=1;            
            
            % gap调整过程展示
            % if ISplotGapCompare,      plotSolutionT(LU,VEH,0,0,0,1,ISplotShowType,'Gap调整成功后展示');         end
            
        end
        
        if flagGap,   break;    end
        
    end
    
    if flagGap,   break;    end
    
end

% 如果本次调整Gap成功，则需要递归, 重新对该bin进行调整
if flagGap && parMulipleGap
    
     [LU] = HGapAdjust(LU,VEH);  %替换成功，需要对LU和Gap的pg重计算
     
     bottomLU = LU(LU.CoordLUBin(:,3)==0, : );  % 底层的托盘 % bottomLU 必要重计算
     pgLU = getPgLU(bottomLU,0);  
     pgVEH = polyshape(pgRectangle(0,0,VEH.LWH(1,1),VEH.LWH(1,2)));
     pgGap = subtract(pgVEH,pgLU);    if pgGap.NumRegions  > 1,  warning('Exsit %d Regions in this pgon', pgGap.NumRegions);  end
     [XBound,YBound]=boundary(pgGap);
     coordGapArray = fliplr(sortrows([YBound,XBound]));  % 先Y后X排序, 继而左右改变顺序
    
end




% 190108 优化
%  1 ：去除整层间隔间隙
gapY=sort(pgGap.Vertices(:,2));
gapY=unique(gapY);
% 当多个pgon时,会出现NaN值, 需要排除
gapY = gapY(~isnan(gapY));

% if ISplotShowGapAdjust && length(gapY)>2
%     figure('name',strjoin({'Gap整层间隔调整过程展示：'}));
%     plotGapPgon(VEH,pgLU,pgGap,pgVEH,coordGapArray,coordX,coordY);
% end

for g=1:length(gapY)-1
    pgRect=polyshape([0 gapY(g);  VEH.LWH(1,1) gapY(g); VEH.LWH(1,1) gapY(g+1); 0 gapY(g+1)]); % 矩阵
    
    if ISplotShowGapAdjust
%         plotGapPgon(VEH,pgLU,pgGap,pgVEH);
                plotGapPgon(VEH,pgLU,pgGap,pgVEH,coordGapArray,coordX,coordY,[],[],[],[],[],pgRect);  %参数13-整层移动
    end
    
    % 如果存在整层间隔且后面有堆垛，就往前移动 bug：因为某些时候托盘lu在车型中间，不贴近边缘; 修复：增加overlaps判定，矩阵Rect和LU不能重叠
    if all(isinterior(pgGap, getBoundaryCentroid(pgRect))) && ~overlaps(pgLU,pgRect)%当pgRect的边界在pgGap中，返回true
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

% gap调整过程不展示，但依旧可能展示多次
if ISplotGapCompare && flagGap
    plotSolutionT(LU,VEH,0,0,0,1,ISplotShowType,'Gap调整成功后展示');         end


% n1=fliplr(sortrows(n))

% LU.CoordLUBin(:,1)
% LU.CoordLUBin(:,2)
% sLU = sortrows(LU,{'LU_Bin'},{'ascend'})
% 1

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

% 作混装间隙图函数
function plotGapPgon(VEH,pgLU,pgGap,pgVEH,coordGapArray,coordX,coordY,pgLU1,pgLU2,pgLU2Rota,flagLU,flagLURota,pgRect)       

        global ISplotPause;        
        % 1 pgVEH-矩形 pgLU-不规则俯视图 浅绿 pgGap-不规则俯视图 浅蓝
        plot(pgVEH,'FaceColor','white','FaceAlpha',0.01);
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
        plot(pgLU,'FaceColor','green','FaceAlpha',0.2)
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
        plot(pgGap,'FaceColor','blue','FaceAlpha',0.2)
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
        
        if nargin == 7        % 1 coordGapArray：所有不规则角点 coordX/coordX：当前角点
        plot(coordGapArray(:,1),coordGapArray(:,2),'.', 'MarkerSize', 8);
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
        plot(coordX,coordY,'.', 'MarkerSize', 20);
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
        end
        
        if nargin == 12
        % 作图不规则角点 
        plot(coordGapArray(:,1),coordGapArray(:,2),'.', 'MarkerSize', 8);
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
        plot(coordX,coordY,'.', 'MarkerSize', 20);
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
           
        % 作图待移动堆垛
        %  pgLU1 - 待移动堆垛-矩形-深绿  
        plot(pgLU1,'FaceColor','green','FaceAlpha',0.5);         %         plot(pgGapNew,'FaceColor','green','FaceAlpha',0.8)
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
        % pgLU2Rota/pgLU2 - 若可移动(外部判定), 则红色表示
        if flagLURota
            plot(pgLU2Rota,'FaceColor','red','FaceAlpha',0.5); axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]); 
        elseif flagLU
            plot(pgLU2,'FaceColor','red','FaceAlpha',0.5); axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]); 
        end                      


        if flagLU || flagLURota
            if ISplotPause>0 ,      pause(ISplotPause*6);  end
        else
            if ISplotPause>0 ,      pause(ISplotPause*1.5);   end
        end
                
        end
        
        if nargin <=12
            hold off;
        end
        
        if nargin == 13
            % 参数13-整层移动时启用
            plot(pgRect,'FaceColor','red','FaceAlpha',0.2)
            axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause*1.5);   end;  hold on;
        end   
        
        
end




%% V1
% % %% 555  子托盘调整重要函数
% % % 均在同一bin内的LU, VEH是一个车
% % function [LU] = HGapAdjust(LU,VEH)
% % global parMulipleGap ISplotShowGapAdjust
% % 
% % % 1 INITILIZE
% % flagGap=0; % 1: 调整 0: 未调整
% % 
% % % 2 Get bottomLU
% % bottomLU = LU(LU.CoordLUBin(:,3)==0, : );  % 底层的托盘 % plotSolutionT(LU,VEH);
% % 
% % % 2 Get pgGap: VEH内的剩余空白区域  
% % pgLU = getPgLU(bottomLU,0);  % plot(pgLU);
% % pgVEH = polyshape(pgRectangle(0,0,VEH.LWH(1,1),VEH.LWH(1,2)));	    
% % pgGap = subtract(pgVEH,pgLU);    if pgGap.NumRegions  > 1,  warning('Exsit %d Regions in this pgon', pgGap.NumRegions);  end
% %                         % TODO pgGap排序? ,gBlanks = sortregions(pgBlanks,'area','ascend');
% % 
% % % 2 Get CoordGap : pgGap 的 坐标 先X后Y 
% % [XBound,YBound]=boundary(pgGap);
% % coordGapArray = fliplr(sortrows([YBound,XBound]));  % 先Y后X排序, 继而左右改变顺序
% % 
% % % 3 循环每个Boundary顶点
% % for iVertex=1:size(coordGapArray,1)
% %     % Get Each GapBound Coord X and Y
% %     % 每个边界顶点的坐标值
% %     coordX = coordGapArray(iVertex,1);
% %     coordY = coordGapArray(iVertex,2);
% %             
% %     % 获取: 满足条件托盘序号 Get idxLUs i.e. Y of LU is larger than coordY of bound vertex 
% %     flagLeqLU = bottomLU.CoordLUBin(:,2) - bottomLU.margin(:,3) > coordY;
% %     idxLUs = find(flagLeqLU);         % 找出比yPoint值 高的且是bottom的LU索引号    
% %     
% %     % 如果没有该类型托盘，继续.
% %     if isempty(idxLUs), continue; end       
% %            
% %     % 节约时间: 如果该Boundary的顶点不能放下最小的LU边(即超出车辆的边界), 则跳出该Boundary到下一个
% %     minLW = min(min(LU.LWH(idxLUs,[1,2])));
% %     if ~isinterior(pgVEH,coordX+minLW,coordY) || ~isinterior(pgVEH,coordX,coordY+minLW) || ~isinterior(pgVEH,coordX+minLW/2, coordY+minLW/2),    continue;    end
% %     
% %     % 循环每个Boundary顶点下每个符合条件的LU    
% %     %   排序: 托盘排序
% %     [~,ff] = sortrows([bottomLU.CoordLUBin(idxLUs,[2,1])]);  idxLUs = idxLUs(ff);
% %     for ii=1:length(idxLUs)
% %         idxLU = idxLUs(ii); % bottomLU的序号
% %         thisLU = bottomLU(idxLU,:);
% %         
% %         % pgLU1 待合并的LU pgGapNew 合并后的新Gap
% %         pgLU1 = getPgLU(thisLU, 0);        
% %         pgGapNew = union(pgGap,pgLU1);
% %         
% %         % pgLU2 待新摆放的LU位置
% %         pgLU2 = getPgLU(thisLU, 0, coordX,coordY);        %         pgLU2 =  polyshape(pgRectangle(coordX,coordY,wLU,lLU));  
% %         flagLU = all(isinterior(pgGapNew, getBoundaryCentroid(pgLU2)));        
% %         
% %         pgIntersect = intersect(pgLU2,subtract(pgVEH,pgGapNew));
% %         if pgIntersect.NumRegions ~= 0, flagLU = 0; end
% %         
% %         % pgLU2Rota 待新摆放的旋转后LU位置
% %         flagLURota=0;
% %         if thisLU.isRota     %bottomLU.isRota(idxLU)
% %             pgLU2Rota =  getPgLU(thisLU, 1, coordX,coordY); %polyshape(pgRectangle(coordX,coordY,lLU,wLU));
% %             flagLURota = all(isinterior(pgGapNew, getBoundaryCentroid(pgLU2Rota))); 
% %             
% %             pgIntersectRota = intersect(pgLU2Rota,subtract(pgVEH,pgGapNew));
% %             if pgIntersectRota.NumRegions ~= 0, flagLURota = 0; end
% %         end        
% %         
% %         % 若都可放下, 二者选择一个，横放优先
% %         if flagLU && flagLURota
% %             if thisLU.LWH(1) > thisLU.LWH(2) %wLU>lLU
% %                 flagLURota=0;
% %             else
% %                 flagLU=0;
% %             end
% %         end
% %             
% %         if ISplotShowGapAdjust
% %             
% %         pausetime = 0.1;   
% %         % plot 某些vertex的尝试过程(明显不可能的点已经排除，若不想排除，注释65行节约时间语句)
% %         plot(pgVEH,'FaceColor','white','FaceAlpha',0.01);
% %         hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %         plot(pgLU,'FaceColor','green','FaceAlpha',0.2)
% %         hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %         plot(pgGap,'FaceColor','blue','FaceAlpha',0.2)
% %         hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %         plot(coordGapArray(:,1),coordGapArray(:,2),'.', 'MarkerSize', 8);
% %         hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %         plot(coordX,coordY,'.', 'MarkerSize', 20);
% %         hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %         plot(pgLU1,'FaceColor','green','FaceAlpha',0.5)
% %         %         plot(pgGapNew,'FaceColor','green','FaceAlpha',0.8)
% %         hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %             plot(pgLU2Rota,'FaceColor','red','FaceAlpha',0.5)
% %             hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %             plot(pgLU2,'FaceColor','red','FaceAlpha',0.5)
% %             if flagLU || flagLURota
% %                 axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime*6);  hold off;
% %             else                
% %                 axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime*1.5);  hold off;
% %             end
% %         
% % %              plot(pgGapNew,'FaceColor','green','FaceAlpha',0.8)
% % %             hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% % 
% % % figure('name',strjoin({'Gap展示：'}));
% % % plotGapPgon(VEH,pgLU,pgGap,pgVEH)%,coordGapArray,coordX,coordY,pgLU1,pgLU2,pgLU2Rota,flagLU,flagLURota,pgGapNew);
% % figure('name',strjoin({'Gap展示：'}));
% % plotGapPgon(VEH,pgLU,pgGap,pgVEH,coordGapArray,coordX,coordY,pgLU1,pgLU2,pgLU2Rota,flagLU,flagLURota);
% % % figure('name',strjoin({'Gap展示：'}));
% % % plotGapPgon(VEH,pgLU,pgGap,pgVEH,coordGapArray,coordX,coordY,pgLU1,pgLU2,pgLU2Rota,flagLU,flagLURota,pgGapNew);
% % 
% %         end
% %         
% %         if flagLU || flagLURota
% %             if ISplotShowGapAdjust
% %                 % plot 仅调整的vertex的过程（作图成功的顶点和调整）
% %                 plot(pgVEH,'FaceColor','white','FaceAlpha',0.01);
% %                 hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %                 plot(pgLU,'FaceColor','green','FaceAlpha',0.2)
% %                 hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %                 plot(pgGap,'FaceColor','blue','FaceAlpha',0.2)
% %                 hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);     pause(pausetime/100);
% %                 plot(coordGapArray(:,1),coordGapArray(:,2),'.', 'MarkerSize', 8);
% %                 hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %                 plot(coordX,coordY,'.', 'MarkerSize', 20);
% %                 hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %                 plot(pgLU1,'FaceColor','green','FaceAlpha',0.5)
% %                 hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %                 if flagLURota
% %                     plot(pgLU2Rota,'FaceColor','red','FaceAlpha',0.5)
% %                 else
% %                     plot(pgLU2,'FaceColor','red','FaceAlpha',0.5)
% %                 end
% %                 axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime*6);  hold off;
% %                 
% %                 plotSolutionT(LU,VEH,0,0,0,1,3,'ID作图Bin');                  %     Bin排序后
% %                 
% %             end
% %             fprintf(1,'       Exsiting 混装间隙 in HBinGapAdjust (do2/do3)...\n');    
% %             
% %             % fLU: 需要调整的LU标记  从botttomLU计算, 但需返回到LU替换, 即非底部LU也要调整
% %             fLU = LU.CoordLUBin(:,1) == thisLU.CoordLUBin(1) & LU.CoordLUBin(:,2) == thisLU.CoordLUBin(2);
% %             
% %             LU.CoordLUBin(fLU,1) = coordX + LU.margin(idxLU,1);
% %             LU.CoordLUBin(fLU,2) = coordY + LU.margin(idxLU,4);
% %             
% %             if flagLURota
% %                 tmpLWH = LU.LWH(fLU,1);
% %                 LU.LWH(fLU,1) = LU.LWH(fLU,2);
% %                 LU.LWH(fLU,2) = tmpLWH;
% %                 LU.Rotaed(fLU)=~LU.Rotaed(fLU);
% %             end
% %             
% %             flagGap=1;            
% %             
% %         end
% %         
% %         if flagGap,   break;    end
% %         
% %     end
% %     
% %     if flagGap,   break;    end
% %     
% % end
% % 
% % % 如果本次调整Gap成功，则需要递归, 重新对该bin进行调整
% % if flagGap && parMulipleGap
% %     
% %      [LU] = HGapAdjust(LU,VEH);  %替换成功，需要对LU和Gap的pg重计算
% %      
% %      bottomLU = LU(LU.CoordLUBin(:,3)==0, : );  % 底层的托盘 % bottomLU 必要重计算
% %      pgLU = getPgLU(bottomLU,0);  
% %      pgVEH = polyshape(pgRectangle(0,0,VEH.LWH(1,1),VEH.LWH(1,2)));
% %      pgGap = subtract(pgVEH,pgLU);    if pgGap.NumRegions  > 1,  warning('Exsit %d Regions in this pgon', pgGap.NumRegions);  end
% %      [XBound,YBound]=boundary(pgGap);
% %     coordGapArray = fliplr(sortrows([YBound,XBound]));  % 先Y后X排序, 继而左右改变顺序
% %     
% % end
% % 
% % % 如果调整Gap不成功, 则做些后处理
% % if ISplotShowGapAdjust
% %     pausetime = 0.0;
% %     % plot 某些vertex的尝试过程
% %     plot(pgVEH,'FaceColor','white','FaceAlpha',0.01);
% %     hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %     plot(pgLU,'FaceColor','green','FaceAlpha',0.2)
% %     hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %     plot(pgGap,'FaceColor','blue','FaceAlpha',0.2)
% %     hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %     plot(coordGapArray(:,1),coordGapArray(:,2),'.', 'MarkerSize', 8);
% %     hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime/100);
% %     axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime*1.5);
% %     hold on;
% % end
% % 
% % 
% % 
% % % 190108 优化
% % %  1 ：去除整层间隔间隙
% % gapY=sort(pgGap.Vertices(:,2));
% % gapY=unique(gapY);
% % % 当多个pgon时,会出现NaN值, 需要排除
% % gapY = gapY(~isnan(gapY));
% % for g=1:length(gapY)-1
% %     pgRect=polyshape([0 gapY(g);  VEH.LWH(1,1) gapY(g); VEH.LWH(1,1) gapY(g+1); 0 gapY(g+1)]); % 矩阵
% %     
% %     if ISplotShowGapAdjust
% %         plot(pgRect,'FaceColor','red','FaceAlpha',0.2)
% %         axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime*1.5);  hold on;
% %     end
% %     
% %     % 如果存在整层间隔且后面有堆垛，就往前移动 bug：因为某些时候托盘lu在车型中间，不贴近边缘; 修复：增加overlaps判定，矩阵Rect和LU不能重叠
% %     if all(isinterior(pgGap, getBoundaryCentroid(pgRect))) && ~overlaps(pgLU,pgRect)%当pgRect的边界在pgGap中，返回true
% %         warning('车内存在整层间隔');
% %         fYLU = LU.CoordLUBin(:,2)>gapY(g);
% %         if any(fYLU)
% %             LU.CoordLUBin(fYLU,2) =  LU.CoordLUBin(fYLU,2) - (gapY(g+1) - gapY(g));  %向前移动一定距离
% %         end
% %     end
% % end
% % 
% % if ISplotShowGapAdjust
% %     hold off;
% % end
% % 
% % 
% % % n1=fliplr(sortrows(n))
% % 
% % % LU.CoordLUBin(:,1)
% % % LU.CoordLUBin(:,2)
% % % sLU = sortrows(LU,{'LU_Bin'},{'ascend'})
% % % 1
% % 
% % end