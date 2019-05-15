%% chkStability Each bin 的四个方面的稳定性判定
%% Form
%    chkHorizontalStability
%    核验bin内是否水平稳定（主要考虑运输途中的稳定性）：即每个托盘的四周是否有相邻托盘或车厢wall支撑

% 相邻定义：前后左右四个方向，在LU的maxAdjDist（100）内的其它托盘与本托盘的交集，若有，则定义为该方向有相邻托盘；否则，无相邻托盘。
% NOTE: 
% ?	左侧和前测：必须有相邻；
% ?	右侧：若该托盘为甩尾托盘/或该托盘右侧任一托盘为甩尾托盘，即使非相邻。
% ?	后侧：无论是否有相邻托盘，均视为相邻。

% subTLU - 属于某个bin内的托盘
% subVeh - 某个车型
function isStablility = chkHorizontalStability(subTLU,subVeh)

        global ISplotShowGapAdjust
        isStablility = 1;
        maxAdjDist = 100; 
        
%         subTLU=sortrows(subTLU,{'LU_Bin'},{'ascend'}); 
        subTLU=sortrows(subTLU,{'CoordLUBin'},{'ascend'}); 
        
            subTLU.leftStab =   zeros(height(subTLU),1);
            subTLU.rightStab = zeros(height(subTLU),1);
            subTLU.frontStab = zeros(height(subTLU),1);
            subTLU.backStab =  zeros(height(subTLU),1);
            
        if istable(subTLU),   [subTLU] = getSturctT(subTLU); end
        if istable(subVeh),   [subVeh] = getSturctT(subVeh); end

        % 循环每个托盘, 进行判定是否stability
        for idxl=1:length(subTLU.ID)
            % 构建polygon, 车辆的宽,高 截面pgVEHLeft + 车辆的长,高 截面pgVEHFront
            pgVEHLeft = polyshape(pgRectangle(0,0,subVeh.LWH(1),subVeh.LWH(3)));
            pgVEHFront = polyshape(pgRectangle(0,0,subVeh.LWH(2),subVeh.LWH(3)));
            
            % 当前托盘的坐标值
            x=subTLU.CoordLUBin(1,idxl);
            y=subTLU.CoordLUBin(2,idxl);
            z=subTLU.CoordLUBin(3,idxl);
            
            % 当前托盘的宽长高 - 对应xyz
            w = subTLU.LWH(1,idxl);
            l = subTLU.LWH(2,idxl);
            d = subTLU.LWH(3,idxl);
            
            % 当前托盘的margin (左上右下)
            m1 = subTLU.margin(1,idxl);
            m2 = subTLU.margin(2,idxl);
            m3 = subTLU.margin(3,idxl);
            m4 = subTLU.margin(4,idxl);
            
            % 构建polygon, 本托盘的宽,高 截面pgLULeft + 本托盘的长,高 截面pgLUFront
            pgLULeft = polyshape(pgRectangle(x,z,w,d));
            pgLUFront = polyshape( pgRectangle(y,z,l,d));
            
            % 坐标值: 本托盘相邻托盘的最大坐标位置 (相邻的定义: 距离本托盘的坐标均在3个margin之内 或 30cm之内 (有时margin为0) )
            YleftLU =  ( y - max(3*m1,maxAdjDist));  % 当前LU左侧Y值
            XfrontLU =  ( x - max(3*m2,maxAdjDist));  % 当前LU前侧x值
            YrightLU =  ( y + l + max(3*m3,maxAdjDist));  % 当前LU右侧Y值
            XbackLU =  ( x + w + max(3*m4,maxAdjDist));  % 当前LU后侧x值
            
            % 坐标值: 所有托盘, 待确定是否相邻的托盘集合的对应坐标值
            Yright = subTLU.CoordLUBin(2,:)+subTLU.LWH(2,:);  % 所有LU的右侧Y坐标值 (Y值+LU的长) % 当前托盘在margin范围内的左侧的托盘（对应托盘的右侧）
            Xfront= subTLU.CoordLUBin(1,:);  % 所有LU的前方X坐标值 (X值+LU的宽)
            Yleft = subTLU.CoordLUBin(2,:);  % 所有LU的左侧Y坐标值 (Y值)
            Xback = subTLU.CoordLUBin(1,:)+subTLU.LWH(1,:);  % 所有LU的后方X坐标值 (X值)
            
            % 判定四个侧面是否均满足要求Horizontal的稳定性要求
            [flagLeftStability,flagRightStability,flagFrontStability,flagBackStability]  = deal(zeros(1,1));
            
            %% 1 判定左侧
            if YleftLU<=0  % 如左侧是车厢左边wall, 则稳定
                flagLeftStability = 1;
                pgLUSet =pgVEHLeft;
            else % 否则, 找出不重叠且位于最大左侧边缘坐标值内的LU集合pgLUSet(不包含自身)
                setidxLU = (Yright >= YleftLU) & (Yright <= y);  setidxLU(idxl) = 0;
                pgLUSet = getPgLUSet(setidxLU,subTLU,1);
                
%                 yc = y - Yright;                yc(yc<0) = Inf;
%                 [a,b] = min(yc);
            end
            % 本托盘和左侧托盘集的交集. 若不空, 则稳定; 若空, 则不稳定
            pgGap = intersect(pgLUSet,pgLULeft);     if ~isempty(pgGap.Vertices),      flagLeftStability = 1;       end
            
            if flagLeftStability==0 && ISplotShowGapAdjust
                plotStabilityPg(pgVEHLeft,pgLULeft,pgLUSet,pgGap,subVeh.LWH,1);    
                subTLUTTT = struct2table(structfun(@(x) x',subTLU,'UniformOutput',false));
                setidxLUTTT = setidxLU;            setidxLUTTT(idxl) = 1;
                plotSolutionT(subTLUTTT(setidxLUTTT,:),subVeh,0,0,0,1,3,'subTLUUUUUUUU'); 
                plotSolutionT(subTLUTTT,subVeh,0,0,0,1,3,'subTLUUUUUUUU'); 
            end
            

            %% 2 判定右侧 如甩尾托盘 不考虑稳定性; 如非甩尾托盘, 但其右侧相邻是甩尾托盘，也不考虑稳定性
            if YrightLU >= subVeh.LWH(2) % 如右侧是车厢右边wall, 则稳定
                flagRightStability = 1;
                pgLUSet =pgVEHLeft;
            else  % 否则, 找出不重叠且位于最大右侧边缘坐标值内的LU集合pgLUSet(不包含自身)
                setidxLU = (Yleft <= YrightLU) & (Yleft >= y + l);           setidxLU(idxl) = 0;
                pgLUSet = getPgLUSet(setidxLU,subTLU,1);
            end
            % 本托盘和右侧托盘集的交集. 若不空, 则稳定; 若空, 则不稳定
            pgGap = intersect(pgLUSet,pgLULeft);     if ~isempty(pgGap.Vertices),      flagRightStability = 1;       end
            if (flagRightStability==0 && subTLU.isShuaiWei(idxl)==1), flagRightStability=1; end % 特殊处理, right方向在甩尾LU不作为unstability.
            if flagRightStability==0 && any(subTLU.isShuaiWei(setidxLU)==1), flagRightStability=1; end % 特殊处理, right方向在相邻任一右侧是甩尾LU不作为unstability.
            

            if flagRightStability==0 && ISplotShowGapAdjust
                plotStabilityPg(pgVEHLeft,pgLULeft,pgLUSet,pgGap,subVeh.LWH,1);           
                subTLUTTT = struct2table(structfun(@(x) x',subTLU,'UniformOutput',false));
                setidxLUTTT = setidxLU;                 setidxLUTTT(idxl) = 1;
                plotSolutionT(subTLUTTT(setidxLUTTT,:),subVeh,0,0,0,1,3,'subTLUUUUUUUU');
                plotSolutionT(subTLUTTT,subVeh,0,0,0,1,3,'subTLUUUUUUUU');
            end
                        
            %% 3 判定前侧
            if XfrontLU<=0 % 如前侧是车厢前边wall, 则稳定
                flagFrontStability = 1;
                pgLUSet =pgVEHFront;
            else  % 否则, 找出不重叠且位于最大前侧边缘坐标值内的LU集合pgLUSet(不包含自身)
                setidxLU = (Xback >= XfrontLU) & (Xback <= x);            setidxLU(idxl) = 0;
                pgLUSet = getPgLUSet(setidxLU,subTLU,2);
            end
            % 本托盘和前侧托盘集的交集. 若不空, 则稳定; 若空, 则不稳定
            pgGap = intersect(pgLUSet,pgLUFront);             if ~isempty(pgGap.Vertices),      flagFrontStability = 1;       end 
            
           if flagFrontStability==0 && ISplotShowGapAdjust
               plotStabilityPg(pgVEHFront,pgLUFront,pgLUSet,pgGap,subVeh.LWH,2);  
               subTLUTTT = struct2table(structfun(@(x) x',subTLU,'UniformOutput',false));
               setidxLUTTT = setidxLU;                setidxLUTTT(idxl) = 1;
               plotSolutionT(subTLUTTT(setidxLUTTT,:),subVeh,0,0,0,1,3,'subTLUUUUUUUU');
               plotSolutionT(subTLUTTT,subVeh,0,0,0,1,3,'subTLUUUUUUUU');
           end
            
                       
            %% 4 判定后侧 Back面向自己的都默认为1,即始终处于back稳定         
            if XbackLU>= subVeh.LWH(1) % 如后侧是车厢后边wall, 则稳定
                flagBackStability = 1;
                pgLUSet =pgVEHFront;
            else  % 否则, 找出不重叠且位于最大后侧边缘坐标值内的LU集合pgLUSet(不包含自身)
                setidxLU = (Xfront <= XbackLU) & (Xfront >= x + w);                 setidxLU(idxl) = 0;
                pgLUSet = getPgLUSet(setidxLU,subTLU,2);
            end
            % 本托盘和后侧托盘集的交集. 若不空, 则稳定; 若空, 则不稳定
            pgGap = intersect(pgLUSet,pgLUFront);             if ~isempty(pgGap.Vertices),      flagBackStability = 1;       end
            if flagBackStability==0, flagBackStability=1; end % 特殊处理, back方向永远不作为unstability.
            
            if flagBackStability==0 && ISplotShowGapAdjust % 不可能发生
                plotStabilityPg(pgVEHFront,pgLUFront,pgLUSet,pgGap,subVeh.LWH,2);                         
            end

            
            %% 记录每个
            subTLU.leftStab(idxl) =   flagLeftStability;
            subTLU.rightStab(idxl) =   flagRightStability;
            subTLU.frontStab(idxl) =   flagFrontStability;
            subTLU.backStab(idxl) =   flagBackStability;
        end
        
        % 如有任何一个LU是不Stability. 则Gap调整失效
        if ~all(subTLU.leftStab) || ~all(subTLU.rightStab) || ~all(subTLU.frontStab) || ~all(subTLU.backStab)
            find(subTLU.leftStab==0);
            find(subTLU.rightStab==0);
            find(subTLU.frontStab==0);
            find(subTLU.backStab==0);
            isStablility = 0;
        end
        
end



% 作混装间隙图函数
function plotStabilityPg(pgVEH,pgLU,pgLUSet,pgGap,LWH,type)       
    if type==1
        a = LWH(1);      
    elseif type==2
        a = LWH(2);  
    else
        error('1');
    end
    plot(pgVEH,'FaceColor','white','FaceAlpha',0.01);         hold on;    axis equal;    grid on;    xlim([0 1.5*a]);    ylim([0 1.2*LWH(3)]);
    plot(pgLU,'FaceColor','green','FaceAlpha',0.2);           hold on;    axis equal;    grid on;    xlim([0 1.5*a]);    ylim([0 1.2*LWH(3)]);
    plot(pgLUSet,'FaceColor','red','FaceAlpha',0.2);            hold on;    axis equal;    grid on;    xlim([0 1.5*a]);    ylim([0 1.2*LWH(3)]);
    plot(pgGap,'FaceColor','red','FaceAlpha',0.9);            hold on;    axis equal;    grid on;    xlim([0 1.5*a]);    ylim([0 1.2*LWH(3)]);
    hold off; 

end
  

% 获取待判断的内点, 由Boundary和Centroid共同构成
function pgLUSet = getPgLUSet(setidxLU,subTLUNewBuff,type);
pgLUSet = polyshape();
if sum(setidxLU) == 0,    return; end
P = [];
setidxLU = find(setidxLU);
for idxl2=1:length(setidxLU)
    fidx2 = setidxLU(idxl2);
    % 采用固定坐标, 若有可变输入参数
    x1=subTLUNewBuff.CoordLUBin(1,fidx2);
    y1=subTLUNewBuff.CoordLUBin(2,fidx2);
    z1=subTLUNewBuff.CoordLUBin(3,fidx2);
    
    w1 = subTLUNewBuff.LWH(1,fidx2);
    l1 = subTLUNewBuff.LWH(2,fidx2);
    d1 = subTLUNewBuff.LWH(3,fidx2);
    
    if type == 1
        P =  [P;pgRectangle(x1,z1,w1,d1);[NaN,NaN]];
    elseif type==2
        P =  [P;pgRectangle(y1,z1,l1,d1);[NaN,NaN]];
    else
        error('1');
    end
end
pgLUSet = polyshape(P);
end