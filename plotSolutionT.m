function [] = plotSolutionT(T,V)
if isstruct(T),  T = struct2table(structfun(@(x) x',T,'UniformOutput',false)); end
if isstruct(V), V = struct2table(structfun(@(x) x',V,'UniformOutput',false)); end

global ISplotPause ISplotShowType

% 0 首次join 获取车辆相关数据
if ~ismember('LU_VehType', V.Properties.VariableNames),
V.LU_VehType = 1:height(V); V = V(:,{'LU_VehType','LWH'}); V.Properties.VariableNames{'LWH'} = 'LWH_V'; end
if ~ismember('LU_VehType', T.Properties.VariableNames),  T.LU_VehType = ones(height(T),1); end
if ~ismember('LWH_V', T.Properties.VariableNames),  T = join(T,V,'Keys','LU_VehType');  end

% 1 获取含 颜色 属性的T  LID/SID/OPID etc   % 作图所需 1 LU的bin序号; 2 LWH 3 COORDLUBIN 4 LID 排序
if ISplotShowType == 1
    tmpT = unique(T(:,{'ID'})); %LID/ID/isNonMixed/isMixTile
elseif ISplotShowType == 2
    tmpT = unique(T(:,{'PID'}));
end

tmpT.LUcolor = 0.8*hsv(height(tmpT));
T = join(T,tmpT);


% 仅作图LU,按给定顺序作图
plotLU = 1;
if plotLU
    subT = T(1:end,:);

    % 5555 构建仅plotLU的LU的坐标系
    XYZ = zeros(height(subT),3);
    XYZ(:,1) = cumsum(subT.LWH(:,1));
    subT.Coord =[0 0 0; XYZ(1:end-1,:)];
    
    for iLU=1:height(subT)
        plotcube(subT.LWH(iLU,:),subT.Coord(iLU,:),0.7, subT.LUcolor(iLU,:));
        axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        view(60,40); %view(111,33);  
    end
end

%%

%%

pos = 1;
% 2 开始作图 - 按Bin划分
if ismember('LU_Bin', T.Properties.VariableNames), 
% 2 T的排序 (LU_Bin递增)
T=sortrows(T,{'LU_Bin'},{'ascend'});    
% 逐个bin作图
nBin = max(T.LU_Bin(:,1)); %bin的个数
figure('name',num2str([nBin]));
for ibin=1:nBin
    subplot(2,ceil((nBin+1)/2),pos); pos = pos+1;    
    subT = T(T.LU_Bin(:,1)==ibin,:); % yxz =  T.LWH_T(f,:); coord = T.CoordLUBin(f,:);  color = T.LUcolor(f,:);  yxzVeh = T.LWH_V(f,:);
    for iLU=1:height(subT)
        plotcube(subT.LWH(iLU,:),subT.CoordLUBin(iLU,:),0.7, subT.LUcolor(iLU,:));        
        axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        view(111,33);    %view(60,40);      
        xlim([0 subT.LWH_V(iLU,1)]);  ylim([0 subT.LWH_V(iLU,2)]); zlim([0 subT.LWH_V(iLU,3)]); % 车辆的长宽高调整到合适的车型
                if ISplotPause>0 ,      pause(ISplotPause);   end
    end
end
end

% 3 开始作图 - 按Strip划分
if ismember('LU_Strip', T.Properties.VariableNames),
    % 2 T的排序 (LU_Strip递增)
    T=sortrows(T,{'LU_Strip'},{'ascend'});
    if ismember('LU_Bin', T.Properties.VariableNames),  subplot(2,ceil((nBin+1)/2),pos); 
    else  subplot(1,1,pos); end
    subT= T;
    for iLU=1:height(subT)
        plotcube(subT.LWH(iLU,:),subT.CoordLUStrip(iLU,:),0.7, subT.LUcolor(iLU,:));
        axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        view(111,33);    %view(60,40);
        xlim([0 subT.LWH_V(iLU,1)]);  zlim([0 subT.LWH_V(iLU,3)]); % 车辆的长宽高调整到合适的车型
        if ISplotPause>0 ,      pause(ISplotPause);   end
    end
end

end

%%
% j = 1;
%         % 当托盘不相连出现时, 中断
%         if j~=1 && tmp ~=m(8,j) 
%                  if ISplotPause>0 ,      pause(ISplotPause);   end
%         end        
%         tmp = m(8,j);
%         j=j+1;

