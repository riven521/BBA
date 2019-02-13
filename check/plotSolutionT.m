function [] = plotSolutionT(T,V,plotLU,plotStrip)
% plotSolutionT ==> 作图LU/STRIP/BIN
%   T: LU; V: Veh.
if isstruct(T),  T = struct2table(structfun(@(x) x',T,'UniformOutput',false)); end
if isstruct(V), V = struct2table(structfun(@(x) x',V,'UniformOutput',false)); end

if nargin < 3
    plotLU = 1;
    plotStrip = 1;
end

global ISplotPause ISplotShowType

%% 0 V和T 首次join 获取车辆相关数据
if ~ismember('LU_VehType', V.Properties.VariableNames)
V.LU_VehType = [1:height(V)]'; 
V = V(:,{'LU_VehType','LWH'}); 
V.Properties.VariableNames{'LWH'} = 'LWH_V'; end

if ~ismember('LU_VehType', T.Properties.VariableNames)
    T.LU_VehType = ones(height(T),1); end

if ~ismember('LWH_V', T.Properties.VariableNames)
    T = join(T,V,'Keys','LU_VehType');  end

%% 仅适用于ITEM的plot todo 后期需要再增加
% V1 仅适用于ITEM的plot 将ITEM视为LU -> plotSolutionT(sItem,Veh)  
% % if iscell(T.LID)
% %     T.LID = cell2mat(T.LID); end

%% 1 获取LUcolor 含 颜色 属性的T  LID/SID/OPID etc   % 作图所需 1 LU的bin序号; 2 LWH 3 COORDLUBIN 4 LID 排序
if ISplotShowType == 1
    if ismember('isShuaiWei', T.Properties.VariableNames)
        tmpT = unique(T(:,{'LID'})); %LID/ID/isNonMixed/isMixTile/isShuaiWei
    else
        tmpT = unique(T(:,{'LID'})); %LID/ID/isNonMixed/isMixTile/isShuaiWei
    end
elseif ISplotShowType == 2
    tmpT = unique(T(:,{'PID'}));
elseif ISplotShowType == 3
    tmpT = unique(T(:,{'ID'}));
elseif ISplotShowType == 4
    tmpT = unique(T(:,{'SID'}));
elseif ISplotShowType == 5
    tmpT = unique(T(:,{'EID'}));
end

tmpT.LUcolor = 0.8*hsv(height(tmpT));
T = join(T,tmpT);

%% 图1：作图LU,按给定LU顺序
if plotLU
    
    figure('name',strjoin({'LU展示：入ITEM时的排序后，合计*个',num2str(height(T))}));
    subT = T(1:end,:);

    % 5555 构建仅plotLU的LU的坐标系
    XYZ = zeros(height(subT),3);
    XYZ(:,1) = cumsum(subT.LWH(:,1));
    subT.Coord =[0 0 0; XYZ(1:end-1,:)];
    
    for iLU=1:height(subT)
        plotcube(subT.LWH(iLU,:),subT.Coord(iLU,:),0.7, subT.LUcolor(iLU,:));
        axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        view(60,40); %view(111,33);  
        if ISplotPause>0 ,      pause(ISplotPause);   end
    end
end

%% 图2：作图Strip划分 （有LU_Strip(CoordLUStrip)即可）
if plotStrip
    if ismember('LU_Strip', T.Properties.VariableNames)
        
        T=sortrows(T,{'LU_Strip'},{'ascend'});     % T的排序 (LU_Strip递增)
        
        figure('name',strjoin({'STRIP展示：先后顺序排序后，合计*个',num2str([ max(T.LU_Strip(:,1));])}));
        subT = T(1:end,:);  % subT = T;
        
        for iLU=1:height(subT)
            plotcube(subT.LWH(iLU,:),subT.CoordLUStrip(iLU,:),0.7, subT.LUcolor(iLU,:));
            axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
            view(111,33);    %view(60,40);
            if ISplotPause>0 ,      pause(ISplotPause);   end
        end
    end
end

%% 图3：作图Bin划分 (有LU_Bin即可)
pos = 1;
if ismember('LU_Bin', T.Properties.VariableNames)
    
    T=sortrows(T,{'LU_Bin'},{'ascend'});  % T的排序 (LU_Bin递增)
    
    % 逐个bin作图
    nBin = max(T.LU_Bin(:,1)); %bin的个数
    figure('name',strjoin({'BIN展示：先后顺序排序后，合计*个',num2str([nBin])}));
    
    for ibin=1:nBin % for ibin=nBin:nBin
            %     if ibin==5
            %         s =T;
            %     end
            
        subplot(2,ceil((nBin+1)/2),pos); pos = pos+1;
        
        subT = T(T.LU_Bin(:,1)==ibin,:);                  % yxz =  T.LWH_T(f,:); coord = T.CoordLUBin(f,:);  color = T.LUcolor(f,:);  yxzVeh = T.LWH_V(f,:);
        
        for iLU=1:height(subT)
            plotcube(subT.LWH(iLU,:),subT.CoordLUBin(iLU,:),0.7, subT.LUcolor(iLU,:));
            axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
            view(111,33);    %view(60,40);
            xlim([0 subT.LWH_V(iLU,1)]);  ylim([0 subT.LWH_V(iLU,2)]); zlim([0 subT.LWH_V(iLU,3)]); % 车辆的长宽高调整到合适的车型
            if ISplotPause>0 ,      pause(ISplotPause);   end
        end
        
    end
end

%% 注释
% % % % % 3 开始作图 - 按Strip划分
% % % if ismember('LU_Strip', T.Properties.VariableNames),
% % %     % 2 T的排序 (LU_Strip递增)
% % %     T=sortrows(T,{'LU_Strip'},{'ascend'});
% % %     if ismember('LU_Bin', T.Properties.VariableNames),  subplot(2,ceil((nBin+1)/2),pos); 
% % %     else  subplot(1,1,pos); end
% % %     subT= T;
% % %     for iLU=1:height(subT)
% % %         plotcube(subT.LWH(iLU,:),subT.CoordLUStrip(iLU,:),0.7, subT.LUcolor(iLU,:));
% % %         axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
% % %         view(111,33);    %view(60,40);
% % %         xlim([0 subT.LWH_V(iLU,1)]);  zlim([0 subT.LWH_V(iLU,3)]); % 车辆的长宽高调整到合适的车型
% % %         if ISplotPause>0 ,      pause(ISplotPause);   end
% % %     end
% % % end

end



%% v1
% % function [] = plotSolutionT1(T,V)
% % if isstruct(T),  T = struct2table(structfun(@(x) x',T,'UniformOutput',false)); end
% % if isstruct(V), V = struct2table(structfun(@(x) x',V,'UniformOutput',false)); end
% % 
% % global ISplotPause ISplotShowType
% % 
% % % 0 首次join 获取车辆相关数据
% % if ~ismember('LU_VehType', V.Properties.VariableNames),
% % V.LU_VehType = [1:height(V)]'; 
% % V = V(:,{'LU_VehType','LWH'}); 
% % V.Properties.VariableNames{'LWH'} = 'LWH_V'; end
% % if ~ismember('LU_VehType', T.Properties.VariableNames),  T.LU_VehType = ones(height(T),1); end
% % if ~ismember('LWH_V', T.Properties.VariableNames),  T = join(T,V,'Keys','LU_VehType');  end
% % 
% % % 1 获取含 颜色 属性的T  LID/SID/OPID etc   % 作图所需 1 LU的bin序号; 2 LWH 3 COORDLUBIN 4 LID 排序
% % if ISplotShowType == 1
% %     tmpT = unique(T(:,{'LID'})); %LID/ID/isNonMixed/isMixTile
% % elseif ISplotShowType == 2
% %     tmpT = unique(T(:,{'PID'}));
% % elseif ISplotShowType == 3
% %     tmpT = unique(T(:,{'ID'}));
% % end
% % 
% % tmpT.LUcolor = 0.8*hsv(height(tmpT));
% % T = join(T,tmpT);
% % 
% % 
% % % 仅作图LU,按给定顺序作图
% % plotLU = 0;
% % if plotLU
% %     subT = T(1:end,:);
% % 
% %     % 5555 构建仅plotLU的LU的坐标系
% %     XYZ = zeros(height(subT),3);
% %     XYZ(:,1) = cumsum(subT.LWH(:,1));
% %     subT.Coord =[0 0 0; XYZ(1:end-1,:)];
% %     
% %     for iLU=1:height(subT)
% %         plotcube(subT.LWH(iLU,:),subT.Coord(iLU,:),0.7, subT.LUcolor(iLU,:));
% %         axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
% %         view(60,40); %view(111,33);  
% %     end
% % end
% % 
% % %%
% % 
% % %%
% % 
% % pos = 1;
% % % 2 开始作图 - 按Bin划分
% % if ismember('LU_Bin', T.Properties.VariableNames), 
% % % 2 T的排序 (LU_Bin递增)
% % T=sortrows(T,{'LU_Bin'},{'ascend'});    
% % % 逐个bin作图
% % nBin = max(T.LU_Bin(:,1)); %bin的个数
% % figure('name',num2str([nBin]));
% % 
% % for ibin=1:nBin
% %     subplot(2,ceil((nBin+1)/2),pos); pos = pos+1;    
% %     subT = T(T.LU_Bin(:,1)==ibin,:); % yxz =  T.LWH_T(f,:); coord = T.CoordLUBin(f,:);  color = T.LUcolor(f,:);  yxzVeh = T.LWH_V(f,:);
% %     for iLU=1:height(subT)
% %         plotcube(subT.LWH(iLU,:),subT.CoordLUBin(iLU,:),0.7, subT.LUcolor(iLU,:));        
% %         axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
% %         view(111,33);    %view(60,40);      
% %         xlim([0 subT.LWH_V(iLU,1)]);  ylim([0 subT.LWH_V(iLU,2)]); zlim([0 subT.LWH_V(iLU,3)]); % 车辆的长宽高调整到合适的车型
% %                 if ISplotPause>0 ,      pause(ISplotPause);   end
% %     end
% % end
% % end
% % 
% % % 3 开始作图 - 按Strip划分
% % if ismember('LU_Strip', T.Properties.VariableNames),
% %     % 2 T的排序 (LU_Strip递增)
% %     T=sortrows(T,{'LU_Strip'},{'ascend'});
% %     if ismember('LU_Bin', T.Properties.VariableNames),  subplot(2,ceil((nBin+1)/2),pos); 
% %     else  subplot(1,1,pos); end
% %     subT= T;
% %     for iLU=1:height(subT)
% %         plotcube(subT.LWH(iLU,:),subT.CoordLUStrip(iLU,:),0.7, subT.LUcolor(iLU,:));
% %         axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
% %         view(111,33);    %view(60,40);
% %         xlim([0 subT.LWH_V(iLU,1)]);  zlim([0 subT.LWH_V(iLU,3)]); % 车辆的长宽高调整到合适的车型
% %         if ISplotPause>0 ,      pause(ISplotPause);   end
% %     end
% % end
% % 
% % end

%%
% j = 1;
%         % 当托盘不相连出现时, 中断
%         if j~=1 && tmp ~=m(8,j) 
%                  if ISplotPause>0 ,      pause(ISplotPause);   end
%         end        
%         tmp = m(8,j);
%         j=j+1;

