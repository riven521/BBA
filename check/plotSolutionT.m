function [] = plotSolutionT(T,V,plotLU,plotItem,plotStrip,plotBin,plotColor,figname,luord,itemord,stripord)
% plotSolutionT ==> 作图LU/STRIP/BIN
% T: table 格式的LU
% V: table 格式的Veh
% plotLU：非0：做LU的图；2：按T.index顺序（输入顺序）作图 3:按T.order（排序后顺序）作图 1：按T的目前顺序作图；
% plotItem：非0：做Item的图； 1：按Item排序后作图
% plotStrip：非0：做Strip的图；1：按Strip的顺序作图
% plotBin：非0：做Bin的图；1：
% example：plotSolutionT(d.LU,d.Veh,1,1,1,1)
%                  plotSolutionT(d.LU,d.Veh,2,0,0,3)

if isstruct(T),  T = struct2table(structfun(@(x) x',T,'UniformOutput',false)); end
if isstruct(V), V = struct2table(structfun(@(x) x',V,'UniformOutput',false)); end

if nargin < 3
    plotLU = 1;  % 1：采用目前T内顺序 2: 采用index顺序
    plotStrip = 1;
    plotBin = 1;
    plotItem = 1;
    
end
if nargin < 8
    plotColor = 3;
    figname ='';
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

ISplotShowType = plotColor;
if ISplotShowType == 1
    tmpT = unique(T(:,{'LID'})); 
elseif ISplotShowType == 2
    tmpT = unique(T(:,{'PID'}));
elseif ISplotShowType == 3
    tmpT = unique(T(:,{'ID'}));
elseif ISplotShowType == 4
    tmpT = unique(T(:,{'SID'}));
elseif ISplotShowType == 5
    tmpT = unique(T(:,{'EID'}));
elseif ISplotShowType == 6
    tmpT = unique(T(:,{'isNonMixed'}));
elseif ISplotShowType == 7
    tmpT = unique(T(:,{'isMixTile'}));
elseif ISplotShowType == 8
    tmpT = unique(T(:,{'isShuaiWei'}));   
end

tmpT.LUcolor = 0.8*hsv(height(tmpT));
T = join(T,tmpT);


%% 图1：作图LU,按给定LU顺序
if plotLU
    figure('name',strjoin({figname,'LU展示：入ITEM时的排序后，合计*个',num2str(height(T))}));
    
    if plotLU==1, subT = T(1:end,:); end                          % LU目前顺序
    if plotLU==2, subT = sortrows(T,'Index'); end           % LU的BBA输入顺序
    
    if plotLU==3, subT = T(T.order,:); end                       %LU排序后顺序  %     if plotLU==3, subT = sortrows(T,'order'); end      %LU排序后顺序
    
    % 5555 构建LU的坐标系Coord X坐标依据LU的宽度L移动 Y为0 Z为LU高度
    XYZ = zeros(height(subT),3);                
    XYZ(:,1) = cumsum(subT.LWH(:,1));
    subT.Coord =[0 0 0; XYZ(1:end-1,:)];
    for iLU=1:height(subT)
        plotcube(subT.LWH(iLU,:), subT.Coord(iLU,:),0.7, subT.LUcolor(iLU,:));
        axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        view(60,40); %view(111,33);  
        if ISplotPause>0 ,      pause(ISplotPause);   end
    end
end

% 实质还是plotLU，同一Item依据顺序，给不同的高度Z值.
if plotItem
    XYZ = zeros(height(T),3);
    
    % 逐个Item作图
    nItem = max(T.LU_Item(:,1)); %Item的个数
    figure('name',strjoin({figname,'Item展示：合计*个',num2str([nItem])}));
    
    if plotItem==1, order = 1:nItem;  end       % 默认顺序作图 从1到n
    if plotItem==2, order = itemord;       end       % 给定顺序作图
    
    % 计算每个Lu在Item的坐标Coord
    sumL = 0;
    for iStrip=1:nItem
        % 从第1个Item开始，计算Lu的x坐标(相同堆垛是一样）
        %fidxItem = T.LU_Item(:,1)==iItem;  %同一堆垛下的逻辑值(多个）
        fidxStrip = T.LU_Item(:,1)==order(iStrip);  %同一堆垛下的逻辑值(多个）
        XYZ(fidxStrip,1) = sumL;
        sumL = sumL + unique(T.LWH(fidxStrip,1));       if numel(sumL) >  1 , error('e'); end
        
        nLU = sum(fidxStrip); % (同一堆垛下的托盘数量)
        sumH = 0;
        for iLU = 1:nLU 
            % 计算Lu的z坐标(相同堆垛高度不一样） （仅1个）
            fidxLU = T.LU_Item(:,2)==iLU;
            XYZ(fidxStrip&fidxLU,3) = sumH;
            sumH = sumH + T.LWH(fidxStrip&fidxLU,3);
        end
    end
    
    % 依据顺序堆垛
    T.XYZ = XYZ;
    T = sortrows(T,{'XYZ'},{'ascend'}); 
    
    for iLU=1:height(T)
        plotcube(T.LWH(iLU,:), T.XYZ(iLU,:),0.7, T.LUcolor(iLU,:));
        axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        view(60,40); %view(111,33);
        if ISplotPause>0 ,      pause(ISplotPause);   end
    end
end

%% 图2：作图Strip划分 v2（有LU_Strip(CoordLUStrip)即可）
if plotStrip
    XYZ = zeros(height(T),3);
    
    % 逐个Strip作图
    nstrip = max(T.LU_Strip(:,1)); %Item的个数
    figure('name',strjoin({figname,'STRIP展示：先后顺序排序后，合计*个',num2str([ nstrip ])}));
        
    if plotStrip==1, order = 1:nstrip;   end    % 默认顺序作图 从1到n                  
    if plotStrip==2, order = stripord;       end     % 给定顺序作图
    
    % 计算每个Lu在Strip的坐标Coord （X值不用外部Coord，因为顺序变化，自己给定顺序）
    sumL = 0;
    for iStrip=1:nstrip
        % 从第1个Strip开始，计算Lu的x坐标(相同堆垛是一样）
        fidxStrip = T.LU_Strip(:,1)==order(iStrip);  %同一strip下的逻辑值(多个）
        XYZ(fidxStrip,2) = sumL;  % Y = LWH'S L
        sumL = sumL + unique(T.LWH(fidxStrip,2));       if numel(sumL) >  1 , error('e'); end
        
        XYZ(fidxStrip,1) = T.CoordLUStrip(fidxStrip,1);  %
        XYZ(fidxStrip,3) = T.CoordLUStrip(fidxStrip,3);
    end
    
    % 依据顺序堆垛    
    T.XYZ = XYZ;
    T.YXZ(:,1) = T.XYZ(:,2);
    T.YXZ(:,2) = T.XYZ(:,1);  
    T = sortrows(T,{'YXZ'},{'ascend'}); 
    for iLU=1:height(T)
        plotcube(T.LWH(iLU,:), T.XYZ(iLU,:),0.7, T.LUcolor(iLU,:));
        axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        view(60,40); %view(111,33);
        if ISplotPause>0 ,      pause(ISplotPause);   end
    end    
    
%     for iLU=1:height(subT)
%         plotcube(subT.LWH(iLU,:),subT.CoordLUStrip(iLU,:),0.7, subT.LUcolor(iLU,:));
%         axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
%         view(111,33);    %view(60,40);
%         if ISplotPause>0 ,      pause(ISplotPause);   end
%     end
end



%% 图3：作图Bin划分 (有LU_Bin即可)
if plotBin
pos = 1;
if ismember('LU_Bin', T.Properties.VariableNames)
    
    T=sortrows(T,{'LU_Bin'},{'ascend'});  % T的排序 (LU_Bin递增)
    
    % 逐个bin作图
    nBin = max(T.LU_Bin(:,1)); %bin的个数
    figure('name',strjoin({figname,'BIN展示：先后顺序排序后，合计*个',num2str([nBin])}));
    
    for ibin=1:nBin % for ibin=nBin:nBin
        
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



%% 图2：作图Strip划分 （有LU_Strip(CoordLUStrip)即可）
% if plotStrip
% %     if ismember('LU_Strip', T.Properties.VariableNames)
% 
%     if plotStrip==1, T=sortrows(T,{'LU_Strip'},{'ascend'});  end                    
% %     if plotStrip==2, T=sortrows(T,{'LU_Strip'},{'ascend'});  end
%             
%         figure('name',strjoin({figname,'STRIP展示：先后顺序排序后，合计*个',num2str([ max(T.LU_Strip(:,1));])}));
%         subT = T(1:end,:);  % subT = T;
%         
%         for iLU=1:height(subT)
%             plotcube(subT.LWH(iLU,:),subT.CoordLUStrip(iLU,:),0.7, subT.LUcolor(iLU,:));
%             axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
%             view(111,33);    %view(60,40);
%             if ISplotPause>0 ,      pause(ISplotPause);   end
%         end
% %     end
% end

