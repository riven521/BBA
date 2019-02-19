function [] = plotSolutionT(T,V,plotLU,plotItem,plotStrip,plotBin,plotColor,figname,luord,itemord,stripord)
% plotSolutionT ==> ��ͼLU/STRIP/BIN
% T: table ��ʽ��LU
% V: table ��ʽ��Veh
% plotLU����0����LU��ͼ��2����T.index˳������˳����ͼ 3:��T.order�������˳����ͼ 1����T��Ŀǰ˳����ͼ��
% plotItem����0����Item��ͼ�� 1����Item�������ͼ
% plotStrip����0����Strip��ͼ��1����Strip��˳����ͼ
% plotBin����0����Bin��ͼ��1��
% example��plotSolutionT(d.LU,d.Veh,1,1,1,1)
%                  plotSolutionT(d.LU,d.Veh,2,0,0,3)

if isstruct(T),  T = struct2table(structfun(@(x) x',T,'UniformOutput',false)); end
if isstruct(V), V = struct2table(structfun(@(x) x',V,'UniformOutput',false)); end

if nargin < 3
    plotLU = 1;  % 1������ĿǰT��˳�� 2: ����index˳��
    plotStrip = 1;
    plotBin = 1;
    plotItem = 1;
    
end
if nargin < 8
    plotColor = 3;
    figname ='';
end

global ISplotPause ISplotShowType

%% 0 V��T �״�join ��ȡ�����������
if ~ismember('LU_VehType', V.Properties.VariableNames)
V.LU_VehType = [1:height(V)]'; 
V = V(:,{'LU_VehType','LWH'}); 
V.Properties.VariableNames{'LWH'} = 'LWH_V'; end

if ~ismember('LU_VehType', T.Properties.VariableNames)
    T.LU_VehType = ones(height(T),1); end

if ~ismember('LWH_V', T.Properties.VariableNames)
    T = join(T,V,'Keys','LU_VehType');  end

%% ��������ITEM��plot todo ������Ҫ������
% V1 ��������ITEM��plot ��ITEM��ΪLU -> plotSolutionT(sItem,Veh)  
% % if iscell(T.LID)
% %     T.LID = cell2mat(T.LID); end

%% 1 ��ȡLUcolor �� ��ɫ ���Ե�T  LID/SID/OPID etc   % ��ͼ���� 1 LU��bin���; 2 LWH 3 COORDLUBIN 4 LID ����

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


%% ͼ1����ͼLU,������LU˳��
if plotLU
    figure('name',strjoin({figname,'LUչʾ����ITEMʱ������󣬺ϼ�*��',num2str(height(T))}));
    
    if plotLU==1, subT = T(1:end,:); end                          % LUĿǰ˳��
    if plotLU==2, subT = sortrows(T,'Index'); end           % LU��BBA����˳��
    
    if plotLU==3, subT = T(T.order,:); end                       %LU�����˳��  %     if plotLU==3, subT = sortrows(T,'order'); end      %LU�����˳��
    
    % 5555 ����LU������ϵCoord X��������LU�Ŀ��L�ƶ� YΪ0 ZΪLU�߶�
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

% ʵ�ʻ���plotLU��ͬһItem����˳�򣬸���ͬ�ĸ߶�Zֵ.
if plotItem
    XYZ = zeros(height(T),3);
    
    % ���Item��ͼ
    nItem = max(T.LU_Item(:,1)); %Item�ĸ���
    figure('name',strjoin({figname,'Itemչʾ���ϼ�*��',num2str([nItem])}));
    
    if plotItem==1, order = 1:nItem;  end       % Ĭ��˳����ͼ ��1��n
    if plotItem==2, order = itemord;       end       % ����˳����ͼ
    
    % ����ÿ��Lu��Item������Coord
    sumL = 0;
    for iStrip=1:nItem
        % �ӵ�1��Item��ʼ������Lu��x����(��ͬ�Ѷ���һ����
        %fidxItem = T.LU_Item(:,1)==iItem;  %ͬһ�Ѷ��µ��߼�ֵ(�����
        fidxStrip = T.LU_Item(:,1)==order(iStrip);  %ͬһ�Ѷ��µ��߼�ֵ(�����
        XYZ(fidxStrip,1) = sumL;
        sumL = sumL + unique(T.LWH(fidxStrip,1));       if numel(sumL) >  1 , error('e'); end
        
        nLU = sum(fidxStrip); % (ͬһ�Ѷ��µ���������)
        sumH = 0;
        for iLU = 1:nLU 
            % ����Lu��z����(��ͬ�Ѷ�߶Ȳ�һ���� ����1����
            fidxLU = T.LU_Item(:,2)==iLU;
            XYZ(fidxStrip&fidxLU,3) = sumH;
            sumH = sumH + T.LWH(fidxStrip&fidxLU,3);
        end
    end
    
    % ����˳��Ѷ�
    T.XYZ = XYZ;
    T = sortrows(T,{'XYZ'},{'ascend'}); 
    
    for iLU=1:height(T)
        plotcube(T.LWH(iLU,:), T.XYZ(iLU,:),0.7, T.LUcolor(iLU,:));
        axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        view(60,40); %view(111,33);
        if ISplotPause>0 ,      pause(ISplotPause);   end
    end
end

%% ͼ2����ͼStrip���� v2����LU_Strip(CoordLUStrip)���ɣ�
if plotStrip
    XYZ = zeros(height(T),3);
    
    % ���Strip��ͼ
    nstrip = max(T.LU_Strip(:,1)); %Item�ĸ���
    figure('name',strjoin({figname,'STRIPչʾ���Ⱥ�˳������󣬺ϼ�*��',num2str([ nstrip ])}));
        
    if plotStrip==1, order = 1:nstrip;   end    % Ĭ��˳����ͼ ��1��n                  
    if plotStrip==2, order = stripord;       end     % ����˳����ͼ
    
    % ����ÿ��Lu��Strip������Coord ��Xֵ�����ⲿCoord����Ϊ˳��仯���Լ�����˳��
    sumL = 0;
    for iStrip=1:nstrip
        % �ӵ�1��Strip��ʼ������Lu��x����(��ͬ�Ѷ���һ����
        fidxStrip = T.LU_Strip(:,1)==order(iStrip);  %ͬһstrip�µ��߼�ֵ(�����
        XYZ(fidxStrip,2) = sumL;  % Y = LWH'S L
        sumL = sumL + unique(T.LWH(fidxStrip,2));       if numel(sumL) >  1 , error('e'); end
        
        XYZ(fidxStrip,1) = T.CoordLUStrip(fidxStrip,1);  %
        XYZ(fidxStrip,3) = T.CoordLUStrip(fidxStrip,3);
    end
    
    % ����˳��Ѷ�    
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



%% ͼ3����ͼBin���� (��LU_Bin����)
if plotBin
pos = 1;
if ismember('LU_Bin', T.Properties.VariableNames)
    
    T=sortrows(T,{'LU_Bin'},{'ascend'});  % T������ (LU_Bin����)
    
    % ���bin��ͼ
    nBin = max(T.LU_Bin(:,1)); %bin�ĸ���
    figure('name',strjoin({figname,'BINչʾ���Ⱥ�˳������󣬺ϼ�*��',num2str([nBin])}));
    
    for ibin=1:nBin % for ibin=nBin:nBin
        
        subplot(2,ceil((nBin+1)/2),pos); pos = pos+1;
        
        subT = T(T.LU_Bin(:,1)==ibin,:);                  % yxz =  T.LWH_T(f,:); coord = T.CoordLUBin(f,:);  color = T.LUcolor(f,:);  yxzVeh = T.LWH_V(f,:);
        
        for iLU=1:height(subT)
            plotcube(subT.LWH(iLU,:),subT.CoordLUBin(iLU,:),0.7, subT.LUcolor(iLU,:));
            axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
            view(111,33);    %view(60,40);
            xlim([0 subT.LWH_V(iLU,1)]);  ylim([0 subT.LWH_V(iLU,2)]); zlim([0 subT.LWH_V(iLU,3)]); % �����ĳ���ߵ��������ʵĳ���
            if ISplotPause>0 ,      pause(ISplotPause);   end
        end
        
    end
end
end

%% ע��
% % % % % 3 ��ʼ��ͼ - ��Strip����
% % % if ismember('LU_Strip', T.Properties.VariableNames),
% % %     % 2 T������ (LU_Strip����)
% % %     T=sortrows(T,{'LU_Strip'},{'ascend'});
% % %     if ismember('LU_Bin', T.Properties.VariableNames),  subplot(2,ceil((nBin+1)/2),pos); 
% % %     else  subplot(1,1,pos); end
% % %     subT= T;
% % %     for iLU=1:height(subT)
% % %         plotcube(subT.LWH(iLU,:),subT.CoordLUStrip(iLU,:),0.7, subT.LUcolor(iLU,:));
% % %         axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
% % %         view(111,33);    %view(60,40);
% % %         xlim([0 subT.LWH_V(iLU,1)]);  zlim([0 subT.LWH_V(iLU,3)]); % �����ĳ���ߵ��������ʵĳ���
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
% % % 0 �״�join ��ȡ�����������
% % if ~ismember('LU_VehType', V.Properties.VariableNames),
% % V.LU_VehType = [1:height(V)]'; 
% % V = V(:,{'LU_VehType','LWH'}); 
% % V.Properties.VariableNames{'LWH'} = 'LWH_V'; end
% % if ~ismember('LU_VehType', T.Properties.VariableNames),  T.LU_VehType = ones(height(T),1); end
% % if ~ismember('LWH_V', T.Properties.VariableNames),  T = join(T,V,'Keys','LU_VehType');  end
% % 
% % % 1 ��ȡ�� ��ɫ ���Ե�T  LID/SID/OPID etc   % ��ͼ���� 1 LU��bin���; 2 LWH 3 COORDLUBIN 4 LID ����
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
% % % ����ͼLU,������˳����ͼ
% % plotLU = 0;
% % if plotLU
% %     subT = T(1:end,:);
% % 
% %     % 5555 ������plotLU��LU������ϵ
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
% % % 2 ��ʼ��ͼ - ��Bin����
% % if ismember('LU_Bin', T.Properties.VariableNames), 
% % % 2 T������ (LU_Bin����)
% % T=sortrows(T,{'LU_Bin'},{'ascend'});    
% % % ���bin��ͼ
% % nBin = max(T.LU_Bin(:,1)); %bin�ĸ���
% % figure('name',num2str([nBin]));
% % 
% % for ibin=1:nBin
% %     subplot(2,ceil((nBin+1)/2),pos); pos = pos+1;    
% %     subT = T(T.LU_Bin(:,1)==ibin,:); % yxz =  T.LWH_T(f,:); coord = T.CoordLUBin(f,:);  color = T.LUcolor(f,:);  yxzVeh = T.LWH_V(f,:);
% %     for iLU=1:height(subT)
% %         plotcube(subT.LWH(iLU,:),subT.CoordLUBin(iLU,:),0.7, subT.LUcolor(iLU,:));        
% %         axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
% %         view(111,33);    %view(60,40);      
% %         xlim([0 subT.LWH_V(iLU,1)]);  ylim([0 subT.LWH_V(iLU,2)]); zlim([0 subT.LWH_V(iLU,3)]); % �����ĳ���ߵ��������ʵĳ���
% %                 if ISplotPause>0 ,      pause(ISplotPause);   end
% %     end
% % end
% % end
% % 
% % % 3 ��ʼ��ͼ - ��Strip����
% % if ismember('LU_Strip', T.Properties.VariableNames),
% %     % 2 T������ (LU_Strip����)
% %     T=sortrows(T,{'LU_Strip'},{'ascend'});
% %     if ismember('LU_Bin', T.Properties.VariableNames),  subplot(2,ceil((nBin+1)/2),pos); 
% %     else  subplot(1,1,pos); end
% %     subT= T;
% %     for iLU=1:height(subT)
% %         plotcube(subT.LWH(iLU,:),subT.CoordLUStrip(iLU,:),0.7, subT.LUcolor(iLU,:));
% %         axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
% %         view(111,33);    %view(60,40);
% %         xlim([0 subT.LWH_V(iLU,1)]);  zlim([0 subT.LWH_V(iLU,3)]); % �����ĳ���ߵ��������ʵĳ���
% %         if ISplotPause>0 ,      pause(ISplotPause);   end
% %     end
% % end
% % 
% % end

%%
% j = 1;
%         % �����̲���������ʱ, �ж�
%         if j~=1 && tmp ~=m(8,j) 
%                  if ISplotPause>0 ,      pause(ISplotPause);   end
%         end        
%         tmp = m(8,j);
%         j=j+1;



%% ͼ2����ͼStrip���� ����LU_Strip(CoordLUStrip)���ɣ�
% if plotStrip
% %     if ismember('LU_Strip', T.Properties.VariableNames)
% 
%     if plotStrip==1, T=sortrows(T,{'LU_Strip'},{'ascend'});  end                    
% %     if plotStrip==2, T=sortrows(T,{'LU_Strip'},{'ascend'});  end
%             
%         figure('name',strjoin({figname,'STRIPչʾ���Ⱥ�˳������󣬺ϼ�*��',num2str([ max(T.LU_Strip(:,1));])}));
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

