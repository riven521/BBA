function [] = plotSolutionT(T,V,plotLU,plotStrip)
% plotSolutionT ==> ��ͼLU/STRIP/BIN
%   T: LU; V: Veh.
if isstruct(T),  T = struct2table(structfun(@(x) x',T,'UniformOutput',false)); end
if isstruct(V), V = struct2table(structfun(@(x) x',V,'UniformOutput',false)); end

if nargin < 3
    plotLU = 1;
    plotStrip = 1;
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

%% ͼ1����ͼLU,������LU˳��
if plotLU
    
    figure('name',strjoin({'LUչʾ����ITEMʱ������󣬺ϼ�*��',num2str(height(T))}));
    subT = T(1:end,:);

    % 5555 ������plotLU��LU������ϵ
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

%% ͼ2����ͼStrip���� ����LU_Strip(CoordLUStrip)���ɣ�
if plotStrip
    if ismember('LU_Strip', T.Properties.VariableNames)
        
        T=sortrows(T,{'LU_Strip'},{'ascend'});     % T������ (LU_Strip����)
        
        figure('name',strjoin({'STRIPչʾ���Ⱥ�˳������󣬺ϼ�*��',num2str([ max(T.LU_Strip(:,1));])}));
        subT = T(1:end,:);  % subT = T;
        
        for iLU=1:height(subT)
            plotcube(subT.LWH(iLU,:),subT.CoordLUStrip(iLU,:),0.7, subT.LUcolor(iLU,:));
            axis equal;         grid on;        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
            view(111,33);    %view(60,40);
            if ISplotPause>0 ,      pause(ISplotPause);   end
        end
    end
end

%% ͼ3����ͼBin���� (��LU_Bin����)
pos = 1;
if ismember('LU_Bin', T.Properties.VariableNames)
    
    T=sortrows(T,{'LU_Bin'},{'ascend'});  % T������ (LU_Bin����)
    
    % ���bin��ͼ
    nBin = max(T.LU_Bin(:,1)); %bin�ĸ���
    figure('name',strjoin({'BINչʾ���Ⱥ�˳������󣬺ϼ�*��',num2str([nBin])}));
    
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
            xlim([0 subT.LWH_V(iLU,1)]);  ylim([0 subT.LWH_V(iLU,2)]); zlim([0 subT.LWH_V(iLU,3)]); % �����ĳ���ߵ��������ʵĳ���
            if ISplotPause>0 ,      pause(ISplotPause);   end
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

