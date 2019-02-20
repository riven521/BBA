function [LU] = HBinGapAdjust(LU,VEH)
% HBinGapAdjust ==> Gap Adjust for Each Bin
 
TLU = getTableLU(LU);
TVEH = getTableLU(VEH);

typeVeh = unique(TLU.LU_Bin(:,1));
numVeh = length(typeVeh);      if numVeh>1, error('���������Ҫ����Gap'); end

% ѭ��ÿ��VEH���ҵ���Ӧ������subTLU���ӳ���subVeh
for idxVeh = 1:numVeh
    subTLU = TLU(TLU.LU_Bin(:,1) == typeVeh(idxVeh), : );    
    subVeh = TVEH(unique(subTLU.LU_VehType), :); %TODO CHECK �Ƿ�LU_VehType���ǳ�������ĵ�һ��?
    
        if height(subVeh)>1,  error('NOT LU in the same Veh type'); end
        
    % 555  �����̵�����Ҫ����
    subTLUNew = HGapAdjust(subTLU,subVeh);   
    
    TLU(TLU.LU_Bin(:,1) == typeVeh(idxVeh), : ) = subTLUNew;
    
end

[LU] = getSturctT(TLU);

end

%% 555  �����̵�����Ҫ����
% ����ͬһbin�ڵ�LU, VEH��һ����
function [LU] = HGapAdjust(LU,VEH)
global parMulipleGap ISplotShowGapAdjust ISplotGapCompare ISplotPause

% 1 INITILIZE
flagGap=0; % 1: ���� 0: δ����

% 2 Get bottomLU
bottomLU = LU(LU.CoordLUBin(:,3)==0, : );  % �ײ������ % plotSolutionT(LU,VEH);

% 2 Get pgGap: VEH�ڵ�ʣ��հ�����  
pgLU = getPgLU(bottomLU,0);  % plot(pgLU);
pgVEH = polyshape(pgRectangle(0,0,VEH.LWH(1,1),VEH.LWH(1,2)));	    
pgGap = subtract(pgVEH,pgLU);    if pgGap.NumRegions  > 1,  warning('Exsit %d Regions in this pgon', pgGap.NumRegions);  end
                        % TODO pgGap����? ,gBlanks = sortregions(pgBlanks,'area','ascend');

% 2 Get CoordGap : pgGap �� ���� ��X��Y 
[XBound,YBound]=boundary(pgGap);
coordGapArray = fliplr(sortrows([YBound,XBound]));  % ��Y��X����, �̶����Ҹı�˳��

% 3 ѭ��ÿ��Boundary����
for iVertex=1:size(coordGapArray,1)
    % Get Each GapBound Coord X and Y
    % ÿ���߽綥�������ֵ
    coordX = coordGapArray(iVertex,1);
    coordY = coordGapArray(iVertex,2);
            
    % ��ȡ: ��������������� Get idxLUs i.e. Y of LU is larger than coordY of bound vertex 
    flagLeqLU = bottomLU.CoordLUBin(:,2) - bottomLU.margin(:,3) > coordY;
    idxLUs = find(flagLeqLU);         % �ҳ���yPointֵ �ߵ�����bottom��LU������    
    
    % ���û�и��������̣�����.
    if isempty(idxLUs), continue; end       
           
    % ��Լʱ��: �����Boundary�Ķ��㲻�ܷ�����С��LU��(�����������ı߽�), ��������Boundary����һ��
    minLW = min(min(LU.LWH(idxLUs,[1,2])));
%     if ~isinterior(pgVEH,coordX+minLW,coordY) || ~isinterior(pgVEH,coordX,coordY+minLW) || ~isinterior(pgVEH,coordX+minLW/2, coordY+minLW/2),    continue;    end
    
    % ѭ��ÿ��Boundary������ÿ������������LU    
    %   ����: ��������
    [~,ff] = sortrows([bottomLU.CoordLUBin(idxLUs,[2,1])]);  idxLUs = idxLUs(ff);
    for ii=1:length(idxLUs)
        idxLU = idxLUs(ii); % bottomLU�����
        thisLU = bottomLU(idxLU,:);
        
        % pgLU1 ���ϲ���LU pgGapNew �ϲ������Gap
        pgLU1 = getPgLU(thisLU, 0);        
        pgGapNew = union(pgGap,pgLU1);
        
        % pgLU2 ���°ڷŵ�LUλ��
        pgLU2 = getPgLU(thisLU, 0, coordX,coordY);        %         pgLU2 =  polyshape(pgRectangle(coordX,coordY,wLU,lLU));  
        flagLU = all(isinterior(pgGapNew, getBoundaryCentroid(pgLU2)));        
        
        pgIntersect = intersect(pgLU2,subtract(pgVEH,pgGapNew));
        if pgIntersect.NumRegions ~= 0, flagLU = 0; end
        
        % pgLU2Rota ���°ڷŵ���ת��LUλ��
        flagLURota=0;
        if thisLU.isRota     %bottomLU.isRota(idxLU)
            pgLU2Rota =  getPgLU(thisLU, 1, coordX,coordY); %polyshape(pgRectangle(coordX,coordY,lLU,wLU));
            flagLURota = all(isinterior(pgGapNew, getBoundaryCentroid(pgLU2Rota))); 
            
            pgIntersectRota = intersect(pgLU2Rota,subtract(pgVEH,pgGapNew));
            if pgIntersectRota.NumRegions ~= 0, flagLURota = 0; end
        end        
        
        % �����ɷ���, ����ѡ��һ�����������
        if flagLU && flagLURota
            if thisLU.LWH(1) > thisLU.LWH(2) %wLU>lLU
                flagLURota=0;
            else
                flagLU=0;
            end
        end
            
        if ISplotShowGapAdjust
            
            plotGapPgon(VEH,pgLU,pgGap,pgVEH,coordGapArray,coordX,coordY,pgLU1,pgLU2,pgLU2Rota,flagLU,flagLURota);
% figure('name',strjoin({'Gapչʾ��'}));
% plotGapPgon(VEH,pgLU,pgGap,pgVEH,coordGapArray,coordX,coordY,pgLU1,pgLU2,pgLU2Rota,flagLU,flagLURota,pgGapNew);
        end
        
        if flagLU || flagLURota
            
            
            if ISplotGapCompare,            plotSolutionT(LU,VEH,0,0,0,1,3,'Gap�����ɹ�չʾ');         end
            
            fprintf(1,'       Exsiting ��װ��϶ in HBinGapAdjust (do2/do3)...\n');    
            
            % fLU: ��Ҫ������LU���  ��botttomLU����, ���践�ص�LU�滻, ���ǵײ�LUҲҪ����
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
            
            if ISplotGapCompare,            plotSolutionT(LU,VEH,0,0,0,1,3,'Gap�����ɹ�չʾ');         end
            
        end
        
        if flagGap,   break;    end
        
    end
    
    if flagGap,   break;    end
    
end

% ������ε���Gap�ɹ�������Ҫ�ݹ�, ���¶Ը�bin���е���
if flagGap && parMulipleGap
    
     [LU] = HGapAdjust(LU,VEH);  %�滻�ɹ�����Ҫ��LU��Gap��pg�ؼ���
     
     bottomLU = LU(LU.CoordLUBin(:,3)==0, : );  % �ײ������ % bottomLU ��Ҫ�ؼ���
     pgLU = getPgLU(bottomLU,0);  
     pgVEH = polyshape(pgRectangle(0,0,VEH.LWH(1,1),VEH.LWH(1,2)));
     pgGap = subtract(pgVEH,pgLU);    if pgGap.NumRegions  > 1,  warning('Exsit %d Regions in this pgon', pgGap.NumRegions);  end
     [XBound,YBound]=boundary(pgGap);
    coordGapArray = fliplr(sortrows([YBound,XBound]));  % ��Y��X����, �̶����Ҹı�˳��
    
end




% 190108 �Ż�
%  1 ��ȥ����������϶
gapY=sort(pgGap.Vertices(:,2));
gapY=unique(gapY);
% �����pgonʱ,�����NaNֵ, ��Ҫ�ų�
gapY = gapY(~isnan(gapY));

% if ISplotShowGapAdjust && length(gapY)>2
%     figure('name',strjoin({'Gap��������������չʾ��'}));
%     plotGapPgon(VEH,pgLU,pgGap,pgVEH,coordGapArray,coordX,coordY);
% end

for g=1:length(gapY)-1
    pgRect=polyshape([0 gapY(g);  VEH.LWH(1,1) gapY(g); VEH.LWH(1,1) gapY(g+1); 0 gapY(g+1)]); % ����
    
    if ISplotShowGapAdjust
%         plotGapPgon(VEH,pgLU,pgGap,pgVEH);
                plotGapPgon(VEH,pgLU,pgGap,pgVEH,coordGapArray,coordX,coordY,[],[],[],[],[],pgRect);
    end
    
    % ��������������Һ����жѶ⣬����ǰ�ƶ� bug����ΪĳЩʱ������lu�ڳ����м䣬��������Ե; �޸�������overlaps�ж�������Rect��LU�����ص�
    if all(isinterior(pgGap, getBoundaryCentroid(pgRect))) && ~overlaps(pgLU,pgRect)%��pgRect�ı߽���pgGap�У�����true
        warning('���ڴ���������');
        fYLU = LU.CoordLUBin(:,2)>gapY(g);
        if any(fYLU)
            LU.CoordLUBin(fYLU,2) =  LU.CoordLUBin(fYLU,2) - (gapY(g+1) - gapY(g));  %��ǰ�ƶ�һ������
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



% ͨ�ú��� �� ��ȡ���̣������Ķ���� 
%   TLU�����̣����ϣ�
%   isRota : �Ƿ���Ҫ��ת���ı䳤��
%   varargin���Ƿ���ù̶����꣨��LU����ʼλ�ã�
function [pgon] = getPgLU(TLU,isRota,varargin)
% polygon of LUs ( Note: only get polyshape whose hight = 0 )
TLU = sortrows(TLU,'CoordLUBin');
P = [];
for idxl=1:height(TLU)
    % ���ù̶�����, ���пɱ��������
    if nargin>2
        x = varargin{1};
        y = varargin{2};
    else
        x=TLU.CoordLUBin(idxl,1)-TLU.margin(idxl,1);
        y=TLU.CoordLUBin(idxl,2)-TLU.margin(idxl,4);
    end
    
    w = TLU.LWH(idxl,1) + TLU.margin(idxl,1 ) + TLU.margin(idxl,2 );
     l = TLU.LWH(idxl,2) + TLU.margin(idxl,3 ) + TLU.margin(idxl,4 );
    
    % ��ת��, ��isRota==1
    if isRota
        P =  [P;pgRectangle(x,y,l,w);[NaN,NaN]];
    else
        P =  [P;pgRectangle(x,y,w,l);[NaN,NaN]];
    end
end
pgon = polyshape(P);
end

% ��ȡ���жϵ��ڵ�, ��Boundary��Centroid��ͬ����
function [P] = getBoundaryCentroid(pgon)
[x,y] = boundary(pgon);
[xcenter,ycenter] = centroid(pgon);
P = [[x;xcenter],[y;ycenter]];
end


function plotGapPgon(VEH,pgLU,pgGap,pgVEH,coordGapArray,coordX,coordY,pgLU1,pgLU2,pgLU2Rota,flagLU,flagLURota,pgRect)       

        global ISplotPause;        
        
        plot(pgVEH,'FaceColor','white','FaceAlpha',0.01);
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
        plot(pgLU,'FaceColor','green','FaceAlpha',0.2)
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
        plot(pgGap,'FaceColor','blue','FaceAlpha',0.2)
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
        
        if nargin == 7        
        plot(coordGapArray(:,1),coordGapArray(:,2),'.', 'MarkerSize', 8);
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
        plot(coordX,coordY,'.', 'MarkerSize', 20);
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
        end
        
        if nargin == 12
            
        plot(pgLU1,'FaceColor','green','FaceAlpha',0.5);         %         plot(pgGapNew,'FaceColor','green','FaceAlpha',0.8)
        hold on;    axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause/100);   end
        
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
            plot(pgRect,'FaceColor','red','FaceAlpha',0.2)
            axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      if ISplotPause>0 ,      pause(ISplotPause*1.5);   end;  hold on;
        end   
        
        
end
        




%% V1
% % %% 555  �����̵�����Ҫ����
% % % ����ͬһbin�ڵ�LU, VEH��һ����
% % function [LU] = HGapAdjust(LU,VEH)
% % global parMulipleGap ISplotShowGapAdjust
% % 
% % % 1 INITILIZE
% % flagGap=0; % 1: ���� 0: δ����
% % 
% % % 2 Get bottomLU
% % bottomLU = LU(LU.CoordLUBin(:,3)==0, : );  % �ײ������ % plotSolutionT(LU,VEH);
% % 
% % % 2 Get pgGap: VEH�ڵ�ʣ��հ�����  
% % pgLU = getPgLU(bottomLU,0);  % plot(pgLU);
% % pgVEH = polyshape(pgRectangle(0,0,VEH.LWH(1,1),VEH.LWH(1,2)));	    
% % pgGap = subtract(pgVEH,pgLU);    if pgGap.NumRegions  > 1,  warning('Exsit %d Regions in this pgon', pgGap.NumRegions);  end
% %                         % TODO pgGap����? ,gBlanks = sortregions(pgBlanks,'area','ascend');
% % 
% % % 2 Get CoordGap : pgGap �� ���� ��X��Y 
% % [XBound,YBound]=boundary(pgGap);
% % coordGapArray = fliplr(sortrows([YBound,XBound]));  % ��Y��X����, �̶����Ҹı�˳��
% % 
% % % 3 ѭ��ÿ��Boundary����
% % for iVertex=1:size(coordGapArray,1)
% %     % Get Each GapBound Coord X and Y
% %     % ÿ���߽綥�������ֵ
% %     coordX = coordGapArray(iVertex,1);
% %     coordY = coordGapArray(iVertex,2);
% %             
% %     % ��ȡ: ��������������� Get idxLUs i.e. Y of LU is larger than coordY of bound vertex 
% %     flagLeqLU = bottomLU.CoordLUBin(:,2) - bottomLU.margin(:,3) > coordY;
% %     idxLUs = find(flagLeqLU);         % �ҳ���yPointֵ �ߵ�����bottom��LU������    
% %     
% %     % ���û�и��������̣�����.
% %     if isempty(idxLUs), continue; end       
% %            
% %     % ��Լʱ��: �����Boundary�Ķ��㲻�ܷ�����С��LU��(�����������ı߽�), ��������Boundary����һ��
% %     minLW = min(min(LU.LWH(idxLUs,[1,2])));
% %     if ~isinterior(pgVEH,coordX+minLW,coordY) || ~isinterior(pgVEH,coordX,coordY+minLW) || ~isinterior(pgVEH,coordX+minLW/2, coordY+minLW/2),    continue;    end
% %     
% %     % ѭ��ÿ��Boundary������ÿ������������LU    
% %     %   ����: ��������
% %     [~,ff] = sortrows([bottomLU.CoordLUBin(idxLUs,[2,1])]);  idxLUs = idxLUs(ff);
% %     for ii=1:length(idxLUs)
% %         idxLU = idxLUs(ii); % bottomLU�����
% %         thisLU = bottomLU(idxLU,:);
% %         
% %         % pgLU1 ���ϲ���LU pgGapNew �ϲ������Gap
% %         pgLU1 = getPgLU(thisLU, 0);        
% %         pgGapNew = union(pgGap,pgLU1);
% %         
% %         % pgLU2 ���°ڷŵ�LUλ��
% %         pgLU2 = getPgLU(thisLU, 0, coordX,coordY);        %         pgLU2 =  polyshape(pgRectangle(coordX,coordY,wLU,lLU));  
% %         flagLU = all(isinterior(pgGapNew, getBoundaryCentroid(pgLU2)));        
% %         
% %         pgIntersect = intersect(pgLU2,subtract(pgVEH,pgGapNew));
% %         if pgIntersect.NumRegions ~= 0, flagLU = 0; end
% %         
% %         % pgLU2Rota ���°ڷŵ���ת��LUλ��
% %         flagLURota=0;
% %         if thisLU.isRota     %bottomLU.isRota(idxLU)
% %             pgLU2Rota =  getPgLU(thisLU, 1, coordX,coordY); %polyshape(pgRectangle(coordX,coordY,lLU,wLU));
% %             flagLURota = all(isinterior(pgGapNew, getBoundaryCentroid(pgLU2Rota))); 
% %             
% %             pgIntersectRota = intersect(pgLU2Rota,subtract(pgVEH,pgGapNew));
% %             if pgIntersectRota.NumRegions ~= 0, flagLURota = 0; end
% %         end        
% %         
% %         % �����ɷ���, ����ѡ��һ�����������
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
% %         % plot ĳЩvertex�ĳ��Թ���(���Բ����ܵĵ��Ѿ��ų����������ų���ע��65�н�Լʱ�����)
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
% % % figure('name',strjoin({'Gapչʾ��'}));
% % % plotGapPgon(VEH,pgLU,pgGap,pgVEH)%,coordGapArray,coordX,coordY,pgLU1,pgLU2,pgLU2Rota,flagLU,flagLURota,pgGapNew);
% % figure('name',strjoin({'Gapչʾ��'}));
% % plotGapPgon(VEH,pgLU,pgGap,pgVEH,coordGapArray,coordX,coordY,pgLU1,pgLU2,pgLU2Rota,flagLU,flagLURota);
% % % figure('name',strjoin({'Gapչʾ��'}));
% % % plotGapPgon(VEH,pgLU,pgGap,pgVEH,coordGapArray,coordX,coordY,pgLU1,pgLU2,pgLU2Rota,flagLU,flagLURota,pgGapNew);
% % 
% %         end
% %         
% %         if flagLU || flagLURota
% %             if ISplotShowGapAdjust
% %                 % plot ��������vertex�Ĺ��̣���ͼ�ɹ��Ķ���͵�����
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
% %                 plotSolutionT(LU,VEH,0,0,0,1,3,'ID��ͼBin');                  %     Bin�����
% %                 
% %             end
% %             fprintf(1,'       Exsiting ��װ��϶ in HBinGapAdjust (do2/do3)...\n');    
% %             
% %             % fLU: ��Ҫ������LU���  ��botttomLU����, ���践�ص�LU�滻, ���ǵײ�LUҲҪ����
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
% % % ������ε���Gap�ɹ�������Ҫ�ݹ�, ���¶Ը�bin���е���
% % if flagGap && parMulipleGap
% %     
% %      [LU] = HGapAdjust(LU,VEH);  %�滻�ɹ�����Ҫ��LU��Gap��pg�ؼ���
% %      
% %      bottomLU = LU(LU.CoordLUBin(:,3)==0, : );  % �ײ������ % bottomLU ��Ҫ�ؼ���
% %      pgLU = getPgLU(bottomLU,0);  
% %      pgVEH = polyshape(pgRectangle(0,0,VEH.LWH(1,1),VEH.LWH(1,2)));
% %      pgGap = subtract(pgVEH,pgLU);    if pgGap.NumRegions  > 1,  warning('Exsit %d Regions in this pgon', pgGap.NumRegions);  end
% %      [XBound,YBound]=boundary(pgGap);
% %     coordGapArray = fliplr(sortrows([YBound,XBound]));  % ��Y��X����, �̶����Ҹı�˳��
% %     
% % end
% % 
% % % �������Gap���ɹ�, ����Щ����
% % if ISplotShowGapAdjust
% %     pausetime = 0.0;
% %     % plot ĳЩvertex�ĳ��Թ���
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
% % % 190108 �Ż�
% % %  1 ��ȥ����������϶
% % gapY=sort(pgGap.Vertices(:,2));
% % gapY=unique(gapY);
% % % �����pgonʱ,�����NaNֵ, ��Ҫ�ų�
% % gapY = gapY(~isnan(gapY));
% % for g=1:length(gapY)-1
% %     pgRect=polyshape([0 gapY(g);  VEH.LWH(1,1) gapY(g); VEH.LWH(1,1) gapY(g+1); 0 gapY(g+1)]); % ����
% %     
% %     if ISplotShowGapAdjust
% %         plot(pgRect,'FaceColor','red','FaceAlpha',0.2)
% %         axis equal;    grid on;    xlim([0 1.5*VEH.LWH(1,1)]);    ylim([0 1.2*VEH.LWH(1,2)]);      pause(pausetime*1.5);  hold on;
% %     end
% %     
% %     % ��������������Һ����жѶ⣬����ǰ�ƶ� bug����ΪĳЩʱ������lu�ڳ����м䣬��������Ե; �޸�������overlaps�ж�������Rect��LU�����ص�
% %     if all(isinterior(pgGap, getBoundaryCentroid(pgRect))) && ~overlaps(pgLU,pgRect)%��pgRect�ı߽���pgGap�У�����true
% %         warning('���ڴ���������');
% %         fYLU = LU.CoordLUBin(:,2)>gapY(g);
% %         if any(fYLU)
% %             LU.CoordLUBin(fYLU,2) =  LU.CoordLUBin(fYLU,2) - (gapY(g+1) - gapY(g));  %��ǰ�ƶ�һ������
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