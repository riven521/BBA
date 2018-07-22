function [Res1_LUBeBinMatrix,Res2_CoordLUBin,Res3_LWHRota,Res4_DrawSeq,da] = ...
    BBA_Main(LUID,LULWH,BINLWH,PARA_whichRotationHori,LUBUFF,BINBUFF)
% ������������ LUID,LULWH,BINLWH,PARA,LUBUFF,BINBUFF
%  LUID - ������(n��) �������� ��ͬ���ֱ���ͬһ����,����Ѷ�
%  LULWH - ����(3��*n��) ���̳���� ��ͼʹ�ø�ֵ
%  BINLWH - ������(3��) ��������� Ŀǰ�����ǵ����� TODO �������ӵ��೵�Ͱ�˳��ʹ��
%  PARA_whichRotationHori - ����(1*1) ����rotation,������rotation�Ĳ���
%  LUBUFF - ������(2��) ���̼䳤����ܼ�϶(2���ĵ��߳����϶) �������̳����=ÿ�����̵�ʵ�ʳ����+���ӵ�buff
%  BINBUFF - ������(3��) �����ĳ���ߵ��ܼ�϶(2���ĵ��߳����϶,1���ĸ߶ȼ�϶) ���ó��ͳ����=���͵�ʵ�ʳ����-BINBUFF
% ����ĸ����� Res1_LUBeBinMatrix,Res2_CoordLUBin,Res3_LWHRota,Res4_DrawSeq
%  Res1_LUBeBinMatrix - ����(2*n��)
%  Res2_CoordLUBin - double����(3*n��)
%  Res3_LWHRota - ����(*n��)
%  Res4_DrawSeq - double����(2*n��)

%% Call ������:
% clear;close all; format long g; format bank; %NOTE ����MATLAB CODE ֧��
% rng('default');rng(1); % NOTE �Ƿ�����ı�־
% LUArray = struct('ID',[],'LWH',[],...
%     'weight',[],'Lbuffer',[],'Wbuffer',[],'Type',[],'Material',[]);
% ItemArray = struct('ID',[],'LWH',[],...
%     'weight',[],'Lbuffer',[],'Wbuffer',[]);
% StripArray = struct('ID',[],'LWH',[],... %ONLY LW
%     'weight',[],'Lbuffer',[],'Wbuffer',[]);
% BinArray = struct('ID',[],'LWH',[],...
%    'Capacity',[],'Lbuffer',[],'Wbuffer',[],'Hbuffer',[]);
% da = struct('LUArray',LUArray,'ItemArray',ItemArray,'StripArray',StripArray,'BinArray',BinArray);
%% 1������ʼ��
% whichStripH 1 best 2 first 3 next; whichBinH 1 best; TODO ��������������ʽ
% whichSortItemOrder 1 ���ߵݼ� 2 ��̱ߵݼ�; 
% whichRotation 1:����rotation 0:��ֹ
% rotation��� 1 1 2 1 0 (1 1 2 1 1 )(1 1 2 1 2) % ��rotation��� 1 1 1 0 0 ��2/3 1 1 0 0��
% whichRotationHori 0:�ڰ���˳��ʱ��FBS_{RG}��ʽ; 1��New/NoNew��Horizon��ʽ 2��New/NoNew��Vertical��ʽ
% ParaArray = struct('whichStripH',1,'whichBinH',1,'whichSortItemOrder',2,...
%     'whichRotation',1,'whichRotationHori',0,'timeLimit',100,'ub0',10);
ParaArray = struct('whichStripH',1,'whichBinH',1,'whichSortItemOrder',2,...
    'whichRotation',1,'whichRotationHori',0,'timeLimit',100,'ub0',10);
%% 2�ṹ��da��ֵ
if nargin ~=0 %���������Ϊ��
    da.LUArray.ID = LUID; 
    da.LUArray.LWH = LULWH;   %LWHREAL ��ʵ�ߴ�    
    da.BinArray.LWH = BINLWH; %LWHREAL ��ʵ�ߴ�    
    % ���Ӽ�϶ -
    da.LUArray.BUFF = LUBUFF; %BUFF ����LU�ļ�϶
    da.BinArray.BUFF = BINBUFF; %BUFF ����BIN�ļ�϶
    % �������� -
    da.BinArray.Weight = 1000;
    da.LUArray.Weight = ones(size(LULWH,1),1);
    % ���������
    ParaArray.whichStripH = PARA_whichRotationHori(1);
    ParaArray.whichBinH = PARA_whichRotationHori(2);
    ParaArray.whichSortItemOrder = PARA_whichRotationHori(3);
    ParaArray.whichRotation = PARA_whichRotationHori(4);
    ParaArray.whichRotationHori = PARA_whichRotationHori(5);
    ParaArray.whichRotationAll = PARA_whichRotationHori(6);
else
%     da.LUArray.ID = [1 1 2 2];
%     da.BinArray.LWH = [5;20;4]';
%     da.LUArray.LWH = [2 2 3 3; 5 5 6 6; 4 4 4 4];
%     da.LUArray.BUFF = [0,0];
%     da.BinArray.BUFF = [0,0,0];
%     da.BinArray.Weight = 1000;
%     da.LUArray.Weight = da.LUArray.LWH(1,:);
    
% %         da.LUArray.ID = [5236934585 4 5236934585 5236934585];
% %     da.BinArray.LWH = [8;10;4];  %[6; 10; 4];
% %     da.LUArray.LWH = [3 2 3 3; 6 4 6 6; 2 2 1 2];
% %     
% %     da.LUArray.BUFF = [0,0]';
% %     da.BinArray.BUFF = [0,0,0]'; 
            
         % �������� -
%          da.BinArray.Weight = 1000;
%     da.LUArray.Weight = da.LUArray.LWH(1,:);
    %% �����������
     n=15;da=getRandDa(n);save('rndDa.mat','da');
     load('rndDa.mat')

    %%
    %     load insLU3.mat;   % load ins.mat;
    %     da.LUArray.ID = LUid;%     da.LUArray.LWH = w;%     da.BinArray.LWH = W;
    %     clear LUid w W;%     error('No input! ');
end
clc;
%% ����Input��������
printstruct(da);
da = GcheckInput(da,ParaArray);
%% ����ʽ: LU��Item���㷨
printstruct(da);
[da] = HLUtoItem(da,ParaArray);
%% ����ʽ��Item��Strip���㷨
printstruct(da);
[da] = HItemToStrip(da,ParaArray);
%% ����װ���� 
da = computeLoadingRate(da);
    function da = computeLoadingRate(da)
        % printstruct(da);
        % ��ʼ��
        nStrip = size(da.StripArray.LW,2);
        da.StripArray.Stripvolume = zeros(1,nStrip);
        da.StripArray.Itemvolume = zeros(1,nStrip);
        da.StripArray.Itemloadingrate = zeros(1,nStrip);
        
        % ����ÿ��strip��װ����
        da.StripArray.Stripvolume = da.StripArray.LW(2,:)*da.BinArray.LWH(1,1);      %ÿ��strip�Ŀ������
        a = da.ItemArray.LWH;
        b = da.ItemArray.itemBeStripMatrix;
        for iStrip =1:nStrip
            %ÿ��strip��װ�����
            da.StripArray.Itemvolume(iStrip)= sum(a(1, (b(1,:)==iStrip)) .* a(2, (b(1,:)==iStrip)));
        end
        %ÿ��strip��װ�ر���
        da.StripArray.Itemloadingrate =  da.StripArray.Itemvolume ./ da.StripArray.Stripvolume;
    end
%% ����ʽ��Strip��Bin���㷨
printstruct(da);
[da] = HStripToBin(da,ParaArray); %todo CHECK CHECK CHECK
%% Item��bin����Ϣ��ȡ:
printstruct(da);
[da] = HItemToBin(da); 
printstruct(da);
%% ����������������LWHRota:��ת���LWH,���Ǽ�ȥbuffer����)
da.LUArray.LWHRota = da.LUArray.LWH;
nLU = numel(da.LUArray.LURotaFlag);
for iLU=1:nLU
    if da.LUArray.LURotaFlag(iLU)
        da.LUArray.LWHRota(1,iLU)=da.LUArray.LWH(2,iLU);
        da.LUArray.LWHRota(2,iLU)=da.LUArray.LWH(1,iLU);
    end
end
printstruct(da);
%% ����������(ԭʼ˳��) ���4������
% ����1 - ��1:Bin��ţ���2����bin��˳��
Res1_LUBeBinMatrix=da.LUArray.LUBeBinMatrix;

% ����2 - LU��Bin�ڵ�����
% ���Ӽ�϶-����CoordLUBinWithBuff����
da.LUArray.CoordLUBinWithBuff = da.LUArray.CoordLUBin + da.LUArray.BUFF./2;
Res2_CoordLUBin=da.LUArray.CoordLUBinWithBuff; %Res2_CoordLUBin��DOUBLE����: Lu��xyzֵ TTTTTTTTTT

% ����3 - LU�ĳ����(��ת��)
% ���Ӽ�϶-�޶�LWHRotaΪ��СBuffer���ʵ�����ݱ���
da.LUArray.LWHRota = da.LUArray.LWHRota - da.LUArray.BUFF;
Res3_LWHRota=da.LUArray.LWHRota;  %Res3_LWHRota��DOUBLE LU�ĳ���ߣ���ת��

% ����4 - ��1��LU��Item��λ�ã���2��LU��ID����
LUBeItemArray=da.LUArray.LUBeItemArray(1,:);
LUID=da.LUArray.ID;
Res4_DrawSeq = [LUBeItemArray; LUID];

% Res5_BinLWH = da.BinArray.LWH; %��ȥBuffer��ʵ�ʿ��õĳ����
%% ����������(����˳��)
% % [~,x]=sort(Res1_LUBeBinMatrix(2,:));
% % Res1_LUBeBinMatrix=Res1_LUBeBinMatrix(:,x)
% % LUBeItemArray=LUBeItemArray(1,x)
% % Res2_CoordLUBin=Res2_CoordLUBin(:,x)
% % Res3_LWHRota=Res3_LWHRota(:,x)
% % fprintf('����������ȫ����� \n');
%% ��ͼ
if nargin == 0
%     printstruct(da);
%     plot3DBPP(da,ParaArray);
% �����޶���Ϊ��ͼʹ��
% plot2DBPP(da,ParaArray);hold on;
da.ItemArray.LWH = da.ItemArray.LWH - da.LUArray.BUFF(:,1:size(da.ItemArray.LWH,2));
% plot2DBPP(da,ParaArray);hold on;
da.ItemArray.CoordItemBin = da.ItemArray.CoordItemBin + da.LUArray.BUFF(:,1:size(da.ItemArray.LWH,2))/2;
plot2DBPP(da,ParaArray);
else
%     % �����޶���Ϊ��ͼʹ��
% plot2DBPP(da,ParaArray);hold on;
da.ItemArray.LWH = da.ItemArray.LWH - da.LUArray.BUFF(:,1:size(da.ItemArray.LWH,2));
% plot2DBPP(da,ParaArray);hold on;
da.ItemArray.CoordItemBin = da.ItemArray.CoordItemBin + da.LUArray.BUFF(:,1:size(da.ItemArray.LWH,2))/2;
plot2DBPP(da,ParaArray);
end

% printstruct(da);
% mcc -W 'java:BBA_Main,Class1,1.0' -T link:lib BBA_Main.m -d '.\new'
%% �����½�
% lb = computerLB(da); fprintf('LB = %d \n', lb); %��ĳ��bin����Ϊ׼
end





%% ********************** ������ts�㷨�Ĵ��� ��ʱ���� ****************

% [ub,x,b] = HnextFit(ItemArray,BinArray);
% [ub,x,b] = HnextFit_origin(ItemArray,BinArray);
% disp(b');
% disp(x);
% fprintf('UB = %d \n', ub);
% [ub,x,b] = TSpack(d,n,w,W,lb,timeLimit, ub0,x,b,whichH);

function [toReturn,x,b] = TSpack(d,n,w,W,lb,timeLimit, ub0,x,b,whichH)
[nb,x,b] = heur(d,n,w,W,x,b,whichH);
ub0 = nb;

if nb == lb
    toReturn = nb;
    fprintf('best lb = %d ', nb);
end

%/* initial (trivial) solution */
cnb = n;
cb = 1:n;
cx = zeros(d,n);
cw = w;

%/* external loop */
D = 1; tt = 0.0;
toReturn = nb;

end

function [nb,x,b] = heur(d,n,w,W,x,b,whichH)
nb = -1;
which = (d-2)*100 + whichH; % which Ϊ0��100
if which == 0
    [nb,x,b]  = HnextFit(n,w,W,x,b,n+1); %/* first heuristic for 2d bin packing */
    %          disp(x);
    %          disp(b);
elseif which == 100
    [nb,x,b]  = HHnextFit(n,w,W,x,b); %/* first heuristic for 3d bin packing */
end
end

    function [ub,px,pb]  = HnextFit(ItemArray,BinArray)
        % Initialize
        d = size(ItemArray.LWH,1)-1;
        n = size(ItemArray.LWH,2);
        nn = n + 1;
        w = ItemArray.LWH(1:d,:);
        W = BinArray.LWH(1:d,:);
        x = zeros(d,n); b = zeros(n,1); bNb = zeros(n,1);
        
        %/* sort the items */
        % sortD = size(w,1);%��ȡ��Ҫ�����ά��
        [~,ord] = sort(w(d,:),'descend');%��w��������,ֻ��Ҫ����˳��ord;����d�����򣨸߶�)
        pw = w(:,ord);
        px = x;
        pb = (b+999); % 0 + 999
        pbNb = bNb;
        
        %/* next fit packing */
        % binLeftArray(1,ub) �� wleft
        % binLeftArray(2,ub) :  hleft
        nBin = n;
        binLeftArray = repmat(W,1,nBin);  %��ʼ
        ub = 1;
        for i=1:n
            %     if (binLeftArray(1,ub) == W(1)) & (binLeftArray(2,ub) == W(2)) %����ǿ�bin
            if pbNb(ub) == 0   %����ǿ�bin
                if (pw(1,i) <= binLeftArray(1,ub)) && (pw(2,i) <= binLeftArray(2,ub)) %�����߶�������
                    px(1,i) = 0; px(2,i) = 0;
                    binLeftArray(1,ub) = binLeftArray(1,ub) - pw(1,i);
                    binLeftArray(2,ub) = binLeftArray(2,ub) - pw(2,i);
                    pbNb(ub) = pbNb(ub) + 1;
                else
                    error('EEE');
                end
            else               %������ǿ�bin
                if pw(1,i) <= binLeftArray(1,ub)  %���i�Ŀ����㵱ǰbin��ʣ���ȣ�ʣ��߶�Ӧ�ò���
                    px(1,i) = W(1) - binLeftArray(1,ub);
                    px(2,i) = W(2) - binLeftArray(2,ub) - pw(2,i);     %�߶�Ϊ????
                    binLeftArray(1,ub) = binLeftArray(1,ub) - pw(1,i);
                    binLeftArray(2,ub) = binLeftArray(2,ub);
                    pbNb(ub) = pbNb(ub) + 1;
                else
                    if pw(2,i)  <= binLeftArray(2,ub)  %���i�ĸ����㵱ǰbin��ʣ��߶�
                        px(1,i) = 0;
                        px(2,i) = binLeftArray(2,ub);
                        binLeftArray(1,ub) = W(1) - pw(1,i);
                        binLeftArray(2,ub) = binLeftArray(2,ub) - pw(2,i);
                        pbNb(ub) = pbNb(ub) + 1;
                    else  %���i�ĸ߲������㵱ǰbin��ʣ��߶�
                        ub = ub + 1;
                        px(1,i) = 0;   px(2,i) = 0;
                        pbNb(ub) = pbNb(ub) + 1;
                        binLeftArray(1,ub) = binLeftArray(1,ub) - pw(1,i);
                        binLeftArray(2,ub) = binLeftArray(2,ub) - pw(2,i);
                    end
                end
            end
            pb(i) = ub-1;
        end
        
        
        %ԭʼ��
        x = px(:,ord);
        b = pb(ord);
        bNb = pbNb(ord);
        ub = ub +1;
    end

    function [ub,px,pb]  = HnextFit_origin(ItemArray,BinArray)
        % Initialize
        d = size(ItemArray.LWH,1);
        if d==3
            d=d-1;
        end
        n = size(ItemArray.LWH,2);
        nn = n + 1;
        w = ItemArray.LWH(1:d,:);
        W = BinArray.LWH(1:d,:);
        x = zeros(d,n); b = zeros(n,1); bNb = zeros(n,1);
        
        %/* sort the items */
        sortD = size(w,1);%��ȡ��Ҫ�����ά��
        [~,ord] = sort(w(sortD,:),'descend');%��w��������,ֻ��Ҫ����˳��ord
        pw = w(:,ord);
        ord;
        
        px = zeros(size(x,1),size(x,2));
        pb = 999*ones(size(b,1),size(b,2));
        %/* next fit packing */
        hleft = W(2) - pw(2,1);
        wleft = W(1);
        ub = 0; hcurr = 0;
        for i=1:n  %�ӵ�һ��item��ʼ����
            if pw(1,i) <= wleft  %���item��w ��wleftС������item����bin���㣺����x1ֵ������wleft��hleft����
                px(1,i) = W(1) - wleft;
                wleft = wleft - pw(1,i);
            else    %��������һ�㰲�š�
                if pw(2,i) <= hleft  %���item��h ��hleftС������bin�߶ȳ��㣬����item����һ������������hleft������hcurr��wleft�� ��������xֵ������wleft
                    hcurr = W(2) - hleft; %������ͬһ�㣬����hcurr���䣬Ҳ����pw(2,1)(��������bin�Ͳ����ˣ�������hcurr)
                    hleft = hleft - pw(2,i);
                else  %����Ų��£�����bin������hcurr������hleft����������ub+1������ﵽnnֵ������������������xֵ0������wleft
                    hcurr = 0;    %�������µ�bin������hcurrΪ0;
                    hleft = W(2) - pw(2,i);
                    if (ub+1 == nn)
                        break;
                    end
                    ub = ub + 1;
                end
                % ���۷����ϲ����bin������x1Ϊ0������wleftΪW(1)-��item�Ŀ�w
                px(1,i) = 0;
                wleft = W(1) - pw(1,i);
            end
            % �˴�ͳһ����x1ֵ�����߶�ֵ��Ϊhcurr��ͳһ����bֵ=ub��
            px(2,i) = hcurr;
            pb(i) = ub;
        end
        
        %ԭʼ��
        x = px(:,ord);
        b = pb(ord);
        ub = ub +1;
    end


    function [ nb ] = HHnextFit(n,w,W,x,b)
        
        nb = 1;
        
    end

    function [ lb ] = computerLB(da)
        sum1 = sum(prod(da.ItemArray.LWH,1));
        % todo �����ж��Ƿ����е�BinArray�����е�bin����ͬ�� ����� �����ִ��
        sum2 = prod(da.BinArray.LWH(:,1));
        lb = ceil(sum1/sum2);
        if lb <=0, error('EEE');end
    end

%% ************************* ������ע�ʹ���  ************************

%% OLD lower�������
% % function [ lb ] = lower(d,n,w,W,whichL)
% % if whichL == 0
% %     lb = 1;
% % elseif whichL == 1
% %      sum1 = sum(prod(w,1));
% %      sum2 = prod(W);
% %      lb = ceil(sum1/sum2);
% % end
% %      if lb <=0, error('EEE');end
% % end

%% �ṹ�������strip�㷨

% % %% function [StripSolutionSort] = HnextFitDH(da)
% % function [StripSolutionSort] = HnextFitDH(da)
% % % ����: da
% % % ���: StripSolutionSort
% % %% ��ȡ������bin,��άitem����
% % % nDim nItem nBin
% % % itemDataMatrix uniBinDataMatrix
% % nDim = size(da.ItemArray.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % itemDataMatrix = da.ItemArray.LWH(1:nDim,:);
% % tmpbinDataMatrix = da.BinArray.LWH(1:nDim,:);
% % uniBinDataMatrix = unique(tmpbinDataMatrix','rows')';
% % nItem = size(itemDataMatrix,2);  nBin = nItem;
% % if size(uniBinDataMatrix,2)==1
% %     fprintf('������ֻ��һ������ ��=%1.0f ��=%1.0f  \n', uniBinDataMatrix);
% %     fprintf('�������� %d ����Ʒ,����ֱ�Ϊ \n',nItem);
% %     fprintf('%1.0f %1.0f \n',itemDataMatrix);
% % else
% %     error('�������ж������,�������� \n');
% % end
% % %% ���StripSolutionSort��ʼ��
% % % StripSolutionSort: stripWidth stripDataMatrix stripBeItemArray
% % % StripSolutionSort: itemOrd itemDataMatrixSort itemCoordMatrixSort itemBeLevelMatrixSort
% % StripSolutionSort.stripWidth = uniBinDataMatrix (1,1); %ֻ��Ҫstrip�Ŀ��,dim1Ϊ���
% % StripSolutionSort.stripDataMatrix = zeros(2,nBin);  %dim2-��(��)��(����ߵļ���)
% % StripSolutionSort.stripDataMatrix(1,:) = StripSolutionSort.stripWidth; %dim1-���ʣ�� ;
% % StripSolutionSort.stripBeItemArray = zeros(1,nBin); %ĳ��strip�������ٸ�item,�����Ų�����
% % 
% %  %/* sort the items */
% % [~,itemOrd] = sort(itemDataMatrix(nDim,:),'descend'); %��w��������,ֻ��Ҫ����˳��ord;����nDim�����򣨳�/�߶�)
% % StripSolutionSort.itemOrd = itemOrd;
% % StripSolutionSort.itemDataMatrixSort = itemDataMatrix(:,itemOrd);
% % StripSolutionSort.itemCoordMatrixSort = zeros(nDim,nItem);
% % StripSolutionSort.itemBeStripMatrixSort = zeros(2,nItem); %dim1:���ڵڼ���level dim2:���ڸ�level�ڼ����ŷ�
% % %% NFѭ��
% % iLevel = 1; iItem = 1;
% % %/* next fit packing */
% % while 1
% %     if iItem > nItem, break; end
% %     % ��ͬ�����µ�ѡ�������ǰitem�Ŀ�<=��ǰstrip�ĵ�ǰlevel�Ŀ�
% %     flag = StripSolutionSort.itemDataMatrixSort(1,iItem) <= StripSolutionSort.stripDataMatrix(1,iLevel);
% %     if ~isempty(flag)
% %         thisLevel = iLevel;
% %         [StripSolutionSort] = insertItemToStrip(thisLevel,iItem,StripSolutionSort);
% %         iItem = iItem + 1;            
% %     else
% %         iLevel = iLevel + 1;% �����Ȳ����㣬��level����
% %     end
% % end
% % printstruct(StripSolutionSort);
% % end
% % 
% % %% function [StripSolutionSort] = HfirstFitDH(da)
% % function [StripSolutionSort] = HfirstFitDH(da)
% % % ����: da
% % % ���: StripSolutionSort
% % %% ��ȡ������bin,��άitem����
% % % nDim nItem nBin
% % % itemDataMatrix uniBinDataMatrix
% % nDim = size(da.ItemArray.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % itemDataMatrix = da.ItemArray.LWH(1:nDim,:);
% % tmpbinDataMatrix = da.BinArray.LWH(1:nDim,:);
% % uniBinDataMatrix = unique(tmpbinDataMatrix','rows')';
% % nItem = size(itemDataMatrix,2);  nBin = nItem;
% % if size(uniBinDataMatrix,2)==1
% %     fprintf('������ֻ��һ������ ��=%1.0f ��=%1.0f  \n', uniBinDataMatrix);
% %     fprintf('�������� %d ����Ʒ,����ֱ�Ϊ \n',nItem);
% %     fprintf('%1.0f %1.0f \n',itemDataMatrix);
% % else
% %     error('�������ж������,�������� \n');
% % end
% % %% ���StripSolutionSort��ʼ��
% % % StripSolutionSort: stripWidth stripDataMatrix stripBeItemArray
% % % StripSolutionSort: itemOrd itemDataMatrixSort itemCoordMatrixSort itemBeLevelMatrixSort
% % StripSolutionSort.stripWidth = uniBinDataMatrix (1,1); %ֻ��Ҫstrip�Ŀ��,dim1Ϊ���
% % StripSolutionSort.stripDataMatrix = zeros(2,nBin);  %dim2-��(��)��(����ߵļ���)
% % StripSolutionSort.stripDataMatrix(1,:) = StripSolutionSort.stripWidth; %dim1-���ʣ�� ;
% % StripSolutionSort.stripBeItemArray = zeros(1,nBin); %ĳ��strip�������ٸ�item,�����Ų�����
% % 
% %  %/* sort the items */
% % [~,itemOrd] = sort(itemDataMatrix(nDim,:),'descend'); %��w��������,ֻ��Ҫ����˳��ord;����nDim�����򣨳�/�߶�)
% % StripSolutionSort.itemOrd = itemOrd;
% % StripSolutionSort.itemDataMatrixSort = itemDataMatrix(:,itemOrd);
% % StripSolutionSort.itemCoordMatrixSort = zeros(nDim,nItem);
% % StripSolutionSort.itemBeStripMatrixSort = zeros(2,nItem); %dim1:���ڵڼ���level dim2:���ڸ�level�ڼ����ŷ�
% % %% FFѭ��
% % iLevel = 1; iItem = 1;
% % %/* next fit packing */
% % while 1
% %     if iItem > nItem, break; end
% %     % ��ͬ�����µ�ѡ�����find����㹻�Ķ��level,�������ڵ�һ�������� Ψһ������thisLevel�Ļ�ȡ
% %     flag = find(StripSolutionSort.stripDataMatrix(1,1:iLevel) >= StripSolutionSort.itemDataMatrixSort(1,iItem));
% %     if ~isempty(flag)
% %         thisLevel = flag(1);
% %         [StripSolutionSort] = insertItemToStrip(thisLevel,iItem,StripSolutionSort);
% %         iItem = iItem + 1;
% %     else
% %         iLevel = iLevel + 1;% �����Ȳ����㣬��level����
% %     end
% % end
% % printstruct(StripSolutionSort);
% % end
% % 
% % %% function [StripSolutionSort] = HbestFitDH(da)
% % function [StripSolutionSort] = HbestFitDH(da)
% % % ����: da
% % % ���: StripSolutionSort
% % %% ��ȡ������bin,��άitem����
% % % nDim nItem nBin
% % % itemDataMatrix uniBinDataMatrix
% % nDim = size(da.ItemArray.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % itemDataMatrix = da.ItemArray.LWH(1:nDim,:);
% % tmpbinDataMatrix = da.BinArray.LWH(1:nDim,:);
% % uniBinDataMatrix = unique(tmpbinDataMatrix','rows')';
% % nItem = size(itemDataMatrix,2);  nBin = nItem;
% % if size(uniBinDataMatrix,2)==1
% %     fprintf('������ֻ��һ������ ��=%1.0f ��=%1.0f  \n', uniBinDataMatrix);
% %     fprintf('�������� %d ����Ʒ,����ֱ�Ϊ \n',nItem);
% %     fprintf('%1.0f %1.0f \n',itemDataMatrix);
% % else
% %     error('�������ж������,�������� \n');
% % end
% % %% ���StripSolutionSort��ʼ��
% % % StripSolutionSort: stripWidth stripDataMatrix stripBeItemArray
% % % StripSolutionSort: itemOrd itemDataMatrixSort itemCoordMatrixSort itemBeLevelMatrixSort
% % StripSolutionSort.stripWidth = uniBinDataMatrix (1,1); %ֻ��Ҫstrip�Ŀ��,dim1Ϊ���
% % StripSolutionSort.stripDataMatrix = zeros(2,nBin);  %dim2-��(��)��(����ߵļ���)
% % StripSolutionSort.stripDataMatrix(1,:) = StripSolutionSort.stripWidth; %dim1-���ʣ�� ;
% % StripSolutionSort.stripBeItemArray = zeros(1,nBin); %ĳ��strip�������ٸ�item,�����Ų�����
% % 
% %  %/* sort the items */
% % [~,itemOrd] = sort(itemDataMatrix(nDim,:),'descend'); %��w��������,ֻ��Ҫ����˳��ord;����nDim�����򣨳�/�߶�)
% % StripSolutionSort.itemOrd = itemOrd;
% % StripSolutionSort.itemDataMatrixSort = itemDataMatrix(:,itemOrd);
% % StripSolutionSort.itemCoordMatrixSort = zeros(nDim,nItem);
% % StripSolutionSort.itemBeStripMatrixSort = zeros(2,nItem); %dim1:���ڵڼ���level dim2:���ڸ�level�ڼ����ŷ�
% % %% FFѭ��
% % iLevel = 1; iItem = 1;
% % %/* next fit packing */
% % while 1
% %     if iItem > nItem, break; end
% %     % ��ͬ�����µ�ѡ����� find����㹻�Ķ��level,����������Сʣ���ȵ� 
% %     flag = find(StripSolutionSort.stripDataMatrix(1,1:iLevel) >= StripSolutionSort.itemDataMatrixSort(1,iItem));
% %     if ~isempty(flag)
% %         % Ψһ��FF������⵽thisLevel�ļ��㣨ѡ��������������С��
% %         tepMin = StripSolutionSort.stripDataMatrix(1,1:iLevel);
% %         tepMin = min(tepMin(flag));
% %         thisLevel = find(StripSolutionSort.stripDataMatrix(1,1:iLevel)==tepMin);
% %         if length(thisLevel)>1
% %             thisLevel = thisLevel(1);
% %         end 
% %         
% %         [StripSolutionSort] = insertItemToStrip(thisLevel,iItem,StripSolutionSort);
% %         iItem = iItem + 1;
% %     else
% %         iLevel = iLevel + 1;% �����Ȳ����㣬��level����
% %     end
% % end
% % printstruct(StripSolutionSort);
% % end
% % 
% % 
% % 
% % 

%% �ǽṹ���HfirstFitDH2�㷨
% % function [stripLeftMatrix,pbelongMatrix,pitemMatrix,pcoordMatrix,ord]  = HfirstFitDH2(ItemArray,BinArray)
% % % ���������ʼ��
% % nDim = size(ItemArray.LWH,1);
% % if nDim ==3, nDim = nDim-1;end
% % nItem = size(ItemArray.LWH,2);
% % nBin = nItem;
% % % nn = n + 1;
% % itemMatrix = ItemArray.LWH(1:nDim,:);
% % binMatrix = BinArray.LWH(1:nDim,:);
% % % ���������ʼ��
% % coordMatrix = zeros(nDim,nItem);
% % stripWidth = binMatrix(1,1); %ֻ��Ҫstrip�Ŀ��,dim1Ϊ���
% % stripLeftMatrix = [stripWidth*(ones(1,nBin));zeros(1,nBin);zeros(1,nBin)];%��ʼ��strip: dim1-���ʣ�� ; dim2-��(��)��(����ߵļ���); dim3-��strip������item����
% % belongMatrix = zeros(2,nItem); %dim1:���ڵڼ���level dim2:���ڸ�level�ڼ����ŷ�
% % 
% % %/* sort the items */
% % [~,ord] = sort(itemMatrix(nDim,:),'descend');%��w��������,ֻ��Ҫ����˳��ord;����d�����򣨸߶�)
% % pitemMatrix = itemMatrix(:,ord);
% % pcoordMatrix = coordMatrix;
% % pbelongMatrix = belongMatrix;
% % iLevel = 1; iItem = 1;  
% % 
% % %%
% % %/* first fit packing */
% % while 1
% %     if iItem > nItem, break; end
% %     % find����㹻�Ķ��level,�������ڵ�һ�������� Ψһ������⵽thisLevel�ļ��� + �����thisLevel���滻
% %     findLevelArray = find(stripLeftMatrix(1,1:iLevel) >= pitemMatrix(1,iItem));
% %     if findLevelArray
% %         thisLevel = findLevelArray(1);
% %         [pcoordMatrix,stripLeftMatrix,pbelongMatrix] = insertItemToStrip2(thisLevel,iItem,pitemMatrix,stripWidth,pcoordMatrix,stripLeftMatrix,pbelongMatrix);
% %         iItem = iItem + 1;        
% %     % �����Ȳ����㣬��level����
% %     else
% %         iLevel = iLevel + 1;
% %     end
% % end
% %      pcoordMatrix
% %      stripLeftMatrix
% %      pbelongMatrix
% % end

%% �ǽṹ���HbesttFitDH2�㷨
% % function [stripLeftMatrix,pbelongMatrix,pitemMatrix,pcoordMatrix,ord]  = HbestFitDH2(ItemArray,BinArray)
% % % ���������ʼ��
% % nDim = size(ItemArray.LWH,1);
% % if nDim ==3, nDim = nDim-1;end
% % nItem = size(ItemArray.LWH,2);
% % nBin = nItem;
% % % nn = n + 1;
% % itemMatrix = ItemArray.LWH(1:nDim,:);
% % binMatrix = BinArray.LWH(1:nDim,:);
% % % ���������ʼ��
% % coordMatrix = zeros(nDim,nItem);
% % stripWidth = binMatrix(1,1); %ֻ��Ҫstrip�Ŀ��,dim1Ϊ���
% % stripLeftMatrix = [stripWidth*(ones(1,nBin));zeros(1,nBin);zeros(1,nBin)];%��ʼ��strip: dim1-���ʣ�� ; dim2-��(��)��(����ߵļ���); dim3-��strip������item����
% % belongMatrix = zeros(2,nItem); %dim1:���ڵڼ���level dim2:���ڸ�level�ڼ����ŷ�
% % 
% % %/* sort the items */
% % [~,ord] = sort(itemMatrix(nDim,:),'descend');%��w��������,ֻ��Ҫ����˳��ord;����d�����򣨸߶�)
% % %         ord = 1:nItem;    % �����Ŀ�Ĳ���items����
% % pitemMatrix = itemMatrix(:,ord);
% % pcoordMatrix = coordMatrix;
% % pbelongMatrix = belongMatrix;
% % iLevel = 1; iItem = 1;  
% % 
% % %%
% % %/* best fit packing */
% % while 1
% %     if iItem > nItem, break; end
% %     % find����㹻�Ķ��level,����������Сʣ���ȵ� 
% %     findLevelArray = find(stripLeftMatrix(1,1:iLevel) >= pitemMatrix(1,iItem));
% %     if findLevelArray
% % %         Ψһ��FF������⵽thisLevel�ļ��㣨ѡ��������������С��
% %         tepMin = stripLeftMatrix(1,1:iLevel);
% %         tepMin = min(tepMin(findLevelArray));
% %         thisLevel = find(stripLeftMatrix(1,1:iLevel)==tepMin);
% %         if length(thisLevel)>1
% %             thisLevel = thisLevel(1);
% %         end
% %         [pcoordMatrix,stripLeftMatrix,pbelongMatrix] = insertItemToStrip2(thisLevel,iItem,pitemMatrix,stripWidth,pcoordMatrix,stripLeftMatrix,pbelongMatrix);
% %         iItem = iItem + 1;
% %         
% %     % �����Ȳ����㣬��level����
% %     else
% %         iLevel = iLevel + 1;
% %     end
% % end
% %      pcoordMatrix
% %      stripLeftMatrix
% %      pbelongMatrix
% % 
% % end

%% �ǽṹ���HbestFitBinDH�㷨
% % function [pbelongItemBinMatrix,pbelongStripBinMatrix,pcoordItemBinMatrix,binLeftMatrix ] = HbestFitBinDH(stripLeftMatrix,pbelongMatrix,pitemMatrix,pcoordMatrix,ItemArray,BinArray)
% % % ���������ʼ��
% % nDim = size(ItemArray.LWH,1);
% % if nDim == 3, nDim = nDim-1;end
% % nItem = size(ItemArray.LWH,2);
% % nBin = nItem;
% % 
% % nStrip = sum(stripLeftMatrix(3,:)>0); %����ʹ�õ�Strip������
% % % nn = n + 1;
% % binMatrix = BinArray.LWH(1:nDim,:);
% % 
% % stripWidth = binMatrix(1,1); %ֻ��Ҫstrip�Ŀ��,dim1Ϊ���
% % % ���������ʼ��
% % pbelongItemBinMatrix = zeros(2,nItem); % dim1:���item��ĳ��bin dim2:����˳��
% % pbelongStripBinMatrix = zeros(2,nBin); % dim1:���strip��ĳ��bin dim2:����˳��
% % pcoordItemBinMatrix = zeros(nDim,nItem); %����ֵ
% % binLeftMatrix = [binMatrix(1,1)*(ones(1,nBin));binMatrix(2,1)*ones(1,nBin);zeros(1,nBin);zeros(1,nBin)];
% % %��ʼ��bin: dim1-bin���ʣ�� ; dim2-bin��(��)��(555ʣ�ࣩ; dim3-��bin������item����; dim4-��bin������strip����;
% % 
% % %/* sort the strips by ��(��) Ĭ���Ѿ��ǰ����˳�� ������������ */
% % 
% % % Best fit 
% % iStrip=1;iBin=1;
% % while 1
% %     if iStrip > nStrip, break; end
% %     findBinArray = find(binLeftMatrix(2,1:iBin) >= stripLeftMatrix(2,iStrip));
% %     if findBinArray
% %         tepMin = binLeftMatrix(2,1:iBin);
% %         tepMin = min(tepMin(findBinArray)); % 555 check
% %         thisBin = find(binLeftMatrix(2,1:iBin)==tepMin);
% %         if length(thisBin)>1
% %             thisBin = thisBin(1);
% %         end
% %         %����strip������Ϣ
% %         pbelongStripBinMatrix(1,iStrip) = thisBin;
% %         binLeftMatrix(4,thisBin) = binLeftMatrix(4,thisBin) + 1; %��bin�µڼ��ΰ���strip
% %         pbelongStripBinMatrix(2,iStrip) = binLeftMatrix(4,thisBin);
% %         
% %         %��ȡ��iStrip�ڵ�item���, ������Item������Ϣ
% %         idxItemStrip = find(pbelongMatrix(1,:)==iStrip);
% %         pbelongItemBinMatrix(1,idxItemStrip) = thisBin;    %�ڼ���bin
% %         
% %         %����bin����Ϣ
% %         binLeftMatrix(1,thisBin) = min(binLeftMatrix(1,thisBin),stripLeftMatrix(1,iStrip)); %����ʣ���ȵ���Сֵ
% %         binLeftMatrix(2,thisBin) = binLeftMatrix(2,thisBin) - stripLeftMatrix(2,iStrip); %����ʣ��߶�
% %         binLeftMatrix(3,thisBin) = binLeftMatrix(3,thisBin) + length(idxItemStrip); %��bin�ºϼƼ���item
% %         
% %         
% %         %����xy������Ϣ x���� yͨ��bin�߶�-binʣ��߶�-����strip�߶�
% %         pcoordItemBinMatrix(1,idxItemStrip) = pcoordMatrix(1,idxItemStrip);
% %         pcoordItemBinMatrix(2,idxItemStrip) = binMatrix(2,1) - (binLeftMatrix(2,thisBin) + stripLeftMatrix(2,iStrip));
% %      
% %         iStrip = iStrip + 1;
% %     else
% %         iBin = iBin + 1;
% %     end    
% % end
% % 
% % % ���Ӹ���pbelongItemBinMatrix�з���˳��Ĳ���
% % for iItem=1:nItem
% %     tmp = find(pbelongItemBinMatrix(1,:)==iItem);
% %     if isempty(tmp),  break;   end
% %     for i=1:length(tmp)
% %         pbelongItemBinMatrix(2,tmp(i)) = pbelongItemBinMatrix(2,tmp(i)) + i;
% %     end
% % end
% % 
% % pbelongStripBinMatrix
% % pbelongItemBinMatrix
% % binLeftMatrix
% % pcoordItemBinMatrix
% % end
