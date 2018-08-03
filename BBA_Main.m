%% BBA_MAIN demo
%% Form
%    [output_LU_Bin,output_CoordLUBin,output_LU_LWH,Res4_DrawSeq] = ... 
%    BBA_Main(LUID,LULWH,BINLWH,PARANOUSE,LUBUFF,BINBUFF,LUWEIGHT,BINWEIGHT,LUISROTA)
%
%% Description
%   BBA_Main.
%
%% Inputs
%   LUID	                (1,n)   �������� ��ͬ���ֱ���ͬһ����,����Ѷ� 
%   LULWH                (3,n)   ���̿��� 
%   BINLWH              (3,1)   ���Ϳ��ߣ�Ŀǰ�����ǵ����� TODO �������ӵ��೵�Ͱ�˳��ʹ�ã�
%   PARANOUSE       (1,1)   ����(1*1) ����rotation,������rotation�Ĳ���
%   LUBUFF               (2,1)   ���̼�����ܼ�϶(2���ĵ��߿���϶) �������̳����=ÿ�����̵�ʵ�ʳ����+���ӵ�buff
%   LUWEIGHT           (1,n)  ��������
%   BINWEIGHT         (1,1)  ��������������
%   LUISROTA            (1,n)  �����Ƿ�������ת
%   LUZWID               (1,n)  �����Ƿ���������
%
%% Outputs
%   output_LU_Bin	(2,n)	 ��1: LU��ĳ��BIN�ڣ���2: LU�ڸ�BIN�ڵİ���˳��
%   output_CoordLUBin      (3,n)    ÿ��LU��X,Y,Z
%   output_LU_LWH          (3,n)    ÿ��LU�Ŀ��ߣ���ת��ģ�ʵ��ֵ��
%   output_LU_Item_ID         (2,n)   ��1: LU��ĳ��ITEM�ڣ���2: LU������ID����
%

function [output_LU_Bin,output_CoordLUBin,output_LU_LWH,output_LU_Item_ID] = ...
    BBA_Main(varargin)
% function [output_LU_Bin,output_CoordLUBin,output_LU_LWH,output_LU_Item_ID] = ...
%     BBA_Main(LUID,LULWH,BINLWH,PARANOUSE,LUBUFF,BINBUFF,LUWEIGHT,BINWEIGHT,LUISROTA)

%% Initialize Data Structure
% clear;close all; format long g; format bank; %NOTE ����MATLAB CODE ֧��
% rng('default');rng(1); % NOTE �Ƿ�����ı�־
close all
if nargin ~= 0
    if nargin <= 6 
        d = DataInitialize( ...
            'LUID', varargin{1},...
            'LULWH',varargin{2}, ...
            'LUBUFF',varargin{5}, ...
            'BINLWH',varargin{3}, ...
            'BINBUFF', varargin{6});
        
    else
        d = DataInitialize( ...
            'LUID', varargin{1},...
            'LULWH',varargin{2}, ...
            'LUBUFF',varargin{5}, ...
            'LUWEIGHT', varargin{7}, ...
            'LUISROTA', varargin{9}, ...
            'BINLWH',varargin{3}, ...
            'BINBUFF', varargin{6}, ...
            'BINWEIGHT', varargin{8});
    end
else
    d = DataInitialize(5); %0 Ĭ��ֵ; >0 �����������n������
end

%% Initialize Parameter
nAlg = 1;
for i = 1:3 %1-3 best first next����
    for j=1:1 %1-2 ����:1 �߶ȣ��������� 2 ��̱�
        for k=1:1 %0-2 Ĭ��0 ������ת 1����ת 2: ����Ϊ�����Ƿ�����Rotation 
            for l=1:2 %0-2 0��ȡ�� ����1-2 RotaHori 1hori 2 vert 555 ��Ų��˻��ݷţ��������ݷź󲻻��ţ��Ų��£���
                for m=1:3 %1-3 best first next����
                % pA nAlg 
                pA(nAlg) = ParameterInitialize( ...
                             'whichStripH', i,...
                             'whichBinH',m, ...
                             'whichSortItemOrder',j, ...
                             'whichRotation',k, ...
                             'whichRotationHori', l);
                 nAlg=nAlg+1;
                end
            end
        end
    end
end
nAlg = nAlg - 1;

%% Simulate - All ALGORITHM

fprintf(1,'\nRunning the simulation...\n');

% Run ALL algorithm configure
for iAlg = 1:nAlg
    dA(iAlg) = RunAlgorithm(d,pA(iAlg));        %��ȡ���н�ṹ��
    
%     plotSolution(dA(iAlg),pA(iAlg));
    
    % �㷨�ж��Ƿ���ͬ�����������ڰڷ�
    flagA(iAlg) =  isAdjacent(dA(iAlg));
end

%  printstruct(dA(1,1),'sortfields',0,'PRINTCONTENTS',1)

%% Simulate - CHOOSE BEST ONE
% 555 �㷨�����ų�bin����ͬ�������̲����ڵĽ�
% % dA = dA(1,logical(flagA));
% % pA = pA(1,logical(flagA));

% �Ӷ���㷨�����ѡ���ӱض�bin�����ڵ����Ž��
if ~isempty(dA)
    [daBest,paBest] = getbestsol(dA,pA); 
else
    error('�����������нⶼ�������̲����ڵ���� \n');
end

%% POST PROCESSING
% Return length(parMax) �� solutions to BBA
if ~isempty(daBest)
    bestOne = 1;
    getReturnBBA(daBest(bestOne)); %���ж��,���ص�һ��
    plotSolution(daBest(bestOne),paBest(bestOne)); % == plotSolution(RunAlgorithm(d,paBest(bestOne)) ,paBest(bestOne));    
else
    error('��������δ�ҳ����Žⷵ��BBA \n');
end

fprintf(1,'Simulation done.\n');

% mcc -W 'java:BBA_Main,Class1,1.0' -T link:lib BBA_Main.m -d '.\new'
% d = rmfield(d, {'Veh', 'LU'});

%END MAIN

      
%% ******* Ƕ�׺���  **********

    function getReturnBBA(daMax) 
        % ����������(ԭʼ˳��) ���4������
        % ����1 - ��1:Bin��ţ���2����bin��˳��
        output_LU_Bin=daMax.LU.LU_Bin;
        
        % ����2 - LU��Bin�ڵ�����
        % ���Ӽ�϶-����CoordLUBinWithBuff����
        daMax.LU.CoordLUBinWithBuff = daMax.LU.CoordLUBin + daMax.LU.BUFF./2;
        output_CoordLUBin=daMax.LU.CoordLUBinWithBuff; %output_CoordLUBin��DOUBLE����: Lu��xyzֵ TTTTTTTTTT
        
        % ����3 - LU�ĳ����(��ת��)
        % ���Ӽ�϶-�޶�LWHΪ��С�����ӦBuffer���ʵ�����ݱ���
        % ������V2      
        daMax.LU.LWHOriRota = daMax.LU.LWH - daMax.LU.BUFF;
        output_LU_LWH=daMax.LU.LWHOriRota;  %output_LU_LWH��DOUBLE LU�ĳ���ߣ���ת��ʵ��ֵ��
            % ������V1
            %         daMax.LU.LWHRota = daMax.LU.LWHRota - daMax.LU.BUFF;
            %         Res3_LWHRota=daMax.LU.LWHRota;  %Res3_LWHRota��DOUBLE LU�ĳ���ߣ���ת��

       
        % ����4 - ��1��LU��Item��λ�ã���2��LU��ID����
        LU_Item=daMax.LU.LU_Item(1,:);
        LUID=daMax.LU.ID;
        output_LU_Item_ID = [LU_Item; LUID];
        
        % Res5_BinLWH = d.Veh.LWH; %��ȥBuffer��ʵ�ʿ��õĳ����
        % ����������(����˳��)
        % % [~,x]=sort(output_LU_Bin(2,:));
        % % output_LU_Bin=output_LU_Bin(:,x)
        % % LU_Item=LU_Item(1,x)
        % % output_CoordLUBin=output_CoordLUBin(:,x)
        % % Res3_LWHRota=Res3_LWHRota(:,x)
        % % fprintf('����������ȫ����� \n');
        end


end



%% ******* �ֲ����� ****************

%%DATAITIALIZE Initialize the DATA data structure.
%
%% Form
%  d = DataInitialize( varargin )
%
%% Description
% Initializes the INPUT data structure using parameter pairs.
%
%% Inputs
% varargin:  ('parameter',value,...)
% varargin:  (value) value's random LU input
%
%% Outputs
%   d	(.) DATA data structure
function d = DataInitialize( varargin )

% Defaults
d.LU.ID = [1 1 2 2];
d.LU.LWH = [2 2 3 3; 5 5 6 6; 4 4 4 4]';
d.Veh.LWH = [5;20;4]';

%     d.LU.ID = [1 1 2 2];
%     d.Veh.LWH = [5;20;4]';
%      d.LU.isRota = [1 1 1 1];
%     d.LU.LWH = [2 2 3 3; 5 5 6 6; 4 4 4 4];
    
% %     d.LU.ID = [5236934585 4 5236934585 5236934585];
% %     d.LU.isRota = [1 1 1];
% %     d.Veh.LWH = [8;10;4];  %[6; 10; 4];
% %     d.LU.LWH = [3 2 3 3; 6 4 6 6; 2 2 1 2];

d.LU.BUFF = [0,0];
d.Veh.BUFF = [0,0,0];
d.Veh.Weight = 1000;

n = length(varargin);

if n == 1 && varargin{1}~=0
    d=getRandDa(varargin{1});
%     save('rndDa.mat','d');
%     load('rndDa.mat');
end

for k = 1 : 2 : length(varargin)
    switch varargin{k}
        case 'LUID'
            d.LU.ID                        = varargin{k+1};
        case 'LULWH'
            d.LU.LWH                    = varargin{k+1};
        case 'LUBUFF'
            d.LU.BUFF                   = varargin{k+1};
        case 'LUWEIGHT'
            d.LU.Weight                = varargin{k+1};            
        case 'LUISROTA'
            d.LU.isRota                = varargin{k+1};            
        case 'BINLWH'
            d.Veh.LWH                  = varargin{k+1};
        case 'BINBUFF'
            d.Veh.BUFF                 = varargin{k+1};
        case 'BINWEIGHT'
            d.Veh.Weight              = varargin{k+1};            
    end
end

% ����ΪĬ�����ֵ
d.LU.Weight = randi([1,10],1,numel(d.LU.ID));
d.LU.isRota = randi([0,1],1,numel(d.LU.ID));

end

%%PARAMETERNITIALIZE Initialize the parameter data structure.
%
%% Form
%  p = ParameterInitialize( varargin )
%
%% Description
% Initializes the algorithm's parameter data structure using parameter pairs.
%
%% Inputs
% varargin:  ('parameter',value,...)
%
% 'whichStripH'                             (1,1) 
% 'whichBinH'                               (1,1) 
% 'whichSortItemOrder'               (1,1)
% 'whichRotation'                         (1,1)
% 'whichRotationHori'                  (1,1)
%
%% Outputs
%   p	(.)  Data structure
function p = ParameterInitialize( varargin )

% Defaults
p.whichStripH = 1;
p.whichBinH = 1;
p.whichSortItemOrder = 1;
p.whichRotation = 1;
p.whichRotationHori = 1;

n = length(varargin);

for k = 1:2:length(varargin)
    switch varargin{k}
        case 'whichStripH'
            p.whichStripH                        = varargin{k+1};
        case 'whichBinH'
            p.whichBinH                         = varargin{k+1};
        case 'whichSortItemOrders'
            p.whichSortItemOrder          = varargin{k+1};
        case 'whichRotation'
            p.whichRotation                 = varargin{k+1};
        case 'whichRotationHori'
            p.whichRotationHori         = varargin{k+1};
    end
end

end


function plotSolution(d,par)
%% ��ͼ
%     printstruct(d);
%     plot3DBPP(d,ParaArray);
% �����޶���Ϊ��ͼʹ��

fields = fieldnames(par);
aField = [];
for idx = 1:length(fields), aField = [aField par.(fields{idx})];   end
figure('name',num2str(aField));
d.Item.LWH = d.Item.LWH - d.LU.BUFF(:,1:size(d.Item.LWH,2));
d.Item.CoordItemBin = d.Item.CoordItemBin + d.LU.BUFF(:,1:size(d.Item.LWH,2))/2;
plot2DBPP(d,par);
end

function [d] = RunAlgorithm(d,p)
        
        %% ����Input��������
        d = GcheckInput(d,p);
        %% ����ʽ: LU��Item���㷨
        printstruct(d);
        [d.LU,d.Item,d.ItemID] = HLUtoItem(d.LU,d.Veh); %Item����ID������򣨵���һ�������仯˳��
        printstruct(d);
        %% �����½�
        lb = computerLB(d.Item,d.Veh);   fprintf('LB = %d \n', lb); %��ĳ��bin����Ϊ׼
        %% ����ʽ��Item��Strip���㷨
%         printstruct(d);
%         printstruct(d.Item);
        [d] = HItemToStrip(d,p);
        %% ����stripװ����
%         printstruct(d);
        d = computeLoadingRateStrip(d);
        function d = computeLoadingRateStrip(d)
            % ��ʼ��
            nStrip = size(d.Strip.LW,2);
            d.Strip.Stripvolume = zeros(1,nStrip);
            d.Strip.StripvolumeLimit = zeros(1,nStrip);
            d.Strip.Itemvolume = zeros(1,nStrip);
            d.Strip.loadingrate = zeros(1,nStrip);
            d.Strip.loadingrateLimit = zeros(1,nStrip);
            
            % ����ÿ��strip��װ����
            %ÿ��strip�Ŀ������ = �߶�*���(�����Ŀ��)
            d.Strip.Stripvolume = d.Strip.LW(2,:)*d.Veh.LWH(1,1);
            %ÿ��strip�����޿������ = �߶�*���(stripʹ�ÿ��=�������-stripʣ����)
            d.Strip.StripvolumeLimit = d.Strip.LW(2,:) .* (d.Veh.LWH(1,1) - d.Strip.LW(1,:));
            a = d.Item.LWH;
            b = d.Item.Item_Strip;
            for iStrip =1:nStrip
                %ÿ��strip��װ�����
                d.Strip.Itemvolume(iStrip)= sum(a(1, (b(1,:)==iStrip)) .* a(2, (b(1,:)==iStrip)));
            end
            %ÿ��strip��װ�ر���
            d.Strip.loadingrate =  d.Strip.Itemvolume ./ d.Strip.Stripvolume;
            %ÿ��strip������װ�ر���
            d.Strip.loadingrateLimit =  d.Strip.Itemvolume ./ d.Strip.StripvolumeLimit;
        end
        %% ��Strip�н���һ���Ҹ�>���Item����ѡ�񲢸�����Ӧ����
         d = modifyStripWithOneItem(d);
        function d = modifyStripWithOneItem(d)
            stripheight = d.Strip.LW(2,:);
            binwidth = d.Veh.LWH(1,1);
            stripleftwidth = d.Strip.LW(1,:);
            stripwidth = ( binwidth - stripleftwidth );
            [tmpset] = find(stripheight > stripwidth);
            if ~isempty(tmpset)
                if isscalar(tmpset) %�Ը�strip�����ڲ�����1��Item����,��������漰CoordItemStrip
                    d.Strip.LW(:,tmpset) = [binwidth-stripheight(tmpset),stripwidth(tmpset)];    %strip�ĳ������
                    %�ڲ�Item��itemRotaFlag���� 
                    idxItem = find(d.Item.Item_Strip(1,:)==tmpset );
                    if isscalar(idxItem)
                        d.Item.itemRotaFlag(idxItem) = ~d.Item.Rotaed(idxItem);
                    end                    
                    %�ڲ�LU��LURotaFlag ���{ ߀δ�� %�ڲ�Item��CoordItemStrip���{                    
                end
            end
        end
        %% ����ʽ��Strip��Bin���㷨
%         printstruct(d);
        [d.Strip,d.Bin]= HStripToBin(d.Strip,d.Veh,p);
        %% Item��bin����Ϣ��ȡ:
%         printstruct(d);
        [d] = HItemToBin(d);
         %% ����binװ����
         % ItemloadingrateLimit - ÿ��bin��Item�������/ÿ��binȥ��ʣ���ߺ�������
         % Itemloadingrate - ÿ��bin��Item�������/ÿ��bin���������
         d = computeLoadingRateBin(d);
        function d = computeLoadingRateBin(d)
            % ��ʼ��
            nBin = size(d.Bin.LW,2);
            d.Bin.Binvolume = zeros(1,nBin);
            d.Bin.Itemvolume = zeros(1,nBin);
            d.Bin.Itemloadingrate = zeros(1,nBin);
            d.Bin.ItemloadingrateLimit = zeros(1,nBin);
            % ����ÿ��Bin��װ����            
            BinWidth = d.Veh.LWH(1,:);
            BinHeight = d.Veh.LWH(2,:);
            BinVolume = BinWidth .* BinHeight;
            %ÿ��Bin�Ŀ������ = �����߶�*�������
            d.Bin.Binvolume = repmat(BinVolume,1,nBin);            
            %ÿ��Bin �����޿������ = ���(binʹ�ÿ��=�������-binʣ����) *�߶�(binʹ�ø߶�=�����߶�-binʣ��߶�)
            d.Bin.BinvolumeLimit = (BinWidth - d.Bin.LW(1,:)) .* (BinHeight - d.Bin.LW(2,:));
            
            a = d.Item.LWH;
            b = d.Item.Item_Bin;
            for iBin =1:nBin
                %ÿ��Bin��װ�����
                d.Bin.Itemvolume(iBin)= sum(a(1, (b(1,:)==iBin)) .* a(2, (b(1,:)==iBin)));
            end
            %ÿ��bin��װ�ر���
            d.Bin.loadingrate =  d.Bin.Itemvolume ./ d.Bin.Binvolume;
            %ÿ��bin������װ�ر���
            d.Bin.loadingrateLimit =  d.Bin.Itemvolume ./ d.Bin.BinvolumeLimit;
        end

end

    %% ************ �ж��Ƿ���ͬ�����������ڰڷ�
    function flag = isAdjacent(d)
        flag = 1;
        printstruct(d);
        % ÿ��bin���ҳ�������ID����Strip�Ƿ�����
        nBin = size(d.Bin.LW,2);
        for iBin = 1:nBin
            t = [d.Item.ID; d.Item.Item_Strip; d.Item.Item_Bin ];
            tiBin = t( : , t(4,:) == iBin );
            nIdType = unique(tiBin(1,:)); %nIdType: ��iBin�ڰ�����LU��ID����
            for iId = 1:nIdType
                tiId = tiBin( : , tiBin(1,:) == iId );
                nIdStrip = unique(tiId(2,:)); %nIdStrip: ��iBin����iID�°�����Strip�����
                % �ж������ķ��뱾ID���͵�Strip����Ƿ�����
                if ~all(diff(sort(nIdStrip))==1)
                    flag = 0;
                end
            end                    
        end
    end
% % %             % ns - ��bin��strip������˳��
% % %             ns = d.Strip.stripBeBinMatrix(2,d.Strip.stripBeBinMatrix(1,:) == iBin);
% % %             % ni - ��bin��item��LU���ͼ�˳��
% % %             d.Item.Item_Bin(1,:) == iBin
% % %             ni = d.Item.ID(d.Item.Item_Bin(1,:) == iBin);
% % %             [a,b] = find(d.Item.ID(d.Item.Item_Bin(1,:) == iBin));
% % %             ni_uni = unique(ni);
% % %             for ini = 1:length(ni_uni)
% % % %                 d.Item.
% % % %                 d.Item.Item_Strip(:,
% % %             end
% % %             nStrip = length(ns);
% % %             % i,j is adjacent strips(levels)
% % %             for iStrip = 1:nStrip
% % %                 for jStrip = (iStrip+1):(nStrip-1)
% % %                 [is] = find(ns==iStrip); %��3��strip�ŵ�1��
% % %                 [js] = find(ns==jStrip); %��1��strip�ŵ�2��
% % %                 LUIDInis = d.Item.ID(1,(d.Item.Item_Strip(1,:)==is))
% % %                 LUIDInjs = d.Item.ID(1,(d.Item.Item_Strip(1,:)==js))
% % %                 
% % %                 end
% % %             end
% % %         end
% % %         

%% **** �㷨ָ��ѡ�����Ž� ****    
function [daMax,parMax] = getbestsol(DaS,Par)

%��ȡ����ָ��Ͷ�Ӧ����
for r=1:length(DaS)
    resLoadingRateStrip(r) = mean(DaS(r).Strip.loadingrate); %strip��װ������� Itemloadingrate ItemloadingrateLimit
    resLoadingRateStripLimit(r) = mean(DaS(r).Strip.loadingrateLimit); %strip��limitװ������� Itemloadingrate ItemloadingrateLimit
    resLoadingRateBinLimit(r) = mean(DaS(r).Bin.loadingrateLimit); %bin��limitװ������� Itemloadingrate ItemloadingrateLimit
    resLoadingRateBin(r) = mean(DaS(r).Bin.loadingrate); %bin��limitװ������� Itemloadingrate ItemloadingrateLimit
%     Par(r);
end

%% �㷨ѡ�����ŵĽ���û�
% maxresBin=max(resLoadingRateBinLimit(1,idxStrip)); %�ҳ�idxStrip�е����bin
% if ~all(ismember(idxBin,idxStrip)),   error('not all member of bin in strip'); end %�����п��ܳ��� 
%% 1 maxresBin�����泵����ƽ��װ����,��Ʒ����һ��,binԽ��,��ֵԽС,��Խ��,�����ֵ�Ǳ���
%% 2 maxresStrip������Strip��ƽ��װ����,��Ʒ����һ��,strip���һ��,�߶�Խ��,��ֵԽС,��Խ��,�����ֵ��һ��ʱ���루��Ϊ���������ÿ���
%% ����ֵ�ڲ�ͬbin�߶�ʱ������Ӱ�죬�Ҹ�ֵ��ʱ����Ϊ���������ܲ�����
%% 3 maxresStripLimit��������Strip��ƽ��װ����,��Ʒ����һ��,strip�ڲ����Խ��,��϶Խ��,ֵԽС,�����ֵ�����Ǳ���
%% ��ֵ��ʱ����Ϊ���������ܺã�strip�ڲ���϶С��������һ��ʱ���ţ����п�����ͬ���̲���һ��
%% 4 maxresBinLimit��������Bin��ƽ��װ����,��Ʒ����һ��?? �������������,���۲�
idxBin=find(resLoadingRateBin==max(resLoadingRateBin));
idxStrip=find(resLoadingRateStrip==max(resLoadingRateStrip));
idxStripLimit=find(resLoadingRateStripLimit==max(resLoadingRateStripLimit));
idxBinLimit=find(resLoadingRateBinLimit==max(resLoadingRateBinLimit));

%% 5 �ҳ�idxStrip��idxBin���ߵĽ���
% % if isempty(intersect(idxBin,idxStrip))
idx =idxBin;
if isempty(idx), error('idxBinΪ�� '); end %���󼸺������ܳ���
idx0 =intersect(idx,idxStripLimit);
if ~isempty(idx0),  
    idx = idx0; 
else
    warning('idx0 is empty');
end
% if isempty(idx), error('idxBin and idxStripLimit �Ľ���Ϊ�� '); end %���󼸺������ܳ���
idx1 = intersect(idx,idxBinLimit);
if ~isempty(idx1),  
%     idx = idx1; 
else
    warning('idx1 is empty');
end
idx2 = intersect(idx,idxStrip);
if ~isempty(idx2),  
%     idx = idx2; 
else
    warning('idx2 is empty');
end

%% ��idxʣ��ķ��ص�������
if ~isempty(idx)
    for tmpidx=1:length(idx)
        daMax(tmpidx) = DaS(idx(tmpidx));
        parMax(tmpidx) = Par(idx(tmpidx));
    end
end

end



%% ********************** ������ts�㷨�Ĵ��� ��ʱ���� ****************

% [ub,x,b] = HnextFit(Item,Veh);
% [ub,x,b] = HnextFit_origin(Item,Veh);
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

    function [ub,px,pb]  = HnextFit(Item,Veh)
        % Initialize
        d = size(Item.LWH,1)-1;
        n = size(Item.LWH,2);
        nn = n + 1;
        w = Item.LWH(1:d,:);
        W = Veh.LWH(1:d,:);
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

    function [ub,px,pb]  = HnextFit_origin(Item,Veh)
        % Initialize
        d = size(Item.LWH,1);
        if d==3
            d=d-1;
        end
        n = size(Item.LWH,2);
        nn = n + 1;
        w = Item.LWH(1:d,:);
        W = Veh.LWH(1:d,:);
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

    function [ lb ] = computerLB(Item,Veh)
        sum1 = sum(prod(Item.LWH,1));        
        % todo �����ж��Ƿ����е�BinArray�����е�bin����ͬ�� ����� �����ִ��
        sum2 = prod(Veh.LWH(:,1));
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

% % %% function [StripSolutionSort] = HnextFitDH(d)
% % function [StripSolutionSort] = HnextFitDH(d)
% % % ����: d
% % % ���: StripSolutionSort
% % %% ��ȡ������bin,��άitem����
% % % nDim nItem nBin
% % % itemDataMatrix uniBinDataMatrix
% % nDim = size(d.Item.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % itemDataMatrix = d.Item.LWH(1:nDim,:);
% % tmpbinDataMatrix = d.Veh.LWH(1:nDim,:);
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
% % %% function [StripSolutionSort] = HfirstFitDH(d)
% % function [StripSolutionSort] = HfirstFitDH(d)
% % % ����: d
% % % ���: StripSolutionSort
% % %% ��ȡ������bin,��άitem����
% % % nDim nItem nBin
% % % itemDataMatrix uniBinDataMatrix
% % nDim = size(d.Item.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % itemDataMatrix = d.Item.LWH(1:nDim,:);
% % tmpbinDataMatrix = d.Veh.LWH(1:nDim,:);
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
% % %% function [StripSolutionSort] = HbestFitDH(d)
% % function [StripSolutionSort] = HbestFitDH(d)
% % % ����: d
% % % ���: StripSolutionSort
% % %% ��ȡ������bin,��άitem����
% % % nDim nItem nBin
% % % itemDataMatrix uniBinDataMatrix
% % nDim = size(d.Item.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % itemDataMatrix = d.Item.LWH(1:nDim,:);
% % tmpbinDataMatrix = d.Veh.LWH(1:nDim,:);
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
% % function [stripLeftMatrix,pbelongMatrix,pitemMatrix,pcoordMatrix,ord]  = HfirstFitDH2(Item,Veh)
% % % ���������ʼ��
% % nDim = size(Item.LWH,1);
% % if nDim ==3, nDim = nDim-1;end
% % nItem = size(Item.LWH,2);
% % nBin = nItem;
% % % nn = n + 1;
% % itemMatrix = Item.LWH(1:nDim,:);
% % binMatrix = Veh.LWH(1:nDim,:);
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
% % function [stripLeftMatrix,pbelongMatrix,pitemMatrix,pcoordMatrix,ord]  = HbestFitDH2(Item,Veh)
% % % ���������ʼ��
% % nDim = size(Item.LWH,1);
% % if nDim ==3, nDim = nDim-1;end
% % nItem = size(Item.LWH,2);
% % nBin = nItem;
% % % nn = n + 1;
% % itemMatrix = Item.LWH(1:nDim,:);
% % binMatrix = Veh.LWH(1:nDim,:);
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
% % function [pbelongItemBinMatrix,pbelongStripBinMatrix,pcoordItemBinMatrix,binLeftMatrix ] = HbestFitBinDH(stripLeftMatrix,pbelongMatrix,pitemMatrix,pcoordMatrix,Item,Veh)
% % % ���������ʼ��
% % nDim = size(Item.LWH,1);
% % if nDim == 3, nDim = nDim-1;end
% % nItem = size(Item.LWH,2);
% % nBin = nItem;
% % 
% % nStrip = sum(stripLeftMatrix(3,:)>0); %����ʹ�õ�Strip������
% % % nn = n + 1;
% % binMatrix = Veh.LWH(1:nDim,:);
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
%% Call ������:
% clear;close all; format long g; format bank; %NOTE ����MATLAB CODE ֧��
% rng('default');rng(1); % NOTE �Ƿ�����ı�־
% LU = struct('ID',[],'LWH',[],...
%     'weight',[],'Lbuffer',[],'Wbuffer',[],'Type',[],'Material',[]);
% Item = struct('ID',[],'LWH',[],...
%     'weight',[],'Lbuffer',[],'Wbuffer',[]);
% Strip = struct('ID',[],'LWH',[],... %ONLY LW
%     'weight',[],'Lbuffer',[],'Wbuffer',[]);
% Veh = struct('ID',[],'LWH',[],...
%    'Capacity',[],'Lbuffer',[],'Wbuffer',[],'Hbuffer',[]);
% d = struct('LU',LU,'Item',Item,'Strip',Strip,'Veh',Veh);

% 1������ʼ��
% whichStripH 1 best 2 first 3 next; whichBinH 1 best; TODO ��������������ʽ
% whichSortItemOrder 1 ���ߵݼ� 2 ��̱ߵݼ�; 
% whichRotation 1:����rotation 0:��ֹ
% rotation��� 1 1 2 1 0 (1 1 2 1 1 )(1 1 2 1 2) % ��rotation��� 1 1 1 0 0 ��2/3 1 1 0 0��
% whichRotationHori 0:�ڰ���˳��ʱ��FBS_{RG}��ʽ; 1��New/NoNew��Horizon��ʽ 2��New/NoNew��Vertical��ʽ
% ParaArray = struct('whichStripH',1,'whichBinH',1,'whichSortItemOrder',2,...
%     'whichRotation',1,'whichRotationHori',0,'timeLimit',100,'ub0',10);
% % ParaArray = struct('whichStripH',1,'whichBinH',1,'whichSortItemOrder',2,...
% %     'whichRotation',1,'whichRotationHori',0,'whichRotationAll',1,'whichRotationBin',1,'timeLimit',100,'ub0',10);
