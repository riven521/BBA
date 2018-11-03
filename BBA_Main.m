%% BBA_MAIN demo
%% Form
%    [output_CoordLUBin,output_LU_LWH,output_LU_Seq] = ... 
%    BBA_Main(LUID,LULWH,VEHID,VEHLWH,varargin)
%        
%% Inputs (varargin) LUID: ������ ; LULID: �����
%   LUID	                (1,n)   �������� ��ͬ���ֱ���ͬһ����,����Ѷ� 
%   LULWH                (3,n)   ���̿���
%   VEHID                 (1,m)  ���ͱ��
%   VEHLWH              (3,m)   ���Ϳ��ߣ����Ƕ೵�ͣ�
%   ------------------------------------------------------
%   LUSID                  (1,n)   ���̹�Ӧ�̱��
%   LUPID                  (1,n)   �����㲿�����
%   LUISROTA            (1,n)  �����Ƿ�������ת
%   LUMARGIN         (1,n)   ���̼�margin(1-4��������)  �������̳����=ÿ�����̵�ʵ�ʳ����+���ӵ�margin
%   LUWEIGHT           (1,n)  ��������
%   BINWEIGHT         (1,m)  ��������������
%   LULID                  (1,n)   �������ͱ��
%% Outputs
%   output_CoordLUBin      (3,n)    ÿ��LU��X,Y,Z
%   output_LU_LWH            (3,n)    ÿ��LU�Ŀ��ߣ���ת��ģ�ʵ��ֵ��
%   output_LU_Seq             (8,n)    ��1: LU��ĳ��BIN�ڣ���2: LU�ڸ�BIN�ڵİ���˳�� ������
%

%%
function [output_CoordLUBin,output_LU_LWH,output_LU_Seq] = ...
    BBA_Main(LUID,LULWH,VEHID,VEHLWH,varargin) %ǰ4������

%% Initialize Global Variable
% clear;close all; format long g; format bank; %NOTE ����MATLAB CODE ֧��
% rng('default');rng(1); % NOTE �Ƿ�����ı�־
%close all;
global ISdiagItem ISshuaiwei ISpingpu ISlastVehType ISreStripToBin ISisNonMixed ISisMixTile ISsItemAdjust ISpingpuAll ISreStripToBinMixed
global ISplotBBA ISplotSolu ISplotEachPingPu ISplotStrip ISplotPause ISplotShowType % plotStrip
global ISisNonMixedLU ISisMixTileLU

% ISsItemAdjust = 1  % ��ʱ���� ��;������
% ISreStripToBinMixed = 1 %��ͷ���ȷ�AllPure����, �ٿ�������LU����������� Ĭ��Ϊ1 Ӧ�ÿ���ɾ���Ĳ���

ISplotBBA = 0
ISplotShowType = 1 % 1 LID 2 isShuaiWei
        % ISplotSolu = 0

        ISplotStrip = 0 % ÿ��Run algorithm ����Strip����ʾ��� ��ϸ��
        ISplotEachPingPu = 0 % ÿ��Main ƽ��ʱ ����Strip����ʾ��� ��ϸ��

ISplotPause = -0.05 % plot���ʱ��

ISdiagItem = 0  % Ĭ��Ϊ 0 �� Ϊ1 ����Щ���ڵ͵ı���ΪItem�߶�����, checkԭ���

% ���滹������, ����Ҫ�� Ŀǰȫ��Ϊ1
ISisNonMixedLU = 1 % 555: ���ȷǻ��LU�γ�ITEM, ͼ�ÿ���� ������ Ĭ��Ϊ 1
ISisMixTileLU = 1      % 555: ���Ȼ��LU�ĵ���ITEM�������γ�ITEM, ͼ�ÿ���� ������ Ĭ��Ϊ 1

ISisNonMixed = 1    % 555: ���ȷǻ��Item�γ�STRIP, ͼ�ÿ���� ������ Ĭ��Ϊ 1
ISisMixTile  = 1    % 555: ���Ȼ��Item�ĵ���Strip�������γ�STRIP, ͼ�ÿ���� ������ Ĭ��Ϊ 1 �����ܳ��ֻ������

ISreStripToBin = 1  % ��ͷ����LU����������� Ĭ��Ϊ1 ����

ISshuaiwei = 1          % 555 : ��Ⱥ͸߶Ȳ���, ˦β   ******  �ò�����Ҫ�������pingpu���ʹ�� ��˦β ƽ���޷�����*******
ISpingpu = 1            % 555 : ��Ⱥ͸߶Ȳ���, �Ҳ���>1, ƽ��. ���������� (����ƽ�̺���ISisNonMixedì��)
ISpingpuAll = 1       %555: ���о�ƽ��, ֻҪ�ó����ŵ���; ���Ų���, ��������˦βƽ������

ISlastVehType = 1   % 555: ���һ���ĵ���, �������޹�, �ݲ�����

%% BUG
%% Initialize Data Structure
if nargin ~= 0
    d = DataInitialize( ...
            'LUID', LUID,...
            'LULWH',LULWH, ...
            'VEHID',VEHID,...
            'VEHLWH',VEHLWH,...
            'LUSID',varargin{1},...
            'LUPID',varargin{2},...
            'LUISROTA',varargin{3},...
            'LUMARGIN',varargin{4},...
            'LUWEIGHT',varargin{5},...
            'VEHWEIGHT',varargin{6},...
            'LULID',varargin{7},...
            'LUINDEX',varargin{8});
else
    n=28; m=2;  % 16��Ҫע�� 250 srng1
    d = DataInitialize(n,m);  %0 Ĭ��ֵ; >0 �����������n������ ����ֱ������BBAʱ����
    
    filename = strcat('GoodIns',num2str(n));
    printstruct(d.Veh);  %��������һ������,�Ѷ��䰴����Ӵ�С����; 
    
     save( strcat( '.\new\', filename), 'd');
%     load .\new\GoodIns200.mat;
end
% printstruct(d);
% TVEHIN = struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false));
% TLUIN.Properties.VariableNames{'PID'} = 'OPID'; TLUIN.Properties.VariableNames{'SID'} = 'OSID';
% s = table2struct(TLUIN,'ToScalar',true)
% t = struct2table(l,'AsArray',true)
%%
% t.ID
% t = [d.LU.ID;d.LU.LWH]
% sortrows(t',[1,4],{'ascend','descend'})

%% û�����Ե���ʱ����
    n = numel(d.LU.Weight);
    % ������ʼ����˳�� ��ʱû����
    if ~isfield(d.LU, 'Index') %��ͬ��ǿ����index����
        d.LU.Index = 1:n;
    end
    
    % ƽ��ʹ������
    if ~isfield(d.LU, 'maxL')
        d.LU.maxL(1,:) =  floor(d.Veh.LWH(1,1)./d.LU.LWH(1,:));
        d.LU.maxL(2,:) =  floor(d.Veh.LWH(2,1)./d.LU.LWH(2,:));
        d.LU.maxL(3,:) =  floor(d.Veh.LWH(3,1)./d.LU.LWH(3,:));   %����ÿ������LU�ĸ߶ȵ�������
    end
    
    if ~isfield(d.LU, 'maxHLayer'),     d.LU.maxHLayer = 10*ones(1,n); end% maximum given height layer
    
%% Initialize Parameter
nAlg = 1;
pA(nAlg) = ParameterInitialize('whichStripH', 3,...
                             'whichBinH',3, ...
                             'whichSortItemOrder',3, ... 
                             'whichRotation',2, ...
                             'whichRotationHori', 1);
                  
% % nAlg = 1;
% % for i = 3:3 %1-3 best first next���� ��Ϊ3: ������ǰ��С��϶���������� ��Ϊһ������, ������Υ������Լ��
% %     for j=3:3 %0-3 ����: 0: Vert��1: Hori; 2:error  3:����϶��С����   Gpreproc �˴����HItemToStrip�����е���Ʒ�ڷ�
% %         for k=2:2 %0-2 Ĭ��0 ������ת 1ȫ������ת 2: ����Ϊ�����Ƿ�����Rotation 
% %             for TLUl=1:1 % ������ :  % 0-2 0��ȡ�� ����1-2 RotaHori 1hori 2 vert 555 ��Ų��˻��ݷţ��������ݷź󲻻��ţ��Ų��£���
% %                 for m=3:3 %1-3 best first next���� ѡ�õ�best fit �Ƿ��λNEXT FIT 1002�ո�Ϊm=3
% %                 % pA nAlg 
% %                 pA(nAlg) = ParameterInitialize( ...
% %                              'whichStripH', i,...
% %                              'whichBinH',m, ...
% %                              'whichSortItemOrder',j, ... 
% %                              'whichRotation',k, ...
% %                              'whichRotationHori', TLUl);
% %                  nAlg=nAlg+1;
% %                 end
% %             end
% %         end
% %     end
% % end
% % nAlg = nAlg - 1;

%% Simulate - All ALGORITHM

fprintf(1,'\nRunning the simulation...\n');

% Run ALL algorithm configure
for iAlg = 1:nAlg
    
    %% 1 �������㷨
    % 1.1 ��ȡd: �����������㷨    
    maind = d; % ��Ҫ���������ݱ���
    
    do = RunAlgorithm(d,pA(iAlg));   %��ȡ���н�ṹ��
    do.LU.LU_VehType = ones(size(do.LU.ID)) * do.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ

    % 1.2 �޶�d�ڵ�LU��Veh��LWH���� % ����֮ǰ���㲻��margin��LU��Item��LWH+Coord.
    [do.LU,do.Item] = updateItemMargin(do.LU,do.Item);
    dA(iAlg)=do;
    
    % 1.3 CHECK ����������LU
    maintLU = struct2table(structfun(@(x) x',maind.LU,'UniformOutput',false));
    tLU = struct2table(structfun(@(x) x',do.LU,'UniformOutput',false));    
    checkLU(maintLU,tLU);
    checktLU(do.LU);
    
%     plotSolutionT(do.LU,do.Veh);
    1
   % plotSolution(do,pA(iAlg)); %��������
    
    %% 2 ���г��͵����㷨,���ı�d ��ȡd1��do1, flaggetSmallVeh : 
    if ISlastVehType
    % 2.1 ��ʼ��
    allidxVehType = length(unique(do.Veh.ID)); %��������������(δ�ų���ͬ����)
    flaggetSmallVeh = 0;
    
                                                                        % ׼���滻d1 %  �Ե�Lu.LWH�ĳ��� -< ֮ǰ�ǿ� (�ѷ���getdinLastVeh��)
                                                                        %     d1.LU.LWH([1,2],:) = flipud(d1.LU.LWH([1,2],:)); 
                                                                        %d1 = getdinLastVeh(do);
    % 2.2 CHANGE ��maind������������
    luIdx = do.LU.LU_Bin(1,:) == max(do.LU.LU_Bin(1,:));
    d1 = getdinThisVeh(maind,luIdx);                    
                    
    % 2.3 ���������Ͷ���2��,�Ž����滻
    while(allidxVehType>=2)
        % 2.1 ��ȡ����Ͳ������㷨 % �����һ����������ǰѭ��; until�ڶ�����; �˴�����
        d1.Veh = structfun(@(x) x(:,allidxVehType), do.Veh,'UniformOutput',false); %�����һ�ֳ��Ϳ�ʼ����
                                            %disp(d1.Veh.LWH)
                                            %TLUIN_LAST = struct2table(structfun(@(x) x',d1.LU,'UniformOutput',false));
        
                                                                                                                            %d1 = RunAlgorithm(d1,pA(iAlg));
        do1 = RunAlgorithm(d1,pA(iAlg));   %������������һ��Bin������lastd�������� 555555555555555555555
                                            %     plotSolution(do1,pA(iAlg));

                            %             do1.LU.LU_VehType = ones(size(do1.LU.ID)) * do1.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ
        do1.LU.LU_VehType = ones(size(do1.LU.ID))*do.Veh.order(allidxVehType); % �������LU_VehType
        [do1.LU,do1.Item] = updateItemMargin(do1.LU,do1.Item);

        % 2.4 �жϸó����Ƿ����
        % ����Veh�ڲ�������ݼ�����,��ȡorder�ĵڸ���Ϊ��ǰ��Ӧ�泵��������
        % �ж�: �Ƿ��Ϊ��allidxVehType(С)���ͺ�,1���������Է���;
        if max(do1.LU.LU_Bin(1,:)) == 1
                                                                                                                                    %do1.LU.LU_VehType = ones(size(do1.LU.ID))*do.Veh.order(allidxVehType); % �������LU_VehType
            flaggetSmallVeh=1;
            break;
        end
        
        % 2.5 ���Ų���,ѡ������� -> allidxVehType�ݼ� do1.Veh�����ֵ
        allidxVehType= allidxVehType-1;
        d1.Veh = [];
    end
    
    %% CHECK 2 ���г��͵����㷨,
    if flaggetSmallVeh==1
        % d1;
        % do1
        t1 = struct2table(structfun(@(x) x',d1.LU,'UniformOutput',false));
        to1 = struct2table(structfun(@(x) x',do1.LU,'UniformOutput',false));
        checkLU(t1,to1);
        checktLU(do1.LU);
    end
    end


    %% 3 ����ƽ���㷨,���ı�d ��ȡd2Array d3Array  do2Array do3Array flagTiled
    if ISpingpu==1
    % 3.1 ��ʼ��5������    
    flagTiledArray = zeros(1,length(do.Bin.Weight));  %1��������ƽ�� 2����˦βƽ��
    
    d2Array(1:length(do.Bin.Weight)) = maind;
    do2Array(1:length(do.Bin.Weight)) = structfun(@(x) [], do, 'UniformOutput', false);
    d3Array(1:length(do.Bin.Weight)) = maind;
    do3Array(1:length(do.Bin.Weight)) = structfun(@(x) [], do, 'UniformOutput', false); %do;
            
    %% 3.2 ����binѭ��ƽ��
    bidx = 1:length(do.Bin.Weight); % NOTE: �޸�ֻ����˦βƽ�̵�ȫ��Bin���뿼��, �Է�˦βƽ�̵Ľ���ƽ���ж� % bidx = find(do.Bin.isTileNeed);
    % ѭ��: ÿ��bin�ֱ���ƽ��(����ƽ�̺�˦βƽ�̶�ѡһ����������ƽ��)
    for i=1:numel(bidx)
        ibin = bidx(i);

        % $1 GET d2 ��ibin�ڴ��㷨���������
        luIdx = do.LU.LU_Bin(1,:) == ibin;
        d2 = getdinThisVeh(maind,luIdx); %�޸ĳɴ�maind��ȡIuIdx������,���Ǵ�������d����ȡ

        %% COMMENT
%         d2.Veh = do.Veh;
%         d2.Veh = rmfield(d2.Veh,{'Volume','order'});
        % 2 �������/һ��strip�ڵ�LU
%         luidx = do.LU.LU_Bin(1,:) == ibin;  %do.LU.LU_Strip(1,:) == istrip

%         d2.LU = structfun(@(x) x(:,luidx),do.LU,'UniformOutput',false);
%         d2.LU.LWH([1,2], d2.LU.Rotaed ) = flipud(d2.LU.LWH([1,2], d2.LU.Rotaed)); %LU.LWH ����ת,��ָ�ԭ��
%         d2.LU.PID = d2.LU.OPID;     d2.LU.SID = d2.LU.OSID;  %  d2.LU.LID = d2.LU.OLID;
%         d2.LU = rmfield(d2.LU,{'Rotaed','order','LU_Item','DOC','LU_Strip',...
%             'LU_Bin','CoordLUBin','CoordLUStrip','LU_VehType','OPID','OSID'});
%         
%         d2.Par = do.Par;
        %% 3.3 �������ȫ��ƽ��(�����Ƿ�˦βƽ��), �۲챾ibin���Ƿ����ȫ��ƽ��,�����,��ȡ��˦βƽ��; ����,����˦βƽ��
        if ISpingpuAll==1
            d3 = d2;
            d3.LU.maxHLayer(:) = 1; %d2��ȫ��LU�Ĳ����趨Ϊ1 55555 ȫ��ƽ�̵���Ҫ����
            d3Array(ibin) = d3;

            % $3.3.1 reRunAlgorithm do3��d3�����Ľ��
            do3 = RunAlgorithm(d3,pA(iAlg));             
            % do3Array(ibin) = do3;
            
            % $3.3.2 ��ȫ��ƽ��û������,������
            if max(do3.LU.LU_Bin(1,:)) == 1
                flagTiledArray(ibin)=1; %1��������ƽ��
                do3.LU.LU_VehType = ones(size(d3.LU.ID)) * do3.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ
                [do3.LU,do3.Item] = updateItemMargin(do3.LU,do3.Item);               %  plot3DBPP(do3,pA(iAlg))
                do3Array(ibin) = do3;
                
%                 plotSolutionT(do3.LU,do3.Veh);
% 1
                % do3�޸ĵ�d�У�����Ŀǰ������do2Array�У�δ��d�ϲ�
                continue;   %continue������������˦βƽ����
            end
        end
        %% 3.4 ������ƽ��ʧ�� ����˦βƽ��

                % 3.4.1 GET do2 ��ibin�ں�������ݣ�˦βƽ��ʹ��
                        stripidx = do.Strip.Strip_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip
                        itemidx = do.Item.Item_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip
                do2.Veh = do.Veh;
                do2.LU = structfun(@(x) x(:,luIdx),do.LU,'UniformOutput',false);
                do2.Bin = structfun(@(x) x(:,ibin),do.Bin,'UniformOutput',false);                       
                do2.Strip = structfun(@(x) x(:,stripidx),do.Strip,'UniformOutput',false);
                do2.Item = structfun(@(x) x(:,itemidx),do.Item,'UniformOutput',false);        
%                plotSolutionT(do2.LU,do2.Veh);
            while do2.Bin.isTileNeed(1) == 1 %do2�ڵ�Bin��Զֻ��1��, ����ƽ�̺��bin����Ҫƽ��,������while�ж�

            % $3.4.2 �޶�d2.LU.maxHLayer (����ibin�����ѡ���ļ���stripƽ��) TODO $4д����Щ����,���ڼ�
            % $4.1 GET luidxPP
            % ѭ���ӱ�ibin�����һ��strip��ʼƽ�� istrip= nbStrip;
            nbStrip = numel(do2.Strip.Weight);
                        if unique(do2.Strip.Strip_Bin(2, :)) ~= nbStrip,    error('��Ԥ�ڴ���');    end
            istrip= nbStrip;
            fi = find(do2.Strip.Strip_Bin(2,:) >= istrip ); % ͬbin��strip��� and ˳��>=istrip
            u=unique(do2.LU.LU_Strip(1,:)); %��ȡStrip��ŵ�Ψһ����ֵ
            luidxPP = ismember(do2.LU.LU_Strip(1,:), u(fi)); %%% fi->u(fi) ��������� ********************* 
                        if ~any(luidxPP),  error('luidxPPȫ��Ϊ��, ������u(fi)��Ӧ��Lu�߼��ж�'); end
        
            % $4.2 �޶�d2.LU.maxHLayer
%              d2.LU
%              maind.LU
            d2.LU.maxHLayer(luidxPP) = min( d2.LU.maxL(3,luidxPP), d2.LU.maxHLayer(luidxPP)) - 1;

            % $4.2 ����ǰluidxPP��ӦLu�Ĳ������Ѿ�Ϊ1��, ����Ҫ���Ӹ����istrip��luidxPP; ���޶�d2.LU.maxHLayer
            % GET ���� d2.LU.maxHLayer(luidxPP)
            while all(d2.LU.maxHLayer(luidxPP)<1)  
               istrip = istrip-1;
               if istrip==0,break;end
                fi = find( do2.Strip.Strip_Bin(2,:) >= istrip ); %fi = find( do2.Strip.Strip_Bin(2,:) == istrip ); 
                luidxPP = ismember(do2.LU.LU_Strip(1,:), u(fi)); %%% fi->u(fi) ��������� ********************* 
                                        if ~any(luidxPP),  error('luidxPPȫ��Ϊ��, ������u(fi)��Ӧ��Lu�߼��ж�'); end
                                        if istrip == 0,  error('��bin������tileneed,��Ԥ�ڴ���');   end
                d2.LU.maxHLayer(luidxPP) = min( d2.LU.maxL(3,luidxPP), d2.LU.maxHLayer(luidxPP)) - 1;
            end
            % �޸�: ������Ļָ�Ϊ1
            d2.LU.maxHLayer(d2.LU.maxHLayer<=1) = 1;

            % $5 reRunAlgorithm
                                                             %    plotSolution(do2,pA(iAlg));

            d2Array(ibin) = d2;
            do2 = RunAlgorithm(d2,pA(iAlg)); 
            
            do2.LU.LU_VehType = ones(size(d2.LU.ID)) * do2.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ
            [do2.LU,do2.Item] = updateItemMargin(do2.LU,do2.Item);

            % do2Array(ibin) = do2; ����ע�ͣ���Ϊ�Ǹ�ѭ��
                                                                if ISplotEachPingPu == 1,     plotSolution(do2,pA(iAlg));       end
            % $6 ����
            if max(do2.LU.LU_Bin(1,:)) == 1
                                                                                                                                %    do2.LU.LU_VehType = ones(size(do2.LU.ID))*do.Veh.order(1); 
                flagTiledArray(ibin)=2; %2����˦βƽ��
                do2.LU.LU_Bin
                do2Array(ibin) = do2;
%                 plotSolutionT(do2.LU,do2.Veh);
                1
                % do2 ���ݲ�����d ����return2bba���޸�
                % do2 ���ݽ���d???? return2bba���޸ģ�����                
            else
                break;  %����˦β�Ų��� ���ټ���˦βƽ���ˣ����˽�����������һ����˦β�ж�
            end
            
            end % END OF WHILE
    end% END OF FOR
    
    %% CHECK
    if any(flagTiledArray)
        for ibin=1:length(flagTiledArray)
            if flagTiledArray(ibin)==1 %����ƽ��
                t3 = struct2table(structfun(@(x) x',d3Array(ibin).LU,'UniformOutput',false));
                to3 = struct2table(structfun(@(x) x',do3Array(ibin).LU,'UniformOutput',false));
                checkLU(t3,to3);
                checktLU(t3);
                checktLU(to3);
                
            end
            if flagTiledArray(ibin)==2 %˦βƽ��
                t2 = struct2table(structfun(@(x) x',d2Array(ibin).LU,'UniformOutput',false));
                to2 = struct2table(structfun(@(x) x',do2Array(ibin).LU,'UniformOutput',false));
                checkLU(t2,to2);
                % checktLU(t2);
                checktLU(to2);
            end
            
            flagTiledArray;
            d2Array;
            do2Array;
            d3Array;
            do3Array;
        end
    end
    end
end

   
%% Simulate - CHOOSE BEST ONE
% 555 �㷨�����жϲ��ų�bin����ͬ�������̲����ڵĽ� TODO ���ݵ�CHECK
%  flagA(iAlg) =  isAdjacent(dA(iAlg));           % �㷨�ж��Ƿ���ͬ�����������ڰڷ� +
% % dA = dA(1,logical(flagA));
% % pA = pA(1,logical(flagA));

% TODO �Ӷ���㷨�����ѡ���ӱض�bin�����ڵ����Ž�� - NOTE: ���õ�����ʱ���迼��
% if isempty(dA), error('�����������нⶼ�������̲����ڵ���� \n'); end
% [daBest,paBest] = getbestsol(dA,pA);  %�ɲ�����d1 -> ������������һ��bin��ϵ����
%%%% Return length(parMax) �� solutions to BBA
% if isempty(daBest), error('��������δ�ҳ����Žⷵ��BBA \n'); end
% bestOne = 1;
                                %%%% dA = do = daBest(1) = daBest(bestOne) %isequal(do,d1)
%% POST PROCESSING


    % daBest(bestOne).LU.OPID
%% 1 ******************��ȡչʾ˳�� T=d.LU����ShowSEQ
T = getTableLU(do);
        checktLU(T) 
        1
%%
% [T_Coord,T_LWH,T_Seq] = getReturnBBA(daBest(bestOne)); %���ж��,���ص�һ�����Ž�
% T_Seq.tblorder
% T_Seq1 = T_Seq;

% T.SID = oD.SID;
% T.PID = oD.PID;
% T.LID = oD.LID;
% Tseq = T(:,{'LU_VehType','BINID','BINSEQ','SID','LID','ITEMID','PID','ShowSEQ','Weight','tblorder'});

% T_Seq1 = finalCheck([T_Coord,T_LWH,T_Seq],TLUIN); %����1�����㣻 ����2��ԭʼ.

%  [output_CoordLUBin,output_LU_LWH,output_LU_Seq]= getReturnBBA1(daBest(bestOne)); %% ���з��ش���

%%
% ��1���������ڳ��ͺ�(���뻻)    ��2���������ڳ����(���,���ܻ�,���ʹ�) ��3�����̳��ڰ���˳��(���뻻) ��4������SID��Ӧ�̱��(�����,���ñ�??)
% ��5������ID�ͺ�LID(�����,���ñ�?) ��6�����̶Ѷ����ITEM(���,���ܻ�,���ʹ�,��;?) ��7�������㲿�����PID(�����,���ñ�?) ������8: չʾ˳��(���뻻)

%% 2 ****************** ��Գ��ͱ仯 do1���� ��ȡ�޶��� output ******************
if ISlastVehType==1 && flaggetSmallVeh == 1 %���е������ҳ����滻�ɹ�
    T1 = getTableLU(do1);
        checktLU(T1) 
        1
    % �滻T�е����һ���Ĳ������� ����T1
    lastVehIdx = max(T{:,'BINID'});
    flaglastLUIdx = T{:,'BINID'}==lastVehIdx;
    
       %% ��Щ��仯??
       % ĳ��bin�ڵ���,��binIDһ������仯; 
       % ��LID/Weight/LWHӦ�ò���仯; SID/PID��仯; ��ΪOPID OSID OID��ԭ��
       % ĳ��bin�ڵ���,��BINSEQ,CoordLUBin,LU_VehTypeһ�������仯 ����bid��binseq����ģ� ��ITEMID�ƺ�û�� �������˰�
       % �ص��Ǹ��������LU_VehType��BINSEQ��LU_VehType �⼸���ض��仯(PID/SID��Ҫ����) % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
    T{flaglastLUIdx,{'CoordLUBin','BINSEQ','LU_VehType'}} = T1{:,{'CoordLUBin','BINSEQ','LU_VehType'}};
    checktLU(T) 
    %%
    
        % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
                % %     T{flaglastLUIdx,{'LU_VehType','BINSEQ','ShowSEQ'}} = ...
                % %         T1{:,{'LU_VehType','BINSEQ','ShowSEQ'}};
                % %     T{flaglastLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}} = ...
                % %         T1{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}}; 
                
%      [TLAST_Coord,TLAST_LWH,TLAST_Seq] = getReturnBBA(do1); %% ���з��ش���
%     TLAST_Seq1 = finalCheck( [TLAST_Coord,TLAST_LWH,TLAST_Seq],TLUIN_LAST); %����1�����㣻 ����2��ԭʼ.
    
%     %����order�ı���,�˴��������һ��bin�����������޸�
%     lastVehIdx = max(T_Seq{:,'BINID'});
%     flaglastLUIdx = T_Seq1{:,'BINID'}==lastVehIdx;
%     
%     % �ص��Ǹ�������ͳ����
%     T_Coord{flaglastLUIdx,:} = TLAST_Coord{:,:};
%     T_LWH{flaglastLUIdx,:} = TLAST_LWH{:,:};
%     % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
%     T_Seq1{flaglastLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}} = ...
%         TLAST_Seq1{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}};
    
%     T_Seq = finalCheck([T_Coord,T_LWH,T_Seq],TLUIN); %����1�����㣻 ����2��ԭʼ.

%     output_CoordLUBin(:,flaglastLUIdx) = TLAST_Coord;
%     output_LU_LWH(:,flaglastLUIdx) = TLAST_LWH;

    % ������ز�����Ҫ�滻: ��1���������ڳ��ͺ� ��3�����̳��ڰ���˳�� ��8: չʾ˳��
        %     output_LU_Seq([1,3,4,5,7],flaglastLUIdx) = output_LU_Seq2([1,3,4,5,7],:); %[1,3,4,5,7]��ʾ���޸�����ļ���
%     output_LU_Seq([1,3,8],flaglastLUIdx) = TLAST_Seq([1,3,8],:); %[1,3,4,5,7,8]��ʾ���޸�����ļ���
%     output_LU_Seq([1,3,4,5,7,8],flaglastLUIdx) = TLAST_Seq([1,3,4,5,7,8],:); %[1,3,4,5,7,8]��ʾ���޸�����ļ���
    
        %     i = [4,5,7]; % ��ЩӦ���ǲ�����
        %     if sum(output_LU_Seq(i,flaglastLUIdx) ~= output_LU_Seq2(i,:) ) >0, error('�����ı���, ����'); end
end

%% 3 ****************** ���ƽ��ѡ�� do2/do3Array���� ��ȡ�޶��� output ******************
if ISpingpu==1
    if ~all(flagTiledArray==0)
        flagTiledArray
        warning('��Ҫƽ��');;end
    for ibin=1:length(do2Array) %do*Array ��������BIN
        if flagTiledArray(ibin)==0  % ��ibinδƽ�� ����ѭ��
            continue;
        end
        if ISlastVehType==1 && flaggetSmallVeh == 1 && ibin == lastVehIdx
            error('���һ���ȸ�������,����Ҫƽ��;');
        end
        
        if flagTiledArray(ibin)==1    % ��ibin����ƽ�̳ɹ�
            T23 = getTableLU(do3Array(ibin));
                    checktLU(T23) ;
        end
        if flagTiledArray(ibin)==2    % ��ibin˦βƽ�̳ɹ�
            T23 = getTableLU(do2Array(ibin));
                    checktLU(T23) ;
            %do2.LU.LU_Bin
        end

       flagTileLUIdx = T{:,'BINID'}==ibin;
       
       %% CHECK 1 �Ƿ��ibin�ڵ�LU�������ϸ����; �Ƿ�T��ѡ�е�flagTileLUIdx����LU���ڸ�ibin��Ҳ���ϸ������
       if ~issorted(sort(T23.BINSEQ),'strictascend') || ~issorted(sort(T.BINSEQ(flagTileLUIdx,:)'),'strictascend')
               T.BINSEQ(flagTileLUIdx,:)'
               T23.BINSEQ'
               sort(T23.BINSEQ)'
               sort(T23.BINID)'
               do2Array(ibin).LU.LU_Bin
               sort(T.BINSEQ(flagTileLUIdx,:))'
               error('1');
       end
       % 2 �Ƿ��ܱ�T��BINID����ͳһbinid���Ƿ��ӱ�T23���Ƿ�ͬһbinid���Ƿ�ibin����T��LU��Ӧ��binid��
       if ~isscalar(unique(sort(T.BINID(flagTileLUIdx,:))))  || ~isscalar(unique(sort(T23.BINID)))  || ibin~= unique(sort(T.BINID(flagTileLUIdx,:)))
            sort(T23.BINID)'
            unique(sort(T23.BINID))
            unique(sort(T.BINID(flagTileLUIdx,:)))
            error('11');
       end
       % 3 �Ƿ�T23����������bin��LU����
       if sum(flagTileLUIdx) ~= height(T23)
           error('111');
       end
       % 4 �Ƿ��ܱ�T�ڳ�ʼOPID�뷵���ֱ���T23��OPID��ͬ;
       a = T.OPID(flagTileLUIdx,:)';
       b = T23.OPID';
       if  ~isequal(a,b)
           error('1');
       end
               
       %% ��Щ��仯??
       % ĳ��bin�ڵ���,��binIDһ������仯; ��bin��LU_VehTypeһ������仯��
       % ��LID/Weight/LWHӦ�ò���仯; SID/PID��仯; ��ΪOPID OSID OID��ԭ�� idExchange����
       % ĳ��bin�ڵ���,��BINSEQ,CoordLUBinһ�������仯 ����bid��binseq����ģ� 
       % ��LU_Itemһ����仯 ���ƺ�û���� �������˰�
       % �ص��Ǹ������� % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
%        sortrows(T.LU_Item)'
%        sortrows(T23.LU_Item)'
%        T.LU_Item
%        T{flagTileLUIdx,{'CoordLUBin','BINSEQ'}} = T23{:,{'CoordLUBin','BINSEQ'}};    
        T{flagTileLUIdx,{'CoordLUBin','BINSEQ','LU_Item'}} = ... %��������LU_Item�����л�,��Ȼ��;����,�����ᱨchecktLU����.
            T23{:,{'CoordLUBin','BINSEQ','LU_Item'}};    
%        T.LU_Item
%        sortrows(T.LU_Item)'
%        checktLU(T) %�����޷�ͨ���Ŀ�����; ��LU_ItemӰ�첻��,������ע�� TODO

            %% �����Ǵ��
            %        T{flagTileLUIdx,{'LU_VehType','BINSEQ','ShowSEQ'}} = ...
%            T23{:,{'LU_VehType','BINSEQ','ShowSEQ'}};
%        T{flagTileLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}} = ...
%            T23{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}};
%        T{flagTileLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ','CoordLUBin','LWH'}} = ...
%            T23{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ','CoordLUBin','LWH'}};


            %     T{flagTileLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ','CoordLUBin','LWH','BINID','ITEMID','Weight'}} = ...
            %         T2{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ','CoordLUBin','LWH','BINID','ITEMID','Weight'}};
            
            % output_LU_Seq=T{:,{'','BINID',',''ITEMID','Weight'}}'
            
            % % %     [TPP_Coord,TPP_LWH,TPP_Seq]= getReturnBBA(do2Array(ibin)); %% ���з��ش���
            % % %     %     [output_CoordLUBin2,output_LU_LWH2,output_LU_Seq2]= getReturnBBA1(do2Array(ibin)); %% ���з��ش���
            % % %     TPP_Seq1 = finalCheck([TPP_Coord,TPP_LWH,TPP_Seq],TLUIN_PP1);
            % % % %     if  sum(flagTiled(:)~=0)>1
            % % % %         sum(flagTiled(:)~=0)
            % % % %         1
            % % % %     end
            % % % %     if flagTiled(ibin)==1
            % % % %         if sum(flagTiled(:)==1)>1
            % % % %             sum(flagTiled(:)==1)
            % % % %             1
            % % % %         end
            % % % %         TPP_Seq = finalCheck([TPP_Coord,TPP_LWH,TPP_Seq],TLUIN_PP1); %����1�����㣻 ����2��ԭʼ.
            % % % %     end
            % % % %     if flagTiled(ibin)==2
            % % % %         if  sum(flagTiled(:)==2)>1
            % % % %             sum(flagTiled(:)==2)
            % % % %             1
            % % % %         end
            % % % %         TPP_Seq = finalCheck([TPP_Coord,TPP_LWH,TPP_Seq],TLUIN_PP2); %����1�����㣻 ����2��ԭʼ.
            % % % %     end
            % % %
            % % %     % �ҳ�ƽ��ibin�����е������߼�ֵ
            % % %     flagTileLUIdx = T_Seq1{:,'BINID'} == ibin;
            % % %
            % % %     % �ص��Ǹ�������ͳ����
            % % %     T_Coord{flagTileLUIdx,:} = TPP_Coord{:,:};
            % % %     T_LWH{flagTileLUIdx,:} = TPP_LWH{:,:};
            % % %     % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
            % % %     T_Seq1{flagTileLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}} = ...
            % % %         TPP_Seq1{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}};
            
            
            %         T_Seq = finalCheck([T_Coord,T_LWH,T_Seq],TLUIN); %����1�����㣻 ����2��ԭʼ.
            
            % ������ز�����Ҫ�滻: ��1���������ڳ��ͺ� ��3�����̳��ڰ���˳�� ��8: չʾ˳��
            %     output_LU_Seq([1,3,8],flagTileLUIdx) = output_LU_Seq2([1,3,8],:);
            %     output_LU_Seq([1,3,4,5,7,8],flagTileLUIdx) = output_LU_Seq2([1,3,4,5,7,8],:); %[1,3,4,5,7,8]��ʾ���޸�����ļ���
            
            %     i = [4,5,7]; % ��ЩӦ���ǲ�����
            %     if sum(output_LU_Seq(i,flagTileLUIdx) ~= output_LU_Seq3(i,:) ) >0,
            %          output_LU_Seq(i,flagTileLUIdx) ~= output_LU_Seq3(i,:) ;         output_LU_Seq(i,flagTileLUIdx);         output_LU_Seq3(i,:) ;
            %         warning('�����ı���, ����'); end
    end
end
% ****************** ��Գ���ѡ�� ��ȡ�޶��� output ******************

T = getShowSeq(T); %����ShowSEQ
%% CHECK 1 ��������

%% tableתΪ�ṹ����ж��Ƿ���������
% lu = table2struct(T,'ToScalar',true)
% lu = (structfun(@(x) x',lu,'UniformOutput',false));
    
    %% ����T������LU_LUinBin���������ж�
%     T2 = sortrows(T,{'BINID','BINSEQ'},{'ascend','ascend'})
    %     x=T(T.BINID==2&T.ID==1,{'ID','LID','PID','H','Weight','CoordLUBin','BINSEQ','ShowSEQ','ITEMID','ITEMSEQ'})
%     x=T2(:,{'ID','LID','PID','H','Weight','X','Y','Z','ITEMID','ITEMSEQ','BINID','BINSEQ','ShowSEQ'})
    % if ITEMID��ͬ ������CoordLUBin�ĳ��������ͬ ��ͬITEMID�����������Ա�
    % ��table��ʽLU����check ��Ҫ������

%        checktLU(T)  %      ���治ͨ��,������LU_Itemδ��ʱ����,�ں���˦βƽ�̺�. TODO

% ����BBA�����ʽ��JAR
[~,T.ttt] = sort(T.tblorder);
T = sortrows(T,'ttt');
% T = sortrows(T,'BINID')
% T.LID
% T.BINID
% T.BINSEQ

output_CoordLUBin=T.CoordLUBin';
output_LU_LWH=T.LWH';
% output_LU_Seq=T{:,{'LU_VehType','BINID','BINSEQ','OSID','LID','ITEMID','OPID','ShowSEQ','Weight'}}'
% ITEMID���岻��, ��
output_LU_Seq=T{:,{'LU_VehType','BINID','BINSEQ','OSID','LID','ITEMID','OPID','ShowSEQ','Weight','Index'}}'; %���ӷ�����10: LuIndex������ǿֻҪ�����־Ϳ���
% output_LU_Seq([2,3,5,8],:)

if ISplotBBA
    plotSolutionBBA(output_CoordLUBin,output_LU_LWH,output_LU_Seq,do); 
    V = struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false));
    plotSolutionT(T,V);
end

    % if  ISplotSolu 
    % %      plotSolution(daBest(bestOne),paBest(bestOne)); %�������� ����plotStrip ��������������ͼ
    % %      if max(do2.LU.LU_Bin(1,:)) == 1
    % %      plotSolution(do2,pA(iAlg));
    % %      end
    % %        if flaggetSmallVeh,   plotSolution(do1,paBest(bestOne));   end %�������� ����plotStrip ��������������ͼ
    % end

    
%% �������
% T.Properties.VariableNames
% % �򵥲鿴ĳ��bin�ڵĲ����ӱ�
% tmpT = sortrows(T,{'ITEMID'})
% x=T(T.BINID==2&T.ID==1,{'ID','LID','PID','H','Weight','CoordLUBin','BINSEQ','ShowSEQ','ITEMID','ITEMSEQ'})
% sortrows(x,{'CoordLUBin'})
% �޳�չʾ˳�� % output_LU_Seq = output_LU_Seq(1:7,:);
% whos

clearvars -except output*
fprintf(1,'Simulation done.\n');

% mcc -W 'java:BBA_Main,Class1,1.0' -T link:lib BBA_Main.m -d '.\new'
% printstruct(do,'sortfields',1,'PRINTCONTENTS',0);    printstruct(do.Veh);
% do = rmfield(do, {'Veh', 'LU'});
% pcode 'H*.m'
end %END MAIN







%% ******* �ֲ����� ****************
% ÿ��RunAlgorithm���ж�d.LU������������Ĳ��
function checkLU(TIN,TOUT)
    TOUT.LWH(TOUT.Rotaed,1) = TIN.LWH(TOUT.Rotaed,1);
    TOUT.LWH(TOUT.Rotaed,2) = TIN.LWH(TOUT.Rotaed,2);
    %% 1.1 LWH CHECK
    if any(TOUT.LWH ~= TIN.LWH),
        any(TOUT.LWH ~= TIN.LWH);
        error('LWH');
    end
    %% 1.2 Weight CHECK
    if any(TOUT.Weight ~= TIN.Weight) || any(TOUT.OPID ~= TIN.PID) || any(TOUT.OSID ~= TIN.SID)
        any(TOUT.Weight ~= TIN.Weight),
        error('other');
    end
end

% ָ��ibin,��ȡ��bin�ڵ�LU,Veh����Ϊ��������,�ص���LU����
function thisd = getdinThisVeh(tmpd,luIdx)
        % 1 Veh��Par
        thisd.Veh = tmpd.Veh;        
        thisd.Par = tmpd.Par;
        thisd.LU = structfun(@(x) x(:,luIdx),tmpd.LU,'UniformOutput',false);

        % 2 bin�ڵ�LU
%         luIdx = tmpd.LU.LU_Bin(1,:) == ibin;    %tmpd.LU.LU_Strip(1,:) == istrip
        
        % thisd.Veh = rmfield(thisd.Veh,{'Volume','order'});
%         thisd.LU.LWH([1,2], thisd.LU.Rotaed ) = flipud(thisd.LU.LWH([1,2], thisd.LU.Rotaed)); %LU.LWH ����ת,��ָ�ԭ��
%         thisd.LU.PID = thisd.LU.OPID;     thisd.LU.SID = thisd.LU.OSID;  %  thisd.LU.LID = thisd.LU.OLID;
%         thisd.LU = rmfield(thisd.LU,{'Rotaed','order','LU_Item','DOC','LU_Strip',...
%             'LU_Bin','CoordLUBin','CoordLUStrip','LU_VehType','OPID','OSID'});
        
        
        
%     % tmpd�е�Bin��������, ����С�Ŀ�ʼ��
%     tmpusedVehIdx = max(tmpd.LU.LU_Bin(1,:)); %tmpusedVehIdx: ���һ��Bin��indexֵ
%     flagusedLUIdx = tmpd.LU.LU_Bin(1,:)==tmpusedVehIdx; % flagused: �ҳ����һ��Bin��Ӧ��LUindexֵ
%     if isSameCol(tmpd.LU)
%         % ��ȡ�����һ��Bin����������
%         lastd.LU = structfun(@(x) x(:,flagusedLUIdx),tmpd.LU,'UniformOutput',false);  %��ȡ���һ�����ڵ�LU
%         lastd.LU.LWH([1,2], lastd.LU.Rotaed ) = flipud(lastd.LU.LWH([1,2], lastd.LU.Rotaed)); %LU.LWH ����ת,��ָ�ԭ��
%         lastd.LU = rmfield(lastd.LU,{'Rotaed','order','LU_Item','DOC','LU_Strip','LU_Bin','CoordLUBin','maxL','CoordLUStrip'}); 
%         lastd.Par = tmpd.Par;
%     else
%         error('����ʹ��structfun');
%     end
end

function lastd = getdinLastVeh(tmpd)
    % tmpd�е�Bin��������, ����С�Ŀ�ʼ��
    tmpusedVehIdx = max(tmpd.LU.LU_Bin(1,:)); %tmpusedVehIdx: ���һ��Bin��indexֵ
    flagusedLUIdx = tmpd.LU.LU_Bin(1,:)==tmpusedVehIdx; % flagused: �ҳ����һ��Bin��Ӧ��LUindexֵ
    if isSameCol(tmpd.LU)
        % ��ȡ�����һ��Bin����������
        lastd.LU = structfun(@(x) x(:,flagusedLUIdx),tmpd.LU,'UniformOutput',false);  %��ȡ���һ�����ڵ�LU
        lastd.LU.LWH([1,2], lastd.LU.Rotaed ) = flipud(lastd.LU.LWH([1,2], lastd.LU.Rotaed)); %LU.LWH ����ת,��ָ�ԭ��
        lastd.LU = rmfield(lastd.LU,{'Rotaed','order','LU_Item','DOC','LU_Strip','LU_Bin','CoordLUBin','maxL','CoordLUStrip'}); 
        lastd.Par = tmpd.Par;
    else
        error('����ʹ��structfun');
    end
end

%% **** �㷨ָ��ѡ�����Ž� ****
function [daMax,parMax] = getbestsol(DaS,Par)

    % �������һ�����н�, ֱ�ӷ���;
    if size(DaS,2)==1    %����dA�ж��ʱ����,Ŀǰ��������,
        daMax = DaS(1); parMax = Par(1);
        return
    end
    
%��ȡ����ָ��Ͷ�Ӧ����
for r=1:length(DaS)
    resLoadingRateBin(r) = mean(DaS(r).Bin.loadingrate); %bin��װ���ʾ�ֵ��� Itemloadingrate ItemloadingrateLimit
    resLoadingRateStripLimit(r) = mean(DaS(r).Strip.loadingrateLimit); %strip��limitװ������� Itemloadingrate ItemloadingrateLimit
    resLoadingRateBinLimit(r) = mean(DaS(r).Bin.loadingrateLimit); %bin��limitװ������� Itemloadingrate ItemloadingrateLimit
    resLoadingRateStrip(r) = mean(DaS(r).Strip.loadingrate); %strip��װ������� Itemloadingrate ItemloadingrateLimit    
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
idxBin=find(resLoadingRateBin==max(resLoadingRateBin)); %ȡ
idxStripLimit=find(resLoadingRateStripLimit==max(resLoadingRateStripLimit));
idxBinLimit=find(resLoadingRateBinLimit==max(resLoadingRateBinLimit));
idxStrip=find(resLoadingRateStrip==max(resLoadingRateStrip));
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
%      idx = idx1; 
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

end % END OF ALL




%% ********************** �����ǻ������õĴ��� ��ʱ���� ****************

        
function plotSolution(do,par)
%% ��ͼ
% V3 margin ��ǰ��RunAlgorithm���к��ִ��:
% plot2DBPP(do,par);
plot3DBPP(do,par);

        % V1 buff version
        % do.Item.LWH = do.Item.LWH - do.LU.buff(:,1:size(do.Item.LWH,2));
        % do.Item.LWH(1,:) = do.Item.LWH(1,:) - ( do.LU.margin(1, 1:size(do.Item.LWH,2) ) + do.LU.margin(2,: )); 
        % do.Item.LWH(2,:) = do.Item.LWH(2,:) - (do.LU.margin(3,: ) + do.LU.margin(4,: )); 
        % do.Item.CoordItemBin = do.Item.CoordItemBin + do.LU.buff(:,1:size(do.Item.LWH,2))/2;

        % V2 margin version
        % ��ͼǰ����LU ITEM��Coord��LW; ����ITEMͬʱ����LU
        % [do.LU,do.Item] = updateItemMargin(do.LU,do.Item);
        

end


    %% ************ �ж��Ƿ���ͬ�����������ڰڷ�
%     function flag = isAdjacent(d)
%         flag = 1;
%         printstruct(d);
%         % ÿ��bin���ҳ�������ID����Strip�Ƿ�����
%         nBin = size(d.Bin.LW,2);
%         for iBin = 1:nBin
%             t = [d.Item.LID; d.Item.Item_Strip; d.Item.Item_Bin ];
%             tiBin = t( : , t(4,:) == iBin );
%             nIdType = unique(tiBin(1,:)); %nIdType: ��iBin�ڰ�����LU��ID����
%             for iId = 1:nIdType
%                 tiId = tiBin( : , tiBin(1,:) == iId );
%                 nIdStrip = unique(tiId(2,:)); %nIdStrip: ��iBin����iID�°�����Strip�����
%                 % �ж������ķ��뱾ID���͵�Strip����Ƿ�����
%                 if ~all(diff(sort(nIdStrip))==1)
%                     flag = 0;
%                 end
%             end                    
%         end
%     end
    
% % %             % ns - ��bin��strip������˳��
% % %             ns = d.Strip.stripBeBinMatrix(2,d.Strip.stripBeBinMatrix(1,:) == iBin);
% % %             % ni - ��bin��item��LU���ͼ�˳��
% % %             d.Item.Item_Bin(1,:) == iBin
% % %             ni = d.Item.LID(d.Item.Item_Bin(1,:) == iBin);
% % %             [a,b] = find(d.Item.LID(d.Item.Item_Bin(1,:) == iBin));
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
% % %                 LUIDInis = d.Item.LID(1,(d.Item.Item_Strip(1,:)==is))
% % %                 LUIDInjs = d.Item.LID(1,(d.Item.Item_Strip(1,:)==js))
% % %                 
% % %                 end
% % %             end
% % %         end
% % %         



%% ********************** ������ts�㷨�Ĵ��� ��ʱ���� ****************

%%

    %% ע��
    
%         while (all(do2.Strip.isHeightFull(fi) == 1) && all(do2.Strip.isWidthFull(fi) == 1))
% %                  || all(d2.LU.maxHLayer(luidxPP)==1)% Ŀǰ���ܶ����һ��strip����, ���������һ��Strip�ڵ�Lu��maxHLayerȫ��Ϊ1
%              
%             istrip = istrip-1;
%             fi = find( do2.Strip.Strip_Bin(2,:) >= istrip ); 
%             luidxPP = ismember(do2.LU.LU_Strip(1,:), u(fi)); %%% fi->u(fi) ��������� ********************* 
%             if ~any(luidxPP),  error('luidxPPȫ��Ϊ��, ������u(fi)��Ӧ��Lu�߼��ж�'); end
%             if istrip == 1,  error('��bin������tileneed,��Ԥ�ڴ���');   end
%             d2.LU.maxHLayer(luidxPP)
%         end
        
        
        % ����LU(luidxPP),�޸���maxHLayer
%         u=unique(do2.LU.LU_Strip(1,:)); %��ȡStrip��ŵ�Ψһ����ֵ
%         luidxPP = ismember(do2.LU.LU_Strip(1,:), u(fi)); %%% fi->u(fi) ��������� *********************
%         if ~any(luidxPP),  error('luidxPPȫ��Ϊ��, ������u(fi)��Ӧ��Lu�߼��ж�'); end
%         do2.LU.LU_Item(1,luidxPP)


%                     d2.LU.tmpHLayer = zeros(size(d2.LU.maxHLayer))
%                     d2.LU.tmpHLayer(luidxPP) = do2.Item.HLayer(do2.LU.LU_Item(1,luidxPP))
%                     dlu=[ d2.LU.maxL(3,luidxPP);
%                      d2.LU.maxHLayer(luidxPP);
%                      d2.LU.tmpHLayer(luidxPP)];
%                      min(dlu)
%         d2.LU.maxHLayer(luidxPP) = min(dlu) - 1;

% % % [ub,x,b] = HnextFit(Item,Veh);
% % % [ub,x,b] = HnextFit_origin(Item,Veh);
% % % disp(b');
% % % disp(x);
% % % fprintf('UB = %d \n', ub);
% % % [ub,x,b] = TSpack(d,n,w,W,lb,timeLimit, ub0,x,b,whichH);
% % 
% % function [toReturn,x,b] = TSpack(d,n,w,W,lb,timeLimit, ub0,x,b,whichH)
% % [nb,x,b] = heur(d,n,w,W,x,b,whichH);
% % ub0 = nb;
% % 
% % if nb == lb
% %     toReturn = nb;
% %     fprintf('best lb = %d ', nb);
% % end
% % 
% % %/* initial (trivial) solution */
% % cnb = n;
% % cb = 1:n;
% % cx = zeros(d,n);
% % cw = w;
% % 
% % %/* external loop */
% % D = 1; tt = 0.0;
% % toReturn = nb;
% % 
% % end
% % 
% % function [nb,x,b] = heur(d,n,w,W,x,b,whichH)
% % nb = -1;
% % which = (d-2)*100 + whichH; % which Ϊ0��100
% % if which == 0
% %     [nb,x,b]  = HnextFit(n,w,W,x,b,n+1); %/* first heuristic for 2d bin packing */
% %     %          disp(x);
% %     %          disp(b);
% % elseif which == 100
% %     [nb,x,b]  = HHnextFit(n,w,W,x,b); %/* first heuristic for 3d bin packing */
% % end
% % end
% % 
% %     function [ub,px,pb]  = HnextFit(Item,Veh)
% %         % Initialize
% %         d = size(Item.LWH,1)-1;
% %         n = size(Item.LWH,2);
% %         nn = n + 1;
% %         w = Item.LWH(1:d,:);
% %         W = Veh.LWH(1:d,:);
% %         x = zeros(d,n); b = zeros(n,1); bNb = zeros(n,1);
% %         
% %         %/* sort the items */
% %         % sortD = size(w,1);%��ȡ��Ҫ�����ά��
% %         [~,ord] = sort(w(d,:),'descend');%��w��������,ֻ��Ҫ����˳��ord;����d�����򣨸߶�)
% %         pw = w(:,ord);
% %         px = x;
% %         pb = (b+999); % 0 + 999
% %         pbNb = bNb;
% %         
% %         %/* next fit packing */
% %         % binLeftArray(1,ub) �� wleft
% %         % binLeftArray(2,ub) :  hleft
% %         nBin = n;
% %         binLeftArray = repmat(W,1,nBin);  %��ʼ
% %         ub = 1;
% %         for i=1:n
% %             %     if (binLeftArray(1,ub) == W(1)) & (binLeftArray(2,ub) == W(2)) %����ǿ�bin
% %             if pbNb(ub) == 0   %����ǿ�bin
% %                 if (pw(1,i) <= binLeftArray(1,ub)) && (pw(2,i) <= binLeftArray(2,ub)) %�����߶�������
% %                     px(1,i) = 0; px(2,i) = 0;
% %                     binLeftArray(1,ub) = binLeftArray(1,ub) - pw(1,i);
% %                     binLeftArray(2,ub) = binLeftArray(2,ub) - pw(2,i);
% %                     pbNb(ub) = pbNb(ub) + 1;
% %                 else
% %                     error('EEE');
% %                 end
% %             else               %������ǿ�bin
% %                 if pw(1,i) <= binLeftArray(1,ub)  %���i�Ŀ����㵱ǰbin��ʣ���ȣ�ʣ��߶�Ӧ�ò���
% %                     px(1,i) = W(1) - binLeftArray(1,ub);
% %                     px(2,i) = W(2) - binLeftArray(2,ub) - pw(2,i);     %�߶�Ϊ????
% %                     binLeftArray(1,ub) = binLeftArray(1,ub) - pw(1,i);
% %                     binLeftArray(2,ub) = binLeftArray(2,ub);
% %                     pbNb(ub) = pbNb(ub) + 1;
% %                 else
% %                     if pw(2,i)  <= binLeftArray(2,ub)  %���i�ĸ����㵱ǰbin��ʣ��߶�
% %                         px(1,i) = 0;
% %                         px(2,i) = binLeftArray(2,ub);
% %                         binLeftArray(1,ub) = W(1) - pw(1,i);
% %                         binLeftArray(2,ub) = binLeftArray(2,ub) - pw(2,i);
% %                         pbNb(ub) = pbNb(ub) + 1;
% %                     else  %���i�ĸ߲������㵱ǰbin��ʣ��߶�
% %                         ub = ub + 1;
% %                         px(1,i) = 0;   px(2,i) = 0;
% %                         pbNb(ub) = pbNb(ub) + 1;
% %                         binLeftArray(1,ub) = binLeftArray(1,ub) - pw(1,i);
% %                         binLeftArray(2,ub) = binLeftArray(2,ub) - pw(2,i);
% %                     end
% %                 end
% %             end
% %             pb(i) = ub-1;
% %         end
% %         
% %         
% %         %ԭʼ��
% %         x = px(:,ord);
% %         b = pb(ord);
% %         bNb = pbNb(ord);
% %         ub = ub +1;
% %     end
% % 
% %     function [ub,px,pb]  = HnextFit_origin(Item,Veh)
% %         % Initialize
% %         d = size(Item.LWH,1);
% %         if d==3
% %             d=d-1;
% %         end
% %         n = size(Item.LWH,2);
% %         nn = n + 1;
% %         w = Item.LWH(1:d,:);
% %         W = Veh.LWH(1:d,:);
% %         x = zeros(d,n); b = zeros(n,1); bNb = zeros(n,1);
% %         
% %         %/* sort the items */
% %         sortD = size(w,1);%��ȡ��Ҫ�����ά��
% %         [~,ord] = sort(w(sortD,:),'descend');%��w��������,ֻ��Ҫ����˳��ord
% %         pw = w(:,ord);
% %         ord;
% %         
% %         px = zeros(size(x,1),size(x,2));
% %         pb = 999*ones(size(b,1),size(b,2));
% %         %/* next fit packing */
% %         hleft = W(2) - pw(2,1);
% %         wleft = W(1);
% %         ub = 0; hcurr = 0;
% %         for i=1:n  %�ӵ�һ��item��ʼ����
% %             if pw(1,i) <= wleft  %���item��w ��wleftС������item����bin���㣺����x1ֵ������wleft��hleft����
% %                 px(1,i) = W(1) - wleft;
% %                 wleft = wleft - pw(1,i);
% %             else    %��������һ�㰲�š�
% %                 if pw(2,i) <= hleft  %���item��h ��hleftС������bin�߶ȳ��㣬����item����һ������������hleft������hcurr��wleft�� ��������xֵ������wleft
% %                     hcurr = W(2) - hleft; %������ͬһ�㣬����hcurr���䣬Ҳ����pw(2,1)(��������bin�Ͳ����ˣ�������hcurr)
% %                     hleft = hleft - pw(2,i);
% %                 else  %����Ų��£�����bin������hcurr������hleft����������ub+1������ﵽnnֵ������������������xֵ0������wleft
% %                     hcurr = 0;    %�������µ�bin������hcurrΪ0;
% %                     hleft = W(2) - pw(2,i);
% %                     if (ub+1 == nn)
% %                         break;
% %                     end
% %                     ub = ub + 1;
% %                 end
% %                 % ���۷����ϲ����bin������x1Ϊ0������wleftΪW(1)-��item�Ŀ�w
% %                 px(1,i) = 0;
% %                 wleft = W(1) - pw(1,i);
% %             end
% %             % �˴�ͳһ����x1ֵ�����߶�ֵ��Ϊhcurr��ͳһ����bֵ=ub��
% %             px(2,i) = hcurr;
% %             pb(i) = ub;
% %         end
% %         
% %         %ԭʼ��
% %         x = px(:,ord);
% %         b = pb(ord);
% %         ub = ub +1;
% %     end
% % 
% % 
% %     function [ nb ] = HHnextFit(n,w,W,x,b)
% %         
% %         nb = 1;
% %         
% %     end



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
