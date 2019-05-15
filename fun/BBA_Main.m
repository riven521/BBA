function [output_CoordLUBin,output_LU_LWH,output_LU_Seq] = ...
    BBA_Main(LUID,LULWH,VEHID,VEHLWH,varargin)  %ǰ4������
% Inputs��LUID: ������ ; LULID: �����
%   LUID	                (1,n)   �������� ��ͬ���ֱ���ͬһ����,����Ѷ�
%   cpuLUVeh����unique������Ǵ�1��ʼ�ϸ����
%   LULWH                (3,n)   ���̿���
%   VEHID                 (1,m)  ���ͱ��
%   VEHLWH              (3,m)   ���Ϳ��ߣ����Ƕ೵�ͣ�
%   ------------------------------------------------------
%   LUSID                   (1,n)   ���̹�Ӧ�̱��
%   cpuLUVeh����unique������Ǵ�1��ʼ�ϸ����
%   LUPID                   (1,n)   �����㲿�����
%   cpuLUVeh����unique������Ǵ�1��ʼ�ϸ����
%   LUISROTA            (1,n)  �����Ƿ�������ת
%   LUMARGIN          (4,n)   ���̼�margin(1-4��������)  �������̳����=ÿ�����̵�ʵ�ʳ����+���ӵ�margin
%   LUWEIGHT           (1,n)  ��������
%   VEHWEIGHT        (1,m)  ��������������
%   LULID                   (1,n)   �������ͱ�� ���޹涨��
%   LUINDEX            (1,n�� ����������-��ǿר��
%   sort������Ǵ�1��ʼ�ϸ����
%   LUEID                  (1,n�� ����EP LOCATION - Milkrun�汾����
%   cpuLUVeh����unique������Ǵ�1��ʼ�ϸ����
% Outputs
%   output_CoordLUBin      (3,n)    ÿ��LU��X,Y,Z
%   output_LU_LWH            (3,n)    ÿ��LU�Ŀ��ߣ���ת��ģ�ʵ��ֵ��
%   output_LU_Seq             (8,n)    ��1: LU��ĳ��BIN�ڣ���2: LU�ڸ�BIN�ڵİ���˳�� ������

%% Initialize Global Variable
% rng('default');rng(1); % NOTE �Ƿ�����ı�־
if nargin < 1
    clear all;
end

close all; format long g; format bank; 
warning('off');

if datetime('now') > datetime(2020,03,03)
    error('There is something wrong in the matlab version used');
end

% ����LOG FILE OUTPUT - ��ʹBBAҲҪ���Ӹ���־
curTime = datestr(now,30);
mkdir('Log');  %Ϊ��־����Ŀ¼
filepath = fullfile(pwd, strcat('.\Log\')); 
filename = strcat(filepath,curTime,'T',num2str(length(LUID)),'.txt');
fid = fopen(filename,'w'); 
fprintf(fid,'LUID = ['); fprintf(fid,' %.0f',LUID); fprintf(fid,'];\n');
fprintf(fid,'LULWH = ['); fprintf(fid,' %.2f',LULWH); fprintf(fid,'];\n');
fprintf(fid,'VEHID = ['); fprintf(fid,' %.0f',VEHID); fprintf(fid,'];\n');
fprintf(fid,'VEHLWH = ['); fprintf(fid,' %.2f',VEHLWH); fprintf(fid,'];\n');
fprintf(fid,'LUSID = ['); fprintf(fid,' %.0f',varargin{1}); fprintf(fid,'];\n');
fprintf(fid,'LUPID = ['); fprintf(fid,' %.0f',varargin{2}); fprintf(fid,'];\n');
fprintf(fid,'LUISROTA = ['); fprintf(fid,' %.0f',varargin{3}); fprintf(fid,'];\n');
fprintf(fid,'LUMARGIN = ['); fprintf(fid,' %.2f',varargin{4}); fprintf(fid,'];\n');
fprintf(fid,'LUWEIGHT = ['); fprintf(fid,' %.4f',varargin{5}); fprintf(fid,'];\n');
fprintf(fid,'VEHWEIGHT = ['); fprintf(fid,' %.2f',varargin{6}); fprintf(fid,'];\n');
fprintf(fid,'LULID = ['); fprintf(fid,' %.0f',varargin{7}); fprintf(fid,'];\n');
fprintf(fid,'BBAID = ['); fprintf(fid,' %.0f',varargin{8}); fprintf(fid,'];\n');
if nargin >= 13
fprintf(fid,'Para/EID = ['); fprintf(fid,' %.0f',varargin{9}); fprintf(fid,'];\n');     %����13��ʱ���ã�����Ҳ��Ӱ��
end

% ȫ�ֱ���0�� ��ͼ����
% ISplotBBA�� �Ƿ�����T��ͼ���ܿ��أ�  ISplotShowGapAdjust: �Ƿ���ʾgap��϶��������ͼ
% ISplotEachPingPuAll/ISplotEachPingPuShuaiWei������ƽ�̺�˦βƽ�̺󻭶Ա�ͼ 
global ISplotBBA ISplotShowGapAdjust ISplotEachPingPuShuaiWei ISplotEachPingPuAll ISplotPauseWait   
% RunAlgorithm��plot����
% ISplotshuaiwei: �Ƿ���ʾ˦β�Ա�ͼ��ISplotStripToBinAgain���Ƿ���ʾ��ͷ���ȶԱ�ͼ��ISplotRunAlgo:�Ƿ���ʾRunAlgo����ͼ
global  ISplotshuaiwei  ISplotStripToBinAgain ISplotGapCompare  ISplotRunAlgo ISplotRunLIS

% ISplotPause: plotSolutionT ����ͣʱ�� 0����ͣ ISplotShowType������ĳ������������ɫ
global ISplotPause ISplotShowType

if isempty(ISplotBBA),  ISplotBBA = 0;   end   % �Ƿ���ʾLU/Strip/Bin�Ľ������������
if isempty(ISplotEachPingPuShuaiWei),  ISplotEachPingPuShuaiWei = 0;   end % ÿ��˦βƽ�̳ɹ���չʾƽ��ǰ��ĶԱ�ͼ
if isempty(ISplotEachPingPuAll),  ISplotEachPingPuAll = 0;   end % ÿ������ƽ�̳ɹ���չʾƽ��ǰ��ĶԱ�ͼ
if isempty(ISplotShowGapAdjust),  ISplotShowGapAdjust = 0;   end  

if isempty(ISplotRunAlgo),  ISplotRunAlgo = 0;   end  
if isempty(ISplotRunLIS),  ISplotRunLIS = 0;   end  
if isempty(ISplotshuaiwei),  ISplotshuaiwei = 0;   end  
if isempty(ISplotStripToBinAgain),  ISplotStripToBinAgain = 0;   end  
if isempty(ISplotGapCompare),  ISplotGapCompare = 0;   end  

if isempty(ISplotPauseWait),  ISplotPauseWait = 0;   end   % �Ƿ�plotsolutinT���ͼֱ�ӵȴ��û���Ӧ
if isempty(ISplotPause),  ISplotPause = 0.0;   end   %-0.05 % plot���ʱ��
if isempty(ISplotShowType),  ISplotShowType = 3;   end   % 1 LID 3 ID 8 ˦β

% ȫ�ֱ���1�� �汾����
% verMilkRun��1 milkrun�汾��������һ��EID���룩 0 ��milkrun��
global verMilkRun

if isempty(verMilkRun),  verMilkRun = 0;   end % 555: Ĭ�ϲ���MilkRun�汾(MilkRun��9�������İ汾)
% MR�Ĳ���9: EP LOCATION����, ��ı������Զ�����
if nargin > 1 && length(varargin) < 9 %'LUEID',varargin{9});
    verMilkRun  = 0;   varargin{9} = ones(1,length(LUID));  
elseif nargin > 1 && length(varargin{9}) ~= length(LUID) %������para���룬��milkrun���� todo fixme �����޸�
    verMilkRun  = 0;   varargin{9} = ones(1,length(LUID));  % 9���������Ϊmilkrun�汾    
elseif nargin > 1 && all(varargin{9}) %% 
    verMilkRun = 1;  % 9���������Ϊmilkrun�汾 ��ҪEID����Ϊ0
elseif nargin > 1 && ~all(varargin{9}) %% && length(varargin{9}) ~= 6
    verMilkRun = 0; varargin{9} = ones(1,length(LUID));  % 9��������� ������9��, ����ֵΪ0��, Ҳ�������汾
else
    error('varargin 9 is wrong');
end

% *********  ���ܿ��� ******* Ĭ�ϲ���

% ȫ�ֱ���2�� ƽ�̿���
%   ISpingpu : �Ƿ�˦βƽ�� ISpingpuall���Ƿ�����ƽ�� ISpingpuShuaiWei���Ƿ�˦βƽ�� ��ѡһ 
%   ��ƽ�̣�ISpingpu���뿪�����˦βƽ�̺�����ƽ�̿�ѡ 
global ISpingpu ISpingpuAll ISpingpuShuaiWei
if verMilkRun == 1
    ISpingpu = 0;         % 555 : ��Ⱥ͸߶Ȳ���, �Ҳ���>1, ƽ��. ���������� (����ƽ�̺���ISisNonMixedì��)
    ISpingpuAll = 0;       %555: ���о�ƽ��, ֻҪ�ó����ŵ���; ���Ų���, ��������˦βƽ������
    ISpingpuShuaiWei = 0; % 555: ˦βƽ��
else
    ISpingpu = 1;          % 555 : ��Ⱥ͸߶Ȳ���, �Ҳ���>1, ƽ��. ���������� (����ƽ�̺���ISisNonMixedì��)
    ISpingpuAll = 1;       %555: ���о�ƽ��, ֻҪ�ó����ŵ���; ���Ų���, ��������˦βƽ������
    ISpingpuShuaiWei = 1; % 555: ˦βƽ��
end

% ȫ�ֱ���3�� ��װ��϶���� ���뿪�Ŷ�
%   parGap �� �Ƿ������϶����������˦βƽ�̻�����ƽ�̳ɹ��ĳ��ڽ��У���parMulipleGap���Ƿ��ڼ�϶���������������ε����������У�
global  parGap parMulipleGap
if isempty(parGap),  parGap = 1;   end                              % �Ƿ������������ļ�϶����
if isempty(parMulipleGap),  parMulipleGap = 1;   end % ����Ϊ1�����ǲ�������ʹ�� �Ƿ������϶�ݹ��ε���

% ȫ�ֱ���3�� ˦β���� ����ͷ���� ���뿪�Ŷ�
global ISshuaiwei ISreStripToBin    % RunAlgorithm��  % ���� + Gpreproc ��V2�汾 �޸�ҵ��3����(��ITEM��������һ��ľ�������)
if isempty(ISshuaiwei),  ISshuaiwei = 1;   end  % 555 : ��Ⱥ͸߶Ȳ���, ˦β   ******  �ò�����Ҫ�������pingpu���ʹ�� ��˦β ˦βƽ���޷�����*******
if isempty(ISreStripToBin),  ISreStripToBin = 1;   end  % ��ͷ����LU����������� Ĭ��Ϊ1 ����
 
% ȫ�ֱ���4�� �ı䳵�Ϳ���
% ISlastVehType�����������һ����Ϊ��С�ĳ���
global ISlastVehType
if isempty(ISlastVehType),  ISlastVehType = 0;   end   % 555: ���һ���ĵ���, �������޹�, �ݲ�����

% *********  ϸ΢�������ܿ��� ******* Ĭ�ϲ���
global ISdiagItem       %cpuItem�� �ж�Item�Ƿ�isHeightFull�߶����㣬��1�����ݶԽ����ж�
ISdiagItem = 0.25;          % Ĭ��Ϊ 0 �� Ϊ1 ����Щ���ڵ͵ı���ΪItem�߶�����, checkԭ���

global parBalance       %cpuStrip�� �ж�Strip�ж���Ѷ��Ƿ�Ϊ�߶Ⱦ��� isHeightBalance
parBalance = 8/30;     % ��ߵ�Item��1/3����������߶Ѷ�Item��1/3��8/30����Ϊ�߶Ⱦ��⣬����֮һ

% ȫ�ֱ���7�� ��ɾ��
% global ISisGpreprocLU1 ISsItemAdjust  ISreStripToBinMixed ISisMixedStrip
% ISisGpreprocLU1 = 1 % ����1; ��ISisNonMixedLUȡֵ�й� ���ISstripbalanceʹ�� 1��ʾ��ͬһˮƽstrip��ITEM������2�����ϵ��ж�ΪISisNonMixedLU=1->�Ѷ����ʹ�� 0 ��ʾ�����ж�
% ISsItemAdjust = 1              % ��ʱ���� ��;������
% ISreStripToBinMixed = 1   %��ͷ���ȷ�AllPure����, �ٿ�������LU����������� Ĭ��Ϊ1 Ӧ�ÿ���ɾ���Ĳ���
% ISisMixedStrip = 1 ����ɾ���� % 1��ʾ����LU.ID�ж��Ƿ��� 0����LU.LID�ж� NOTE: ����STRIP  NOTE: ITEM��ΪID�滻LID

% ȫ�ֱ���8�� ��ע��
% global ISisNonMixedLU ISisMixTileLU         %LU2ITEM
% global ISisNonMixed ISisMixTile                 %ITEM2STRIP
          
% ISisNonMixedLU = 1; % 555: ���ȷǻ��LU�γ�ITEM, ͼ�ÿ���� ������ Ĭ��Ϊ 1 ��1 LU�����γ����� 0 �ض��з��������� LU��������
% ISisMixTileLU = 1;       % 555: ���Ȼ��LU�ĵ���ITEM�������γ�ITEM, ͼ�ÿ���� ������ Ĭ��Ϊ 1 ��1 ��isNonMixed=0ʱ, ���������Ӧ��LU��ֵΪ1�����LU��������ITEM֪ʶ��

% ISisNonMixed = 1;    % 555: ���ȷǻ��Item�γ�STRIP, ͼ�ÿ���� ������ Ĭ��Ϊ 1
% ISisMixTile  = 1;         % 555: ���Ȼ��Item�ĵ���Strip�������γ�STRIP, ͼ�ÿ���� ������ Ĭ��Ϊ 1 �����ܳ��ֻ������

% global ISmaxLayer
% ISmaxLayer = 1

% global ISstripbalance  ISplotstripbalance
% ISstripbalance = 1;     % ���ø߶Ⱦ��⿪�� 555������ͼ�κÿ�, �Ѷ����ʹ�� ͬһStrip�ǻ���Ҹ߶Ȳ�������LU��������ֵ>1ʱ���� ����������ӦLU���������ݼ�; ���޷��� % ISstripbalance =0;  %��������ʱ,����ISisGpreprocLU1=0; ��ʾLUʹ���߶ȶ�,�����ISisNonMixedLU=0.
%  ISplotstripbalance = 1;  

% RunAlgorithm��
% global   ISplotStrip  % plotStrip
% ISplotStrip = 0;             % ÿ��Run algorithm ����Strip����ʾ��� ��ϸ�� �����滻ΪplotSolutionT


%% Initialize Parameter Variable
nAlg = 1;
pA(nAlg) = InitializeParameter('whichStripH', 3,...
                             'whichBinH',3, ...
                             'whichSortItemOrder',3, ... 
                             'whichRotation',2, ...
                             'whichRotationHori', 1);
                         
%% Initialize Data Structure  - d
if nargin < 1 % Randome Generate
    n=32; m=1;                                          % 16��Ҫע�� 250 srng1
    d = InitializeRandomData(n,m);        % 0 Ĭ��ֵ; >0 �����������n������ ����ֱ������BBAʱ����
%     save( strcat( '.\useless\', strcat('GoodIns',num2str(n))), 'd');          %     load .\new\GoodIns200.mat;
else
    d = InitializeInputData( ...
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
            'LUINDEX',varargin{8},...
            'LUEID',varargin{9}); 
end

%% 1��ȷ��LU/Veh�ṹ���������ͬ �� �Ǹ��ľ���
if ~isSameCol(d.LU) || ~isSameCol(d.Veh)
    warning('LU��Veh��������ͬ');
end
structfun(@(x) validateattributes(x,{'numeric'},{'nonnegative','2d'}), d.LU, 'UniformOutput', false);
structfun(@(x) validateattributes(x,{'numeric'},{'nonnegative','2d'}), d.Veh, 'UniformOutput', false);

%% 2:  �ڴ˻����Ͻ����������е����������޸ģ�LWH, ID�ŵȣ�
fprintf(1,'\nRunging Checking and preprocessing Input Data ...\n');

d = chkInput(d);

if ~isSameCol(d.LU) || ~isSameCol(d.Veh)
    error('LU��Veh��������ͬ');
end

%% 3���������󽫽�����������ݽ���check Gpreproc=cpuLUVeh d����Ҫ������

[d.LU,d.Veh] = cpuVehLU(d.LU,d.Veh);

% d = chkInput(d); % �Ͻ���check����Ϊ��������ݣ��ص�ԭʼ״��

% ************* ���������Ҵ����������� ***********
maind = d; % ��Ҫ���������ݱ�����
%��NOTE : ********** maind��LWH����Ϊ����margin�� ********��
%��NOTE : ********** maind��Rotaed����Ϊfalse�� ********��
%��NOTE : ********** maind��LU���ڽ���runalgorithm֮ǰ�ٶ������margin********��

if verMilkRun == 0 && ~isscalar(unique(d.LU.SID))
    error('ĿǰΪ��MilkRun�汾,�������ж�ҹ�Ӧ�̱��');
end
if verMilkRun == 1 && (isscalar(unique(d.LU.SID)) && isscalar(unique(d.LU.EID)))
    error('ĿǰΪMilkRun�汾,������ֻ�е�����Ӧ�̱��/����EP LOCATION���');
end            % if verMilkRun == 1 && (~isscalar(unique(d.LU.SID)) && ~isscalar(unique(d.LU.EID))),  error('ĿǰΪMilkRun���԰汾,ֻ��������Ӧ�̱�� �� ���EP LOCATION���,����ͬʱ����');  end
      
%% Simulate - All ALGORITHM

for iAlg = 1:nAlg
    %% 1 �������㷨 ����d ��� do
    fprintf(1,'\nRunning RunAlgorithm with d ...\n');
    
    % 1.1 �޶�LU��LWH���� ����margin ����Rotaed
    d = maind;
    
    [d.LU] = setLULWHwithbuff(d.LU, d.Veh);

    % 1.2 ������RunAlgorithm
    do = RunAlgorithm(d,pA(iAlg));   %��ȡ��ṹ��do ( ����d���̶�LU����͹̶�LU��С�͹̶�margin��LWH��   % ��ɾ - ����RunAlgorithm������ do.LU.LU_VehType = ones(size(do.LU.ID)) * do.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ

    % 1.3 �޶�LU��Item��LWH/Coord���� ɾ��margin ����rotated
    [do.LU,do.Item] = setLCwithoutbuff(do.LU,do.Item);
    
    % 1.4 CHECK ����������LU
    chkLUnewold(maind.LU,do.LU); % bug���� : chkLUnewold(d.LU,do.LU);  %�������ݶԱ� todo �����������ݶԱ� Ԥ���㷨�ڲ�����
    
    chktLU(do.LU); %����LU/Item/Strip���� �Ƿ��������������أ��Ƿ�Z�߶���ItemSEQһ�£��Ƿ�ITEMID��XY����һ��; todo :���ƽ���˶Ժ���
    
    % 1.5 % dA(iAlg)=do; % ���ڶ��㷨�汾ʹ��    

    ISplotRunAlgo = 0;  %     ISplotRunLIS = 0;
    %% 2 ���г��͵����㷨,���ı�d ��ȡd1��do1, flaggetSmallVeh
    if ISlastVehType == 1
        fprintf(1,'\nRunning HBinChange with output do1 ...\n');
        [flaggetSmallVeh,do1] = HBinChange(maind,do,pA(iAlg));  
    end

    %% 3 ���жѶ�ƽ�̣�������ƽ�̺�˦βƽ�̣��㷨,���ı�d ����maind 
    % d2Array(˦βƽ�̣� d3Array������ƽ�̣� ���������ݣ�ÿ����һ��bin
    % ��ȡ do2Array do3Array  ��������ݣ�ÿ��������һ��bin
    % ��ȡ flagTiledArray �߼��ж���1 ��bin˦β�ɹ� 0 ���ɹ�
    if ISpingpu==1
        fprintf(1,'\nRunning HBinpingpu with output do2 and do3 ...\n');
        [flagTiledArray,do2Array,do3Array] = HBinpingpu(maind,do,pA(iAlg));            
    end
    
    %% 4 ���м�϶�����㷨 ���do3
    if parGap == 1
        fprintf(1,'\nRunning HBinGapAdjust with output do2 and do3 ...\n');
        if ISpingpu==1 && any(flagTiledArray ~= 0) %��������˦βƽ�̺�����ƽ�̳ɹ��ĳ����У�todo ���Ÿ�����bin
            idxArray = find(flagTiledArray~=0);
            for i=1:length(idxArray)
                if flagTiledArray(idxArray(i)) == 1
                    fprintf(1,'       Running ��װ��϶_����ƽ�� (do3)...\n');    
                    dd = do3Array(idxArray(i)); 
                elseif flagTiledArray(idxArray(i)) == 2
                    fprintf(1,'       Running ��װ��϶_˦βƽ�� (do2)...\n');    
                    dd = do2Array(idxArray(i)); 
                end
                
                [dd.LU] = HBinGapAdjust(dd.LU, dd.Veh);
                
                if flagTiledArray(idxArray(i)) == 1
                    do3Array(idxArray(i)) = dd;
                elseif flagTiledArray(idxArray(i)) == 2
                    do2Array(idxArray(i)) = dd;
                end                
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
% ����do do1 do2 do3 flagTiledArray flaggetSmallVeh �Ƚ������ ���ϵ�do/T��
% ����table��ʽ ���㴦�� - �������ȫ����Tȡֵ (T�����������1��n����)

fprintf(1,'\nRunning HBinCombine with output T...\n');

if ISlastVehType && ISpingpu
    T = HBinCombine(do,flaggetSmallVeh,do1,flagTiledArray,do2Array,do3Array);  
elseif ISlastVehType && ~ISpingpu
    T = HBinCombine(do,flaggetSmallVeh,do1);
elseif ~ISlastVehType && ISpingpu
    T = HBinCombine(do,[],[],flagTiledArray,do2Array,do3Array);
elseif ~ISlastVehType && ~ISpingpu
    T = HBinCombine(do);
end

%% return����Ĳ���
output_CoordLUBin=T.CoordLUBin';
output_LU_LWH=T.LWH';

% V2 get output_LU_Seq
% ��1���������ڳ��ͺ�(���뻻)    ��2���������ڳ����(���,���ܻ�,���ʹ�) ��3�����̳��ڰ���˳��(���뻻) ��4������SID��Ӧ�̱��(�����,���ñ�??)
% ��5������ID�ͺ�LID(�����,���ñ�?) ��6�����̶Ѷ����ITEM(���,���ܻ�,���ʹ�,��;?) ��7�������㲿�����PID(�����,���ñ�?) ������8: չʾ˳��(���뻻)

if verMilkRun == 1 
output_LU_Seq=T{:,{'LU_VehType','BINID','BINSEQ','OSID','LID','ITEMID','OPID','ShowSEQ','Weight','Index','OEID'}}'; %���ӷ�����10: LuIndex������ǿֻҪ�����־Ϳ���
else % ����ȥOPID OSID OEID �� ��ʵ��;Ҳ���� %��;�����
output_LU_Seq=T{:,{'LU_VehType','BINID','BINSEQ','OSID','LID','ITEMID','OPID','ShowSEQ','Weight','Index'}}'; %���ӷ�����10: LuIndex������ǿֻҪ�����־Ϳ���
end
% V1�� % output_LU_Seq=T{:,{'LU_VehType','BINID','BINSEQ','OSID','LID','ITEMID','OPID','ShowSEQ','Weight'}}'; % ITEMID���岻��    output_LU_Seq([2,3,5,8],:)

if ISplotBBA
    V = struct2table(structfun(@(x) x',do.Veh,'UniformOutput',false));%V = struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false));
    plotSolutionT(T,V,0,0,0,1,ISplotShowType,'����BIN'); 
end

%% �������
% T.Properties.VariableNames
% % �򵥲鿴ĳ��bin�ڵĲ����ӱ�
% tmpT = sortrows(T,{'ITEMID'})
% x=T(T.BINID==2&T.ID==1,{'ID','LID','PID','H','Weight','CoordLUBin','BINSEQ','ShowSEQ','ITEMID','ITEMSEQ'})
% sortrows(x,{'CoordLUBin'})
% �޳�չʾ˳�� % output_LU_Seq = output_LU_Seq(1:7,:);
% whos
fprintf(fid,'���гɹ�\n');
fclose(fid);  %�������

clearvars -except output*
fprintf(1,'\nSimulation done.\n');

% mcc -W 'java:BBA_Main,Class1,1.0' -T link:lib BBA_Main.m -d '.\new'
% mcc -W 'java:BBA_MR_Main,Class1,1.0' -T link:lib BBA_MR_Main.m -d '.\new'
% mlintrpt(fullfile(pwd,'fun'),'dir') �������
% mlintrpt(fullfile(pwd,'fun'),'dir') todo/fixme����
% updateContents('F:\Users\HHBS\Documents\555GitHub\ROS\fun') 

% printstruct(do,'sortfields',1,'PRINTCONTENTS',0);    printstruct(do.Veh);
% do = rmfield(do, {'Veh', 'LU'});
% pcode 'H*.m'.
% lu = table2struct(T,'ToScalar',true)
% lu = (structfun(@(x) x',lu,'UniformOutput',false));

% �������
% printstruct(d);
% printstruct(d.LU,'sortfields',1,'PRINTCONTENTS',0)
% TVEHIN = struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false));
% % TLUIN = struct2table(structfun(@(x) x',d.LU,'UniformOutput',false));
% TLUIN.Properties.VariableNames{'PID'} = 'OPID'; 
% TLUIN.Properties.VariableNames{'SID'} = 'OSID';
% s = table2struct(TLUIN,'ToScalar',true)
% t = struct2table(l,'AsArray',true)

                                                        %         [match, er1, er2] = comp_struct(ti,d.Strip,1);
                                                        %         [match, er1, er2] = comp_struct(tl,d.Bin,1);
                                                        %         list_struct(er1)
%%
% t = [d.LU.ID;d.LU.LWH]
% sortrows(t',[1,4],{'ascend','descend'})
% ȱ��field������
% % %     if ~isfield(d.LU, 'isRota'),  d.LU.isRota = ones(1,n); end %��IDһ��
% ��������
% if ~isfield(d, 'Par') || ~isfield(d.Par, 'H')
%     d.Par.maxHeavey = 100; end

end %END MAIN




%% ********************** �����ǻ������õĴ��� ��ʱ���� ****************



%% dd: ��Ҫ����Gap��Bin
% function [flagGap,T23LU] = getMixedGap1(dd)
% % Initilize output
% flagGap=0;
% T23LU = struct2table(structfun(@(x) x',dd.LU,'UniformOutput',false));
% 
% % find and sort 'fidx' (fidx: ��Gap��Strip��Index)
% fidx=find(dd.Strip.isGapBalance==0);  %fidx: ���ڻ�װ��϶��strip�����
% [~,ff]=sort(dd.Strip.Strip_Bin(2,fidx)); %
% fidx=fidx(ff); %����С������������Ľ���˳�����, ����bin����12��ǰ,11�ں�
% 
% for i=1:length(fidx)
%     flagGap = 0;
%     idxs = fidx(i);  % strip12 �ҵ���Ӧ��LU LU_Strip ��12�� ����LWH������CoordLUBin CoordLUStrip
%     
%     %% 1: pgLUinStrip: strip��Ӧ��pg �� pgLUinStrip: ��strip��lu��Ӧ��pg -> pgBlanksinStrip:�������pg
%     %pgLUinStrip: strip�ڵ�LU�Ķ����,�����ǲ�����Ļ��Ƕ����
%     pgLUinStrip = pgStripLU(idxs,T23LU);    % pgLUinStrip.NumRegions
%     
%     %pgStrip: ��strip�Ķ����,һ���Ǹ�����
%     stripWidth = dd.Veh.LWH(1, unique(T23LU.LU_VehType));     %strip���
%     stripHeight = dd.Strip.LW(2,idxs);                                       %strip�߶�
%     [x,y]=boundary(pgLUinStrip);
%     pgStrip = polyshape(pgRectangle(min(x),min(y),stripWidth,stripHeight));  if size(pgStrip.Vertices,1) ~= 4, error('strip���Ǿ���'); end
%     
%     %pgBlanksinStrip: strip�ڵ�ʣ��հ�����,ÿ����������Ǿ���
%     pgBlanksinStrip = regions(subtract(pgStrip,pgLUinStrip)); % ��ȡ��strip�ڵ�ʣ���������    
%     %pgBlanksinStrip����, ��Ŀǰ: ���С����ǰ�棩 todo: yС����ǰ��
%     pgBlanksinStrip = sortregions(pgBlanksinStrip,'area','ascend');
%     %                 plot(pgBlanksinStrip)
%     %                 hold on
%     %                 plot(pgLUinStrip)
%     
%     %% 2: idxLUArray���ҵ�idxs��һ��strip�ڵ�LU,����LU�Ĵ�С�ж��ܷ����
%     % Find nextStrip
%     nextStrip = dd.Strip.Strip_Bin(2,idxs) + 1;                                 %��strip��bin�ڵĽ���˳�� + 1������һ��strip�����
%     nextStripIdx = find(dd.Strip.Strip_Bin(2,:) == nextStrip);
%     
%     % Get and sort idxLUArray by area of LU in nextStrip  %����LU��˳�򣬴Ӵ��������С���
%     idxLUArray = find(T23LU.LU_Strip(:,1)==nextStripIdx);                   %�ڻ�װ��϶��һ��strip�ڰ�����LU��� 
%     
%     idxLUareaArray = T23LU.LWH(idxLUArray,1).*T23LU.LWH(idxLUArray,2);    % areas
%     [~, ord] = sort(idxLUareaArray,'descend');
%     idxLUArray=idxLUArray(ord);
%     
%     % TODO �ҳ��ܷŵ��µ�region������ĳ�ֹ��򣨺���������.
%     % ���LU���Է���pgBlandsinStrip
%     for j=1:length(idxLUArray)
%         % �ӵ�j��LU��LUӦ�����򣬴������ģ���ʼ, ���Է��뵽�����stripȥ
%         l = idxLUArray(j);
%         
%         for i = 1:length(pgBlanksinStrip)
%             % 2.1 get pgBlank
%             pgBlank= pgBlanksinStrip(i);
%             
%             % 2.2.1 get pgLU
%             % Find pgBlank's �������� x and y) 
%             [originx,originy] = boundary(pgBlank);
%             originLUx=min(originx);
%             originLUy=min(originy);
%             
%             % 2.2.2 get pgLU ����ת ����pgBlank's Region
%             w=T23LU.LWH(l,1);
%             h=T23LU.LWH(l,2);
%             pgLU = polyshape(pgRectangle(originLUx,originLUy,w,h));
%             
%             % 2.2.3 get pgLURota ��ת ����pgBlank's Region
%             if T23LU.isRota(l)
%                 wRota = T23LU.LWH(l,2);
%                 hRota = T23LU.LWH(l,1);
%                 pgLURota = polyshape(pgRectangle(originLUx,originLUy,wRota,hRota));
%             end
%             
%             % 2.3 if pgLU/pgLURota can be assigned into pgBlank
%             [xLU,yLU] = boundary(pgLU);
%             [xLURota,yLURota] = boundary(pgLURota);
%             xLU
%             yLU
%             %% �ж��Ƿ���Է��µ���Ҫ����
%             % 555: if all of pgLU's boundary points belong to pgBlank
%             flagLUx =  all(isinterior(pgBlank,xLU([1,4]),yLU([1,4])));
%             flagLU =  all(isinterior(pgBlank,xLU,yLU));
%             flagLURotax =  all(isinterior(pgBlank,xLURota([1,4]),yLURota([1,4])));
%             flagLURota =  all(isinterior(pgBlank,xLURota,yLURota));
%             
%             % 3 ��������滻,
%             if flagLU || flagLUx
% %                 [x,y]=boundary(pgBlank);
% %                 originLUx = min(xLU)
% %                 originLUy = min(yLU)
%                 %TODO : update ��l��LU��strip���/����˳��.CoordLUStrip�����
%                 T23LU.CoordLUBin(l,1) = originLUx;
%                 T23LU.CoordLUBin(l,2) = originLUy;
%                 T23LU.LWH(l,1)=w;
%                 T23LU.LWH(l,2)=h;
%                 flagGap=1;
%                 break
%             end
%             
%             if flagLURota || flagLURotax
% %                 [x,y]=boundary(pgBlank);
% %                 originLUx = min(xLU);
% %                 originLUy = min(yLU);
%                 %TODO : update ��l��LU��strip���/����˳��.CoordLUStrip�����
%                 T23LU.CoordLUBin(l,1) = originLUx;
%                 T23LU.CoordLUBin(l,2) = originLUy;
%                 T23LU.LWH(l,1)=h;
%                 T23LU.LWH(l,2)=w;
%                 T23LU.Rotaed(l)=~T23LU.Rotaed(l);
%                 flagGap=1;
%                 break
%             end
%             
%             %     area(pgLU)                         area(pgBlank)
%             
%             %                         figure
%             %                         plot(pgLU)
%             %                         hold on
%             %                         plot(pgBlank)
%             
%             %                     area(pg12)
%             %                     area(pg1)
%             
%         end
%         if flagGap,   break;    end
%     end
%     if flagGap,   break;    end
%     
% end
% end






% % function lastd = getdinLastVeh(tmpd)
% %     % tmpd�е�Bin��������, ����С�Ŀ�ʼ��
% %     tmpusedVehIdx = max(tmpd.LU.LU_Bin(1,:)); %tmpusedVehIdx: ���һ��Bin��indexֵ
% %     flagusedLUIdx = tmpd.LU.LU_Bin(1,:)==tmpusedVehIdx; % flagused: �ҳ����һ��Bin��Ӧ��LUindexֵ
% %     if isSameCol(tmpd.LU)
% %         % ��ȡ�����һ��Bin����������
% %         lastd.LU = structfun(@(x) x(:,flagusedLUIdx),tmpd.LU,'UniformOutput',false);  %��ȡ���һ�����ڵ�LU
% %         lastd.LU.LWH([1,2], lastd.LU.Rotaed ) = flipud(lastd.LU.LWH([1,2], lastd.LU.Rotaed)); %LU.LWH ����ת,��ָ�ԭ��
% %         lastd.LU = rmfield(lastd.LU,{'Rotaed','order','LU_Item','DOC','LU_Strip','LU_Bin','CoordLUBin','maxL','CoordLUStrip'}); 
% %         lastd.Par = tmpd.Par;
% %     else
% %         error('����ʹ��structfun');
% %     end
% % end

%% **** �㷨ָ��ѡ�����Ž� ****
% % function [daMax,parMax] = getbestsol(DaS,Par)
% % 
% %     % �������һ�����н�, ֱ�ӷ���;
% %     if size(DaS,2)==1    %����dA�ж��ʱ����,Ŀǰ��������,
% %         daMax = DaS(1); parMax = Par(1);
% %         return
% %     end
% %     
% % %��ȡ����ָ��Ͷ�Ӧ����
% % for r=1:length(DaS)
% %     resLoadingRateBin(r) = mean(DaS(r).Bin.loadingrate); %bin��װ���ʾ�ֵ��� Itemloadingrate ItemloadingrateLimit
% %     resLoadingRateStripLimit(r) = mean(DaS(r).Strip.loadingrateLimit); %strip��limitװ������� Itemloadingrate ItemloadingrateLimit
% %     resLoadingRateBinLimit(r) = mean(DaS(r).Bin.loadingrateLimit); %bin��limitװ������� Itemloadingrate ItemloadingrateLimit
% %     resLoadingRateStrip(r) = mean(DaS(r).Strip.loadingrate); %strip��װ������� Itemloadingrate ItemloadingrateLimit    
% % %     Par(r);
% % end
% % 
% % %% �㷨ѡ�����ŵĽ���û�
% % % maxresBin=max(resLoadingRateBinLimit(1,idxStrip)); %�ҳ�idxStrip�е����bin
% % % if ~all(ismember(idxBin,idxStrip)),   error('not all member of bin in strip'); end %�����п��ܳ��� 
% % %% 1 maxresBin�����泵����ƽ��װ����,��Ʒ����һ��,binԽ��,��ֵԽС,��Խ��,�����ֵ�Ǳ���
% % %% 2 maxresStrip������Strip��ƽ��װ����,��Ʒ����һ��,strip���һ��,�߶�Խ��,��ֵԽС,��Խ��,�����ֵ��һ��ʱ���루��Ϊ���������ÿ���
% % %% ����ֵ�ڲ�ͬbin�߶�ʱ������Ӱ�죬�Ҹ�ֵ��ʱ����Ϊ���������ܲ�����
% % %% 3 maxresStripLimit��������Strip��ƽ��װ����,��Ʒ����һ��,strip�ڲ����Խ��,��϶Խ��,ֵԽС,�����ֵ�����Ǳ���
% % %% ��ֵ��ʱ����Ϊ���������ܺã�strip�ڲ���϶С��������һ��ʱ���ţ����п�����ͬ���̲���һ��
% % %% 4 maxresBinLimit��������Bin��ƽ��װ����,��Ʒ����һ��?? �������������,���۲�
% % idxBin=find(resLoadingRateBin==max(resLoadingRateBin)); %ȡ
% % idxStripLimit=find(resLoadingRateStripLimit==max(resLoadingRateStripLimit));
% % idxBinLimit=find(resLoadingRateBinLimit==max(resLoadingRateBinLimit));
% % idxStrip=find(resLoadingRateStrip==max(resLoadingRateStrip));
% % %% 5 �ҳ�idxStrip��idxBin���ߵĽ���
% % % % if isempty(intersect(idxBin,idxStrip))
% % idx =idxBin;
% % if isempty(idx), error('idxBinΪ�� '); end %���󼸺������ܳ���
% % idx0 =intersect(idx,idxStripLimit);
% % if ~isempty(idx0),
% %     idx = idx0; 
% % else
% %     warning('idx0 is empty');
% % end
% % % if isempty(idx), error('idxBin and idxStripLimit �Ľ���Ϊ�� '); end %���󼸺������ܳ���
% % idx1 = intersect(idx,idxBinLimit);
% % if ~isempty(idx1),  
% % %      idx = idx1; 
% % else
% %     warning('idx1 is empty');
% % end
% % idx2 = intersect(idx,idxStrip);
% % if ~isempty(idx2),  
% % %     idx = idx2; 
% % else
% %     warning('idx2 is empty');
% % end
% % 
% % %% ��idxʣ��ķ��ص�������
% % if ~isempty(idx)
% %     for tmpidx=1:length(idx)
% %         daMax(tmpidx) = DaS(idx(tmpidx));
% %         parMax(tmpidx) = Par(idx(tmpidx));
% %     end
% % end
% % 
% % end % END OF ALL
% % 
% % 
% % 
% % 
% % 
% %         
% % function plotSolution(do,par)
% % %% ��ͼ
% % % V3 margin ��ǰ��RunAlgorithm���к��ִ��:
% % % plot2DBPP(do,par);
% % plot3DBPP(do,par);
% % 
% %         % V1 buff version
% %         % do.Item.LWH = do.Item.LWH - do.LU.buff(:,1:size(do.Item.LWH,2));
% %         % do.Item.LWH(1,:) = do.Item.LWH(1,:) - ( do.LU.margin(1, 1:size(do.Item.LWH,2) ) + do.LU.margin(2,: )); 
% %         % do.Item.LWH(2,:) = do.Item.LWH(2,:) - (do.LU.margin(3,: ) + do.LU.margin(4,: )); 
% %         % do.Item.CoordItemBin = do.Item.CoordItemBin + do.LU.buff(:,1:size(do.Item.LWH,2))/2;
% % 
% %         % V2 margin version
% %         % ��ͼǰ����LU ITEM��Coord��LW; ����ITEMͬʱ����LU
% %         % [do.LU,do.Item] = setLCwithoutbuff(do.LU,do.Item);
% %         
% % 
% % end


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
% %         tepMin = min(tepMin(findBinArray)); % 555 chk
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
