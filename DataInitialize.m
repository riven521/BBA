%% V2 - ����Random����
function d = DataInitialize( varargin )
% DATAITIALIZE Initialize the DATA data structure using parameter pairs.
% Inputs
%     varargin:  ('parameter',value,...)
% Outputs
%     d (.) structure : d.LU d.Veh

% ���ذ�˳��, Ҳ����ȫ��������
for k = 1 : 2 : length(varargin)
    switch varargin{k}
        case 'LUID' %LU�Ѷ��ж�
            d.LU.ID                        = varargin{k+1};
        case 'LULID' %LU�����ж�, չʾ��
            d.LU.LID                      = varargin{k+1};
        case 'LULWH'
            d.LU.LWH                    = varargin{k+1};
        case 'LUMARGIN'
            d.LU.margin                 = varargin{k+1};
        case 'LUWEIGHT'
            d.LU.Weight                = varargin{k+1};
        case 'LUISROTA'
            d.LU.isRota                = varargin{k+1};
        case 'LUSID'
            d.LU.SID                     = varargin{k+1};
        case 'LUPID'
            d.LU.PID                     = varargin{k+1};
        case 'LUINDEX'              %LUINDEX ��ǿר�÷���ϵ��
            d.LU.Index                  = varargin{k+1};
        case 'LUEID'                    %LUEID MR��LU��EID����
            d.LU.EID                     = varargin{k+1};
        case 'VEHID'
            d.Veh.ID                    = varargin{k+1};
        case 'VEHLWH'
            d.Veh.LWH                  = varargin{k+1};
        case 'VEHWEIGHT'
            d.Veh.Weight              = varargin{k+1};
        otherwise
            error('����δ֪����');
    end
end

% %     % ��������
    if ~isfield(d, 'Par') || ~isfield(d.Par, 'H')
        d.Par.maxHeavey = 100; end % judge isHeavey

end















%% V1 - �汾1 �����������
% % function d = DataInitialize( varargin )
% % % DATAITIALIZE Initialize the DATA data structure using parameter pairs.
% % % Inputs
% % %     varargin:  ('parameter',value,...)
% % % Outputs
% % %     d (.) DATA data structure : LU Strip Bin
% %
% % % 1 �������
% % if isnumeric(varargin{1})
% %
% %     if varargin{1}~=0
% %         % �������������
% %         d=getRandDa(varargin{1},varargin{2});         %   save('rndDa.mat','d');   load('rndDa.mat');
% %     else
% %
% % % % %         % TODO ���������, ����Defaultsֵ
% % % % %         d.LU.ID = [111 111 222 111]; %LU ID
% % % % %         d.LU.isRota =[0 0 1 0]; %��IDһ��
% % % % %         d.LU.maxL = [3 3 2 3]; % maximum layer in any vehicle  ��IDһ��
% % % % %         d.LU.yID = [0 0 1 0]; % palce in any type of vehicle ��IDһ��
% % % % %         % d.LU.xID = zeros(size(d.LU.ID)); % ��IDһ��
% % % % %         d.LU.margin = zeros(4,length(d.LU.ID)); %L R F B (4,n) ��IDһ��
% % % % %
% % % % %         d.LU.PID = [100 100 100 200];  %part ID
% % % % %         d.LU.SID = [2 3 3 2];   %supplier ID
% % % % %
% % % % %         d.LU.Weight = [10 10 10 10];
% % % % %         d.LU.isH = zeros(size(d.LU.ID)); %1 isHeavy
% % % % %
% % % % %
% % % % %         d.LU.LWH = [2 2 3 2; 5 5 6 5; 1 1 1 1];
% % % % %         d.LU.buff =  zeros(size(d.LU.LWH)); %�Ժ�����
% % % % %
% % % % %         % ��: ��������; ��: ���岻ͬ
% % % % % %         d.Veh.ID = [1, 2, 3, 4];
% % % % % %         d.Veh.Weight = [1000, 800, 500,300];
% % % % % %
% % % % % %         d.Veh.LWH = [5 4 3 3; 20 15 10 10; 4 4 2 2];
% % % % % %         d.Veh.buff = zeros(size(d.Veh.LWH)); %�Ժ�����
% % % % % %
% % % % % %         d.Veh.yID = [8 10 10 10; 20 15 0 0];   %�Ծ���ֵΪ׼  0��ʾ����   % ?*m
% % % % % %         % d.Veh.xID = d.Veh.LWH(1,:);   %
% % % % %
% % % % %         d.Veh.ID = [1];
% % % % %         d.Veh.Weight = [1000];
% % % % %
% % % % %         d.Veh.LWH = [5;20;4];
% % % % %         d.Veh.buff = zeros(size(d.Veh.LWH)); %�Ժ�����
% % % % %
% % % % %         d.Veh.yID = [20];   %�Ծ���ֵΪ׼  0��ʾ����   % ?*m
% % % % %         % d.Veh.xID = d.Veh.LWH(1,:);   %
% % % % %
% % % % %         % ��������
% % % % %         d.Par.H = 100; % judge isHeavey
% %
% %     end
% % % 2 test����
% % else
% %     % ���ذ�˳��, Ҳ����ȫ��������
% %     for k = 1 : 2 : length(varargin)
% %         switch varargin{k}
% %             case 'LUID' %LU�Ѷ��ж�
% %                 d.LU.ID                        = varargin{k+1};
% %             case 'LULID' %LU�����ж�, չʾ��
% %                 d.LU.LID                      = varargin{k+1};
% %             case 'LULWH'
% %                 d.LU.LWH                    = varargin{k+1};
% %             case 'LUMARGIN'
% %                 d.LU.margin                 = varargin{k+1};
% %             case 'LUWEIGHT'
% %                 d.LU.Weight                = varargin{k+1};
% %             case 'LUISROTA'
% %                 d.LU.isRota                = varargin{k+1};
% %             case 'LUSID'
% %                 d.LU.SID                     = varargin{k+1};
% %             case 'LUPID'
% %                 d.LU.PID                     = varargin{k+1};
% %             case 'LUINDEX'      %LUINDEX ��ǿר�÷���ϵ��
% %                 d.LU.Index              = varargin{k+1};
% %             case 'LUEID'            %LUEID MR��LU��EID����
% %                 d.LU.EID              = varargin{k+1};
% %             case 'VEHID'
% %                 d.Veh.ID                    = varargin{k+1};
% %             case 'VEHLWH'
% %                 d.Veh.LWH                  = varargin{k+1};
% %             case 'VEHWEIGHT'
% %                 d.Veh.Weight              = varargin{k+1};
% %                         %             case 'VEHBUFF'
% %                         %                 d.Veh.buff                 = varargin{k+1};
% %             otherwise
% %                 error('����δ֪����');
% %         end
% %     end
% %     printstruct(d.LU,'sortfields',1,'PRINTCONTENTS',0)
% %
% %
% %
% %     % ������initCheck
% %     [d.LU, d.Veh] = initCheck(d.LU,d.Veh);
% %
% % % % %     % TODO �Ƿ���Ҫ? ������������
% % % % %     n = length(d.LU.ID);
% % % % %     m = length(d.Veh.ID);
% % % % %
% % % % %     if ~isfield(d.LU, 'isRota'),  d.LU.isRota = ones(1,n); end %��IDһ��
% % % % %     if ~isfield(d.LU, 'maxL'),     d.LU.maxL = ones(3,n); end% maximum layer in any vehicle  ��IDһ��
% % % % %     if ~isfield(d.LU, 'yID'),       d.LU.yID = zeros(1,n);  end% palce in any type of vehicle ��IDһ��
% % % % %     if ~isfield(d.LU, 'xID'),     % d.LU.xID = zeros(size(d.LU.ID)); % ��IDһ��
% % % % %     if ~isfield(d.LU, 'margin'),     d.LU.margin = zeros(4,n); end %L R F B (4,n) ��IDһ��
% % % % %
% % % % %     if ~isfield(d.LU, 'PID'),     d.LU.PID = ones(1,n);   end%part ID
% % % % %     if ~isfield(d.LU, 'SID'),     d.LU.SID = ones(1,n);   end%supplier ID
% % % % %     if ~isfield(d.LU, 'UID'),     d.LU.UID = ones(1,n);    end%unloading ID
% % % % %
% % % % %     if ~isfield(d.LU, 'Weight'),     d.LU.Weight = 10*ones(1,n);  end
% % % % %     if ~isfield(d.LU, 'isH'),           d.LU.isH = zeros(1,n);         end        %1 isHeavy
% % % % %
% % % % %     if ~isfield(d.LU, 'buff'),           d.LU.buff = zeros(size(d.LU.LWH));  end%�Ժ�����
% % % % %
% % % % %
% % % % %     % ��: ��������; ��: ���岻ͬ
% % % % %     if ~isfield(d.Veh, 'Weight'),           d.Veh.Weight = 1000*ones(1,m);  end
% % % % %     if ~isfield(d.Veh, 'buff'),                d.Veh.buff = zeros(size(d.Veh.LWH));  end     %�Ժ�����
% % % % %     if ~isfield(d.Veh, 'yID')
% % % % %         for i=1:m
% % % % %             d.Veh.yID(1,i) = d.Veh.LWH(2,i);   %�Ծ���ֵΪ׼  0��ʾ����   % ?*m
% % % % %         end
% % % % %     end
% % % % %     % d.Veh.xID = d.Veh.LWH(1,:);
% % % % %
% %
% %     % ��������
% %     if ~isfield(d, 'Par') || ~isfield(d.Par, 'H')
% %         d.Par.maxHeavey = 100; end % judge isHeavey
% %
% % end
% % end