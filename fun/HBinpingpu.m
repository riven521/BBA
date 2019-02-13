function [flagTiledArray,do2Array,do3Array] = HBinpingpu(maind,do,p)
% HBinpingpu ==> ����������maind�����������do������p ��ȡ����ƽ��do3Array/˦βƽ��do2Array���������
%   ������
%   d2inIbin=d3inIbin��ibin�ڵ��������ݣ�����maind
%   do2inIbin=do3inIbin��ibin�ڵ�������ݣ�����do  
        
        global ISpingpuShuaiWei ISpingpuAll ISshuaiwei
        global ISplotEachPingPuAll ISplotEachPingPuShuaiWei 
        

        % 3.1 ��ʼ��5���������
        nbin = length(do.Bin.Weight);

        flagTiledArray = zeros(1,nbin);  %1��������ƽ�� 2����˦βƽ��
        d2Array(1:nbin) = maind;    % fixme �޸�maind��ʼ�� ����������  ˦βƽ��
        do2Array(1:nbin) = structfun(@(x) [], do, 'UniformOutput', false);
        d3Array(1:nbin) = maind;    %����ƽ��
        do3Array(1:nbin) = structfun(@(x) [], do, 'UniformOutput', false); %do;

        %% 3.2 ����binѭ��ƽ��
        bidx = 1:nbin; % NOTE: �޸�ֻ����˦βƽ�̵�ȫ��Bin���뿼��, �Է�˦βƽ�̵Ľ���ƽ���ж� % bidx = find(do.Bin.isTileNeed);
        
        % ѭ��: ÿ��bin�ֱ���ƽ��(����ƽ�̺�˦βƽ��ֻ�ܶ�ѡһ����������ƽ��)
        for i=1:numel(bidx)
            ibin = bidx(i);
            
            % 1 GET d2inIbin/do2inIbin(d3inIbin/do3inIbin)
            d2inIbin = getPartDinThisBin(do, ibin, maind);  % V2  % V1 %    luIdx = do.LU.LU_Bin(1,:) == ibin;   %      d2inIbin = getPartDinThisBin(maind,luIdx,do); %maind��ȡ����margin %�޸ĳɴ�maind��ȡIuIdx������,���Ǵ�������do����ȡ %         d2inIbin = getdinThisVeh(do,luIdx)
            do2inIbin = getPartDinThisBin(do,ibin);             % v2  % v1  %    luIdx = do.LU.LU_Bin(1,:) == ibin; %                 stripidx = do.Strip.Strip_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip %                 itemidx = do.Item.Item_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip %                 do2inIbin.Veh = do.Veh; %                 do2inIbin.Bin = structfun(@(x) x(:,ibin),do.Bin,'UniformOutput',false); %                 do2inIbin.LU = structfun(@(x) x(:,luIdx),do.LU,'UniformOutput',false);                %                 do2inIbin.Strip = structfun(@(x) x(:,stripidx),do.Strip,'UniformOutput',false); %                 do2inIbin.Item = structfun(@(x) x(:,itemidx),do.Item,'UniformOutput',false);

            d3inIbin = d2inIbin;
            do3inIbin = do2inIbin;
            
            %% COMMENT
            %         d2inIbin.Veh = do.Veh;
            %         d2inIbin.Veh = rmfield(d2inIbin.Veh,{'Volume','order'});
            % 2 �������/һ��strip�ڵ�LU
            %         luidx = do.LU.LU_Bin(1,:) == ibin;  %do.LU.LU_Strip(1,:) == istrip
            
            %         d2inIbin.LU = structfun(@(x) x(:,luidx),do.LU,'UniformOutput',false);
            %         d2inIbin.LU.LWH([1,2], d2inIbin.LU.Rotaed ) = flipud(d2inIbin.LU.LWH([1,2], d2inIbin.LU.Rotaed)); %LU.LWH ����ת,��ָ�ԭ��
            %         d2inIbin.LU.PID = d2inIbin.LU.OPID;     d2inIbin.LU.SID = d2inIbin.LU.OSID;  %  d2inIbin.LU.LID = d2inIbin.LU.OLID;
            %         d2inIbin.LU = rmfield(d2inIbin.LU,{'Rotaed','order','LU_Item','DOC','LU_Strip',...
            %             'LU_Bin','CoordLUBin','CoordLUStrip','LU_VehType','OPID','OSID'});
            %
            %         d2inIbin.Par = do.Par;
            %% 2 �������ȫ��ƽ��(�����Ƿ�˦βƽ��), �۲챾ibin���Ƿ����ȫ��ƽ��,�����,��ȡ��˦βƽ��; ����,����˦βƽ��
            if ISpingpuAll==1
                
                % 2.1 ȷ��minmaxLayer��
                % 1 ����case��������һ�� �� һ����ȫһ����ID�ķ���1-2-3�ĵ���ƽ��
                if ibin == numel(bidx) && range(d3inIbin.LU.ID) == 0 
                    minmaxLayer = min(max(d3inIbin.LU.maxHLayer), max(d3inIbin.LU.maxL(3,:))); %maxHLayer��ָ������; maxL(3,:)������������;
                else % 2һ��case������ƽ�̼�����һ��
                    minmaxLayer = 1;
                end
                
                iLayer = 1;
                
                while 1
                    if iLayer > minmaxLayer, break; end
                    
                    % update d3inIbin -> ��ֵ��d3array
                    d3inIbin.LU.maxHLayer(:) = iLayer;   %d3��ȫ��LU�Ĳ����趨Ϊ1  55555 ȫ��ƽ�̵���Ҫ����                    
                    d3Array(ibin) = d3inIbin;
                    
                    iLayer=iLayer+1;
                    
                    % 2.2  reRunAlgorithm do3��d3�����Ľ��
                    
                    d3inIbin.LU = setLULWHwithbuff(d3inIbin.LU, d3inIbin.Veh);
                    do3 = RunAlgorithm(d3inIbin,p);   % do3Array(ibin) = do3; (����ƽ�̳ɹ��ķ���do3Array�������ķŲ�������ν���ò���
                    [do3.LU,do3.Item] = setLCwithoutbuff(do3.LU,do3.Item);               %  margin�滻����
                    
                    % 2.3   ��ȫ��ƽ��û������,������
                    if max(do3.LU.LU_Bin(1,:)) == 1
                        fprintf(1,'       Exsiting ����ƽ�� in HBinpingpu (do3)...\n');
                        flagTiledArray(ibin)=1;   %1��������ƽ��                       %do3.LU.LU_VehType = ones(size(d3inIbin.LU.ID)) * do3.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ                        
                        
                        do3Array(ibin) = do3;
                                   
                        if ISplotEachPingPuAll == 1  % plot ˦βƽ��ǰ��˦βƽ�̺��binͼ
                            plotSolutionT(do3inIbin.LU,do3inIbin.Veh, 0, 0);
                            plotSolutionT(do3.LU,do3.Veh, 0, 0);   % plot����ƽ�̺��bin % do3 �޸ĵ�d�У�����Ŀǰ������do2Array�У�δ��d�ϲ�
                        end
                       
                        break;    % ����ƽ�̳ɹ�������while��ѭ����һ��bin
                    end
                    
                end % END WHILE
                
            end % END ISpingpuAll
            
            %% �ж���������ƽ�̳ɹ���������һ��bin��������˦βƽ��
            if flagTiledArray(ibin)==1
                continue;
            end
            
            %% 3 ������ƽ��ʧ�� ����˦βƽ�̣�ǰ�᣺��˦β������ ��Ϊ����
            if ISpingpuShuaiWei==1 && ISshuaiwei==1
                             
                %% ��ȡluidxPP �������� d2inIbin.LU.maxHLayer(luidxPP)
                %            plotSolutionT(do2inIbin.LU,do2inIbin.Veh);  % ˦βƽ��ǰ�� �۲�
                while do2inIbin.Bin.isTileNeed(1) == 1 %do2inIbin����ƽ�̺��bin����Ҫƽ��,������while�ж�
                    % $3.4.2 �޶�d2.LU.maxHLayer (����ibin�����ѡ���ļ���stripƽ��) TODO $4д����Щ����,���ڼ�
                    % $4.1 GET luidxPP �� ĳ��strip��Ӧ��LU�߼�ֵ
                    % ѭ���ӱ�ibin�����һ��strip��ʼƽ�� istrip= nbStrip;
                    nbStrip = numel(do2inIbin.Strip.Weight);  
                        if unique(do2inIbin.Strip.Strip_Bin(2, :)) ~= nbStrip,    error('��Ԥ�ڴ���');    end
                    istrip= nbStrip;
                    fi = find(do2inIbin.Strip.Strip_Bin(2,:) >= istrip ); % ͬbin��strip��� and ˳��>=istrip
                    u=unique(do2inIbin.LU.LU_Strip(1,:)); %��ȡStrip��ŵ�Ψһ����ֵ
                    luidxPP = ismember(do2inIbin.LU.LU_Strip(1,:), u(fi)); %%% fi->u(fi) ��������� *********************
                    if ~any(luidxPP),  error('luidxPPȫ��Ϊ��, ������u(fi)��Ӧ��Lu�߼��ж�'); end
                    
                    % $4.2 �޶�d2.LU.maxHLayer           d2inIbin.LU    maind.LU
                    d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP)) - 1;
                    
                    % $4.2 ����ǰluidxPP��ӦLu�Ĳ������Ѿ�Ϊ1��, ����Ҫ���Ӹ����istrip��luidxPP; ���޶�d2.LU.maxHLayer
                    % GET ���� d2inIbin.LU.maxHLayer(luidxPP) ����luidxPP�Ĳ���>=2��
                    while all(d2inIbin.LU.maxHLayer(luidxPP)<1)
                        istrip = istrip-1;
                        if istrip==0,break;end
                        fi = find( do2inIbin.Strip.Strip_Bin(2,:) >= istrip ); %fi = find( do2inIbin.Strip.Strip_Bin(2,:) == istrip );
                        luidxPP = ismember(do2inIbin.LU.LU_Strip(1,:), u(fi)); %%% fi->u(fi) ��������� *********************
                        if ~any(luidxPP),  error('luidxPPȫ��Ϊ��, ������u(fi)��Ӧ��Lu�߼��ж�'); end
                        if istrip == 0,  error('��bin������tileneed,��Ԥ�ڴ���');   end
                        d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP)) - 1;
                    end
                    % �޸�: ������Ļָ�Ϊ1
                    d2inIbin.LU.maxHLayer(d2inIbin.LU.maxHLayer<=1) = 1;
                    
                    %% �������㷨������ $5 reRunAlgorithm
                    d2Array(ibin) = d2inIbin;
                    
                    [d2inIbin.LU] = setLULWHwithbuff(d2inIbin.LU, d2inIbin.Veh);                    
                    do2 = RunAlgorithm(d2inIbin,p);        %  do2 = RunAlgorithmPP(d2inIbin,p);  %do2.LU.LU_VehType = ones(size(d2inIbin.LU.ID)) * do2.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ                   
                    [do2.LU,do2.Item] = setLCwithoutbuff(do2.LU,do2.Item);  % do2Array(ibin) = do2; ����ע�ͣ���Ϊ�Ǹ�ѭ��
                    
                    % plot
                    
                    
                    % $6 ���� ˦βƽ��һ�β�������Ҫѭ���ٽ���˦βƽ�� ֱ������һ������
                    if max(do2.LU.LU_Bin(1,:)) == 1  %    do2.LU.LU_VehType = ones(size(do2.LU.ID))*do.Veh.order(1);
                        fprintf(1,'       Exsiting ˦βƽ�� in HBinpingpu (do2)...\n');
                        flagTiledArray(ibin)=2; %2����˦βƽ��
                        
                        do2Array(ibin) = do2;                        % todo -> do2 ���ݲ�����d ����return2bba���޸�    % do2 ���ݽ���d???? return2bba���޸ģ�����
                        
                        if ISplotEachPingPuShuaiWei == 1  % plot ˦βƽ��ǰ��˦βƽ�̺��binͼ
                            plotSolutionT(do2inIbin.LU,do2inIbin.Veh, 0, 0);
                            plotSolutionT(do2.LU,do2.Veh, 0, 0);  % plot˦βƽ�̺��bin
                        end
                        
                    else
                        break;  %����˦β�Ų��� ���ټ���˦βƽ���ˣ����˽�����������һ����˦β�ж�
                    end                    
                
                end % END OF do2 WHILE ����˦β
                
            end % END OF ISpingpuShuaiwei
            
        end% END OF FOR
        
        %% chk
        if any(flagTiledArray)
            for ibin=1:nbin
                if flagTiledArray(ibin)==1 %����ƽ��
                    %                 t3 = struct2table(structfun(@(x) x',d3Array(ibin).LU,'UniformOutput',false));
                    %                 to3 = struct2table(structfun(@(x) x',do3Array(ibin).LU,'UniformOutput',false));
                    chkLUnewold(d3Array(ibin).LU,do3Array(ibin).LU);
                    chktLU(d3Array(ibin).LU);
                    chktLU(do3Array(ibin).LU);
                    
                end
                if flagTiledArray(ibin)==2 %˦βƽ��
                    %                 t2 = struct2table(structfun(@(x) x',d2Array(ibin).LU,'UniformOutput',false));
                    %                 to2 = struct2table(structfun(@(x) x',do2Array(ibin).LU,'UniformOutput',false));
                    chkLUnewold(d2Array(ibin).LU,do2Array(ibin).LU);
                    chktLU(d2Array(ibin).LU);
                    chktLU(do2Array(ibin).LU);
                end
                
%                 flagTiledArray;
%                 d2Array;
%                 do2Array;
%                 d3Array;
%                 do3Array;
            end
        end

end


