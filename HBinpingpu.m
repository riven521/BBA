function [flagTiledArray,do2Array,do3Array] = HBinpingpu(maind,do,p)
global ISpingpuAll ISplotEachPingPu

        % 3.1 ��ʼ��5������    
        flagTiledArray = zeros(1,length(do.Bin.Weight));  %1��������ƽ�� 2����˦βƽ��

        d2Array(1:length(do.Bin.Weight)) = maind; % fixme �޸�maind��ʼ�� ����������
        do2Array(1:length(do.Bin.Weight)) = structfun(@(x) [], do, 'UniformOutput', false);
        d3Array(1:length(do.Bin.Weight)) = maind;
        do3Array(1:length(do.Bin.Weight)) = structfun(@(x) [], do, 'UniformOutput', false); %do;
        
        
        
        %% 3.2 ����binѭ��ƽ��
        bidx = 1:length(do.Bin.Weight); % NOTE: �޸�ֻ����˦βƽ�̵�ȫ��Bin���뿼��, �Է�˦βƽ�̵Ľ���ƽ���ж� % bidx = find(do.Bin.isTileNeed);
        % ѭ��: ÿ��bin�ֱ���ƽ��(����ƽ�̺�˦βƽ��ֻ�ܶ�ѡһ����������ƽ��)
        for i=1:numel(bidx)
            ibin = bidx(i);
            
            % $1 GET d2/d3 ��ibin�ڴ��㷨�������������(������maind)
            luIdx = do.LU.LU_Bin(1,:) == ibin;
            d2 = getPartdinThisBin(maind,luIdx); %�޸ĳɴ�maind��ȡIuIdx������,���Ǵ�������d����ȡ %         d2 = getdinThisVeh(do,luIdx)
            d3 = d2;
            
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
                if ibin == numel(bidx) && range(d3.LU.ID) == 0 % ���һ�� �� һ����ȫһ����ID�ķ���1-2-3�ĵ���ƽ��
                    minmaxLayer = min(max(d3.LU.maxHLayer), max(d3.LU.maxL(3,:))); %maxHLayer��ָ������; maxL(3,:)������������;
                else
                    minmaxLayer = 1;
                end
                iLayer = 1;
                while 1
                    if iLayer > minmaxLayer, break; end
                    d3.LU.maxHLayer(:) = iLayer; %d2��ȫ��LU�Ĳ����趨Ϊ1 55555 ȫ��ƽ�̵���Ҫ����
                    iLayer=iLayer+1;
                    
                    % $3.3.1 reRunAlgorithm do3��d3�����Ľ��
                    d3Array(ibin) = d3;
                    [d3.LU] = setLULWHwithbuff(d3.LU, d3.Veh);
                    do3 = RunAlgorithm(d3,p);   % do3Array(ibin) = do3;
                    
                    % $3.3.2 ��ȫ��ƽ��û������,������
                    if max(do3.LU.LU_Bin(1,:)) == 1
                        flagTiledArray(ibin)=1; %1��������ƽ��
                        %do3.LU.LU_VehType = ones(size(d3.LU.ID)) * do3.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ
                        [do3.LU,do3.Item] = setLCwithoutbuff(do3.LU,do3.Item);               %  plot3DBPP(do3,p)
                        do3Array(ibin) = do3;
                        %                 plotSolutionT(do3.LU,do3.Veh);
                        % do3 �޸ĵ�d�У�����Ŀǰ������do2Array�У�δ��d�ϲ�
                        break;    %  continue;   %�����������˦βƽ���� ����while����continue��
                    end
                end
            end
            if flagTiledArray(ibin)==1
                continue;
            end
            %% 3.4 ������ƽ��ʧ�� ����˦βƽ�� ��Ϊ����
            
            % 3.4.1 GET do2 ��ibin�ں�������ݣ�˦βƽ��ʹ��
            stripidx = do.Strip.Strip_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip
            itemidx = do.Item.Item_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip
            do2.Veh = do.Veh;
            do2.LU = structfun(@(x) x(:,luIdx),do.LU,'UniformOutput',false);
            do2.Bin = structfun(@(x) x(:,ibin),do.Bin,'UniformOutput',false);
            do2.Strip = structfun(@(x) x(:,stripidx),do.Strip,'UniformOutput',false);
            do2.Item = structfun(@(x) x(:,itemidx),do.Item,'UniformOutput',false);
            
            %% ��ȡluidxPP �������� d2.LU.maxHLayer(luidxPP)
            % 555 plotSolutionT(do2.LU,do2.Veh);  % ˦βƽ��ǰ�� �۲�
            while do2.Bin.isTileNeed(1) == 1 %do2�ڵ�Bin��Զֻ��1��, ����ƽ�̺��bin����Ҫƽ��,������while�ж�
                % $3.4.2 �޶�d2.LU.maxHLayer (����ibin�����ѡ���ļ���stripƽ��) TODO $4д����Щ����,���ڼ�
                % $4.1 GET luidxPP �� ĳ��strip��Ӧ��LU�߼�ֵ
                % ѭ���ӱ�ibin�����һ��strip��ʼƽ�� istrip= nbStrip;
                nbStrip = numel(do2.Strip.Weight);
                if unique(do2.Strip.Strip_Bin(2, :)) ~= nbStrip,    error('��Ԥ�ڴ���');    end
                istrip= nbStrip;
                fi = find(do2.Strip.Strip_Bin(2,:) >= istrip ); % ͬbin��strip��� and ˳��>=istrip
                u=unique(do2.LU.LU_Strip(1,:)); %��ȡStrip��ŵ�Ψһ����ֵ
                luidxPP = ismember(do2.LU.LU_Strip(1,:), u(fi)); %%% fi->u(fi) ��������� *********************
                if ~any(luidxPP),  error('luidxPPȫ��Ϊ��, ������u(fi)��Ӧ��Lu�߼��ж�'); end
                
                % $4.2 �޶�d2.LU.maxHLayer           d2.LU    maind.LU
                d2.LU.maxHLayer(luidxPP) = min( d2.LU.maxL(3,luidxPP), d2.LU.maxHLayer(luidxPP)) - 1;
                
                % $4.2 ����ǰluidxPP��ӦLu�Ĳ������Ѿ�Ϊ1��, ����Ҫ���Ӹ����istrip��luidxPP; ���޶�d2.LU.maxHLayer
                % GET ���� d2.LU.maxHLayer(luidxPP) ����luidxPP�Ĳ���>=2��
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
                
                %% �������㷨������ $5 reRunAlgorithm
                d2Array(ibin) = d2;                
                [d2.LU] = setLULWHwithbuff(d2.LU, d2.Veh);
                
                do2 = RunAlgorithm(d2,p);        %  do2 = RunAlgorithmPP(d2,p);  %do2.LU.LU_VehType = ones(size(d2.LU.ID)) * do2.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ                
                
                [do2.LU,do2.Item] = setLCwithoutbuff(do2.LU,do2.Item);  % do2Array(ibin) = do2; ����ע�ͣ���Ϊ�Ǹ�ѭ��
                
                % plot
                if ISplotEachPingPu == 1,     plotSolution(do2,p);       end
                
                % $6 ����
                if max(do2.LU.LU_Bin(1,:)) == 1  %    do2.LU.LU_VehType = ones(size(do2.LU.ID))*do.Veh.order(1);
                    
                    flagTiledArray(ibin)=2; %2����˦βƽ��
                    do2Array(ibin) = do2;
                    %                  plotSolutionT(do2.LU,do2.Veh);
                    %                  pause(0.2)
                    
                    % todo -> do2 ���ݲ�����d ����return2bba���޸�    % do2 ���ݽ���d???? return2bba���޸ģ�����
                else
                    break;  %����˦β�Ų��� ���ټ���˦βƽ���ˣ����˽�����������һ����˦β�ж�
                end
                
            end % END OF WHILE
        end% END OF FOR
        
        %% chk
        if any(flagTiledArray)
            for ibin=1:length(flagTiledArray)
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
                
                flagTiledArray;
                d2Array;
                do2Array;
                d3Array;
                do3Array;
            end
        end

end


