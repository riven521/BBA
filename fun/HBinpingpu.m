%% V2 
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
            doinIbin = getPartDinThisBin(do,ibin);             % v2  % v1  %    luIdx = do.LU.LU_Bin(1,:) == ibin; %                 stripidx = do.Strip.Strip_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip %                 itemidx = do.Item.Item_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip %                 do2inIbin.Veh = do.Veh; %                 do2inIbin.Bin = structfun(@(x) x(:,ibin),do.Bin,'UniformOutput',false); %                 do2inIbin.LU = structfun(@(x) x(:,luIdx),do.LU,'UniformOutput',false);                %                 do2inIbin.Strip = structfun(@(x) x(:,stripidx),do.Strip,'UniformOutput',false); %                 do2inIbin.Item = structfun(@(x) x(:,itemidx),do.Item,'UniformOutput',false);

            d3inIbin = d2inIbin;  
            
            %% 2 �������ȫ��ƽ��-������(�����Ƿ�˦βƽ��), �۲챾ibin���Ƿ����ȫ��ƽ��,�����,��ȡ��˦βƽ��; ����,����˦βƽ��
            if ISpingpuAll==1
                
                % 2.1 ȷ��minmaxLayer��
                % 1 ����case��������һ�� �� һ����ȫһ����ID�ķ���1-2-3�ĵ���ƽ��
                if ibin == numel(bidx) && range(d3inIbin.LU.ID) == 0 
                    % minmaxLayer = min(max(d3inIbin.LU.maxHLayer), max(d3inIbin.LU.maxL(3,:))); %maxHLayer��ָ������; maxL(3,:)������������;
                    minmaxLayer = min(max(d3inIbin.LU.maxHLayer), max(d3inIbin.LU.maxL(3,:)));  % fixme ��ȡ����ɾ�� maxL���ж�, �����Ǹ���
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
                    do3inIbin = RunAlgorithm(d3inIbin,p);   % do3Array(ibin) = do3; (����ƽ�̳ɹ��ķ���do3Array�������ķŲ�������ν���ò���
                    [do3inIbin.LU,do3inIbin.Item] = setLCwithoutbuff(do3inIbin.LU,do3inIbin.Item);               %  margin�滻����
                    
                    % 2.3   ��ȫ��ƽ��û������,������
                    if max(do3inIbin.LU.LU_Bin(1,:)) == 1
                        fprintf(1,'       Exsiting ����ƽ�� in HBinpingpu (do3)...\n');
                        
                        flagTiledArray(ibin)=1;   %1��������ƽ��                       %do3.LU.LU_VehType = ones(size(d3inIbin.LU.ID)) * do3.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ                        
                        
                        do3Array(ibin) = do3inIbin;
                                   
                        if ISplotEachPingPuAll  % plot ����ƽ��ǰ������ƽ�̺��binͼ
%                             plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'����ƽ��ǰ Bin'); 
                            plotSolutionT(do3inIbin.LU,do3inIbin.Veh,  0, 0, 0 , 1 ,3,'����ƽ�̺� Bin');   % plot����ƽ�̺��bin % do3 �޸ĵ�d�У�����Ŀǰ������do2Array�У�δ��d�ϲ�
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
                
                % ˦βƽ��: ����1:�ҳ���Ҫ˦βƽ�̵�bin,��bin.isTileNeed==1��bin
                %                ����2:�ӳ�β���ҳ���Ӧ�߶ȿ�Ȳ���Strip,��������߶Ѷ�ĶѶ����,����reLoading,ֱ���ռ䲻��������˦β��
                
                do2inIbin = doinIbin;  % ÿ������Ҫ���¿�,�տ�ʼû������,������������Runalgorithm,û��Ҫ,�͸�ֵ֮ǰmaind��do����Ϳ���
                if ~isscalar(do2inIbin.Bin.isTileNeed), error('1'); end
                                               
                while 1   
                    
                    %do2inIbinѭ���ж�,ֱ��Run��û����Ҫ˦βƽ�̻�˦βƽ�̿ռ䲻��(��Ҫ����1��bin)
                    if do2inIbin.Bin.isTileNeed ~= 1, break; end                    
                    if ~isscalar(unique(do2inIbin.LU.LU_Bin(1,:))), break; end  %if max(do2inIbin.LU.LU_Bin(1,:)) ~= 1, break; end
                    
                    nStrip = numel(do2inIbin.Strip.Weight);      if unique(do2inIbin.Strip.Strip_Bin(2, :)) ~= nStrip,    error('��Ԥ�ڴ���');    end
                                        
                    % ���ϸ���:d2inIbin��LU.maxHLayer �Խ���reLoading
                    freLoading = 0;                   
                
                    for istrip = nStrip : -1 : 1  % �ӳ�β������ʼ

                        fstrip = find(do2inIbin.Strip.Strip_Bin(2,:) == istrip); %real strip ���
                        i1 = ~do2inIbin.Strip.isWidthFull(fstrip); %��Ȳ���
                        i2 = ~do2inIbin.Strip.isHeightFull(fstrip) && ~do2inIbin.Strip.isHeightBalance(fstrip); %�߶Ȳ����ҷǾ���
                        
                        i3 = do2inIbin.Strip.isHeightFull(fstrip) &&  ~do2inIbin.Strip.isHeightBalance(fstrip); %�߶����� , ���߶ȷǾ���
                        if i3
                            error('impossible');
                        end
                        if ~(i1 || i2 ) % �������������һ��, �Ϳ��н��н��Ͳ����㷨
                            continue;
                        end

                        % ���е��� �� istrip����������
                        luidxPP = getLuIdx(istrip,do2inIbin);                        
                        [curLuLayer,curLuItemHeight] = getLuLayer(luidxPP,do2inIbin);
                        
                        if all(curLuLayer(luidxPP) <=1) % �������Ϊ1,����
                            continue;
                        else
                            
                            % ��ȡҪ���Ͳ�����idx: 1:���Ѷ� 2:����>1 3:��ߵ�
                            idx = (curLuLayer > 1 & luidxPP);  
                            idx = curLuItemHeight == max(curLuItemHeight(idx));  %�����,�ҳ���߲�Ѷ� , ����ָ������
                            
                            if ~isscalar(unique(do2inIbin.LU.LU_Item(1,idx)))
                                uniItem = unique(do2inIbin.LU.LU_Item(1,idx));
                                idx = do2inIbin.LU.LU_Item(1,:) == uniItem(1);                % ������,ѡ���һ����ͬ�߶ȶѶ⽵�Ͳ���                
                                warning('����Ѷ�߶�һ��,�߶Ⱦ���,������ֻ�����ڿ�Ȳ�����,Ҳ���ܳ�����3���Ѷ�,2��һ����'); 
                            end
                            
                            d2inIbin.LU.maxHLayer(idx) = curLuLayer(idx) - 1;  if any(d2inIbin.LU.maxHLayer<1), error('1'); end % ���Ŀǰ����һ��,�Ҳ�����=0.
                            freLoading = 1;
                        end
           
                        if freLoading
                            break; 
                        end
                        
                    end

                    % ѭ������: û��i���Ͳ���,����reloading
                    if ~freLoading
                        break;
                    end
                    
                    %% �������㷨������ $5 reRunAlgorithm
                    tmpd2 = d2inIbin;
                    [tmpd2.LU] = setLULWHwithbuff(tmpd2.LU, tmpd2.Veh);          
                    do2inIbin = RunAlgorithm(tmpd2,p);        %  do2 = RunAlgorithmPP(d2inIbin,p);  %do2.LU.LU_VehType = ones(size(d2inIbin.LU.ID)) * do2.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ                   
                    [do2inIbin.LU,do2inIbin.Item] = setLCwithoutbuff(do2inIbin.LU,do2inIbin.Item);  % do2Array(ibin) = do2; ����ע�ͣ���Ϊ�Ǹ�ѭ��

                    % $6 ���� ˦βƽ��һ�β�������Ҫѭ���ٽ���˦βƽ�� ֱ������һ������
                    if isscalar(unique(do2inIbin.LU.LU_Bin(1,:)))  % ������LU��һ������,�Ե� % if max(do2inIbin.LU.LU_Bin(1,:)) == 1  %    do2.LU.LU_VehType = ones(size(do2.LU.ID))*do.Veh.order(1);
                        fprintf(1,'       Exsiting ˦βƽ�� in HBinpingpu (do2)...\n');
                        
                        flagTiledArray(ibin)=2; %2����˦βƽ��
              
                        do2Array(ibin) = do2inIbin;                        % todo -> do2 ���ݲ�����d ����return2bba���޸�    % do2 ���ݽ���d???? return2bba���޸ģ�����
                        
                        if ISplotEachPingPuShuaiWei  % plot ˦βƽ��ǰ��˦βƽ�̺��binͼ
%                             plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'˦βƽ��ǰ Bin');
                            plotSolutionT(do2inIbin.LU,do2inIbin.Veh, 0, 0, 0 , 1 ,3,'˦βƽ�̺� Bin');  % plot˦βƽ�̺��bin
                        end
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
%                     chkLUnewold(d2Array(ibin).LU,do2Array(ibin).LU);
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

%% �ֲ�����
 
                                                            
%% ����istrip��ӦLu���߼�ֵ
function LuIdx = getLuIdx(istrip,do)

if istrip == 0,  error('��bin������tileneed,��Ԥ�ڴ���');   end

u=unique(do.LU.LU_Strip(1,:));                                % ��ȡStrip��ŵ�Ψһ����ֵ

fi = find(do.Strip.Strip_Bin(2,:) == istrip );              % ͬbin��strip��� and ˳��>=istrip fi = find(do.Strip.Strip_Bin(2,:) >= istrip );          

LuIdx = ismember(do.LU.LU_Strip(1,:), u(fi));          %% fi->u(fi) ��������� *********************

if ~any(LuIdx),  error('LuIdxȫ��Ϊ��, ������u(fi)��Ӧ��Lu�߼��ж�'); end


end

%% ����LuIdx��ӦLu�Ĳ���/�߶ȵ�
function   [curLuLayer,curLuItemHeight] = getLuLayer(LuIdx,do)

    [curLuLayer,curLuItemHeight] = deal(zeros(size(LuIdx)));
      
    fidx = find(LuIdx);    
    for i=1:length(fidx)
        curLuLayer(fidx(i)) = sum(do.LU.LU_Item(1,:) == do.LU.LU_Item(1, fidx(i))); % ��idx�ڱ�bin�ڶ�Ӧ�Ѷ�Ĳ���
        curLuItemHeight(fidx(i)) = sum(do.LU.LWH(3, (do.LU.LU_Item(1,:) == do.LU.LU_Item(1, fidx(i))))); % ��idx�ڱ�bin�ڶ�Ӧ�Ѷ�Ĳ���
    end 
end


%%  V1 HBinpingpu V2: ɾ��ע��
% % function [flagTiledArray,do2Array,do3Array] = HBinpingpu(maind,do,p)
% % % HBinpingpu ==> ����������maind�����������do������p ��ȡ����ƽ��do3Array/˦βƽ��do2Array���������
% % %   ������
% % %   d2inIbin=d3inIbin��ibin�ڵ��������ݣ�����maind
% % %   do2inIbin=do3inIbin��ibin�ڵ�������ݣ�����do  
% %         
% %         global ISpingpuShuaiWei ISpingpuAll ISshuaiwei
% %         global ISplotEachPingPuAll ISplotEachPingPuShuaiWei 
% %         
% %         % 3.1 ��ʼ��5���������
% %         nbin = length(do.Bin.Weight);
% % 
% %         flagTiledArray = zeros(1,nbin);  %1��������ƽ�� 2����˦βƽ��
% %         d2Array(1:nbin) = maind;    % fixme �޸�maind��ʼ�� ����������  ˦βƽ��
% %         do2Array(1:nbin) = structfun(@(x) [], do, 'UniformOutput', false);
% %         d3Array(1:nbin) = maind;    %����ƽ��
% %         do3Array(1:nbin) = structfun(@(x) [], do, 'UniformOutput', false); %do;
% % 
% %         %% 3.2 ����binѭ��ƽ��
% %         bidx = 1:nbin; % NOTE: �޸�ֻ����˦βƽ�̵�ȫ��Bin���뿼��, �Է�˦βƽ�̵Ľ���ƽ���ж� % bidx = find(do.Bin.isTileNeed);
% %         
% %         % ѭ��: ÿ��bin�ֱ���ƽ��(����ƽ�̺�˦βƽ��ֻ�ܶ�ѡһ����������ƽ��)
% %         for i=1:numel(bidx)
% %             ibin = bidx(i);
% %             
% %             % 1 GET d2inIbin/do2inIbin(d3inIbin/do3inIbin)
% %             d2inIbin = getPartDinThisBin(do, ibin, maind);  % V2  % V1 %    luIdx = do.LU.LU_Bin(1,:) == ibin;   %      d2inIbin = getPartDinThisBin(maind,luIdx,do); %maind��ȡ����margin %�޸ĳɴ�maind��ȡIuIdx������,���Ǵ�������do����ȡ %         d2inIbin = getdinThisVeh(do,luIdx)
% %             doinIbin = getPartDinThisBin(do,ibin);             % v2  % v1  %    luIdx = do.LU.LU_Bin(1,:) == ibin; %                 stripidx = do.Strip.Strip_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip %                 itemidx = do.Item.Item_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip %                 do2inIbin.Veh = do.Veh; %                 do2inIbin.Bin = structfun(@(x) x(:,ibin),do.Bin,'UniformOutput',false); %                 do2inIbin.LU = structfun(@(x) x(:,luIdx),do.LU,'UniformOutput',false);                %                 do2inIbin.Strip = structfun(@(x) x(:,stripidx),do.Strip,'UniformOutput',false); %                 do2inIbin.Item = structfun(@(x) x(:,itemidx),do.Item,'UniformOutput',false);
% % 
% %             d3inIbin = d2inIbin;  
% %             
% %             %% 2 �������ȫ��ƽ��-������(�����Ƿ�˦βƽ��), �۲챾ibin���Ƿ����ȫ��ƽ��,�����,��ȡ��˦βƽ��; ����,����˦βƽ��
% %             if ISpingpuAll==1
% %                 
% %                 % 2.1 ȷ��minmaxLayer��
% %                 % 1 ����case��������һ�� �� һ����ȫһ����ID�ķ���1-2-3�ĵ���ƽ��
% %                 if ibin == numel(bidx) && range(d3inIbin.LU.ID) == 0 
% %                     % minmaxLayer = min(max(d3inIbin.LU.maxHLayer), max(d3inIbin.LU.maxL(3,:))); %maxHLayer��ָ������; maxL(3,:)������������;
% %                     minmaxLayer = min(max(d3inIbin.LU.maxHLayer), max(d3inIbin.LU.maxL(3,:)));  % fixme ��ȡ����ɾ�� maxL���ж�, �����Ǹ���
% %                 else % 2һ��case������ƽ�̼�����һ��
% %                     minmaxLayer = 1;
% %                 end
% %                 
% %                 iLayer = 1;
% %                 
% %                 while 1
% %                     if iLayer > minmaxLayer, break; end
% %                     
% %                     % update d3inIbin -> ��ֵ��d3array
% %                     d3inIbin.LU.maxHLayer(:) = iLayer;   %d3��ȫ��LU�Ĳ����趨Ϊ1  55555 ȫ��ƽ�̵���Ҫ����                    
% %                     d3Array(ibin) = d3inIbin;
% %                     
% %                     iLayer=iLayer+1;
% %                     
% %                     % 2.2  reRunAlgorithm do3��d3�����Ľ��                    
% %                     d3inIbin.LU = setLULWHwithbuff(d3inIbin.LU, d3inIbin.Veh);
% %                     do3inIbin = RunAlgorithm(d3inIbin,p);   % do3Array(ibin) = do3; (����ƽ�̳ɹ��ķ���do3Array�������ķŲ�������ν���ò���
% %                     [do3inIbin.LU,do3inIbin.Item] = setLCwithoutbuff(do3inIbin.LU,do3inIbin.Item);               %  margin�滻����
% %                     
% %                     % 2.3   ��ȫ��ƽ��û������,������
% %                     if max(do3inIbin.LU.LU_Bin(1,:)) == 1
% %                         fprintf(1,'       Exsiting ����ƽ�� in HBinpingpu (do3)...\n');
% %                         
% %                         flagTiledArray(ibin)=1;   %1��������ƽ��                       %do3.LU.LU_VehType = ones(size(d3inIbin.LU.ID)) * do3.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ                        
% %                         
% %                         do3Array(ibin) = do3inIbin;
% %                                    
% %                         if ISplotEachPingPuAll  % plot ����ƽ��ǰ������ƽ�̺��binͼ
% %                             plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'����ƽ��ǰ Bin'); 
% %                             plotSolutionT(do3inIbin.LU,do3inIbin.Veh,  0, 0, 0 , 1 ,3,'����ƽ�̺� Bin');   % plot����ƽ�̺��bin % do3 �޸ĵ�d�У�����Ŀǰ������do2Array�У�δ��d�ϲ�
% %                         end
% %                        
% %                         break;    % ����ƽ�̳ɹ�������while��ѭ����һ��bin
% %                     end
% %                     
% %                 end % END WHILE
% %                 
% %             end % END ISpingpuAll
% %             
% % %% V1 ISpingpuAll
% % % %          if ISpingpuAll==1
% % % %                 
% % % %                 % 2.1 ȷ��minmaxLayer��
% % % %                 % 1 ����case��������һ�� �� һ����ȫһ����ID�ķ���1-2-3�ĵ���ƽ��
% % % %                 if ibin == numel(bidx) && range(d3inIbin.LU.ID) == 0 
% % % %                     minmaxLayer = min(max(d3inIbin.LU.maxHLayer), max(d3inIbin.LU.maxL(3,:))); %maxHLayer��ָ������; maxL(3,:)������������;
% % % %                 else % 2һ��case������ƽ�̼�����һ��
% % % %                     minmaxLayer = 1;
% % % %                 end
% % % %                 
% % % %                 iLayer = 1;
% % % %                 
% % % %                 while 1
% % % %                     if iLayer > minmaxLayer, break; end
% % % %                     
% % % %                     % update d3inIbin -> ��ֵ��d3array
% % % %                     d3inIbin.LU.maxHLayer(:) = iLayer;   %d3��ȫ��LU�Ĳ����趨Ϊ1  55555 ȫ��ƽ�̵���Ҫ����                    
% % % %                     d3Array(ibin) = d3inIbin;
% % % %                     
% % % %                     iLayer=iLayer+1;
% % % %                     
% % % %                     % 2.2  reRunAlgorithm do3��d3�����Ľ��                    
% % % %                     d3inIbin.LU = setLULWHwithbuff(d3inIbin.LU, d3inIbin.Veh);
% % % %                     do3inIbin = RunAlgorithm(d3inIbin,p);   % do3Array(ibin) = do3; (����ƽ�̳ɹ��ķ���do3Array�������ķŲ�������ν���ò���
% % % %                     [do3inIbin.LU,do3inIbin.Item] = setLCwithoutbuff(do3inIbin.LU,do3inIbin.Item);               %  margin�滻����
% % % %                     
% % % %                     % 2.3   ��ȫ��ƽ��û������,������
% % % %                     if max(do3inIbin.LU.LU_Bin(1,:)) == 1
% % % %                         fprintf(1,'       Exsiting ����ƽ�� in HBinpingpu (do3)...\n');
% % % %                         
% % % %                         flagTiledArray(ibin)=1;   %1��������ƽ��                       %do3.LU.LU_VehType = ones(size(d3inIbin.LU.ID)) * do3.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ                        
% % % %                         
% % % %                         do3Array(ibin) = do3inIbin;
% % % %                                    
% % % %                         if ISplotEachPingPuAll  % plot ����ƽ��ǰ������ƽ�̺��binͼ
% % % %                             plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'����ƽ��ǰ Bin'); 
% % % %                             plotSolutionT(do3inIbin.LU,do3inIbin.Veh,  0, 0, 0 , 1 ,3,'����ƽ�̺� Bin');   % plot����ƽ�̺��bin % do3 �޸ĵ�d�У�����Ŀǰ������do2Array�У�δ��d�ϲ�
% % % %                         end
% % % %                        
% % % %                         break;    % ����ƽ�̳ɹ�������while��ѭ����һ��bin
% % % %                     end
% % % %                     
% % % %                 end % END WHILE
% % % %                 
% % % %             end % END ISpingpuAll
% %             
% %                 %% �ж���������ƽ�̳ɹ���������һ��bin��������˦βƽ��
% %             if flagTiledArray(ibin)==1
% %                 continue;
% %             end
% %             
% %             %% 3 ������ƽ��ʧ�� ����˦βƽ�̣�ǰ�᣺��˦β������ ��Ϊ����  
% %             if ISpingpuShuaiWei==1 && ISshuaiwei==1
% %                 
% %                 % ˦βƽ��: ����1:�ҳ���Ҫ˦βƽ�̵�bin,��bin.isTileNeed==1��bin
% %                 %                ����2:�ӳ�β���ҳ���Ӧ�߶ȿ�Ȳ���Strip,��������߶Ѷ�ĶѶ����,����reLoading,ֱ���ռ䲻��
% %                 
% %                 do2inIbin = doinIbin;  % ÿ������Ҫ���¿�,�տ�ʼû������,������������Runalgorithm,û��Ҫ,�͸�ֵ֮ǰmaind��do����Ϳ���
% %                 if ~isscalar(do2inIbin.Bin.isTileNeed), error('1'); end
% %                 
% %                 plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'˦βƽ��ǰ Bin');
% %                 plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,8,'˦βƽ��ǰ Bin');
% %                                
% %                 while 1   
% %                     
% %                     %do2inIbinѭ���ж�,ֱ��Run��û����Ҫ˦βƽ�̻�˦βƽ�̿ռ䲻��(��Ҫ����1��bin)
% %                     if do2inIbin.Bin.isTileNeed ~= 1, break; end                    
% %                     if ~isscalar(unique(do2inIbin.LU.LU_Bin(1,:))), break; end  %if max(do2inIbin.LU.LU_Bin(1,:)) ~= 1, break; end
% %                     
% %                     nStrip = numel(do2inIbin.Strip.Weight);      if unique(do2inIbin.Strip.Strip_Bin(2, :)) ~= nStrip,    error('��Ԥ�ڴ���');    end
% %                                         
% %                     % ���ϸ���:d2inIbin��LU.maxHLayer �Խ���reLoading
% %                     freLoading = 0;
% %                     
% % %                     uniStrip = unique(do2inIbin.LU.LU_Strip(1,:));
% %                     
% %                     for istrip = nStrip : -1 : 1  % �ӳ�β������ʼ
% % 
% % %                         do2inIbin.LU.LU_Strip(1,:) == istrip
% % %                         do2inIbin.Strip.Strip_Bin
% %                         fstrip = find(do2inIbin.Strip.Strip_Bin(2,:) == istrip)
% %                         
% %                         
% %                         i1 = ~do2inIbin.Strip.isWidthFull(fstrip); %��Ȳ���
% %                         i2 = ~do2inIbin.Strip.isHeightFull(fstrip) && ~do2inIbin.Strip.isHeightBalance(fstrip); %�߶Ȳ����ҷǾ���
% %                         
% %                         if ~(i1 || i2) % �������һ��������, �ͼ�����Ҫ˦βƽ��
% %                             continue;
% %                         end
% % 
% %                         % ���е��� �� istrip����������
% %                         luidxPP = getLuIdx(istrip,do2inIbin);                        
% %                         [curLuLayer,curLuItemHeight] = getLuLayer(luidxPP,do2inIbin);
% %                         
% %                         if all(curLuLayer(luidxPP) <=1)
% %                             continue;
% %                         else                            
% %                             idx = (curLuLayer > 1 & luidxPP);  % �϶���curLuLayer>=2
% %                             % ������stripbalance��
% %                             idx2 = curLuItemHeight == max(curLuItemHeight(idx));  %�����,�ҳ���߲�,����ָ������
% %                             
% %                             
% %                                     if ~isscalar(unique(do2inIbin.LU.LU_Item(1,idx2))), 
% %                                         uniItem = unique(do2inIbin.LU.LU_Item(1,idx2))
% %                                         idx2 = do2inIbin.LU.LU_Item(1,:) == uniItem(1)
% %                                         
% %                                         warning('����Ѷ�߶�һ��,�߶Ⱦ���,������ֻ�����ڿ�Ȳ�����'); end
% %                                     
% %                             d2inIbin.LU.maxHLayer(idx2) = curLuLayer(idx2) - 1; % ���Ŀǰ����һ��,�Ҳ�����=0.
% %                                                     if any(d2inIbin.LU.maxHLayer<1), error('1'); end
% %                             freLoading = 1;
% %                         end
% %                         
% % %                         d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP));                                 
% % %                         d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxHLayer(luidxPP));                        
% % %                         d2inIbin.LU.maxHLayer(luidxPP) 
% % %                         curLuLayer(luidxPP)
% % 
% % %                         while any(d2inIbin.LU.maxHLayer(luidxPP) > 1)  % ����ǰluidxPP�ڵ����� ���ڷ�ƽ�̵�LU, ����
% %                             
% %                             % ��ǰLu,�ҶѶ����1,����ͬ�Ѷ���ߵ�
% % %                             a = curLuItemHeight == max(curLuItemHeight(luidxPP)) & luidxPP & curLuLayer > 1;
% %                             
% % %                              m2 = max(curLuLayer(luidxPP)) %��߶Ѷ�Ĳ���
% % %                              idx2 = curLuLayer(:) == m2 %��߶Ѷ�Ĳ�������index��
% %                              
% % %                              m = max(d2inIbin.LU.maxHLayer(idx2))
% % %                             idx = d2inIbin.LU.maxHLayer(:) == m  %ֻ�����ֵӦ�ÿ���,����ֵ����maxHlayer��Ϊ��ֵtoo
% % %                             a = idx'&luidxPP
% %                             
% % %                             if sum(a) <=0 , 
% % %                                 continue;
% % %                             end
% % %                             d2inIbin.LU.maxHLayer(a) = d2inIbin.LU.maxHLayer(a) - 1;
% %                             
% % %                             d2inIbin.LU.maxHLayer(luidxPP) = d2inIbin.LU.maxHLayer(luidxPP)  - 1;
% % %                             d2inIbin.LU.maxHLayer(d2inIbin.LU.maxHLayer<=1) = 1;  % �޸�: ������Ļָ�Ϊ1
% %                             
% % %                             d2inIbin.LU.maxHLayer(luidxPP) 
% %                             
% % %                             if any(d2inIbin.LU.maxHLayer(luidxPP) < curLuLayer(luidxPP))  % �����º������� < ��ǰLU��Ӧ����
% % %                                 d2inIbin.LU.maxHLayer(luidxPP)
% % %                                 curLuLayer(luidxPP)
% % %                                 freLoading = 1;
% % %                                 break;
% % %                             else
% % %                                 continue;
% % %                             end
% % 
% % %                         end
% %                         
% %                         if freLoading
% %                             break; 
% %                         end
% %                         
% %                     end
% % %%                    
% % %                     for istrip = nStrip : -1 : 1
% % %                         
% % %                         luidxPP = getLuIdx(istrip,do2inIbin);
% % %                         
% % %                         curLuLayer = getLuLayer(luidxPP,do2inIbin);
% % %                         
% % % %                         d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP));                                 
% % % %                         d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxHLayer(luidxPP));                        
% % %                         d2inIbin.LU.maxHLayer(luidxPP) 
% % %                         
% % %                         if any(d2inIbin.LU.maxHLayer(luidxPP) > 1)  % ����ǰluidxPP�ڵ����� ���ڷ�ƽ�̵�LU, ����
% % %                             d2inIbin.LU.maxHLayer(luidxPP) = d2inIbin.LU.maxHLayer(luidxPP)  - 1;
% % %                             d2inIbin.LU.maxHLayer(d2inIbin.LU.maxHLayer<=1) = 1;  % �޸�: ������Ļָ�Ϊ1
% % %                             
% % %                             d2inIbin.LU.maxHLayer(luidxPP) 
% % %                             
% % %                             if any(d2inIbin.LU.maxHLayer(luidxPP) < curLuLayer(luidxPP))  % �����º������� < ��ǰLU��Ӧ����
% % %                                 d2inIbin.LU.maxHLayer(luidxPP)
% % %                                 curLuLayer(luidxPP)
% % %                                 fre = 1;
% % %                                 break;
% % %                             else
% % %                                 continue;
% % %                             end
% % %                             
% % %                         else
% % %                             continue; 
% % %                         end
% % %                         
% % %                     end
% %                     
% %                     
% % %% �ϰ汾����                    
% % %                     istrip= nStrip;
% % %                     
% % %                     luidxPP = getLuIdx(istrip,do2inIbin);
% % %                                         
% % %                     d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP)) - 1;
% % %                     
% % %                     while all(d2inIbin.LU.maxHLayer(luidxPP)<1)
% % %                         
% % %                         istrip = istrip-1;
% % %                         
% % %                         if istrip==0,  break;  end
% % %                         
% % %                         luidxPP = getLuIdx(istrip,do2inIbin);
% % %                         
% % %                        d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP)) - 1;
% % %                         
% % %                     end
% % %                     
% % %                     d2inIbin.LU.maxHLayer(d2inIbin.LU.maxHLayer<=1) = 1;  % �޸�: ������Ļָ�Ϊ1
% %                     if ~freLoading
% %                         break;
% %                     end
% %                     %% �������㷨������ $5 reRunAlgorithm
% %                     tmpd2 = d2inIbin;
% %                     d2inIbin.LU.maxHLayer(luidxPP) 
% %                     [tmpd2.LU] = setLULWHwithbuff(tmpd2.LU, tmpd2.Veh);          
% %                     do2inIbin = RunAlgorithm(tmpd2,p);        %  do2 = RunAlgorithmPP(d2inIbin,p);  %do2.LU.LU_VehType = ones(size(d2inIbin.LU.ID)) * do2.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ                   
% %                     [do2inIbin.LU,do2inIbin.Item] = setLCwithoutbuff(do2inIbin.LU,do2inIbin.Item);  % do2Array(ibin) = do2; ����ע�ͣ���Ϊ�Ǹ�ѭ��
% %                     do2inIbin.LU.maxHLayer(luidxPP) 
% % 
% %                     % $6 ���� ˦βƽ��һ�β�������Ҫѭ���ٽ���˦βƽ�� ֱ������һ������
% %                     if isscalar(unique(do2inIbin.LU.LU_Bin(1,:)))  % ������LU��һ������,�Ե�
% %                     % if max(do2inIbin.LU.LU_Bin(1,:)) == 1  %    do2.LU.LU_VehType = ones(size(do2.LU.ID))*do.Veh.order(1);
% %                         fprintf(1,'       Exsiting ˦βƽ�� in HBinpingpu (do2)...\n');
% %                         
% %                         flagTiledArray(ibin)=2; %2����˦βƽ��
% %               
% %                         do2Array(ibin) = do2inIbin;                        % todo -> do2 ���ݲ�����d ����return2bba���޸�    % do2 ���ݽ���d???? return2bba���޸ģ�����
% %                         
% %                         if ISplotEachPingPuShuaiWei  % plot ˦βƽ��ǰ��˦βƽ�̺��binͼ
% % %                             plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'˦βƽ��ǰ Bin');
% %                             plotSolutionT(do2inIbin.LU,do2inIbin.Veh, 0, 0, 0 , 1 ,3,'˦βƽ�̺� Bin');  % plot˦βƽ�̺��bin
% %                         end
% %                     end                    
% %                 
% %                 end % END OF do2 WHILE ����˦β
% %                 
% %                 
% %                 
% % % %                 %% ��ȡluidxPP �������� d2inIbin.LU.maxHLayer(luidxPP)
% % % %                 while do2inIbin.Bin.isTileNeed == 1 %do2inIbin����ƽ�̺��bin����Ҫƽ��,������while�ж�
% % % %                     % $3.4.2 �޶�d2.LU.maxHLayer (����ibin�����ѡ���ļ���stripƽ��) TODO $4д����Щ����,���ڼ�
% % % %                     % $4.1 ��do2inIbin��ȡ luidxPP �� ĳ��strip��Ӧ��LU�߼�ֵ,
% % % %                     % ѭ���ӱ�ibin�����һ��strip��ʼƽ�� istrip= nbStrip;
% % % %                     nStrip = numel(do2inIbin.Strip.Weight);  if unique(do2inIbin.Strip.Strip_Bin(2, :)) ~= nStrip,    error('��Ԥ�ڴ���');    end
% % % %                     
% % % %                     % ��istrip�����һ��������ʼ
% % % %                     istrip= nStrip;
% % % %                     
% % % %                     luidxPP = getLuIdx(istrip,do2inIbin);
% % % %                                         
% % % %                     % $4.2 �޶�d2.LU.maxHLayer           d2inIbin.LU    maind.LU
% % % %                     % ����bin�ڵ� ĳЩ���� (luidxPP) ���������ڷŲ��� - 1 (����һ��,���ܴ���ͬһStrip������ͬcase,TODO)
% % % % %                     d2inIbin.LU.maxHLayer(luidxPP)
% % % % %                     d2inIbin.LU.maxL(3,luidxPP)
% % % %                     d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP)) - 1;
% % % %                     
% % % %                     % $4.2 ����ǰluidxPP��ӦLu�Ĳ������Ѿ�Ϊ1��, ����Ҫ���Ӹ����istrip��luidxPP; ���޶�d2.LU.maxHLayer
% % % %                     % GET ���� d2inIbin.LU.maxHLayer(luidxPP) ����luidxPP�Ĳ���>=2��
% % % %                     while all(d2inIbin.LU.maxHLayer(luidxPP)<1)
% % % %                         
% % % %                         istrip = istrip-1;
% % % %                         
% % % %                         if istrip==0,  break;  end
% % % %                         
% % % %                         luidxPP = getLuIdx(istrip,do2inIbin);
% % % %                         
% % % %                        d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP)) - 1;
% % % %                         
% % % %                     end
% % % %                     % �޸�: ������Ļָ�Ϊ1
% % % %                     d2inIbin.LU.maxHLayer(d2inIbin.LU.maxHLayer<=1) = 1;
% % % %                     
% % % %                     %% �������㷨������ $5 reRunAlgorithm
% % % % %                     d2Array(ibin) = do2inIbin;
% % % % %                     d2inIbin.LU.Rotaed(:)
% % % % %                     d2inIbin.LU.Rotaed(:) = 0;
% % % % %                     d2inIbin.LU.Rotaed(:)
% % % % d2inIbin111 = d2inIbin
% % % %                     [d2inIbin111.LU] = setLULWHwithbuff(d2inIbin111.LU, d2inIbin111.Veh);          
% % % % %                     d2inIbin.LU.Rotaed(:)
% % % %                      d2inIbin.LU.LWH
% % % %                     do2inIbin = RunAlgorithm(d2inIbin111,p);        %  do2 = RunAlgorithmPP(d2inIbin,p);  %do2.LU.LU_VehType = ones(size(d2inIbin.LU.ID)) * do2.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ                   
% % % %                     [do2inIbin.LU,do2inIbin.Item] = setLCwithoutbuff(do2inIbin.LU,do2inIbin.Item);  % do2Array(ibin) = do2; ����ע�ͣ���Ϊ�Ǹ�ѭ��
% % % % %                     [d2inIbin.LU,d2inIbin.Item] = setLCwithoutbuff(d2inIbin.LU,d2inIbin.Item);
% % % %                     
% % % % %                     d2inIbin.LU.LWH = do2inIbin.LU.LWH
% % % % %                     d2inIbin.LU.CoordLUBin = do2inIbin.LU.CoordLUBin
% % % % %                     d2inIbin.Item.LWH = do2inIbin.Item.LWH
% % % % %                     d2inIbin.Item.CoordItemBin = do2inIbin.Item.CoordItemBin
% % % %                     
% % % %                     plotSolutionT(do2inIbin.LU,do2inIbin.Veh, 0, 0, 0 , 1 ,3,'˦βƽ�̺���ܶ��bin Bin');
% % % %                     % plot
% % % %                     
% % % %                     
% % % %                     % $6 ���� ˦βƽ��һ�β�������Ҫѭ���ٽ���˦βƽ�� ֱ������һ������
% % % %                     if max(do2inIbin.LU.LU_Bin(1,:)) == 1  %    do2.LU.LU_VehType = ones(size(do2.LU.ID))*do.Veh.order(1);
% % % %                         fprintf(1,'       Exsiting ˦βƽ�� in HBinpingpu (do2)...\n');
% % % %                         flagTiledArray(ibin)=2; %2����˦βƽ��
% % % %                         
% % % % %                         [do2inIbin.LU,do2inIbin.Item] = setLCwithoutbuff(do2inIbin.LU,do2inIbin.Item); 
% % % %                         
% % % %                         do2Array(ibin) = do2inIbin;                        % todo -> do2 ���ݲ�����d ����return2bba���޸�    % do2 ���ݽ���d???? return2bba���޸ģ�����
% % % %                         
% % % %                         if ISplotEachPingPuShuaiWei  % plot ˦βƽ��ǰ��˦βƽ�̺��binͼ
% % % %                             plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'˦βƽ��ǰ Bin');
% % % %                             plotSolutionT(do2inIbin.LU,do2inIbin.Veh, 0, 0, 0 , 1 ,3,'˦βƽ�̺� Bin');  % plot˦βƽ�̺��bin
% % % %                         end
% % % %                         
% % % %                     else
% % % %                         break;  %����˦β�Ų��� ���ټ���˦βƽ���ˣ����˽�����������һ����˦β�ж�
% % % %                     end                    
% % % %                 
% % % %                 end % END OF do2 WHILE ����˦β
% %                 
% %             end % END OF ISpingpuShuaiwei
% %             
% %         end% END OF FOR
% %         
% %         %% chk
% %         if any(flagTiledArray)
% %             for ibin=1:nbin
% %                 if flagTiledArray(ibin)==1 %����ƽ��
% %                     %                 t3 = struct2table(structfun(@(x) x',d3Array(ibin).LU,'UniformOutput',false));
% %                     %                 to3 = struct2table(structfun(@(x) x',do3Array(ibin).LU,'UniformOutput',false));
% %                     chkLUnewold(d3Array(ibin).LU,do3Array(ibin).LU);
% %                     chktLU(d3Array(ibin).LU);
% %                     chktLU(do3Array(ibin).LU);
% %                     
% %                 end
% %                 if flagTiledArray(ibin)==2 %˦βƽ��
% %                     %                 t2 = struct2table(structfun(@(x) x',d2Array(ibin).LU,'UniformOutput',false));
% %                     %                 to2 = struct2table(structfun(@(x) x',do2Array(ibin).LU,'UniformOutput',false));
% % %                     chkLUnewold(d2Array(ibin).LU,do2Array(ibin).LU);
% %                     chktLU(d2Array(ibin).LU);
% %                     chktLU(do2Array(ibin).LU);
% %                 end
% %                 
% % %                 flagTiledArray;
% % %                 d2Array;
% % %                 do2Array;
% % %                 d3Array;
% % %                 do3Array;
% %             end
% %         end
% % 
% % end