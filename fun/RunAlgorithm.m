function [d] = RunAlgorithm(d,p)
        global  ISshuaiwei ISreStripToBin  %����ִ��HStripToBin����
        global  ISplotshuaiwei  ISplotStripToBinAgain   ISplotRunAlgo  ISplotRunLIS

        %% Ԥ����:����Input��������
        fprintf(1,'   Running RunAlgorithm  ...\n');
        
        %         d = GcheckInput(d);    %���Բ���   ��Ϊmain��һ�¼���
                
        % ����Ԥ�����ص㣺��ȡLU.Rotaed,�����Ƿ����� ����cpuLU���� ���������LU������  ������֮ǰLWH�Ѿ�ʱ����margin�ģ�        
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        [d.LU, d.Veh] = cpuVehLU(d.LU, d.Veh);            %������>??  in case �������

        %% ����ʽ: cpuLU 
        % ����LU��isNonMixed��LU��isMixedTile ���� ��;�ƺ�����, Ŀ��Ϊ������Itemʱ,�����һ��,������ĺ���
         
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % ADD: isNonmixed/ismixedtile
        %         [d.LU] = cpuLU(d.LU,d.Veh);     % REMOVE ��ΪisNonmixed/ismixedtile ���ٱ�Ҫ
        
        %% ����ʽ: HLUtoItem LU��Item���㷨 
        %fprintf(1,'     Running HLUtoItem...\n');
        
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % ADD: order/LU_Item
        % Item:
        % ADD: isRota/Rotaed/HLayer/LWH/Weight
        [d.LU,d.Item] = HLUtoItem(d.LU,d.Veh);          %Item����ID������򣨵���һ�������仯˳��
        
        %% ����ʽ: repairItems �޸�LU_Item���Ƿ���������(�޸�˳�����,�����ʱ��δ����)
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % REPAIR: LU_Item
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        d.LU.LU_Item(2,:)= repairItems(d.LU);
        
        % �˲��Ƿ����������ص����
        chktLU(d.LU)
    
        %% ����ʽ: cpuItem - ����Item�� MixOrder isHeightFull nbItem ��������
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % ADD: /isHeightFull/MixOrder/nbItem/IDs
        [d.Item] = cpuItem(d.Item,d.LU,d.Veh);        % printstruct(d,'sortfields',1,'PRINTCONTENTS',0);
        %        printstruct(d.LU,'sortfields',1,'PRINTCONTENTS',1);

        %% ����ʽ: HItemToStrip �� Item��Strip���㷨        
        % 1 Item����: %1: SID ; 2: MixOrder; һ�����濪ʼ: 
        %                       3: Longth/Height; 4:Width; 5: LID; (3,4,5,����һ��) 6: Height
        % 2 ��Item����˳������»�ɵ�Strip
        
        %fprintf(1,'     Running HItemToStrip...\n');     
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/MixOrder/nbItem/IDs
        % ADD: itemorder/Item_Strip/CoordItemStrip
        % Strip:
        % ADD: LW/Weight
        [d.Item,d.Strip] = HItemToStrip(d.LU,d.Item,d.Veh,p);     %   printstruct(d);   %  printstruct(d.Item);
      
        %% ����ʽ: cpuLU2Strip �� ����LU_Strip, CoordLUStrip  
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % ADD: LU_Strip, CoordLUStrip
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/MixOrder/nbItem/IDs
        % itemorder/Item_Strip/CoordItemStrip
        % Strip:
        % LW/Weight        
        [d.LU] = cpuLU2Strip(d.LU,d.Item);     %   printstruct(d);   %  printstruct(d.Item);
        
        %% ����ʽ: cpuStrip �����ؼ�strip���Լ���         
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % LU_Strip, CoordLUStrip
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/MixOrder/nbItem/IDs
        % itemorder/Item_Strip/CoordItemStrip
        % Strip:
        % LW/Weight
        % ADD: isMixedIDs/nbIDs/isHeightFull/isHeightBalance/isWidthFull/nbItem/nbLU/nbLULID
        % ADD: isAllPured/isSingleItem/mlmdHeight/IDs/~Limit
        % ADD: nLUID / nLULID
         [d.Strip] = cpuStrip(d.Strip,d.Item,d.LU,d.Veh); %          [d.Strip,d.LU] = cpuStrip(d.Strip,d.Item,d.LU,d.Veh);

        %% ����ʽ��HStripToBin - Strip��Bin���㷨
        % ********* 1 Strip����: % 1: SID 2: priorityofLID
        %fprintf(1,'     Running HStripToBin...\n');
        
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % LU_Strip, CoordLUStrip
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/MixOrder/nbItem/IDs
        % itemorder/Item_Strip/CoordItemStrip
        % Strip:
        % LW/Weight
        % isMixedIDs/nbIDs/isHeightFull/isHeightBalance/isWidthFull/nbItem/nbLU/nbLULID
        % isAllPured/isSingleItem/mlmdHeight/IDs/~Limit
        % ADD: striporder/strip_bin
        % Bin:
        % ADD: LW/Weight
        [d.Strip,d.Bin] = HStripToBin(d.Strip,d.Veh,p);
        
         %% ����ʽ��HItemToBin - % ����LU��Bin������and˳������ݼ��� 
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % LU_Strip, LU.CoordLUStrip
        % ADD: LU_Bin/CoordLUBin/ lbrt_A / nIDBin / nLIDBin 
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/MixOrder/nbItem/IDs
        % itemorder/Item_Strip/CoordItemStrip
        % ADD: Item_bin/CoordItemBin
        % Strip:
        % LW/Weight
        % isMixedIDs/nbIDs/isHeightFull/isHeightBalance/isWidthFull/nbItem/nbLU/nbLULID
        % isAllPured/isSingleItem/mlmdHeight/IDs/~Limit
        % striporder/strip_bin
        % ADD: nLUID / nLULID����ʱ����)  / nLUIDBin / nLULIDBin  
        % Bin:
        % LW/Weight
        [d.LU,d.Item,d.Strip] = HItemToBin(d.LU,d.Item,d.Strip); % ����LU��Bin������and˳��   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin

        %% ����ʽ��cpuBin - ����isTileNed �Ƿ���Ҫ˦βƽ�̵�binָ��
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % LU_Strip, LU.CoordLUStrip
        % LU_Bin/CoordLUBin/ lbrt_A
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/MixOrder/nbItem/IDs
        % itemorder/Item_Strip/CoordItemStrip
        % Item_bin/CoordItemBin
        % Strip:
        % LW/Weight
        % isMixedIDs/nbIDs/isHeightFull/isHeightBalance/isWidthFull/nbItem/nbLU/nbLULID
        % isAllPured/isSingleItem/mlmdHeight/IDs/~Limit
        % striporder/strip_bin
        % Bin:
        % LW/Weight        
        % ADD: isTileNed/IDs
        [d.Bin] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %����Bin��isTileNeed ���߶ȣ���Ȳ�����
        
        
        if ISplotRunAlgo && ISplotRunLIS
            %             plotSolutionT(d.LU,d.Veh,1,0,0,0,3,'ԭ˳��LU'); %LU��ǰ˳��
                        plotSolutionT(d.LU,d.Veh,3,0,0,0,3,'�����LU');     %LU�������ΪLu.order��ʱ��֪
                        plotSolutionT(d.LU,d.Veh,0,1,0,0,3,'����ǰItem',[],[],[]);   %Item����ǰ����ΪItem->Stripǰ������LU_Item��û��
                        plotSolutionT(d.LU,d.Veh,0,2,0,0,3,'�����Item',[],d.Item.itemorder,[]);
                        plotSolutionT(d.LU,d.Veh,0,0,1,0,3,'����ǰStrip');     %Strip����ǰ����ΪStrip->Binǰ������LU_Strip��û��
            plotSolutionT(d.LU,d.Veh,0,0,2,0,1,'�����Strip',[],[],d.Strip.striporder);     %Strip�����
            plotSolutionT(d.LU,d.Veh,0,0,0,1,1,'ID��ͼBin');                  %     Bin�����
            %             plotSolutionT(d.LU,d.Veh,0,0,0,1,1,'LID��ͼBin');            % Bin����� 
        end
                                                    
        %% ����ʽ������ͷ ����2: ÿ��ʣ��stripȫ���ڱȽ��� better than ����1 �й���? % �ر��ǶԵ����������Ժ�
        if ISreStripToBin==1
                fprintf(1,'     Running HStripToBinAgain...\n');                
                % LU:
                % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
                % Roated/maxHLayer/nbID/nbLID
                % isNonmixed/ismixedtile
                % order/LU_Item
                % LU_Strip, LU.CoordLUStrip
                % LU_Bin/CoordLUBin/ lbrt_A
                % Item:
                % isRota/Rotaed/HLayer/LWH/Weight
                % isHeightFull/MixOrder/nbItem/IDs
                % itemorder/Item_Strip/CoordItemStrip
                % Item_bin/CoordItemBin
                % Strip:
                % LW/Weight
                % isMixedIDs/nbIDs/isHeightFull/isHeightBalance/isWidthFull/nbItem/nbLU/nbLULID
                % isAllPured/isSingleItem/mlmdHeight/IDs/~Limit
                % striporder/strip_bin
                % Bin:
                % LW/Weight
                % isTileNed/IDs
                [d.Strip,d.Bin,TFHStripToBinAgain] = HStripToBinAgain(d.Bin,d.Strip,d.Item,d.LU,d.Veh,p); 
                
                % sum(~d.Strip.isHeightBalance & d.Strip.isHeightFull)
                
                [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip); % ����LU��Bin������and˳��   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin

                [d.Bin] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %����Bin��������� % ����isTileNeed

                if TFHStripToBinAgain &&  ISplotStripToBinAgain && ISplotRunAlgo
                      plotSolutionT(d.LU,d.Veh,0,0,0,1,1,'����ͷ��Bin'); % Bin�����
                end
        end
                                                        
        %% ����ʽ������Strip��˦β�Ż� *********** �޸� Strip.Strip_Bin �ڶ���
        if ISshuaiwei==1   
                fprintf(1,'     Running HStripSW...\n');

                % ����strip��strip_bin˳��, strip�Ƿ�˦β,strip˦β˳�� 
                [d.Strip.Strip_Bin(2,:),  d.Strip.isShuaiWei, d.LU.isShuaiWei, TFHStripSW] = HStripSW(d.Strip,d.LU);

                % ����LU/Item ��Bin������and˳��LU_Bin  NOTE: Lu_Strip and Item_Strip ����˦βӰ��
                [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip);    %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin
       
                [d.Bin] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %����Bin��isTileNeed

                if TFHStripSW && ISplotshuaiwei && ISplotRunAlgo
                        plotSolutionT(d.LU,d.Veh,0,0,0,1,1,'˦β��Bin'); % Bin����� 
%                     plotSolutionT(d.LU,d.Veh,0,0,0,1,8,'˦β��Bin'); % Bin����� 
                end                
        end

        %% ����ʽ�����һ������ֵ����
       d.LU.LU_VehType = ones(size(d.LU.ID)) * d.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ
       fprintf(1,'   Running RunAlgorithm  DONE ...\n');                

end

%         printstruct(d.Item,'sortfields',1,'PRINTCONTENTS',0)
%         printstruct(d.Bin,'sortfields',1,'PRINTCONTENTS',1)
%         printstruct(d,'sortfields',1,'PRINTCONTENTS',0)

% cpuLU2Strip ����LU.LU_Strip, LU.CoordLUStrip
function [LU] = cpuLU2Strip(LU,Item)

nbLU = size(LU.LWH,2);
LU.LU_Strip = zeros(2,nbLU);
LU.CoordLUStrip = zeros(3,nbLU);

% ����LU_Strip
for iLU=1:nbLU
    % ����LU_Strip��һ��
    iItem = LU.LU_Item(1,iLU);   %iLU���ڵڼ���Item, Item���ڵڼ���Strip,��Lu���ڵڼ���Strip
    LU.LU_Strip(1,iLU)= Item.Item_Strip(1,iItem);
    % ����LU_Strip�ڶ���
    fiItem = find(Item.Item_Strip(1,:) == Item.Item_Strip(1,iItem) & Item.Item_Strip(2,:) < Item.Item_Strip(2,iItem));
    nbLUfiItem = sum(ismember(LU.LU_Item(1,:),fiItem));
    LU.LU_Strip(2,iLU) = nbLUfiItem+LU.LU_Item(2,iLU); % ����Strip˳��: ͬһStrip����ǰ�������nbLUfiItem + ��iLU��Item��˳��
    % ����LU.CoordLUStrip
    LU.CoordLUStrip(1,iLU) = Item.CoordItemStrip(1,iItem);
    LU.CoordLUStrip(2,iLU) = Item.CoordItemStrip(2,iItem);
    % fLU: ��iLUͬ��iItem �� ˳�����ڱ�iLU; ����Ϊ��, ��Ӱ��.
    fLU = LU.LU_Item(1,:) == iItem & LU.LU_Item(2,:) < LU.LU_Item(2,iLU);
    LU.CoordLUStrip(3,iLU) = sum(LU.LWH(3,fLU));
end

end



%% ����ʽ�����⣺�Ѷ�������� ��cpuStrip֮��
        % ********* 1 5555 �Ѷ�������� *********** �޸� d.LU.maxHLayer(luidxPP) 
%          if ISstripbalance==0
%              [d.Strip,d.Item,d.LU,TFStripBalance] = HStripBalance(d.Strip,d.Item,d.LU,d.Veh,p);  end
%          if ISstripbalance && TFStripBalance && ISplotstripbalance && ISplotRunAlgo
%              plotSolutionT(d.LU,d.Veh,0,0,1,0,3,'�����Strip'); %Strip:HStripBalance��
%          end     % V1:  if ISplotStrip==1,      figure(222);      plot3DStrip(d.LU,d.Item,d.Veh,'LU');        end    %  plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % ����LU.CoordLUBin %   figure(222);     plot3DStrip(d.LU,d.Item,d.Veh,'LU');         %    igure(111);          plot3DStrip(d.LU,d.Item,d.Veh,'Item'); 

       %% ����: PID/SID�ȷ��� ��Ϊ����˳�򣬴��಻����ֵ��ԭʼ��¼�л�ȡ
%        d.LU.SID = d.LU.OSID;
%        d.LU.PID = d.LU.OPID;

        %% ��ӡ���
%     printOut(d.Bin,d.Strip,d.Item,d.LU,d.Veh); %����,��ʱע��
% %         function printOut(Bin,Strip,Item,LU,Veh)
% %             nBin = size(Bin.LW,2);
% %             for iBin = 1: nBin
% %                 [~,ibin] = find(Item.Item_Bin(1,:)==iBin); %��iBin�µ�item������
% %                 idxSeq = Item.Item_Bin(2,ibin); %��iBin��item����˳��Seq
% %                 fprintf('bin �Ŀ�+��+��Ϊ: ' );
% %                 fprintf(' %d  ',Veh.LWH);
% %                 fprintf('\n');
% %                 fprintf('bin %d ��ʣ���+ʣ�೤Ϊ:  ',iBin);fprintf('\n');
% %                 fprintf('( %d ) ',Bin.LW(:,iBin));fprintf('\n');
% %                 fprintf('\n');
% % 
% %                 fprintf('bin %d ���� original item ������{˳��}(����)[��ת��־]{����}Ϊ  \n  ',iBin);
% %                 fprintf('%d ',ibin);fprintf('\n');
% %                 fprintf('{%d} ',idxSeq);fprintf('\n');
% %                 fprintf(' (%d %d %d) ', Item.LWH(1:3,ibin));fprintf('\n');
% %                 fprintf(' [%d]     ', Item.Rotaed(:,ibin));fprintf('\n');
% %                 fprintf(' {%d %d %d} ', Item.CoordItemBin(:,ibin));fprintf('\n');
% %                 fprintf('\n');
% % 
% %                 [~,idxLU] = find(LU.LU_Bin(1,:)==iBin); %��iBin�µ�item������
% %                 fprintf('bin %d ���� original LU ������{˳��}[item���](����)[��ת��־]{����}Ϊ  \n  ',iBin);
% %                 idxLUSeq = LU.LU_Bin(2,idxLU); %��iBin��item����˳��Seq
% %                 idxLUItem = LU.LU_Item(1,idxLU);
% %                 fprintf('%d ',idxLU);fprintf('\n');
% %                 fprintf('{%d} ',idxLUSeq);fprintf('\n');
% %                 fprintf('[%d] ',idxLUItem);fprintf('\n');
% %                 fprintf(' (%d %d %d) ', LU.LWH(1:3,idxLU));fprintf('\n');
% %                 fprintf(' [%d]     ', LU.Rotaed(:,idxLU));fprintf('\n');
% %                 fprintf(' {%d %d %d} ', LU.CoordLUBin(:,idxLU));fprintf('\n');
% %                 fprintf('\n');
% % 
% %                 % ������˳��չʾ
% %                 % %     [~,x]=sort(LU.LU_Bin(2,idxLU));
% %                 % %     idxLUSeq = idxLUSeq(x); %��iBin��item����˳��Seq
% %                 % %     idxLUItem = idxLUItem(x);
% %                 % %     fprintf('%d ',idxLU);fprintf('\n');
% %                 % %     fprintf('{%d} ',idxLUSeq);fprintf('\n');
% %                 % %     fprintf('[%d] ',idxLUItem);fprintf('\n');
% %                 % %     fprintf(' (%d %d %d) ', LU.LWH(1:nDim,idxLU(x)));fprintf('\n');
% %                 % %     fprintf(' [%d]     ', LU.LURotaFlag(:,idxLU(x)));fprintf('\n');
% %                 % %     fprintf(' {%d %d %d} ', LU.CoordLUBin(:,idxLU(x)));fprintf('\n');
% %                 % %     fprintf('\n');
% %             end
% %         end

%% ������Ҫע��
        % %         %  *******  ����д���,�ͱ���; 
        %         ti = d.Item;  tl = d.LU;
        %         [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip);
        %         [match, er1, er2] = comp_struct(ti,d.Item,1);
        %         list_struct(er1)        
        %         if isequal(ti,d.Item)~=1, warning('eeeee'); end
        %         if isequal(tl,d.LU)~=1, warning('eeeee'); end
        
                        %         pgon = getPolyshape(d.LU.LWH);    maxX = sum(d.LU.LWH(1,:))+10;    maxY = max(max(d.LU.LWH'))+10;  maxV = max(maxX,maxY);
                %                    plot(pgon);        axis equal;   axis ([0 maxX 0 maxY]);   
                %                    plot3Dshape(d.LU.LWH);
                
                                        %  pgon = getPolyshape(d.Item.LWH);    % figure; plot(pgon);  axis equal;  axis ([0 maxX 0 maxY]);
                                        
                                        % �����½�
                                        %         lb = computerLB(d.Item,d.Veh);   fprintf('LB = %d \n', lb); %��ĳ��bin����Ϊ׼
                                        
                                                % ��Strip�н���һ���Ҹ�>���Item����ѡ�񲢸�����Ӧ����
                                                %         d = modifyStripWithOneItem(d);
                                                %         function d = modifyStripWithOneItem(d)
                                                %             stripheight = d.Strip.LW(2,:);
                                                %             binwidth = d.Veh.LWH(1,1);
                                                %             stripleftwidth = d.Strip.LW(1,:);
                                                %             stripwidth = ( binwidth - stripleftwidth );
                                                %             [tmpset] = find(stripheight > stripwidth);
                                                %             if ~isempty(tmpset)
                                                %                 if isscalar(tmpset) %�Ը�strip�����ڲ�����1��Item����,��������漰CoordItemStrip
                                                %                     d.Strip.LW(:,tmpset) = [binwidth-stripheight(tmpset),stripwidth(tmpset)];    %strip�ĳ������
                                                %                     %�ڲ�Item��itemRotaFlag���� 
                                                %                     idxItem = find(d.Item.Item_Strip(1,:)==tmpset );
                                                %                     if isscalar(idxItem)
                                                %                         d.Item.itemRotaFlag(idxItem) = ~d.Item.Rotaed(idxItem);
                                                %                     end                    
                                                %                     %�ڲ�LU��LURotaFlag ���{ ߀δ�� %�ڲ�Item��CoordItemStrip���{                    
                                                %                 end
                                                %             end
                                                %         end
                                                
                                                                % ********** ƽ���Ż� **********
                                                                % %         if any(d.Strip.isWidthFull)
                                                                % %             fw = d.Strip.isWidthFull == 0;
                                                                % %             fm = d.Strip.isMixed == 0;
                                                                % %             %             fi = find(fw);
                                                                % %             nbbin = numel(d.Bin.Weight);
                                                                % %             %    for i=1:sum(fi)
                                                                % %             for i=1:nbbin
                                                                % %                 f2 = d.Strip.Strip_Bin(1,:) == i;
                                                                % %                 f3 = f2 & fw & fm;   %����bin i��, ����isWidthFull��strip            
                                                                % % %                 f2 = d.Strip.Strip_Bin(1,:) == d.Strip.Strip_Bin(1,fi(i)); % �뱾fi(i)����ͬһbin�ڵ�isWidthFull�ĸ���
                                                                % % %                 f3 = f2 & f;                
                                                                % %                 if sum(f3) > 1                    
                                                                % %                     error('�뱾fi(i)����ͬһbin�ڵ�isWidthFull�ĸ���>1');
                                                                % %                 end
                                                                % %                 fi3 = find(f3);
                                                                % %                 
                                                                % %                 d.Strip
                                                                % %                 d.Strip.LW
                                                                % %                 
                                                                % %                 % ��Ҫƽ�̵�strip���ڵ�item��lu��logicalֵ
                                                                % %                 fitem =d.Item.Item_Strip(1,:) == fi3;
                                                                % %                 flu =d.LU.LU_Strip(1,:) == fi3;
                                                                % %                 d.LU
                                                                % %                 d.Item
                                                                % %                 1
                                                                % %                 % �޸�LU.CoordLUBin��ֵ
                                                                % %                 % �޸�LU.LU_Item�仯�Ķ��ֵ
                                                                % %                 % �޸�Item.LWH��ֵ
                                                                % %                 % �޸�Item.isFull1��ֵ,Weight,isHeightFull,isWeightFine,Layer,PID,LID,SID
                                                                % %                 % CoordItemBin, CoordItemStrip, Item_Strip(���µ�����仯)
                                                                % %                 
                                                                % %                 %             d.Strip
                                                                % %                 d.Strip.isMixed
                                                                % %                 d.Strip.isAllPured
                                                                % %                 d.Strip.maxHeight
                                                                % %             end    
                                                                % %         1
                                                                % %         end

                                                                

%% �����½�
% % % % %     function [ lb ] = computerLB(Item,Veh)
% % % % %         sum1 = sum(prod(Item.LWH,1));        
% % % % %         % todo �����ж��Ƿ����е�BinArray�����е�bin����ͬ�� ����� �����ִ��
% % % % %         sum2 = prod(Veh.LWH(:,1));
% % % % %         lb = ceil(sum1/sum2);
% % % % %         if lb <=0, error('EEE');end
% % % % %     end
    
    
%%     d.LU = doc(d.LU,d.Item,d.Strip,d.Bin);
% % % %     function LU = doc(LU,Item,Strip,Bin)
% % % %         
% % % %         LU.DOC=[LU.PID;LU.ID;LU.SID;zeros(size(LU.ID));zeros(size(LU.ID));...
% % % %             LU.LU_Item;LU.LU_Strip;LU.LU_Bin];
% % % %         
% % % %         nItem = size(Item.LWH,2);
% % % %         for iItem=1:nItem
% % % %             tmp = LU.DOC([1,2,3], LU.DOC(6,:) == iItem);
% % % %             Item.PID1(:,iItem) = num2cell(unique(tmp(1,:))',1); %unique(tmp(1,:))';
% % % %             Item.LID1(:,iItem) = num2cell(unique(tmp(2,:))',1);
% % % %             Item.SID1(:,iItem) =num2cell(unique(tmp(3,:))',1);
% % % %         end
% % % %         
% % % %         nStrip = size(Strip.LW,2);
% % % %         for iStrip=1:nStrip
% % % %             tmp = LU.DOC([1,2,3], LU.DOC(8,:) == iStrip);
% % % %             Strip.PID1(:,iStrip) = num2cell(unique(tmp(1,:))',1);
% % % %             Strip.LID1(:,iStrip) = num2cell(unique(tmp(2,:))',1);
% % % %             Strip.SID1(:,iStrip) = num2cell(unique(tmp(3,:))',1);
% % % %         end
% % % %         
% % % %         nBin = size(Bin.LW,2);
% % % %         for iBin=1:nBin
% % % %             tmp = LU.DOC([1,2,3], LU.DOC(10,:) == iBin);
% % % %             Bin.PID1(:,iBin) = num2cell(unique(tmp(1,:))',1);
% % % %             Bin.LID1(:,iBin) = num2cell(unique(tmp(2,:))',1);
% % % %             Bin.SID1(:,iBin) = num2cell(unique(tmp(3,:))',1);
% % % %         end
% % % %         
% % % %         LU.DOC
% % % %         LU.DOC([1,2,3,8],:)
% % % %     end
