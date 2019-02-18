function [d] = RunAlgorithm(d,p)
        global ISshuaiwei ISreStripToBin ISstripbalance %����ִ��HStripToBin����
        global  ISplotshuaiwei  ISplotStripToBinAgain  ISplotstripbalance ISplotRunAlgo ISplotStrip

        %% Ԥ����:����Input��������
        fprintf(1,'   Running RunAlgorithm  ...\n');
        %         d = GcheckInput(d);    %���Բ���   ��Ϊmain��һ�¼���          
        
        % ����Ԥ�����ص㣺��ȡLU.Rotaed,�����Ƿ����� ����cpuLU���� ���������LU������  ������֮ǰLWH�Ѿ�ʱ����margin�ģ�        
        
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        [d.LU, d.Veh] = cpuVehLU(d.LU, d.Veh);            %������ in case �������
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        

        %% ����ʽ: LU��Item���㷨    
        % ����LU��isNonMixed��LU��isMixedTile ���� ��;�ƺ�����, Ŀ��Ϊ������Itemʱ,�����һ��,������ĺ���
         
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % ADD: isNonmixed/ismixedtile
        [d.LU] = cpuLU(d.LU,d.Veh);     
        
        
        %fprintf(1,'     Running HLUtoItem...\n');
        
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % ADD: order/LU_Item
        % Item:
        % ADD: isRota/Rotaed/HLayer/LWH/Weight
        [d.LU,d.Item] = HLUtoItem(d.LU,d.Veh);   %Item����ID������򣨵���һ�������仯˳��
        
                                %printstruct(d.LU,'sortfields',1,'PRINTCONTENTS',1);
%             plotSolutionT(d.LU,d.Veh,1,0,0,0,3,'ԭ˳��LU'); %LUԭʼ
%             plotSolutionT(d.LU,d.Veh,3,0,0,0,3,'�����LU'); %LU�����

        % SECTION 2 ����LU_Item���Ƿ���������(�޸�˳�����,�����ʱ��δ����)
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
    
        % ����Item�� isNonMixed isMixedTile isHeightFull ��������
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % ADD: /isHeightFull/isNonMixed/isMixedTile/nbItem/IDs
        [d.Item] = cpuItem(d.Item,d.LU,d.Veh);        % printstruct(d,'sortfields',1,'PRINTCONTENTS',0);
        %        printstruct(d.LU,'sortfields',1,'PRINTCONTENTS',1);

        %% ����ʽ: Item��Strip���㷨        
        %  *******  *******  *******
        % 1 Item����: %1: SID ; 2: isNonMixed; һ�����濪ʼ: 
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
        % isHeightFull/isNonMixed/isMixedTile/nbItem/IDs
        % ADD: itemorder/Item_Strip/CoordItemStrip
        % Strip:
        % ADD: LW/Weight
        [d.Item,d.Strip] = HItemToStrip(d.LU,d.Item,d.Veh,p);     %   printstruct(d);   %  printstruct(d.Item);

        % Item.H HLyaer Weight isHeightFull LU_Item �ȱ仯
        % SECTION 3 ����LU_Strip���Ƿ����Ҳ��� (�޸�˳�����,�����ʱ��δ����)
        %         d.LU.LU_Item(2,:)= repairStrips(d.LU); % �޸�LU_Strip��˳������
        % �˲��Ƿ����������ص����
        %         chktLU(d.LU)

      
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % ADD: .LU_Strip, LU.CoordLUStrip
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/isNonMixed/isMixedTile/nbItem/IDs
        % itemorder/Item_Strip/CoordItemStrip
        % Strip:
        % LW/Weight
        % ADD: isMixedIDs/nbIDs/isHeightFull/isHeightBalance/isWidthFull/nbItem/nbLU/nbLULID
        % ADD: isAllPured/isSingleItem/mlmdHeight/IDs/~Limit
        [d.Strip,d.LU] = cpuStrip(d.Strip,d.Item,d.LU,d.Veh);

        if ISplotRunAlgo & 0
%             plotSolutionT(d.LU,d.Veh,1,0,0,0,3,'ԭ˳��LU'); %LUԭʼ
            plotSolutionT(d.LU,d.Veh,3,0,0,0,3,'�����LU'); %LU�����
            plotSolutionT(d.LU,d.Veh,0,1,0,0,3,'�����Item'); %Item�����
            plotSolutionT(d.LU,d.Veh,0,0,1,0,3,'�����Strip'); %Strip�����
        end
        % ISplotRunAlgo = 1
        %% ����ʽ�����⣺�Ѷ��������
        % ********* 1 5555 �Ѷ�������� *********** �޸� d.LU.maxHLayer(luidxPP) 
%          if ISstripbalance==0
%              fprintf(1,'     Running HStripBalance...\n');
%              [d.Strip,d.Item,d.LU,TFStripBalance] = HStripBalance(d.Strip,d.Item,d.LU,d.Veh,p);  end
%          if ISstripbalance && TFStripBalance && ISplotstripbalance && ISplotRunAlgo
%              plotSolutionT(d.LU,d.Veh,0,0,1,0,3,'�����Strip'); %Strip:HStripBalance��
%          end     % V1:  if ISplotStrip==1,      figure(222);      plot3DStrip(d.LU,d.Item,d.Veh,'LU');        end    %  plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % ����LU.CoordLUBin %   figure(222);     plot3DStrip(d.LU,d.Item,d.Veh,'LU');         %    igure(111);          plot3DStrip(d.LU,d.Item,d.Veh,'Item'); 

        %% ����ʽ��Strip��Bin���㷨        
        %% ********* 1 Strip����: % 1: SID 2: priorityofLID
        %fprintf(1,'     Running HStripToBin...\n');
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % ADD: .LU_Strip, LU.CoordLUStrip
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/isNonMixed/isMixedTile/nbItem/IDs
        % itemorder/Item_Strip/CoordItemStrip
        % Strip:
        % LW/Weight
        % isMixedIDs/nbIDs/isHeightFull/isHeightBalance/isWidthFull/nbItem/nbLU/nbLULID
        % isAllPured/isSingleItem/mlmdHeight/IDs/~Limit
        % ADD: striporder/strip_bin
        % Bin:
        % ADD: LW/Weight
        [d.Strip,d.Bin] = HStripToBin(d.Strip,d.Veh,p);
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % LU_Strip, LU.CoordLUStrip
        % ADD: LU_Bin/CoordLUBin/ lbrt_A
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/isNonMixed/isMixedTile/nbItem/IDs
        % itemorder/Item_Strip/CoordItemStrip
        % ADD: Item_bin/CoordItemBin
        % Strip:
        % LW/Weight
        % isMixedIDs/nbIDs/isHeightFull/isHeightBalance/isWidthFull/nbItem/nbLU/nbLULID
        % isAllPured/isSingleItem/mlmdHeight/IDs/~Limit
        % striporder/strip_bin
        % Bin:
        % LW/Weight
        [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip); % ����LU��Bin������and˳��   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin

        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % LU_Strip, LU.CoordLUStrip
        % LU_Bin/CoordLUBin/ lbrt_A
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/isNonMixed/isMixedTile/nbItem/IDs
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
        [d.Bin] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %����Bin��������� % ����isTileNeed

        list_struct(d.LU)
        list_struct(d.Bin)        
        
        if ISplotRunAlgo
            plotSolutionT(d.LU,d.Veh,0,0,0,1,3,'�����Bin'); % Bin�����
        end         % if ISplotStrip==1,      plot3DBPP(d,p);      end    % igure(111);   plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % ����LU.CoordLUBin 
                                                    
        %% ���⣺2 ����ͷ����2: ÿ��ʣ��stripȫ���ڱȽ��� better than ����1 �й���
        % �ر��ǶԵ����������Ժ�
                                                        %         ti = d.Strip;  tl = d.Bin;
        if ISreStripToBin==1
                fprintf(1,'     Running HStripToBinAgain...\n');
                [d.Strip,d.Bin,TFHStripToBinAgain] = HStripToBinAgain(d.Bin,d.Strip,d.Item,d.LU,d.Veh,p); 
                
                [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip); % ����LU��Bin������and˳��   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin

                [d.Bin] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %����Bin��������� % ����isTileNeed

                if TFHStripToBinAgain &&  ISplotStripToBinAgain && ISplotRunAlgo
%                     plotSolutionT(d.LU,d.Veh,0,0,0,1,3,'����ͷ��Bin'); % Bin�����
                end %if ISplotStrip==1,      plot3DBPP(d,p);      end    % igure(111);   plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % ����LU.CoordLUBin
        end     % ����ͷ����1: ÿ��Bin��strip�Ƚ���      % [d.Strip,d.Bin]= HreStripToEachBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh,p);
        
        %         close all;  plotSolutionT(d.LU,d.Veh);
    
                                                        %         [match, er1, er2] = comp_struct(ti,d.Strip,1);
                                                        %         [match, er1, er2] = comp_struct(tl,d.Bin,1);
                                                        %         list_struct(er1)
        
         
                                                        
        %% ���⣺3 ����Strip��˦β�Ż� *********** �޸� Strip.Strip_Bin �ڶ���
        if ISshuaiwei==1   
                fprintf(1,'     Running HStripSW...\n');


                % ����strip��strip_bin˳��, strip�Ƿ�˦β,strip˦β˳�� 
%                 [d.Strip,d.LU.isShuaiWei,TFHStripSW] = HStripSW(d.Strip,d.LU);       
                [d.Strip.Strip_Bin(2,:),  d.Strip.isShuaiWei, d.LU.isShuaiWei, TFHStripSW] = HStripSW(d.Strip,d.LU);  

                % ����LU/Item ��Bin������and˳��LU_Bin  NOTE: Lu_Strip and Item_Strip ����˦βӰ��
                [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip);    %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin
       
%                 plotSolutionT(d.LU,d.Veh,0,0,0,1,3,'˦β��Bin');
%                 plotSolutionT(d.LU,d.Veh,0,0,0,1,8,'˦β��Bin');
                
                [d.Bin] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %����Bin��isTileNeed

                if TFHStripSW && ISplotshuaiwei && ISplotRunAlgo
                    plotSolutionT(d.LU,d.Veh,0,0,0,1,3,'˦β��Bin'); % Bin����� 
                end %if ISplotStrip==1,      plot3DBPP(d,p);      end    % igure(111);   plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % ����LU.CoordLUBin
                
        end

       %% ��ֵ����
       d.LU.LU_VehType = ones(size(d.LU.ID)) * d.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ
       fprintf(1,'   Running RunAlgorithm  DONE ...\n');
                %         printstruct(d.Item,'sortfields',1,'PRINTCONTENTS',0)  
                %         printstruct(d.Bin,'sortfields',1,'PRINTCONTENTS',1)  
                %         printstruct(d,'sortfields',1,'PRINTCONTENTS',0)  

end

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
