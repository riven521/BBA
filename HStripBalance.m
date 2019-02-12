%% ����: Strip����
function   [Strip,Item,LU] = HStripBalance(Strip,Item,LU,Veh,p)
        flagBalanced = 0;
        while 1
            % �ҵ�strip����ϵ��Ҳ�ƽ���
            fstrip = Strip.isMixed==0 & Strip.isHeightBalance==0;
            if ~any(fstrip) || flagBalanced==1, break; end
            idxstr = find(fstrip);
            nbstrip = numel(idxstr);
            
            for idx=1:nbstrip
                luidxPP = ismember(LU.LU_Strip(1,:), idxstr(idx));
                
                % FINDOUT ��strip��ӦItem��LU����,ͨ��Item
                tmpItems = LU.LU_Item(1,luidxPP);
                [~,~,ic] = unique(tmpItems);
                nbLuArray= accumarray(ic,1);
                
                if length(unique(nbLuArray)) >2
                                %                     nbLuArray
                                %                     diff(nbLuArray)
                                                    %          plotSolutionT(LU,d.Veh);  %???
                                                    %          Ϊ�λ���nbLuArrayΪ���� 
                                                    %          Ϊ�λ���Strip��ֻ��1��Item�Ĳ����? ��ҪCHECK
                    error('����unblanced�ĶѶ�С��2��Ψһ������2����Ԥ�ڴ���');        
                end
                if isscalar(nbLuArray) 
                    warning('nbLuArray����1��ֵ,�������');
                end
                
               
                if  isscalar(nbLuArray)                         %��stripû��balance�ı�Ҫ ֻ��1��Item
                    flagBalanced = 1;
                    continue;
                elseif abs(diff(nbLuArray)) < 2          %��stripû��balance�ı�Ҫ �������2��
                    flagBalanced = 1;
                    continue;
                else
                    % ֱ�ӽ���Ŀǰ���Item�ѷŵĲ�����������-1
                    LU.maxHLayer(luidxPP) = max(nbLuArray);
                end

                LU.maxHLayer(luidxPP) = min( max(LU.maxL(3,luidxPP)), max(LU.maxHLayer(luidxPP))) - 1;

                % �����Ͳ�����0 �� ���͵�ԭ��LU����С����, �����ܴ�����, ������strip�ĸ߶Ⱦ���
                if any(LU.maxHLayer(luidxPP)<=0) || any(LU.maxHLayer(luidxPP) <= min(nbLuArray))
                    LU.maxHLayer(LU.maxHLayer(luidxPP)<=0) = 1;
                    flagBalanced = 1;
                    continue;
                end
                
                % 555 �Խ��Ͳ�����LU���»�ȡ�µ�Strip
                %   [LU,Veh] = cpuLUVeh(LU,Veh);
                [LU] = cpuLU(LU,Veh);
                [LU,Item] = HLUtoItem(LU,Veh);          %Item����ID������򣨵���һ�������仯˳��
                [Item,LU] = cpuItem(Item,LU,Veh);        % printstruct(d,'sortfields',1,'PRINTCONTENTS',0);
                [LU,Item,Strip] = HItemToStrip(LU,Item,Veh,p);
                [Strip,LU] = cpuStrip(Strip,Item,LU,Veh);
                % fstrip = Strip.isMixed==0 & Strip.isHeightBalance==0;
                flagBalanced = 0;
            end
        end           
end
