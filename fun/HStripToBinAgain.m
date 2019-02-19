%% ����1: ����ͷ����2: V3: 
    % s1 :HStripToBin����ʱ:Ϊ����Strip�ų��Ѱ���Bin���ʣ��Strip; cpuStripnbItemΪȫ��Strip����
function [Strip,Bin,TF] = HStripToBinAgain(Bin,Strip,Item,LU,Veh,p)
        % Ŀ��: ��������LU���������Ϊ��С���԰ڷų�ͷ;
        % ����: ****** ÿ��Bin�����ų�ǰ���Ѱ���Strip���ʣ��Strip, �������򲢷ֱ�ִ��HStripToBin�㷨. ******
        nbBin = max(Strip.Strip_Bin(1,:));
        TF = false;
        if nbBin>1
            fprintf(1,'       Exsiting ����ͷ in HStripToBinAgain (nBin>1)...\n');
            ibin=2;
            while 1
                % 1 ��ȡf and fidx: �ų��׸�bin���ʣ��Strip�߼�ֵ and ����ֵ
                Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) >= ibin ) = -1; % ���з���bin��bin��Ÿ�ֵ-1, ��һ��bin��Զ���ᶯ
                % �ҳ���Ҫ������������strip,lu,item�ȵ�flag
                fStrip = Strip.Strip_Bin(1,:) == -1; 
                fidx=find(fStrip);    if ~any(fStrip), break; end
                
                % V1 ����ԭ����LU.LU_Bin�ᷢ���仯, ���˴����ǲ���
                %                 fItem = Item.Item_Bin(1,:) >= ibin ; %
                %                 fLu = LU.LU_Bin(1,:) >= ibin ;    

                % V2 �޸������LU.LU_Bin�仯���Ż�
                Item.Item_Bin(1, Item.Item_Bin(1,:) >= ibin |  Item.Item_Bin(1,:) == 0) = -1;
                fItem = Item.Item_Bin(1,:) == -1;
                                
                LU.LU_Bin(1, LU.LU_Bin(1,:) >= ibin |  LU.LU_Bin(1,:) == 0) = -1; 
                fLu = LU.LU_Bin(1,:) == -1;

                % 2 ��ȡs1,i1,l1; �ֱ�����λ��ʣ��bin�ڵı��
                s1 = rmfield(Strip,{'striporder','Strip_Bin','nbLU','nbItem'});   s1.f = fStrip; %Strip
                i1  = Item;  i1.f = fItem;
                l1 = LU;      l1.f = fLu;

                % 3 ����cpuStripnbItem����ȫ��,������flag���
                [s1.nbItem,s1.nbLU, s1.nbLULID] = cpuStripnbItem(s1,i1,l1);           
                
%                 [s1.isMixed,s1.isMixedSID,s1.isMixedEID] = isMixedStrip(s1);  % ������Ҳû���ã�
                
                % 4 ��ȡ��s2��b2 and ����ִ������ʽS2B (����HStripToBinֻ�ܸ�ʣ�ಿ�ֵ� x(:,f) )
                s1 = structfun(@(x) x(:,fStrip),s1,'UniformOutput',false);    % �ų��Ѿ�ȷ����bin��ʣ��Strip
                [s2,b2]= HStripToBin(s1,Veh,p);
        
                % 5 ��ȡfs: ʣ��Strip�ڷź���׸�bin�ڵ� so �Ҳ�==1
                % �滻ԭʼStrip_Bin��ֵ
                fs = s2.Strip_Bin(1,:)==1; %�ҳ���i1��bin���߼�ֵ                
                Strip.Strip_Bin(1,fidx(fs)) = ibin;
                Strip.Strip_Bin(2,fidx(fs)) = s2.Strip_Bin(2,fs);
                
                % 6 Bin���¸�ֵ���
                Bin.Weight(ibin) = b2.Weight(1);
                Bin.LW(:,ibin) = b2.LW(1);

                % 7 ��Ҫ���Ӷ�LU_Bin�ļ���,
                [LU,Item] = HItemToBin(LU,Item,Strip); % ����LU��Bin������and˳��   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin
                [Bin] = cpuBin(Bin,Strip,Item,LU,Veh);  %����Bin��������� % ����isTileNeed
                
                % 8 ibin����
                ibin = ibin+1;
                TF = true;
            end
        end
end

                                                                
%% ����1: ����ͷ����2: V1: s1 :ȫ��Ϊ����Strip�ų��Ѱ���Bin���ʣ��Strip;;cpuStripnbItemΪ����Strip����
% % function [Strip,Bin] = HStripToBinAgain(Bin,Strip,Item,LU,Veh,p)
% %         % DONE: ����ͷ����2: 
% %         % Ŀ��: �����LU���������Ϊ��С���԰ڷų�ͷ;
% %         % ����: ÿ��Bin�����ų�ǰ��Strip���ʣ��Strip�������򲢷ֱ�ִ��HStripToBin�㷨.
% %         nbBin = max(Strip.Strip_Bin(1,:));
% %         if nbBin>1
% %             ibin=2;
% %             while 1
% %                 % 1 ��ȡf and fidx: �ų��׸�bin���ʣ��Strip�߼�ֵ and ����ֵ
% %                 Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) >= ibin ) = -1; % ���з���bin��bin��Ÿ�ֵ-1
% %                 % �ҳ�δ��ȷ������Strip
% %                 f = Strip.Strip_Bin(1,:) == -1; 
% %                 fidx = find(f);
% %                 if ~any(f), break; end
% %                 
% %                 % 2 ��ȡs1 : ����f ��ȡ��ʣ��Strip�ṹ��
% %                 s1 = rmfield(Strip,{'striporder','Strip_Bin'});
% %                 s1 = structfun(@(x) x(:,f),s1,'UniformOutput',false);
% %                 
% %                 % 3 ��ȡi1 : ����ʣ��Strip��ȡ��ʣ��Item�ṹ��
% %                 i1  = Item;
% %                 f2 = zeros(1,length(i1.Weight));
% %                 iIdxs= find(f);
% %                 for j=1:length(iIdxs)
% %                     f2 = f2+(i1.Item_Strip(1,:) == iIdxs(j));
% %                 end                
% %                 i1 = structfun(@(x) x(:,logical(f2)), i1, 'UniformOutput',false);
% %                 
% %                 % 4 ��ȡ��s2��b2 : ����s1��i1, ���¼���Strip.nbItem and ����ִ������ʽS2B
% %                 [s1] = cpuStripnbItem(s1,i1,LU);
% %                 %[s1,~] = cpuStrip(s1,i1,LU,Veh);  % �˺�������������Ӻ���, �Ƿ��б�Ҫ����ȫ��, ����ֻ���������Ӻ��� TODO
% %                 
% %                 [s2,b2]= HStripToBin(s1,Veh,LU,p);
% %                 
% %                 % 5 ��ȡfs: ʣ��Strip�ڷź���׸�bin�ڵ� so �Ҳ�==1
% %                 % �滻ԭʼStrip_Bin��ֵ
% %                 fs = s2.Strip_Bin(1,:)==1;
% %                 Strip.Strip_Bin(1,fidx(fs)) = ibin;
% %                 Strip.Strip_Bin(2,fidx(fs)) = s2.Strip_Bin(2,fs);
% %   
% %                 % 6 ��ֵ���
% %                 Bin.Weight(ibin) = b2.Weight(1);
% %                 Bin.LW(:,ibin) = b2.LW(1);
% %                 ibin = ibin+1;                
% %             end
% %         end
% % end


%% %% ����1: ����ͷ����2: V2: s1 :HStripToBin����ʱ:Ϊ����Strip�ų��Ѱ���Bin���ʣ��Strip;cpuStripnbItemΪȫ��Strip����
% % function [Strip,Bin] = HStripToBinAgain(Bin,Strip,Item,LU,Veh,p)
% %         % DONE: ����ͷ����2: 
% %         % Ŀ��: ��������LU���������Ϊ��С���԰ڷų�ͷ;
% %         % ����: ÿ��Bin�����ų�ǰ���Ѱ���Strip���ʣ��Strip, �������򲢷ֱ�ִ��HStripToBin�㷨.
% %         nbBin = max(Strip.Strip_Bin(1,:));
% %         if nbBin>1
% %             ibin=2;
% %             while 1
% %                 % 1 ��ȡf and fidx: �ų��׸�bin���ʣ��Strip�߼�ֵ and ����ֵ
% %                 Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) == 1 ) = -1; 
% %                 Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) >= ibin ) = 0; % ���з���bin��bin��Ÿ�ֵ-1
% %                 % �ҳ�δ��ȷ������Strip
% %                 f = Strip.Strip_Bin(1,:) == 0; 
% %                 fidx=find(f);
% %                 if ~any(f), break; end
% %                 
% %                 % 2 ��ȡs1 : �ų�����field,����f���½ṹ��
% %                 s1 = rmfield(Strip,{'striporder','Strip_Bin'});
% %                 s1.f = f;
% %                 
% %                 % 3 ��ȡi1 : ����ʣ��Strip��ȡ��ʣ��Item�ṹ��
% %                 i1  = Item;
% %                 % 4 ��ȡ��s2��b2 : ����s1��i1, ���¼���Strip.nbItem and ����ִ������ʽS2B
% %                 [s1] = cpuStripnbItem(s1,i1,LU);
% %                 %[s1,~] = cpuStrip(s1,i1,LU,Veh);  % �˺�������������Ӻ���, �Ƿ��б�Ҫ����ȫ��, ����ֻ���������Ӻ��� TODO
% %                 
% %                 s1 = structfun(@(x) x(:,f),s1,'UniformOutput',false); %����HStripToBinֻ�ܸ����ֵ� x(:,f)
% %                 [s2,b2]= HStripToBin(s1,Veh,LU,p);
% %                 
% %                 % 5 ��ȡfs: ʣ��Strip�ڷź���׸�bin�ڵ� so �Ҳ�==1
% %                 % �滻ԭʼStrip_Bin��ֵ
% %                 fs = s2.Strip_Bin(1,:)==1; %�ҳ���i1��bin���߼�ֵ                
% %                Strip.Strip_Bin(1,fidx(fs)) = ibin;
% %                Strip.Strip_Bin(2,fidx(fs)) = s2.Strip_Bin(2,fs);
% % 
% %                 % 6 ��ֵ���
% %                 Bin.Weight(ibin) = b2.Weight(1);
% %                 Bin.LW(:,ibin) = b2.LW(1);
% %                 ibin = ibin+1;                
% %             end
% %             Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) == -1 ) = 1;             
% %             % sortrows(Strip.Strip_Bin',[1,2],{'ascend','ascend'})'
% %         end
% % end

        
%% ����2: ����ͷ����1: 
        % Ŀ��: �����LU���������Ϊ��С���԰ڷų�ͷ;
        % ����: ÿ��Bin���ֱ𵫸�Bin��Stripִ��HStripToBin�㷨.
% % function [Strip,Bin] = HreStripToEachBin(Bin,Strip,Item,LU,Veh,p)
% %        % ��Ժ��������������򲢰��ŵ���Bin��
% %         nbBin = max(Strip.Strip_Bin(1,:));
% %         if nbBin>100
% %             for ibin=2:nbBin
% %                 f = Strip.Strip_Bin(1,:) == ibin;
% %                 
% %                 s1 = rmfield(Strip,{'striporder','Strip_Bin'});
% %                 s1 = structfun(@(x) x(:,f),s1,'UniformOutput',false);
% %                 i1  = Item;  
% %                 
% %                 iIdxs= find(f);
% %                 f2 = zeros(1,length(i1.Weight));
% %                 for j=1:length(iIdxs)
% %                     f2 = f2+(i1.Item_Strip(1,:) == iIdxs(j));
% %                 end
% %                 
% %                 i1 = structfun(@(x) x(:,logical(f2)), i1, 'UniformOutput',false);
% %                 
% %                 % 1 ���¼���Strip.nbItem and ����ִ������ʽS2B
% %                 [s1] = cpuStripnbItem(s1,i1,Veh);
% %                 
% %                 [s2,b2]= HStripToBin(s1,Veh,LU,p);
% %                 
% %                 % 2 �滻��ԭʼStrip_Bin��ibin�Ľ���˳��
% %                 Strip.Strip_Bin(2,f) = s2.Strip_Bin(2,:);
% %                 
% %                 % 3 �������
% %                 if any(s2.Strip_Bin(1,:)~=1), error('���·����bin�Ų���'); end
% %                 if size(b2.Weight,2) ~=1, error('��������ͬ'); end
% %                 if size(b2.LW,2) ~=1, error('������ͬ'); end
% %                 if Bin.Weight(ibin) ~= b2.Weight, error('��������ͬ'); end
% %                 if Bin.LW(:,ibin) ~= b2.LW, error('������ͬ'); end                
% %                 Bin.Weight(ibin) = b2.Weight;
% %                 Bin.LW(:,ibin) = b2.LW;
% %                 
% %             end
% %         end
% % end
        