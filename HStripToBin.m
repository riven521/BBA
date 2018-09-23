function [Strip,Bin]= HStripToBin(Strip,Veh,LU,p)
% ��Ҫ����:Strip����Bin�� %  ����:�����(row);  ����:��������(coloum);
% Input ---  Strip/Veh:  
% Output --- Strip: 
% Output --- Bin: 

%% ��ʼ��
sz = size(Strip.LW);
nStrip = sz(2);

wVeh  = Veh.LWH(1,1); 
lVeh  = Veh.LWH(2,1); 

%% Strip���� 555 (���ȷ����ͬ�������͵����ڰڷţ����� TODO FIX ME)
    % ��ȡStrip��˳��(�ص���Strip�߶ȵݼ����򣨵���������strip�߶�һ���ģ�) %
    [Strip.striporder] = getStriporder(Strip);  % Strip������ʽ �߶�/���ȵݼ�
    % ��ȡ��order������Strip: sStrip    
    if isSameCol(Strip)
        sStrip = structfun(@(x) x(:,Strip.striporder),Strip,'UniformOutput',false);
    else
        error('����ʹ��structfun');
    end
    
   
%% LU->Item->Strip->Binת�� 
% ��ȡstripBeBinMatrixSort: ÿ�������strip���ĸ�bin��  �Լ�˳��
% ��ȡLWBin:  �����ɵ�Bin��ʣ�೤��
Bin.LW = zeros(2,nStrip);    %��ʼ��bin: dim1-bin���ʣ�� ; dim2-bin��(��)��(555ʣ�ࣩ;
Bin.LW(1,:) = wVeh;
Bin.LW(2,:) = lVeh;
Bin.Weight = zeros(1,nStrip); % ��ʼ��ֵ

    % ��ʼ������nItem��
%     Bin.LID = zeros(numel(unique(LU.ID)),nStrip);
%     Bin.PID = zeros(numel(unique(LU.PID)),nStrip);
%     Bin.SID = zeros(numel(unique(LU.SID)),nStrip);
%     Bin.UID = zeros(numel(unique(LU.UID)),nStrip);
    
tmpBin_Strip = zeros(1,nStrip);    % ÿ��Bin�ڵ�Strip���� ���ڲ���
% sStrip����
sStrip.Strip_Bin = zeros(2,nStrip); % dim1:��� strip��ĳ��bin dim2:����˳�� 555

% LWStripSort = sStrip.LW; %��sorted���
% StripWeightSort = sStrip.Weight;

% 55 ��ȡthisBin - ��ǰstripҪ�����bin���
% ѭ����bin�а���strip,���̶�strip,�仯ѡ��ͬbin(thisBin)
% ע�ͣ���ȡ FLAG        �ɷ��µ�ǰiStrip������һ��bin�ļ��� 
% ע�ͣ���ȡ thisBin   ��FLAG���ҵ���������Ǹ�thisBin, ��ִ�� insert����

iStrip=1; iBin=1;
while 1
    if iStrip > nStrip, break; end
    
    [thisBin,iBin] = getThisBin(iBin, iStrip, sStrip, Veh, Bin, p);    % ��ȡBin��
    
    [Bin,sStrip.Strip_Bin,tmpBin_Strip] = insertStripToBin(iStrip, thisBin, sStrip, Bin, sStrip.Strip_Bin, tmpBin_Strip);
        
    iStrip = iStrip + 1;
end

% plot2DBin();


    
% ���� ����ֵ��da
% ��ȡStrip_Bin: ÿ��strip���ĸ�bin��  �Լ�˳��
                % Strip.Strip_Bin( : , Strip.striporder) = sStrip.Strip_Bin;
% Strip�ڲ�����,sStrip����order�仯����
if isSameCol(sStrip)
    Strip = reorderStruct(Strip.striporder, sStrip);
else
    error('����ʹ��structfun');
end



% �ɻ�ϵ�LU.DOC����LU_BIN, ����BIN�ڰ�����PID,LID,SID������ 1808����
    nbLU = size(LU.LWH,2);
    LU.LU_Bin = [zeros(1,nbLU);zeros(1,nbLU)];
    for iLU=1:nbLU
         theStrip = LU.LU_Strip(1,iLU); %iLU���ڵڼ���Item
         LU.LU_Bin(1,iLU)= Strip.Strip_Bin(1,theStrip);
    end

    LU.DOC=[LU.DOC; LU.LU_Bin];
    nBin = size(Bin.LW,2);
    for iBin=1:nBin
        tmp = LU.DOC([1,2,3], LU.DOC(10,:) == iBin);
        Bin.PID2(:,iBin) = num2cell(unique(tmp(1,:))',1);
        Bin.LID2(:,iBin) = num2cell(unique(tmp(2,:))',1);
        Bin.SID2(:,iBin) = num2cell(unique(tmp(3,:))',1);
    end
    
    
% ��ȡBin: ȥ��δʹ�õ�Bin ע��Bin�ṹ��ı仯
if isSameCol(Bin)
    Bin = structfun(@(x) x( : , Bin.Weight(1,:)>0 ), Bin, 'UniformOutput', false);
else
    error('����ʹ��structfun');
end

% ITEM�����ж��Ƿ��������ص��ж�Item.isWeightFine
Strip.Strip_Bin
Strip.isFull
Strip.isSingleItem

% ����ڵ���Item��strip��case �� ����strip�в�����strip
if any(Strip.isSingleItem | ~Strip.isFull )
    [~,bsingle] = find(Strip.isSingleItem == 1);
    [~,bnotfull] = find(Strip.isFull == 0);
    b = unique([bnotfull bsingle],'stable');    % ���ڷų�β��Ҫ������unique�����
   for i=1:length(b)
%          Strip.Strip_Bin
        Strip = repairStripPlace(Strip,b(i)); % Strip.Strip_Bin 
%          Strip.Strip_Bin
    end
end

    function Strip = repairStripPlace(Strip,stripidx)
        % �ҵ�stripidx��Ӧ��BIN�µ�����Strip�������߼�ֵ
        flagIdx = Strip.Strip_Bin(1,:) == Strip.Strip_Bin(1,stripidx);    % �������ڱ�Bin�ڵ��߼��ж�
        flagBigIdx = Strip.Strip_Bin(2,:) > Strip.Strip_Bin(2,stripidx);  % ���аڷ�˳������stripidx���߼��ж�
        
        % �������ڱ�Bin�� & �Ұڷ�˳������stripidx ��˳���1, ����ǰ�ڷ�
        Strip.Strip_Bin(2,flagBigIdx & flagIdx)  = Strip.Strip_Bin(2,flagBigIdx & flagIdx)  - 1;
        Strip.Strip_Bin(2,stripidx) = sum(flagIdx); % ��ǰstripidx�ڷŵ���β, ��˳�����õ����
        
        % [~,maxSeq]=max(Strip.Strip_Bin(2,Strip.Strip_Bin(1,:) == Strip.Strip_Bin(1,stripidx) ));
    end

   
%% ����script
% �����Ҫ���:���ÿ��item������ ԭʼ LU���z
printscript();
    

%% Ƕ�׺���


    function printscript()
        % �����Ҫ���:��ô�1��ʼÿ��bin����������
        % Strip.stripBeBinMatrix
        for iBin = 1:max(Strip.Strip_Bin(1,:))
            [~,idx] = find(Strip.Strip_Bin(1,:)==iBin); %��iBin�µ�strip������
            idxSeq = Strip.Strip_Bin(2,idx); %��iBin��strip����˳��Seq
            fprintf('bin �Ŀ�+��Ϊ: ' );
            fprintf(' %d  ', Veh.LWH(:,1));
            fprintf('\n');
            fprintf('bin %d ��ʣ���+ʣ�೤Ϊ:  ',iBin);
            fprintf('( %d ) ',Bin.LW(:,iBin));
            fprintf('\n');
            fprintf('bin %d ���� original strip ������{˳��}(����)Ϊ  \n  ',iBin);
            fprintf('%d ',idx);fprintf('\n');
            fprintf('{%d} ',idxSeq);fprintf('\n');
            fprintf('( %d ) ', Strip.LW(1:2,idx));fprintf('\n');
            fprintf('\n');
        end
    end

    % δ��ɺ��� TODO
    function plot2DBin()
    % ��ʼ��
            % ��ʼ��
        Bin.LW
        sStrip.LW
        tmpBin_Strip
        sStrip.Strip_Bin
        sStrip
            %% ��ʼ��
        nThisItem = size(d.Item.LWH,2);
        nIDType = unique(d.Item.LID);
        nColors = hsv(length(nIDType)); %��ͬ����LU���費ͬ��ɫ        
%         tmpUniqueBin = unique(Veh.LWH(1:2,:)','rows')';
        %         wBin = tmpUniqueBin(1);
%         hBin = tmpUniqueBin(2);     
        wBin = Veh.LWH(1,1);
        hBin = Veh.LWH(2,1);
   
    
        nUsedBin = sum(sStrip.Strip_Bin(2,:)>0);

%         %% ��ͼ
        % 1 �������� ���ΪnUsedBin+1��bin�� �����ߣ���Ϊbin��
        figure();
        DrawRectangle([wBin*(nUsedBin+1)/2 hBin/2 wBin*(nUsedBin+1) hBin 0],'--');
        hold on;
        % 2 ���bin ��ͼ
        iterWidth=0;    %ÿ��bin��ǰ1��bin���Ҳ� ��Ϊ���ӱ���
        for iBin = 1:nUsedBin
            % �ҳ���ǰiBin����Ʒ����
            idxDrawStrip = find(sStrip.Strip_Bin(1,:)==iBin);
            % ������ ����û��Strip��bin�ڵ�Coord���˺�����ͣ
        end
        % ���strip��ͼ
        
    end

end


    function order = getStriporder(Strip)
%         tmpLWStrip = Strip.LW(1:2,:);
%         [~,order] = sort(tmpLWStrip(2,:),'descend');  %��strip��������,ֻ��Ҫ����˳��ord;����nDim=2�����򣨳�/�߶�)
       
%��SID����: SID������˳������,���С����ǰ��; 
% �ص�����ͬһSTRIP�����SID: �����С�ĵ�����SID��ǰ, ������ų���SID���ں�, �̶������ͷǻ�ϵģ�
% ���������������ϵĻ�ϣ���ʵ���Ҳ���٣����������ʾ����
SIDorder = getOrderofSID(Strip.SID); %SIDһ���Ǵ�1-n�Ĺ���
if ~issorted(SIDorder), error('SIDδ����С������������'); end
if any(diff(SIDorder)>1), error('SIDδ����,���ж�����'); end

%��LID����: ���ڰڷŵ���Ҫԭ�� 555555555555555555555555 
% LID��ָ��˳��, ����SID����ȫ��һ��,�ٰ�LID��С��������,��ʵû������(��SID/LID����ͬһITEM),��󿴸߶�
% S = [Strip.SID; Strip.LID; Strip.LW(1:2,:); Strip.loadingrateLimit;Strip.loadingrate];
IDorder = getOrderofLID(SIDorder, Strip.isSingleItem, Strip.isAllPured, Strip.nbLID, Strip.isFull,Strip.isMixed, Strip.LID, Strip.LW(1:2,:), Strip.loadingrateLimit, Strip.loadingrate);

% 555�����䣺ͬһSID��,���������ظ���LID
s=[SIDorder;IDorder];
for i=min(SIDorder):max(SIDorder)
    si = s(2,s(1,:)==i);
    if numel(unique(si)) ~= numel(si)
        error('eeeee');
    end
end

tmpSort = [SIDorder; IDorder];     % Strip.LW(1:2,:); Strip.loadingrateLimit;Strip.loadingrate];
[~,order] = sortrows(tmpSort',[1,2],{'ascend','ascend'}); %[~,order] = sortrows(tmpSort',[1,2],{'ascend','ascend'});

        
% LIDorder = getOrderofLID([Strip.SID;Strip.LID]); 
% Strip.LID;
% LIDorder = getOrderofID(Strip.LID); %��LID������: ֻ��һ�ֵ�LID���ȼ���, �����������LID��ϵ�2��STRIP��
% LIDorder = ones(1,length(SIDorder)); 

        % ����Ӧ��SID/LID����
%         zs
%         zl
%         tmpSort = [SIDorder; LIDorder; Strip.LW(1:2,:); Strip.loadingrateLimit;Strip.loadingrate];
%         [~,order] = sortrows(tmpSort',[1,2],{'ascend','ascend'}); %[~,order] = sortrows(tmpSort',[1,2],{'ascend','ascend'}); 
%         order = LIDorder'
%         [~,order] = sortrows(tmpSort',[1],{'descend'});

%         [~,order] = sortrows(tmpSort',[1],{'ascend'}); %��strip��������;����nDim=2�����򣨳�/�߶�)���ٿ�strip�ڲ�loadingrateLimit
%         tmpSort = [zorder; Strip.LW(1:2,:); Strip.loadingrateLimit;Strip.loadingrate];
%         [~,order] = sortrows(tmpSort',[1 3 4 5 ],{'ascend','descend','descend','descend'}); %��strip��������;����nDim=2�����򣨳�/�߶�)���ٿ�strip�ڲ�loadingrateLimit
       
%         tmpLWH = [tmpIDItem; tmpLWHItem]; %��������ITEM��ID����һ���γ���ʱ����
%         [~,order] = sortrows(tmpLWH',[3 1 ],{'descend','descend'}); %���߶�,ID(��ͬ�߶�ʱ)�ݼ�����
%         tmpLWH = [tmpLWHItem; min(tmpLWHItem(1:2,:))]; %����������̱ߵ��������γ���ʱ����tmpLWH
%         [~,order] = sortrows(tmpLWH',[3 2 1],{'descend','descend','descend'}); %way1 ����̱�,�߶�,��ȵݼ�����
      
        if ~isrow(order), order=order'; end
        
                    %  t = Strip.SID;
                    %  ss = sum(t);  %ÿ��STRIP�ڰ�����SID����
                    %  
                    %   for i=1:size(t,1)
                    %      if sum(ss( find(t(i,:))  ) > 1) > 1
                    %          error('��ͬһ��SID��2��������STRIP����');
                    %      end
                    %  end
                    % if  any(sum(t)>2)
                    %     error('��ͬһ��strip����3�������ϸ�SID');
                    % end
                    % 
                    %  z = zeros(1,size(t,2)); 
                    %  k=1;
                    %  for i=1:numel(ss)
                    %     if ss(i) ==1
                    %         if i>1 && ss(i-1) ==1 && find(t(:, i)==1) ~= find(t(:, i-1)==1) %�жϵ�ǰ��ǰһ��strip�Ƿ�����ͬ��SID
                    %             k=k+1;
                    %         end
                    %         z(i) = k;
                    %     elseif ss(i) >1 %ֻҪ����STRIP����2�������ϵ�STRIPʱ������˳��        
                    %         k=k+1;
                    %         z(i)=k;
                    %         k=k+1;
                    %     end
                    %  end

                                         %  ss = sum(tSID);
                                        %  torder = zeros(1,length(sorder));
                                        %  k=1;
                                        %  for i=1:length(sorder)
                                        %      tt = tSID(sorder(i),:)
                                        %      to
                                        % %      tSID(sorder(i)) = [];
                                        %      other = sorder;
                                        %      other(sorder(i)) = [];
                                        %      tti = tSID(sorder(other),:);
                                        %      torder(tt==1) = k;
                                        %      k = k+1;
                                        %      f = tti==1 & tt==1 %����SID�У�����������Ҳ��,����Ϊk+1
                                        %      if any(f)
                                        %          torder(f) = k; 
                                        %          k = k+1;
                                        %      end
                                        %      1
                                        %  end
                                        % 
                                        % Strip.SID
                                        % [a,b,~]=find(Strip.SID==1)
                                        % [a,b,~]=find(Strip.SID(:,:)==1)
                                        % Strip.LID
    end

        function [thisBin,iBin]  = getThisBin( iBin, iStrip, sStrip, Veh, Bin, p)        
        if p.whichBinH == 1 % 1 bestfit
            % ����: Ѱ�� bin��ʣ��߶� >= ��strip�ĸ߶� &�� bin��ʣ������ >= ��strip������ (�����е���Сֵ)
            flag = find(Bin.LW(2,1:iBin) >= sStrip.LW(2,iStrip)  & ...
                Veh.Weight(1) - Bin.Weight(1 : iBin) >= sStrip.Weight(iStrip) ); %
            if isempty(flag)
                iBin = iBin + 1; % ����߶Ȳ����㣬��bin����
                [thisBin,iBin]  = getThisBin( iBin, iStrip, sStrip, Veh, Bin, p);                
            else
                tepBins = Bin.LW(2,1 : iBin); %��ȡ�����Ѱ��Ż��°��ŵ�bin��ʣ��߶�����tepBins
                tepMin = min(tepBins(flag)); % 555 check �ҳ�bin���ܷ�istrip�Ҹ߶���СֵtepMin��TODO �Ƿ�����������
                thisBin = find(Bin.LW(2,1:iBin)==tepMin); %�ҵ���ֵtepMin��Ӧ���Ǹ�/Щbin���
                if ~all(ismember(thisBin,flag)),      error('Not all thisBin belongs to flag ');        end
                if length(thisBin)>1
                    thisBin = thisBin(1);
                end
            end
        elseif p.whichBinH == 2 % 1 firstfit
            flag = find(Bin.LW(2,1:iBin) >= sStrip.LW(2,iStrip)  & ...
                Veh.Weight(1) - Bin.Weight(1 : iBin) >= sStrip.Weight(iStrip) );            
            if isempty(flag)
                iBin = iBin + 1; 
                [thisBin,iBin]  = getThisBin( iBin, iStrip, sStrip, Veh, Bin, p);    
            else
                thisBin = flag(1);
                if ~all(ismember(thisBin,flag)),     error('Not all thisBin belongs to flag ');       end
            end
        elseif p.whichBinH == 3 % 1 nextfit
            flaged = find(Bin.LW(2, iBin) >= sStrip.LW(2,iStrip) & ...
                Veh.Weight(1) - Bin.Weight(iBin) >= sStrip.Weight(iStrip) );            
            if  isempty(flaged)  %ע����֮ǰ~flag������
                iBin = iBin + 1; 
                [thisBin,iBin]  = getThisBin( iBin, iStrip, sStrip, Veh, Bin, p);    
            else
                if  isempty(flaged) ,   error(' �����ܵĴ��� ');      end
                thisBin = iBin; % ��ǰbinһ���ŵ���
            end
        else
            error('�����������');
        end
        end

        
    
    function [Bin,Strip_Bin,Bin_Strip] = insertStripToBin(iStrip, thisBin,sStrip,Bin,Strip_Bin,Bin_Strip)
%         binBeStripArray=binBeStripArray;stripBeBinMatrixSort=stripBeBinMatrixSort;Bin.LW=Bin.LW;

        % 1 ����Bin��ط�Sort����        
        %  1.1 ����strip����bin����Ϣ ��stripBeBinMatrixSort
        Bin_Strip(thisBin) = Bin_Strip(thisBin) + 1; %��bin�µڼ��ΰ���strip
        
        %  1.2 ���±�bin��ʣ�೤��ʣ��ߣ�Bin.LW
        Bin.LW(1,thisBin) = min(Bin.LW(1,thisBin),sStrip.LW(1,iStrip)); %����binʣ���ȵ���Сֵ
        Bin.LW(2,thisBin) = Bin.LW(2,thisBin) - sStrip.LW(2,iStrip);    %����binʣ��߶�
            
        %  1.3 ���±�bin��Ӧ��BinWeight: 
        Bin.Weight(thisBin) =  Bin.Weight(thisBin) + sStrip.Weight(iStrip);
        
        % 2 ����Strip���Sort����
        %  2.1 ����stripBeBinMatrixSort
        Strip_Bin(1,iStrip) = thisBin;
        Strip_Bin(2,iStrip) = Bin_Strip(thisBin);
   
            % ����bIN�а���ID�����
%             Bin.LID(:,thisBin) = Bin.LID(:,thisBin) + sStrip.LID(:,iStrip); % ��ֵΪ���ִ���
%             Bin.LID(Bin.LID>0) = 1; % ��ֵ��Ϊ�������
%             Bin.PID(:,thisBin) = Bin.PID(:,thisBin) + sStrip.PID(:,iStrip); % ��ֵΪ���ִ���
%             Bin.PID(Bin.PID>0) = 1; % ��ֵ��Ϊ�������
%             Bin.SID(:,thisBin) = Bin.SID(:,thisBin) + sStrip.SID(:,iStrip); % ��ֵΪ���ִ���
%             Bin.SID(Bin.SID>0) = 1; % ��ֵ��Ϊ�������
%             Bin.UID(:,thisBin) = Bin.UID(:,thisBin) + sStrip.UID(:,iStrip); % ��ֵΪ���ִ���
%             Bin.UID(Bin.UID>0) = 1; % ��ֵ��Ϊ�������
            
            
       %% ����ŵ�ItemToBin�ڼ���
    end
