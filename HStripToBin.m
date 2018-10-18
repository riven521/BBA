function [Strip,Bin]= HStripToBin(Strip,Veh,p)
% ��Ҫ����:Strip����Bin�� %  ����:�����(row);  ����:��������(coloum);
% Input ---  Strip/Veh:  
% Output --- Strip: 
% Output --- Bin: 

%% ��ʼ��
sz = size(Strip.LW);
nStrip = sz(2);

wVeh  = Veh.LWH(1,1); 
lVeh  = Veh.LWH(2,1); 

%% Strip���� 555 (���ȷ����ͬ�������͵����ڰڷ� DONE)
    % ��ȡStrip��˳��(�ص���Strip�߶ȵݼ����򣨵���������strip�߶�һ���ģ�) %
    Strip
    Strip.nbLU
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

% %  plot2DBin();

% Strip�ڲ�����,sStrip����order�仯����
if isSameCol(sStrip)
    Strip = reorderStruct(Strip.striporder, sStrip);
else
    error('����ʹ��structfun');
end

% ��ȡBin: ȥ��δʹ�õ�Bin ע��Bin�ṹ��ı仯
if isSameCol(Bin)
    Bin = structfun(@(x) x( : , Bin.Weight(1,:)>0 ), Bin, 'UniformOutput', false);
else
    error('����ʹ��structfun');
end

%% ����script �����Ҫ���:���ÿ��item������ ԭʼ LU���z
% printscript();

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

%% �ֲ�����

function order = getStriporder(Strip)
%��SID����: SID������˳������,���С����ǰ��; 
% �ص�����ͬһSTRIP�����SID: �����С�ĵ�����SID��ǰ, ������ų���SID���ں�, �̶������ͷǻ�ϵģ�
% ���������������ϵĻ�ϣ���ʵ���Ҳ���٣����������ʾ����
        % Strip.SID
SIDorder = getOrderofSID(Strip.SID); %SIDһ���Ǵ�1-n�Ĺ���
if ~issorted(SIDorder), error('SIDδ����С������������'); end
if any(diff(SIDorder)>1), error('SIDδ����,���ж�����'); end

%��LID����: ���ڰڷŵ���Ҫԭ�� 555555555555555555555555 
% LID��ָ��˳��, ����SID����ȫ��һ��,�ٰ�LID��С��������,��ʵû������(��SID/LID����ͬһITEM),��󿴸߶�
IDorder = getOrderofLID(SIDorder, Strip.isSingleItem, Strip.isAllPured, Strip.nbItem,Strip.nbLU,Strip.nbLULID,Strip.isHeightFull,...
                                        Strip.isMixed, Strip.LID, Strip.LW(1:2,:), Strip.loadingrateLimit, Strip.loadingrate);

% 555 ������䣺ͬһSID��,���������ظ���LID
s=[SIDorder;IDorder];
for i=min(SIDorder):max(SIDorder)
    si = s(2,s(1,:)==i);
    if numel(unique(si)) ~= numel(si),     error('ͬһSID��, ���ظ���LID');   end
end

% ����Priority��������������
tmpSort = [SIDorder; IDorder];
[~,order] = sortrows(tmpSort',[1,2],{'ascend','ascend'});
if ~isrow(order), order=order'; end

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
