function [da]= HStripToBin(da,ParaArray)
% ��Ҫ����:Strip����Bin�� %  ����:�����(row);  ����:��������(coloum);
% Input ---  Strip:  
% Output --- Strip: 
% Output --- Bin: 

%% ��ʼ��
% nDim stripά��(2) 
nDim = size(da.StripArray.LW,1);  
nStrip = size(da.StripArray.LW,2); %����ʹ�õ�Strip������ 
nBin = nStrip;
uniBinDataMatrix = unique((da.BinArray.LWH(1:nDim,:))','rows')';

% % printstruct(da);

%% Strip���� 555 (���ȷ����ͬ�������͵����ڰڷţ����� TODO FIX ME)
% ��ȡstriporder = da.StripArray.striporder
% ��ȡsortedStripArray 
    % getStriporder - ��ȡStrip��˳��(�ص���Strip�߶ȵݼ����򣨵���������strip�߶�һ���ģ�) %
    % Strip������ʽ �߶�/���ȵݼ�
    da.StripArray.striporder = getStriporder(); 
    % getSortedStrip - ��ȡ��order������Strip :sortedStripArray
    sortedStripArray = getSortedStrip(da.StripArray.striporder);
%     printstruct(da) ;printstruct(sortedStripArray)    
    
%% LU->Item->Strip->Binת�� 
% ��ȡstripBeBinMatrixSort: ÿ�������strip���ĸ�bin��  �Լ�˳��
% ��ȡLWBin:  �����ɵ�Bin��ʣ�೤��
LWBin = zeros(nDim,nBin);    %��ʼ��bin: dim1-bin���ʣ�� ; dim2-bin��(��)��(555ʣ�ࣩ;
LWBin(1,:) = uniBinDataMatrix(1);
LWBin(2,:) = uniBinDataMatrix(2);
binBeStripArray = zeros(1,nBin);    % ÿ��Bin�ڵ�Strip���� ���ڲ���
stripBeBinMatrixSort = zeros(2,nStrip); % dim1:��� strip��ĳ��bin dim2:����˳�� 555

LWStripSort = sortedStripArray.LW; %��sorted���

% 55 ��ȡthisBin - ��ǰstripҪ�����bin���
% ѭ����bin�а���strip,���̶�strip,�仯ѡ��ͬbin(thisBin)
% ע�ͣ���ȡ FLAG        �ɷ��µ�ǰiStrip������һ��bin�ļ��� 
% ע�ͣ���ȡ thisBin   ��FLAG���ҵ���������Ǹ�thisBin, ��ִ�� insert����

iStrip=1; iBin=1;
while 1
    if iStrip > nStrip, break; end
    
    thisBin = getThisBin();    % ��ȡBin��
    
    insertStripToBin(); %����strip��Bin
    
    iStrip = iStrip + 1;
end

% ���� ����ֵ��da
% ��ȡstripBeBinMatrix: ÿ��strip���ĸ�bin��  �Լ�˳��
% ��ȡLWBin:  �����ɵ�bin��ʣ�೤��
% ��ȡstriporder: strip������
stripBeBinMatrix=stripBeBinMatrixSort;
stripBeBinMatrix( : , da.StripArray.striporder) = stripBeBinMatrixSort;
da.StripArray.stripBeBinMatrix = stripBeBinMatrix;

LWBin = LWBin(:,LWBin(2,:)~=uniBinDataMatrix(2));  % LWBin = LWBin(:,LWBin(2,:)~=uniBinDataMatrix(2));
da.BinSArray.LW = LWBin; % ȥ��δʹ�õ�Strip

%% ����script
% �����Ҫ���:���ÿ��item������ ԭʼ LU���
printscript();
    

%% Ƕ�׺���
    function order = getStriporder()
        tmpLWStrip = da.StripArray.LW(1:nDim,:);
        
        [~,order] = sort(tmpLWStrip(nDim,:),'descend');  %��strip��������,ֻ��Ҫ����˳��ord;����nDim=2�����򣨳�/�߶�)
        
%         tmpLWH = [tmpIDItem; tmpLWHItem]; %��������ITEM��ID����һ���γ���ʱ����
%         [~,order] = sortrows(tmpLWH',[3 1 ],{'descend','descend'}); %���߶�,ID(��ͬ�߶�ʱ)�ݼ�����
%         tmpLWH = [tmpLWHItem; min(tmpLWHItem(1:nDim,:))]; %����������̱ߵ��������γ���ʱ����tmpLWH
%         [~,order] = sortrows(tmpLWH',[3 2 1],{'descend','descend','descend'}); %way1 ����̱�,�߶�,��ȵݼ�����
      
        %         if ~isrow(order), order=order'; end
    end

    function strip = getSortedStrip(order)
        strip = structfun(@(x) x(:,order),da.StripArray,'UniformOutput',false);
    end

    function thisBin = getThisBin()        
        if ParaArray.whichBinH == 1 % 1 bestfit
            % ����: Ѱ�� bin��ʣ��߶� >= ��strip�ĸ߶� �����е���Сֵ
            flag = find(LWBin(2,1:iBin) >= LWStripSort(2,iStrip), 1);
            
            if isempty(flag)
                iBin = iBin + 1;
                thisBin = getThisBin();
                
            else
                tepMin = LWBin(2,1:iBin);
                tepMin = min(tepMin(flag)); % 555 check �ҳ�bin���ܷ�istrip�Ҹ߶���СֵtepMin
                thisBin = find(LWBin(2,1:iBin)==tepMin); %�ҵ���ֵtepMin��Ӧ��bin���
                if length(thisBin)>1
                    thisBin = thisBin(1);
                end
            end
        elseif ParaArray.whichBinH == 2 % 1 firstfit
            
        elseif ParaArray.whichBinH == 3 % 1 nextfit
            
        else
            error('�����������');
        end
    end

    function insertStripToBin()
        % 1 ����strip����bin����Ϣ (stripBeBinMatrixSort)
        binBeStripArray=binBeStripArray;stripBeBinMatrixSort=stripBeBinMatrixSort;LWBin=LWBin;
        binBeStripArray(thisBin) = binBeStripArray(thisBin) + 1; %��bin�µڼ��ΰ���strip
        stripBeBinMatrixSort(1,iStrip) = thisBin;
        stripBeBinMatrixSort(2,iStrip) = binBeStripArray(thisBin);
        
        % 2 ��ȡ��iStrip�ڵ�item���, ������Item������Ϣ
%         idxItemStrip = find(ItemArray.itemBeStripMatrixSort(1,:)==iStrip);
%         itemBeBinMatrixSort(1,idxItemStrip) = thisBin;    %�ڼ���bin

         % 3 ����LWBin
            LWBin(1,thisBin) = min(LWBin(1,thisBin),LWStripSort(1,iStrip)); %����binʣ���ȵ���Сֵ
            LWBin(2,thisBin) = LWBin(2,thisBin) - LWStripSort(2,iStrip);    %����binʣ��߶�
%        binBeItemArray(thisBin) = binBeItemArray(thisBin) + length(idxItemStrip);                          %��bin�ºϼƼ���item
            
            %����xy������Ϣ x���� yͨ��bin�߶�-binʣ��߶�-����strip�߶�
%             CoordItemBinSort(1,idxItemStrip) = CoordItemBinSort(1,idxItemStrip);
%             CoordItemBinSort(2,idxItemStrip) = uniBinDataMatrix(2,1) - (LWBin(2,thisBin) + LWStripSort(2,iStrip));
   

% %         if ParaArray.whichRotation == 1
% %             %����bin����Ϣ
% %             LWBin(1,thisBin) = min(LWBin(1,thisBin),LWStripSort(1,iStrip)); %����binʣ���ȵ���Сֵ
% %             LWBin(2,thisBin) = LWBin(2,thisBin) - LWStripSort(2,iStrip);    %����binʣ��߶�
% % %             binBeItemArray(thisBin) = binBeItemArray(thisBin) + length(idxItemStrip);                          %��bin�ºϼƼ���item
% %             
% %             %����xy������Ϣ x���� yͨ��bin�߶�-binʣ��߶�-����strip�߶�
% % %             CoordItemBinSort(1,idxItemStrip) = CoordItemBinSort(1,idxItemStrip);
% % %             CoordItemBinSort(2,idxItemStrip) = uniBinDataMatrix(2,1) - (LWBin(2,thisBin) + LWStripSort(2,iStrip));
% %         else
% %             %����bin����Ϣ
% %             LWBin(1,thisBin) = min(LWBin(1,thisBin),LWStripSort(1,iStrip)); %����binʣ���ȵ���Сֵ
% %             LWBin(2,thisBin) = LWBin(2,thisBin) - LWStripSort(2,iStrip);    %����binʣ��߶�
% % %             binBeItemArray(thisBin) = binBeItemArray(thisBin) + length(idxItemStrip);                          %��bin�ºϼƼ���item
% %             
% %             %����xy������Ϣ x���� yͨ��bin�߶�-binʣ��߶�-����strip�߶�
% % %          CoordItemBinSort(1,idxItemStrip) = ItemArray.itemCoordMatrixSort(1,idxItemStrip);        %
% % %          CoordItemBinSort(2,idxItemStrip) = uniBinDataMatrix(2,1) - (LWBin(2,thisBin) + LWStripSort(2,iStrip));
% %         end        
    end

    function printscript()
        % �����Ҫ���:��ô�1��ʼÿ��bin����������
        % da.StripArray.stripBeBinMatrix
        for iBin = 1:max(da.StripArray.stripBeBinMatrix(1,:))
            [~,idx] = find(da.StripArray.stripBeBinMatrix(1,:)==iBin); %��iBin�µ�strip������
            idxSeq = da.StripArray.stripBeBinMatrix(2,idx); %��iBin��strip����˳��Seq
            fprintf('bin �Ŀ�+��Ϊ: ' );
            fprintf(' %d  ',uniBinDataMatrix);
            fprintf('\n');
            fprintf('bin %d ��ʣ���+ʣ�೤Ϊ:  ',iBin);
            fprintf('( %d ) ',da.BinSArray.LW(:,iBin));
            fprintf('\n');
            fprintf('bin %d ���� original strip ������{˳��}(����)Ϊ  \n  ',iBin);
            fprintf('%d ',idx);fprintf('\n');
            fprintf('{%d} ',idxSeq);fprintf('\n');
            fprintf('( %d ) ', da.StripArray.LW(1:nDim,idx));fprintf('\n');
            fprintf('\n');
        end
    end

end
