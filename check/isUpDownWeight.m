function [TF] = isUpDownWeight(T)
%isUpDownWeight �ж���ͬһItem��, �Ƿ���������
%   return a logical scalar

TF = false;

if ~any(strcmp('ITEMSEQ', T.Properties.VariableNames)) || ~any(strcmp('Weight', T.Properties.VariableNames)) 
    error('��������');
end
    
T = sortrows(T,'ITEMSEQ');

if ~issorted(T.Weight,'descend')
    warning('��ͬITEM,���������ǵݼ�');  % T.isWeightUpDown(flagLUIdx) = 1;
    TF = true;
end

end
