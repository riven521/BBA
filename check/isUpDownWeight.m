function [TF] = isUpDownWeight(T)
%isUpDownWeight 判定在同一Item下, 是否上轻下重
%   return a logical scalar

TF = false;

if ~any(strcmp('ITEMSEQ', T.Properties.VariableNames)) || ~any(strcmp('Weight', T.Properties.VariableNames)) 
    error('输入有误');
end
    
T = sortrows(T,'ITEMSEQ');

if ~issorted(T.Weight,'descend')
    warning('相同ITEM,但重量不是递减');  % T.isWeightUpDown(flagLUIdx) = 1;
    TF = true;
end

end
