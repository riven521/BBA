% ������߶Ե�,����Rotaed����; 
function [LWHRotaed]  =getRotaedLWH(LWH,Rotaed, Buff)
        % ����ת�Ĳ���������buff������ת�ĳ����滻������buff��
        LWHRotaed = LWH;
        LWHRotaed(1, Rotaed) = LWH(2, Rotaed) - Buff(2, Rotaed) +Buff(1, Rotaed);
        LWHRotaed(2, Rotaed) = LWH(1, Rotaed) - Buff(1, Rotaed) +Buff(2, Rotaed); 
      
end