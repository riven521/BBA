% 将长宽高对调,依据Rotaed变量; 
function [LWHRotaed]  =getRotaedLWH(LWH,Rotaed, Buff)
        % 不旋转的不动（包括buff）；旋转的长宽替换（包括buff）
        LWHRotaed = LWH;
        LWHRotaed(1, Rotaed) = LWH(2, Rotaed) - Buff(2, Rotaed) +Buff(1, Rotaed);
        LWHRotaed(2, Rotaed) = LWH(1, Rotaed) - Buff(1, Rotaed) +Buff(2, Rotaed); 
      
end