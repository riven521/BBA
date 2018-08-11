% V2 : buff: 每个托盘的margin . 将长宽高对调,依据Rotaed变量; 
% function [LWHRotaed]  =getRotaedLWH(LWH,Rotaed, buff)
%         % 不旋转的不动（包括buff）；旋转的长宽替换（包括buff）
%         LWHRotaed = LWH;
%         LWHRotaed(1, Rotaed) = LWH(2, Rotaed) - (buff(3, Rotaed) +buff(4, Rotaed)) + (buff(1, Rotaed) +buff(2, Rotaed));
%         LWHRotaed(2, Rotaed) = LWH(1, Rotaed) - (buff(1, Rotaed) +buff(2, Rotaed)) + (buff(3, Rotaed) +buff(4, Rotaed)); 
% end

% V1 : buff: 托盘间的间隙 . 将长宽高对调,依据Rotaed变量; 
function [LWHRotaed]  =getRotaedLWH(LWH,Rotaed, Buff)
        % 不旋转的不动（包括buff）；旋转的长宽替换（包括buff）
        if nargin < 3
            Buff = zeros(size(LWH));
        end
        LWHRotaed = LWH;
        LWHRotaed(1, Rotaed) = LWH(2, Rotaed) - Buff(2, Rotaed) +Buff(1, Rotaed);
        LWHRotaed(2, Rotaed) = LWH(1, Rotaed) - Buff(1, Rotaed) +Buff(2, Rotaed); 
      
end