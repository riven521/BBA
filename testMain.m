function tests = testMain
    tests = functiontests(localfunctions);
end

%% runtests('testMain.m') 
% https://ww2.mathworks.cn/help/matlab/ref/runtests.html

function testMixedGap(testCase)
% testMixedGap
LUID=[1,2,3,4];
LULWH=[760.0,1955.0,550.0;
                760.0,1140.0,990.0;
                760.0,1955.0,550.0;
                820.0,1220.0,1000.0;];
VEHID=[1];
VEHLWH=[2400.0,12500.0,2400.0;];
LUSID=[1,1,1,1];
LUPID=[1,2,3,4];
LUISROTA=[1,1,1,1];
LUMARGIN=[25.0,25.0,25.0,25.0;25.0,25.0,25.0,25.0;25.0,25.0,25.0,25.0;25.0,25.0,25.0,25.0;];
LUWEIGHT=[2987.392,2919.183,1831.5,269.683];
VEHWEIGHT=[9999999.0];
LULID=[1,3,4,2];
BBAID=[1,2,3,4];

warning('off')
actSolution = ...
        BBA_Main(LULID,LULWH,VEHID,VEHLWH,...
        LUSID, LUPID, LUISROTA, LUMARGIN, LUWEIGHT, VEHWEIGHT,LUID,BBAID);
warning('on')
expLength = length(LUID);

verifyLength(testCase,actSolution,expLength)

end

% % function testBBAmain(testCase)
% % % 
% % LUID=[1,2,3,4,5,6,7,8,9];
% % LULWH=[820.0,1220.0,1000.0;1140.0,1710.0,990.0;760.0,1140.0,990.0;820.0,1220.0,1000.0;1140.0,1710.0,990.0;1140.0,1710.0,990.0;760.0,1140.0,990.0;1140.0,1710.0,990.0;1140.0,1710.0,990.0;];
% % VEHID=[1];
% % VEHLWH=[2400.0,12500.0,2400.0;];
% % % 
% % LUSID=[1,1,1,1,1,1,1,1,1];
% % LUPID=[1,2,3,4,5,6,7,8,9];
% % LUISROTA=[1,1,1,1,1,1,1,1,1];
% % LUMARGIN=[25.0,25.0,25.0,25.0;25.0,25.0,25.0,25.0;25.0,25.0,25.0,25.0;25.0,25.0,25.0,25.0;25.0,25.0,25.0,25.0;25.0,25.0,25.0,25.0;25.0,25.0,25.0,25.0;25.0,25.0,25.0,25.0;25.0,25.0,25.0,25.0;];
% % LUWEIGHT=[778.36,1815.0,246.16,140.14,8521.57,537.45,18533.06,1431.36,182.5];
% % VEHWEIGHT=[9999999.0];
% % LULID=[4,3,1,4,3,3,5,3,3];
% % BBAID=[1,2,3,4,5,6,7,8,9];
% % LUINDEX = BBAID;
% % actSolution = ...
% %         BBA_Main(LULID,LULWH,VEHID,VEHLWH,...
% %         LUSID, LUPID, LUISROTA, LUMARGIN, LUWEIGHT, VEHWEIGHT,LUID,LUINDEX);
% % expLength = [9];
% % verifyLength(testCase,actSolution,expLength)
% % end


