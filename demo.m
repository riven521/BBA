% DEMO FILE
function demo(varargin)
% Example :  
%   demo(2, 1, 1)
 
% clear all
close all
clc;
fprintf('\n BBA 3D Loading Plan TOOLBOX v108 2019 - For Matlab \n');

% Include dependencies
addpath('./check'); % check
addpath('./fun'); 
% addpath('./test'); 
addpath('./utility'); 

global nDemo
%190424 gap�������ƻ�horizontal stability ������gap  % 190421 : gap������˦β��Ų��ǰ�� ����gap
%190420����ͬYֵ��gap��������������X������ ����gap
if nargin < 1
    nDemo = 190420; %190508; %190415/190420(gap) /190424 /190421 / 190508 190510
else
    nDemo = varargin{1}; 
end


global ISplotBBA ISplotShowGapAdjust ISplotEachPingPuShuaiWei ISplotEachPingPuAll ISplotshuaiwei  ISplotStripToBinAgain  ISplotRunAlgo ISplotRunLIS ISplotPause ISplotShowType  ISplotGapCompare ISplotPauseWait
global  parGap ISshuaiwei ISreStripToBin

% ��ͼ����
if isempty(ISplotEachPingPuShuaiWei),  ISplotEachPingPuShuaiWei =0 ;   end % 555 �Ƚ���Ҫչʾ���� ÿ��˦βƽ�̳ɹ���չʾƽ��ǰ��ĶԱ�ͼ
if isempty(ISplotEachPingPuAll),  ISplotEachPingPuAll = 0;   end                        % 555 �Ƚ���Ҫչʾ���� ÿ������ƽ�̳ɹ���չʾƽ��ǰ��ĶԱ�ͼ

if isempty(ISplotRunLIS),  ISplotRunLIS = 0;   end                                    % 3 �����Ҫ����չʾ���� 5 ��Ҫ, �����̵��Ѷ⵽strip��չʾ,��ͼ
if isempty(ISplotshuaiwei),  ISplotshuaiwei =0 ;   end                              % 555 �Ƚ���Ҫչʾ���� - ˦β������
if isempty(ISplotStripToBinAgain),  ISplotStripToBinAgain = 0;   end           % ����ͷ Ĭ��1

if isempty(ISplotGapCompare),  ISplotGapCompare = 1;   end                      % 555 �Ƚ���Ҫչʾ���� - ÿ��gap��������ʾ
if isempty(ISplotShowGapAdjust),  ISplotShowGapAdjust = 1;   end                    % 3 �����Ҫ ��϶��������չʾ���� ��gap��������ԭ��

% ��۲�����
if isempty(ISplotPause),  ISplotPause = 0.0;   end   % 0.5 % plot���ʱ��
if isempty(ISplotShowType),  ISplotShowType = 3;   end   % 1 LID 3 ID 8 ˦β
if isempty(ISplotPauseWait),  ISplotPauseWait = 0;   end   % �Ƿ�plotsolutinT���ͼֱ�ӵȴ��û���Ӧ

% ������������ͼ����
if isempty(ISplotBBA),  ISplotBBA = 1;   end   %  555 �Ƚ���Ҫչʾ���������ս��չʾ
if isempty(ISplotRunAlgo),  ISplotRunAlgo = 1;   end                       % Ĭ��1

% �������������ܲ���
if isempty(parGap),  parGap = 1;   end                              % �Ƿ������������ļ�϶����
if isempty(ISshuaiwei),  ISshuaiwei = 1;   end  % 555 : ��Ⱥ͸߶Ȳ���, ˦β   ******  �ò�����Ҫ�������pingpu���ʹ�� ��˦β ˦βƽ���޷�����*******
if isempty(ISreStripToBin),  ISreStripToBin = 1;   end  % ��ͷ����LU����������� Ĭ��Ϊ1 ����


if nargin >= 2
    ISshuaiwei = varargin{2}; 
end

if nargin >= 3
    if varargin{3}
        [ISplotshuaiwei,ISplotGapCompare,ISplotEachPingPuAll,ISplotEachPingPuShuaiWei,ISplotRunAlgo] = deal(1);
    else
        [ISplotshuaiwei,ISplotGapCompare,ISplotEachPingPuAll,ISplotEachPingPuShuaiWei,ISplotRunAlgo] = deal(0);
    end
end


%% ������
runtests('testBBA.m')                                                 % ������Ŀ¼���ļ�

% results = runtests('test/testBBA.m');                               % ������Ŀ¼���ļ�
% results = runtests(pwd,'IncludeSubfolders',true);       % ���Ժ���Ŀ¼�����в����ļ�

% testBBA

%%


% % %%
% % % Select a feature selection method from the list
% % listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','rfe','L0','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};
% % 
% % [ methodID ] = readInput( listFS );
% % selection_method = listFS{methodID}; % 4Selected
% % 
% % % Load the data and select features for classification
% % load fisheriris
% % X = meas; clear meas
% % % Extract the Setosa class
% % Y = nominal(ismember(species,'setosa')); clear species
% % 
% % % Randomly partitions observations into a training set and a test
% % % set using stratified holdout
% % P = cvpartition(Y,'Holdout',0.20);
% % 
% % X_train = double( X(P.training,:) );
% % Y_train = (double( Y(P.training) )-1)*2-1; % labels: neg_class -1, pos_class +1
% % X_train = [X_train,rand(120,4)];
% % X_test = double( X(P.test,:) );
% % Y_test = (double( Y(P.test) )-1)*2-1; % labels: neg_class -1, pos_class +1
% % X_test = [X_test,rand(30,4)];
% % % number of features
% % numF = size(X_train,2);
% % 
% % 
% % % feature Selection on training data
% % switch lower(selection_method)
% %     case 'inffs'
% %         % Infinite Feature Selection 2015 updated 2016
% %         alpha = 0.5;    % default, it should be cross-validated.
% %         sup = 1;        % Supervised or Not
% %         [ranking, w] = infFS( X_train , Y_train, alpha , sup , 0 );
% %         
% %     case 'ilfs'
% %         % Infinite Latent Feature Selection - ICCV 2017
% %         [ranking, weights] = ILFS(X_train, Y_train , 6, 0 );
% %         
% %     case 'fsasl'
% %         options.lambda1 = 1;
% %         options.LassoType = 'SLEP';
% %         options.SLEPrFlag = 1;
% %         options.SLEPreg = 0.01;
% %         options.LARSk = 5;
% %         options.LARSratio = 2;
% %         nClass=2;
% %         [W, S, A, objHistory] = FSASL(X_train', nClass, options);
% %         [v,ranking]=sort(abs(W(:,1))+abs(W(:,2)),'descend');
% %     case 'lasso'
% %         lambda = 25;
% %         B = lasso(X_train,Y_train);
% %         [v,ranking]=sort(B(:,lambda),'descend');
% %         
% %     case 'ufsol'
% %         para.p0 = 'sample';
% %         para.p1 = 1e6;
% %         para.p2 = 1e2;
% %         nClass = 2;
% %         [~,~,ranking,~] = UFSwithOL(X_train',nClass,para) ;
% %         
% %     case 'dgufs'
% %         
% %         S = dist(X_train');
% %         S = -S./max(max(S)); % it's a similarity
% %         nClass = 2;
% %         alpha = 0.5;
% %         beta = 0.9;
% %         nSel = 2;
% %         [Y,L,V,Label] = DGUFS(X_train',nClass,S,alpha,beta,nSel);
% %         [v,ranking]=sort(Y(:,1)+Y(:,2),'descend');
% %         
% %         
% %     case 'mrmr'
% %         ranking = mRMR(X_train, Y_train, numF);
% %         
% %     case 'relieff'
% %         [ranking, w] = reliefF( X_train, Y_train, 20);
% %         
% %     case 'mutinffs'
% %         [ ranking , w] = mutInfFS( X_train, Y_train, numF );
% %         
% %     case 'fsv'
% %         [ ranking , w] = fsvFS( X_train, Y_train, numF );
% %         
% %     case 'laplacian'
% %         W = dist(X_train');
% %         W = -W./max(max(W)); % it's a similarity
% %         [lscores] = LaplacianScore(X_train, W);
% %         [junk, ranking] = sort(-lscores);
% %         
% %     case 'mcfs'
% %         % MCFS: Unsupervised Feature Selection for Multi-Cluster Data
% %         options = [];
% %         options.k = 5; %For unsupervised feature selection, you should tune
% %         %this parameter k, the default k is 5.
% %         options.nUseEigenfunction = 4;  %You should tune this parameter.
% %         [FeaIndex,~] = MCFS_p(X_train,numF,options);
% %         ranking = FeaIndex{1};
% %         
% %     case 'rfe'
% %         ranking = spider_wrapper(X_train,Y_train,numF,lower(selection_method));
% %         
% %     case 'l0'
% %         ranking = spider_wrapper(X_train,Y_train,numF,lower(selection_method));
% %         
% %     case 'fisher'
% %         ranking = spider_wrapper(X_train,Y_train,numF,lower(selection_method));
% %         
% %         
% %     case 'ecfs'
% %         % Features Selection via Eigenvector Centrality 2016
% %         alpha = 0.5; % default, it should be cross-validated.
% %         ranking = ECFS( X_train, Y_train, alpha )  ;
% %         
% %     case 'udfs'
% %         % Regularized Discriminative Feature Selection for Unsupervised Learning
% %         nClass = 2;
% %         ranking = UDFS(X_train , nClass );
% %         
% %     case 'cfs'
% %         % BASELINE - Sort features according to pairwise correlations
% %         ranking = cfs(X_train);
% %         
% %     case 'llcfs'
% %         % Feature Selection and Kernel Learning for Local Learning-Based Clustering
% %         ranking = llcfs( X_train );
% %         
% %     otherwise
% %         disp('Unknown method.')
% % end
% % 
% % k = 2; % select the first 2 features
% % 
% % % Use a linear support vector machine classifier
% % svmStruct = svmtrain(X_train(:,ranking(1:k)),Y_train,'showplot',true);
% % C = svmclassify(svmStruct,X_test(:,ranking(1:k)),'showplot',true);
% % err_rate = sum(Y_test~= C)/P.TestSize; % mis-classification rate
% % conMat = confusionmat(Y_test,C); % the confusion matrix
% % 
% % disp('X_train size')
% % size(X_train)
% % 
% % disp('Y_train size')
% % size(Y_train)
% % 
% % disp('X_test size')
% % size(X_test)
% % 
% % disp('Y_test size')
% % size(Y_test)
% % 
% % ranking
% % fprintf('\nMethod %s (Linear-SVMs): Accuracy: %.2f%%, Error-Rate: %.2f \n',selection_method,100*(1-err_rate),err_rate);
% % 
% % 
% % 
% % % MathWorks Licence
% % % Copyright (c) 2016-2017, Giorgio Roffo
% % % All rights reserved.
% % %
% % % Redistribution and use in source and binary forms, with or without
% % % modification, are permitted provided that the following conditions are
% % % met:
% % %
% % %     * Redistributions of source code must retain the above copyright
% % %       notice, this list of conditions and the following disclaimer.
% % %     * Redistributions in binary form must reproduce the above copyright
% % %       notice, this list of conditions and the following disclaimer in
% % %       the documentation and/or other materials provided with the distribution
% % %     * Neither the name of the University of Verona nor the names
% % %       of its contributors may be used to endorse or promote products derived
% % %       from this software without specific prior written permission.
% % %
% % % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% % % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% % % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% % % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% % % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% % % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% % % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% % % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% % % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% % % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% % % POSSIBILITY OF SUCH DAMAGE.
