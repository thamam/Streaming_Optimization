function Out = wavesCD(tK,params,opts)
global DEBUG_MODE
%% Debug Variables
% opts=WCD_opts;
% DEBUG_MODE = false;
%% Figure maximal distance of propogation - number of frames
W= ceil(params.duration/opts.sectionLength);  % number of sections
opts.W = W;
[params] = initializeSectionParams(params,opts,tK);
%% Initialize
% alphaM=zeros(W,d);
sec_nIndtMat = params.sec_nIndtMat;
tK_sec_ind = params.tK_sec_ind;
k=1;
% alphaValM = cell(opts.MaxIter,1); alphaValM{1}=ones(size(params.nInd));
alphaValM = zeros(opts.MaxIter,params.d);
alphaValM(1,:)= 1;
% alpha_Kp1 = alphaValM(1,:);
alpha_diff_norm = zeros(opts.MaxIter,1);
alpha_diff_norm(1) = 1;
activeFramesMat =  genWavesMat( W,opts.MaxIter,opts.ActivationMode );
%% Iterative Solution
% activeFramesMat = ones(opts.MaxIter,4 );
nInd = params.nInd;
indMatList = extractIndList(params.secSupportMat);
astrixSpaceLine =[ repmat('*',[1,300])];
fprintf(' %s \n Starting Optimization for %s \n %s \n', astrixSpaceLine, opts.ActivationMode,astrixSpaceLine );
while (k<= opts.MaxIter  && alpha_diff_norm(k) > params.iterThs)% not converge
    %% Select Active frames for comming iteration according to activation matrix
    activeFramesMask  = activeFramesMat(k,:);
    activeFramesList = find(activeFramesMask==1);
    temp_alphaValM = zeros(1,params.d);
    for w=activeFramesList
        % Extract Alpha to be included in current iteration opt. Problem
        iterAlpha = intersect(indMatList(w,:),nInd);
        %         iterAlphaMask = ismember(nInd,iterAlpha);
        % Sort which alpha are fixed
        iterFixedAlpha = iterAlpha(iterAlpha<sec_nIndtMat(w,1) | iterAlpha>sec_nIndtMat(w,2));
        iterFixedAlphaMask = (iterAlpha<sec_nIndtMat(w,1) | iterAlpha>sec_nIndtMat(w,2)); 
        % Extract values of Alpha that are fixed
        iterFixedValues = alphaValM(k,iterFixedAlpha);
        % Solve
        alphaIndUpdate = sec_nIndtMat(w,1):sec_nIndtMat(w,2);
        alphaIndUpdateMask = ismember(nInd,alphaIndUpdate);
        %%       
        CVX_RES =  CVX_GCD4(tK(tK_sec_ind(w,1):tK_sec_ind(w,2)),params,w,...
            iterAlpha,iterFixedAlpha,iterFixedValues,iterFixedAlphaMask);
        temp_alphaValM( alphaIndUpdateMask)= CVX_RES;%(~iterFixedAlphaMask);
    end
     % Solution at iteration (k+1) is  alpha(k) with updates in actiuve cooridnates  
    alpha_Kp1 = ~logical(temp_alphaValM).*alphaValM(k,:)+temp_alphaValM;
    alphaValM(k+1,:) = alpha_Kp1;
    k=k+1; %advance counter
    if isfield(opts,'alpha_hat_cvx')
        alpha_diff_norm(k) = norm(alphaValM(k,:) - opts.alpha_hat_cvx.')/params.d.^2;
    else
        alpha_diff_norm(k) = norm(alphaValM(k,:) - alphaValM(k-1,:))/params.d.^2;
    end
    
    fprintf('itr(k):%2i, Error:%4.2e , Active Frames:[ %s ]\n',k-1, alpha_diff_norm(k),num2str(activeFramesList) );

end
%% Assign Output values
Out.alphaValM = alphaValM;
Out.SectionIterationCount=k;
Out.tK_sec_ind=tK_sec_ind;
Out.alpha_diff_norm = alpha_diff_norm;
Out.activeFramesMat = activeFramesMat;
Out.ActivationMode = opts.ActivationMode;
