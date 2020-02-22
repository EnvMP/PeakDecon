function [data,Res]=PeakDeconMain(inputStr, indExport)
%% ========================================================
% [data,Res]=PeakDeconMain(inputStr)
% ========================================================
% Subcode for PeakDecon (ver 0.1)
% ========================================================
% Output
% - data: a structure containig input data
%           Folders:    (matrix) Folder locations where data were located.
%           x:          (vector) x-axis data
%           nSample:    (scalar) the total number of samples
%           Filenames:  (matrix) The file names of all samples
%           X:          (matrix) m x n data matrix.
%                       Row is data corresponding to x axis (n data points)
%                       Column indicates each sample (m samples)
%           groupLabel: (matrix) Labels for group
%                       (extracted from subfolder names)
%           groupInd:   (vector) numeric indices for group (from 1 to n)
%           groupName:  ([]) Empty variable that will be used
%                       for postprocessing (i.e., guiPostprocess)
%           XNorm:      (matrix) Normalized X
% - Res: a structure containg results
%           components: (vector) component numbers that used for simulation
%           splitInd:   (cell) each cell contains vector indices for sample
%                       selection for split half analysis
%                       e.g., Res.splitInd{3}{2} contains the second split
%                       for component 3
%           cosSimilarity:
%                       (vetor) cosine similarity results
%           W:          (cell) each cell contains W for each omponent
%           H:          (cell) each cell contains H for each omponent
% ========================================================
% Input
% - inputStr: a structure variable containing model parameters, including
%           folder:     the data file location (character vector)
%           nComp:      No. of components (numeric vector)
%           norm:       1 -> Normalize each data by its maximum value
%           initGuess:  a numeric vector for initial guess
%                       1 -> Gauss Kernel PCA initialization
%                       0 -> Random initialization
%           HALS:       a numeric vector to decide HALS implementation
%                       1 -> Run HALS after the preceding initialization
%                       0 -> No run of HALS
%           model:      a numeric vector for the selection of NMF algorithm
%                       1 -> ANLS-BPP
%                       2 -> ANLS-ASGIVENS
%                       3 -> ANLS-ASGROUP
%                       4 -> ALS
%                       5 -> HALS
%                       6 -> MU
%           alpha:      Regularization parameters for W (a numeric vector)
%           beta:       Regularization parameters for H (a numeric vector)
%           splitYes:   a numeric vector to decide split half analysis run
%                       1 -> Run split half analysis
%                       0 -> No split half analysis
%           splitType:  1 -> Sorting data based on the Euclidean distance
%                            followed by alternation
%                       2 -> Sorting data based on the Spearman rank
%                            followed by alternation
%                       3 -> Alteration [1 3 ...] vs. [2 4 ...]
% indExport: Index for selecting mode of code run: NMF run vs Export
%                       1 -> Export
%                       Otherwise -> Run NMF algorihtm
% ========================================================
% == Version history ==
% 2/21/2020: ver. 0.1
% ========================================================
% Minkyu Park

global data Res method

%% Parameters that were not incoporated in the input structure
tolHALS=1e-3;       % Tolerance value for HALS algorithm for initialization
minIterHALS=20;     % Minimum iteration for HALS for initialization
maxIterHALS=500;    % Maximum iteration for HALS for initialization

tolModel=1e-5;      % Tolerance value for the main NMF algorithm selected
minIterModel=20;    % Minimum iteration for the main NMF algorithm selected
maxIterModel=1000;  % Maximum iteration for the main NMF algorithm selected

randDefault=1;      % Use default random for reproducibility
                    
if indExport~=1
    %% Load and check folder names for data loading
    fileNames	= ls([inputStr.folder '\*.csv']);
    folderNames = ls([inputStr.folder '\*.']);
    folderNames(1:2, :)=[];  
    
    % To clear up subfolders that do not contain CSV data files
    if ~isempty(folderNames)
        indExist=zeros(1,size(folderNames,1));
        for i=1:size(folderNames,1)
            fileNamesSub{i}=ls([inputStr.folder '\' deblank(folderNames(i,:)) '\*.csv']);
            if ~isempty(fileNamesSub{i})
                indExist(i)=1;
            end
        end
    elseif ~isempty(fileNames)
        indExist=1;
        data=struct('Folders', inputStr.folder);
        folderNames = '.';
        fileNamesSub{1}=fileNames;
        groupLabel=folderNames;
    end
    
    % Conditions to check data files are appropriately input.
    if sum(indExist)==0
        msgbox('No CSV data files or subfolders with data files exist in the selected folder location.')
        return
    elseif and(~isempty(fileNames), ~strcmp(folderNames, '.'))
        msgbox(['Both the designated folder and its subfolders seem to have data files (CSV).' newline ...
            'Place your data files either in subfolders only or in the designated folder only'])
        return
    elseif and(sum(indExist)~=0,~strcmp(folderNames,'.'))
        folderNames(indExist~=1,:)=[];
        fileNamesSub(indExist~=1)=[];
        data=struct('Folders', [repmat([inputStr.folder '\'], size(folderNames,1),1) folderNames]);
        groupLabel=folderNames;
    end

    dataAll=[];         % Variable to store all data points
    groupInd=[];        % Indices of group.
    
    %% Load data files
    for i=1:size(folderNames,1)

        for j=1:size(fileNamesSub{i},1)

            fileName=[inputStr.folder '\' deblank(folderNames(i,:)) '\' deblank(fileNamesSub{i}(j,:))];

            dataEach = csvread(fileName,1,0);

            if and(i==1, j==1)
                if and(size(dataEach,1)>1000, inputStr.initGuess==1)
                    opts.Interpreter = 'default';
                    opts.Default = 'Yes';
                    selBox = questdlg(['Data points are >1000 (n=' num2str(size(dataEach,1)) ')' newline ...
                        'This can cause exremely long calculation time' newline ...
                        'due to a large size of matrices for Gaussian kernel PCA.' newline ...
                        'Use of guiDataReduction is recommended.' newline ...
                        'Do you still want to run with Gaussian kernel PCA?'],'Data reduction?',...
                        'Yes','No',opts);
                    if strcmp(selBox, 'No')
                        data.X=[];
                        Res=[];
                        return;
                    end
                end
                indBackslash=regexpi(inputStr.folder, '\');
                fileNameDiary=[inputStr.folder '\log_' inputStr.folder(indBackslash(end-1)+1:indBackslash(end)-1) '.txt'];
                warning off
                delete(fileNameDiary)
                        % Delete the existing diary file if there is.
                warning on
                diary(fileNameDiary)
                        % To save a log txt file
                fprintf('<strong>========================================================</strong>\n')
                fprintf('<strong>= Load data files</strong>\n')
                fprintf('<strong>========================================================</strong>\n')
                fprintf('.............\n')

                x = dataEach(:,1);
            end

            fprintf('%s is read and processed.\n', ['.' fileName(length(inputStr.folder)+1:end)])
            dataAll=[dataAll dataEach(:,2)];
            groupInd=[groupInd i];
        end
    end

    %% Create a structure of data
    nComp=inputStr.nComp;
    sqText={'st', 'nd'};

    for i=1:length(fileNamesSub)
        [nRow(i), nCol(i)]=size(fileNamesSub{i});
    end
    
    nRow = [0 nRow];
    fileNamesAll = repmat(' ', sum(nRow), max(nCol));
    for i=1:numel(fileNamesSub)
        fileNamesAll(sum(nRow(1:i))+1:sum(nRow(1:i+1)), 1:size(fileNamesSub{i},2))=fileNamesSub{i};
    end

    data = setfield(data, 'x', x);
    data = setfield(data, 'nSample', size(dataAll,2));
    data = setfield(data, 'Filenames', fileNamesAll);
    data = setfield(data, 'X', dataAll);
    data = setfield(data, 'groupLabel', groupLabel);
    data = setfield(data, 'groupInd', groupInd);
    data = setfield(data, 'groupName', []);
    Res = struct('components', nComp);

    fprintf('Data files loaded.\n')
    fprintf('<strong>Data structure was created.</strong>\n\n')

    %% NMF & Split half analysis

    % Normalization
    max_val=max(data.X);
    data.XNorm(:,:)=data.X./repmat(max_val,size(data.X,1),1);

    if inputStr.norm==1
        X=data.XNorm;
    else
        X=data.X;
    end

    if inputStr.splitYes==1
        jVal=1:2;
    else
        jVal=1;
    end
    
    % Generate a cell with component labels (C1, ... Cn) for legends of a plot
    plotLegend=[repmat('C', max(nComp),1) num2str([1:max(nComp)]')];
    plotLegend=mat2cell(plotLegend,ones(size(plotLegend,1),1),size(plotLegend,2));

    hFigComp=figure;
    for i=nComp
        fprintf('<strong>========================================================</strong>\n')
        fprintf('<strong>= NMF simulation for ALL samples: %i component. N=%i</strong>\n', i, data.nSample)
        fprintf('<strong>========================================================</strong>\n')
        
        indSumCs=0; % To count the number of cosine similarity combinations
        for j=jVal

            regW = [inputStr.alpha 0]; regH = [0 inputStr.beta];


            if inputStr.initGuess ==1
                kpca=gkpca(X,i);
                init_W=kpca.coeff(:,1:i);
            else
                if randDefault==1
                    rng default
                end
                init_W=rand(size(X,1),i);
            end
            if randDefault==1
                rng default
            end
            init_H=rand(i,size(X,2));
            
            % Selection of main NMF algorithm
            switch inputStr.model
                case 1
                    method='anls_bpp';
                case 2
                    method='anls_asgivens';
                case 3
                    method='anls_asgroup';
                case 4
                    method='als';
                case 5
                    method='hals';
                case 6
                    method='mu';
            end

            % Selection of method to calculate distance for split half analysis
            distCrt=0;
            switch inputStr.splitType
                case 1
                    distType='Euclidean'; 
                case 2
                    distType='Spearman';
                case 3
                    distCrt=1;
                    distType='Alternating';
            end

            if j==1
                % NMF run for all data files
                if inputStr.initGuess==1
                    fprintf('<strong>Initial guess: Gaussian Kernel PCA</strong>\n')
                    if inputStr.HALS==1
                        fprintf('<strong>HALS with the Gaussian Kernel PCA</strong> result as initial W was implemented.\n')
                    end
                else
                    fprintf('<strong>Initial guess: Randomized values</strong>\n')
                    if inputStr.HALS==1
                        fprintf('<strong>HALS with randomized initial W</strong> was implemented.')
                    end
                end

                fprintf('<strong>Nonnegative matrix factorization</strong> with <strong>%s</strong> is running for ALL samples\n', upper(method))
                fprintf('Regularization parameters <strong>alpha=%1.2f, beta=%1.2f</strong> were used.\n', inputStr.alpha, inputStr.beta)

                if inputStr.HALS==1
                    % Initialization with HALS
                    [W_all{i},H_all{i},iters_all{i}, HIS_all{i}] = nmf(X,i,'verbose', 1, 'method','hals', 'TOL', tolHALS, 'MIN_ITER', minIterHALS, 'MAX_ITER', maxIterHALS ,...
                        'REG_W', regW, 'REG_H', regH, 'INIT', struct('W',init_W,'H',init_H));   
                else
                    % No HALS initialization
                    W_all{i}=init_W;
                    H_all{i}=init_H;
                end
                [W_all{i},H_all{i},iters_all{i}, HIS_all{i}] = nmf(X,i,'verbose', 1, 'method',method, 'TOL', tolModel, 'MIN_ITER', minIterModel, 'MAX_ITER', maxIterModel ,...
                    'REG_W', regW, 'REG_H', regH, 'INIT', struct('W',W_all{i},'H',H_all{i}));
                if HIS_all{i}.final.iterations==maxIterModel
                    fprintf('Maximum iteration (n=%i) was reached. The model was ', maxIterModel); 
                else
                    fprintf('Model was <strong>Converged</strong>.\n')
                end

                MwAll{i}=sum(W_all{i}.*repmat(data.x,1,size(W_all{i},2)))./sum(W_all{i});
                MnAll{i}=sum(W_all{i})./sum(W_all{i}./repmat(data.x,1,size(W_all{i},2)));

                [~,WInd]=sort(MwAll{i}, 'Ascend');
                W_all{i}=W_all{i}(:,WInd);
                H_all{i}=H_all{i}(WInd,:);
                MwAll{i}=MwAll{i}(WInd);                    % Molecular weight
                MnAll{i}=MnAll{i}(WInd);                    % Number

                poly_disp_all=MwAll{j}./MnAll{j};



                if distCrt==1
                    distInd=1:size(X,2);
                else
                    distMat= squareform(pdist(H_all{i}'./repmat(median(H_all{i}',2),1,size(H_all{i}',2)),distType));
                    [~,distInd] = sort(distMat(1,:));
                end

                if inputStr.splitYes==1
                    subplot(3,numel(nComp),find(nComp==i))
                else
                    nColSubplot=mod(numel(nComp),3);
                    if or(nColSubplot==0, ceil(numel(nComp)/3)>1)
                        nColSubplot=3;
                    end
                    subplot(ceil(numel(nComp)/3),nColSubplot,find(nComp==i))     
                end
                plot(data.x, W_all{i})
                title(['C' num2str(i) ': All'], 'fontweight', 'bold', 'fontsize', 12)
                if i==max(nComp)
                    hLegend=legend(plotLegend(1:max(nComp)));
                end
                fprintf('All samples were processed.\n')


            end

            if inputStr.splitYes==1
                XNew=X;
                XNew(:,distInd(j:2:end))=[];
                distIndCell{i}{j}=distInd(j:2:end);

                if j==1
                    fprintf('<strong>--------------------------------------------------------</strong>\n')
                    fprintf('<strong>[[[[ Split Half Analysis: %s ]]]]</strong>\n', distType)
                    fprintf('<strong>--------------------------------------------------------</strong>\n')
                end
                fprintf('<strong>--------------------------------------------------------</strong>\n')
                fprintf('<strong>Split Half Analysis: %i%s, N=%i </strong>\n', j, sqText{j}, size(XNew,2))
                fprintf('<strong>--------------------------------------------------------</strong>\n')
                fprintf('Split half analysis is in progress......\n')

                if inputStr.initGuess ==1
                    kpca=gkpca(XNew,i);
                    init_W=kpca.coeff(:,1:i);
                else
                    if randDefault==1
                        rng default
                    end
                    init_W=rand(size(XNew,1),i);
                end
                if randDefault==1
                    rng default
                end
                init_H=rand(i,size(XNew,2));


                if inputStr.HALS==1
                    [W{i}{j},H{i}{j},iters{i}{j}, HIS{i}{j}] = nmf(XNew,i,'verbose', 1, 'method', 'hals', 'TOL', tolHALS, 'MIN_ITER', minIterHALS, 'MAX_ITER', maxIterHALS ,...
                        'REG_W', regW, 'REG_H', regH, 'INIT', struct('W',init_W,'H',init_H));
                else
                    W{i}{j}=init_W;
                    H{i}{j}=init_H;
                end
                [W{i}{j},H{i}{j},iters{i}{j}, HIS{i}{j}] = nmf(XNew,i,'verbose', 1, 'method',method, 'TOL', tolModel, 'MIN_ITER', minIterModel, 'MAX_ITER', maxIterModel ,...
                    'REG_W', regW, 'REG_H', regH, 'INIT', struct('W',W{i}{j},'H',H{i}{j}));
                if HIS{i}{j}.final.iterations==maxIterModel
                    fprintf('Maximum iteration (n=%i) was reached. The model was ', maxIterModel); 
                else
                    fprintf('Model was <strong>Converged</strong>.\n')
                end
                
                clear Mw Mn
                Mw{i}{j}=sum(W{i}{j}.*repmat(data.x,1,size(W{i}{j},2)))./sum(W{i}{j});
                Mn{i}{j}=sum(W{i}{j})./sum(W{i}{j}./repmat(data.x,1,size(W{i}{j},2)));

                clear WInd
                
                NormW_all=W_all{i}./repmat(mean(W_all{i}),size(W_all{i},1),1);
                NormW=W{i}{j}./repmat(mean(W{i}{j}(:,i)),1,size(W{i}{j},2),1);
                
                for k=1:i
                    for l=1:i
                        CSMat(k,l)=cosSim(NormW_all(:,k), ...
                            NormW(:,l)) ...
                            *(1-abs(MwAll{i}(k)-Mw{i}{j}(l))/MwAll{i}(k));
                    end
                end
                
                aaa=0;
                bbb=zeros(1,i);
                if sum(isnan(CSMat))==0
                   
                    while aaa<i
                        [~,indMax]=max(CSMat(:));
                        if mod(indMax,i)==0
                            mVal=i;
                        else
                            mVal=mod(indMax,i);
                        end
                        nVal=ceil(indMax/i);

                        if bbb(ceil(indMax/i))==0
                            aaa=aaa+1;
                            bbb(nVal)=mVal;
                            CSMat(mVal,:)=0;
                            CSMat(:, nVal)=0;
                        else
                            CSMat(indMax)=0;
                        end

                    end
                end

                [~,WInd]=sort(bbb);

                W{i}{j}=W{i}{j}(:,WInd);
                H{i}{j}=H{i}{j}(WInd,:);
                Mw{i}{j}=Mw{i}{j}(WInd);
                Mn{i}{j}=Mn{i}{j}(WInd);
                
                poly_disp=Mw{i}{j}./Mn{i}{j};

                % Cosine similiarity calculation
                for k=1:i
                    Cs(k)   = cosSim(W_all{i}(:,k),W{i}{j}(:,k));  
                end
                SumCs(j)=sum(Cs);

                subplot(3,numel(nComp),j*numel(nComp)+find(nComp==i))
                plot(data.x, W{i}{j})
                title(['C' num2str(i) ': Split ' num2str(j)], 'fontweight', 'bold', 'fontsize', 12)
                fprintf('Split half analysis was completed.\n')
                

            end
        end
        

        
        SumCs=0;
        W{i}{j+1}=W_all{i}; % Put W_all in W to calculate all combinations for cosine similarity calculation
        combos=nchoosek([jVal max(jVal)+1],2);
        
        if inputStr.splitYes==1
            for j=1:size(combos,1)
                for k=1:i
                    SumCs=SumCs+cosSim(W{i}{combos(j,1)}(:,k), ...
                        W{i}{combos(j,2)}(:,k));
                end
            end
            W{i}{j+1}=[];
            CsFinal(i)=SumCs/size(combos,1)/i;
            fprintf('Cosine Similiarity value: %1.2f\n', CsFinal(i))
        end
        fprintf('<strong>%i component</strong> simulation was <strong>completed</strong>.\n',i)
        fprintf('<strong>========================================================</strong>\n\n\n')
    end

    if inputStr.splitYes==1
        Res = setfield(Res, 'splitInd', distIndCell);
    else
        Res = setfield(Res, 'splitInd', {NaN});
    end
    
    if inputStr.norm==1
        for i=nComp
            H_all{i}=H_all{i}.*repmat(max_val,i,1);
        end
    end

    if inputStr.splitYes==1
        hFigCos=figure;
        bar(CsFinal(nComp))
        set(gca, 'xticklabel', nComp, ...
            'ylim', [0 1])
        title('Average Cosine Similarity', 'fontweight', 'bold', 'fontsize', 12)
        xlabel('No. of components', 'fontsize', 12)
        ylabel('Similarity index', 'fontsize', 12)
        Res = setfield(Res, 'cosSimilarity', CsFinal);
    end
    tElapsed=toc;
    fprintf('Total calculation time = %1.1f seconds\n ', tElapsed)
    msgbox('Data process is completed!')
    diary off
    fprintf('Log was saved in %s\n', ...
        [inputStr.folder '\log_' inputStr.folder(indBackslash(end-1)+1:indBackslash(end)-1) '.txt'])
    Res.W=W_all;
    Res.H=H_all;
    if inputStr.splitYes==1
        Res.WSplit=W;
        Res.HSplit=H;
    end

else
    indBackslash=regexpi(inputStr.folder, '\');
    if inputStr.HALS==1
        fileName=[inputStr.folder '\' inputStr.folder(indBackslash(end)+1:end) ' HALS-' upper(method) '.mat'];
    else
        fileName=[inputStr.folder '\' inputStr.folder(indBackslash(end)+1:end) ' ' upper(method) '.mat'];
    end
    save(fileName, 'data', 'Res');
    msgbox(['Data was stored in ' fileName])
end

    function kpca = gkpca(X, nComp)
        Kmat = exp(-pdist2(X,X, 'Euclidean').^2/2./1.^2);        
        nSample=size(Kmat,1);
        sumK = sum(Kmat,2);
        Hmat = ones(nSample,nSample)/nSample;
        Kmat = Kmat - Kmat*Hmat - Hmat*Kmat + Hmat*Kmat*Hmat;
        [V, D] = eig(Kmat);
        [~,indEig]=sort(diag(D),'descend');
        V=V(:,indEig);
        normV=sqrt(sum(V.^2));
        V=V./repmat(normV,size(V,1),1);
        V=V(:,1:nComp);
        kpca.coeff=Kmat*V;

    end

    function Cs=cosSim(x,y)
        xy   = dot(x,y);
        normx   = norm(x);
        normy   = norm(y);
        normxy = normx*normy;
        Cs   = xy/normxy;
    end
end