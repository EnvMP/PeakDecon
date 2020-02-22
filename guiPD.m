function guiPD()
%% ========================================================
%  PeakDecon ver 0.1
%  GUI-based Peak Deconvolution with NMF algorithm
% ========================================================
% function guiPD()  
% ========================================================
%   This software was intested to deconvolute potentially
% overlapped chromatograms of size exclusion chromagoraphy.
% Regardless of the developer's intention, this software
% may be used for deconvoluting any 2D signals composed of
% multiple components.
%   This software package also includes
%     - guiPre.m
%           : To separate a XLSX files into CSV files corresponding
%             to each sample.
%           : To reduce the size of data
%     - guiPost.m
%           : To separate postprocess data (mat file)
% ========================================================
% Note for selecting folders to load data files:
%     - Software provids two options of loading data files:
%       1) Loading a folder containing data files
%       2) Loading a folder containing multiple subfolders
%          that contain data files.
%   In the second case, the root folder of subfolders should
% not contain any CSV files 
% ========================================================
% == Version history ==
% 2/21/2020: ver. 0.1
%               - The first version
% ========================================================
% Coded by Minkyu Park

clear all; close all hidden; clc
global indExport Res data inputStr

indExport=0;        % Index for export. indExport=1 -> Enable export

% Get the screen size to put the software panel at the center of the screen
ss = get(0,'screensize');
wHorSize=500;
wVerSize=450;
posCenter=[ss(3)/2-wHorSize/2 ss(4)/2-wVerSize/2];

% Create main figure
h_main=figure;
set(gcf, 'pos', [posCenter wHorSize wVerSize], ...
    'Numbertitle', 'off', 'name', 'PeakDecon: Peak Deconvolution GUI. ver. 0.1');
set(h_main, 'handlevisibility', 'off');
                    % Set handle visibility off to facilitate use 
                    % of 'close all' command to close all the figuresexcept
                    % for this main software window
                    
% UI Button Group 1: Load a folder location
h_group1=uibuttongroup(h_main, 'Title', 'Load a folder', ...
    'Fontsize', 11, 'pos', [0.05 0.88 0.9 0.11], ...
    'ForegroundColor', 'b');
h_group1_text1=uicontrol(h_group1, 'Style', 'text', 'String', 'Folder:', ...
    'Pos', [5 0 50 25], 'Fontsize', 10);
h_group1_text2=uicontrol(h_group1, 'Style', 'edit', 'String', '', ...
    'Pos', [60 3 330 25], 'Fontsize', 10);
h_group1_but1=uicontrol(h_group1, 'Style', 'pushbutton', 'String', 'Load', 'Fontsize', 10, ...
    'pos', [395 3 50 25]);
h_group1_but1.Callback=@getLoc;

% UI Button Group 2: No. of components
h_group2=uibuttongroup(h_main, 'Title', 'Number of components', ...
    'Fontsize', 11, 'pos', [0.05 0.77 0.9 0.1], ...
    'ForegroundColor', 'b');
h_group2_text1=uicontrol(h_group2, 'Style', 'text', 'String', 'No of Components:', ...
    'Pos', [5 0 135 25], 'Fontsize', 10);
h_group2_text2=uicontrol(h_group2, 'Style', 'edit', 'String', '2:9', ...
    'Pos', [145 3 50 25], 'Fontsize', 10);
h_group2_text3=uicontrol(h_group2, 'Style', 'text', 'String', 'e.g., For 2,4,5 & 6 comp, type "2 4:6"', ...
    'Pos', [195 0 250 25], 'Fontsize', 10);

% UI Button Group 3: Data normalization
h_group3=uibuttongroup(h_main, 'Title', 'Data Normalization parameters', ...
    'Fontsize', 11, 'pos', [0.05 0.66 0.9 0.11], ...
    'ForegroundColor', 'b');
h_group3_text1=uicontrol(h_group3, 'Style', 'text', 'String', 'Data normalization with maximum', ...
    'Pos', [5 0 200 25], 'Fontsize', 10);
h_group3_checkbox1=uicontrol(h_group3, 'Style','checkbox','String','Yes', ...
    'Pos', [230 0 50 25],'HandleVisibility','off', 'Fontsize', 10, ...
    'Value', false);

% UI Button Group 4: Initialization parameters
h_group4=uibuttongroup(h_main, 'Title', 'Initialization parameters', ...
    'Fontsize', 11, 'pos', [0.05 0.49 0.9 0.16], ...
    'ForegroundColor', 'b');
h_group4_text1=uicontrol(h_group4, 'Style', 'text', 'String', 'Initial Guess for W matrix:', ...
    'Pos', [5 27 160 25], 'Fontsize', 10);
h_group4_radio1=uicontrol(h_group4, 'Style','radiobutton','String','Gauss Kernel PCA', ...
    'Pos', [170 30 150 25],'HandleVisibility','off', 'Fontsize', 10, ...
    'value', 0);
h_group4_radio2=uicontrol(h_group4, 'Style','radiobutton','String','Random', ...
    'Pos', [330 30 120 25],'HandleVisibility','off', 'Fontsize', 10, ...
    'value', 1);

h_group4_text2=uicontrol(h_group4, 'Style', 'text', 'String', 'Subsequent initialization with HALS', ...
    'Pos', [0 7 220 25], 'Fontsize', 10);
h_group4_checkbox1=uicontrol(h_group4, 'Style','checkbox','String','Yes', ...
    'Pos', [230 10 135 25],'HandleVisibility','off', 'Fontsize', 10, ...
    'Value', true);

% UI Button Group 5: Nonnegative Matrix Factorization (NMF) run parameters
h_group5=uibuttongroup(h_main, 'Title', 'Nonnegative Matrix Factorization (NMF) run parameters', ...
    'Fontsize', 11, 'pos', [0.05 0.29 0.9 0.19], ...
    'ForegroundColor', 'b');
h_group5_text1=uicontrol(h_group5, 'Style', 'text', 'String', 'Model:', ...
    'Pos', [5 44 45 25], 'Fontsize', 10);
h_group5_radio{1}=uicontrol(h_group5, 'Style','radiobutton','String','ANLS-BPP', ...
    'Pos', [70 47 100 25],'HandleVisibility','off', 'Fontsize', 10);
h_group5_radio{2}=uicontrol(h_group5, 'Style','radiobutton','String','ANLS-ASGIVENS', ...
    'Pos', [170 47 130 25],'HandleVisibility','off', 'Fontsize', 10);
h_group5_radio{3}=uicontrol(h_group5, 'Style','radiobutton','String','ANLS-ASGROUP', ...
    'Pos', [300 47 130 25],'HandleVisibility','off', 'Fontsize', 10);
h_group5_radio{4}=uicontrol(h_group5, 'Style','radiobutton','String','ALS', ...
    'Pos', [70 28 130 25],'HandleVisibility','off', 'Fontsize', 10);
h_group5_radio{5}=uicontrol(h_group5, 'Style','radiobutton','String','HALS', ...
    'Pos', [170 28 130 25],'HandleVisibility','off', 'Fontsize', 10);
h_group5_radio{6}=uicontrol(h_group5, 'Style','radiobutton','String','MU', ...
    'Pos', [300 28 130 25],'HandleVisibility','off', 'Fontsize', 10);
h_group5_text2=uicontrol(h_group5, 'Style', 'text', 'String', 'Model parameters :', ...
    'Pos', [5 2 120 25], 'Fontsize', 10);
h_group5_text2=uicontrol(h_group5, 'Style', 'text', 'String', 'alpha:', ...
    'Pos', [130 2 50 25], 'Fontsize', 10);
h_group5_text3=uicontrol(h_group5, 'Style', 'edit', 'String', '0', ...
    'Pos', [180 5 30 25], 'Fontsize', 10);
h_group5_text4=uicontrol(h_group5, 'Style', 'text', 'String', 'beta:', ...
    'Pos', [260 2 50 25], 'Fontsize', 10);
h_group5_text5=uicontrol(h_group5, 'Style', 'edit', 'String', '0', ...
    'Pos', [310 5 30 25], 'Fontsize', 10);

% UI Button Group 6: Split Half Analysis
h_group6=uibuttongroup(h_main, 'Title', 'Split Half Analysis parameters', ...
    'Fontsize', 11, 'pos', [0.05 0.07 0.9 0.23], ...
    'ForegroundColor', 'b');
h_group6_text1=uicontrol(h_group6, 'Style', 'text', 'String', 'Run split half analysis?', ...
    'Pos', [0 62 145 25], 'Fontsize', 10);
h_group6_checkbox1=uicontrol(h_group6, 'Style','checkbox','String','Yes', ...
    'Pos', [230 65 135 25],'HandleVisibility','off', 'Fontsize', 10, ...
    'Value', true);
h_group6_text1=uicontrol(h_group6, 'Style', 'text', 'String', 'Type: ', ...
    'Pos', [50 40 50 25], 'Fontsize', 10);
h_group6_radio{1}=uicontrol(h_group6, 'Style', 'radiobutton', 'String', 'Euclidean distance sort followed by alternation', ...
    'Pos', [100 42 320 25], 'Fontsize', 10);
h_group6_radio{2}=uicontrol(h_group6, 'Style', 'radiobutton', 'String', 'Spearman rank sort followed by alternation', ...
    'Pos', [100 22 320 25], 'Fontsize', 10);
h_group6_radio{3}=uicontrol(h_group6, 'Style', 'radiobutton', 'String', 'Alteration [1:3: ...] vs [2:4: ...]', ...
    'Pos', [100 2 320 25], 'Fontsize', 10);

% UI Button Group 7: Push buttons for 'Run' and 'Export (Generate NMFRes.mat)'
h_group7_but1 = uicontrol(h_main, 'Style', 'pushbutton', 'String', 'Start', 'Fontsize', 10, 'FontWeight', 'bold', ...
    'pos', [100 4 100 25]);
h_group7_but1.Callback = @callbackRun;
h_group7_but2=uicontrol(h_main, 'Style', 'pushbutton', 'String', 'Export (Generate mat file)', 'Fontsize', 10, 'FontWeight', 'bold', ...
    'pos', [220 4 200 25]);
h_group7_but2.Callback = @callbackReport;

%% Functions
    % Function to open dialogue for directory selection
    function getLoc(h_main, event)
        dir_loc=uigetdir('title', ['Select folder where data are stored ' ...
            'or folder containg subfolders with data files.']);
        set(h_group1_text2, 'String', dir_loc)
    end

    % Function to run subfunction for NMF algorithm
    function callbackRun(h_main, event)
        if isempty(get(h_group1_text2, 'string'))
            msgbox('Folder containing data files was not selected.')
        else
            inputStr.folder=get(h_group1_text2, 'string');
            inputStr.nComp=str2num(get(h_group2_text2, 'string'));
            inputStr.norm=get(h_group3_checkbox1, 'value');
            inputStr.initGuess=get(h_group4_radio1, 'value');
            inputStr.HALS=get(h_group4_checkbox1, 'value');
            for i=1:6       
                if get(h_group5_radio{i}, 'value')==1
                    inputStr.model=i;
                end
            end
            inputStr.alpha=str2num(get(h_group5_text3, 'string'));
            inputStr.beta=str2num(get(h_group5_text5, 'string'));
            inputStr.splitYes=get(h_group6_checkbox1, 'value');
            for i=1:3
                if get(h_group6_radio{i}, 'value')==1
                    inputStr.splitType=i;
                end
            end

            clc
            close all
                    % Close all the windows and clear screen 
                    % before running NMF algorithm
            tic
            indExport=0;
            [data,Res]=PeakDeconMain(inputStr, indExport);
            if isempty(data.X)
                clear Res data
            else
                assignin('base', 'Res', Res);
                assignin('base', 'data', data);
                indExport=1;
            end
            
        end
    end

    % Function to export processed data
    function callbackReport(h_main, event)
        if indExport==1
            [data,Res]=PeakDeconMain(inputStr, indExport);
        else
            msgbox('No data present for export.')
        end
    end
 
end