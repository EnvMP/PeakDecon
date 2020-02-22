function guiPre()
%% ========================================================
%  GUI for PeakDecon Pre-processor (ver 0.1)
%  packaged with PeakDecon(ver. 0,1)
% ========================================================
% function guiDataReduction()
% ========================================================
% This software was intested to process data prior to
% use of the peak deconvolution software.
% This includes mainly two features:
%   1) Separation of a XLSX file containing data of
%       multiple samples into separate CSV files
% 	2) Reduce data size to speed up deconvolution process
% ========================================================
% == Version history ==
% 2/21/2020: ver. 0.1
% ========================================================
% Minkyu Park

clear all; close all hidden; clc
global dataGraph listFold fileLocCSV 

% Get the screen size to put the software panel at the center of the screen
ss = get(0,'screensize');
wHorSize=500;
wVerSize=450;
posCenter=[ss(3)/2-wHorSize/2 ss(4)/2-wVerSize/2];

% List of No. fold for data reduction;
listFoldMatPreset=[1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 200 500 1000];

% Create main figure
h_main=figure;
set(gcf, 'pos', [posCenter wHorSize wVerSize], ...
    'Numbertitle', 'off', 'name', 'PeakDecon Pre-processor ver. 0.1');
set(h_main, 'handlevisibility', 'off');

% UI Button Group 1: Data File Separation
h_group1=uibuttongroup(h_main, 'Title', '1. Data File Separation (a XLSX -> multiple CSVs)', ...
    'Fontsize', 11, 'pos', [0.05 0.77 0.9 0.22], ...
    'ForegroundColor', 'b');
h_group1_text1=uicontrol(h_group1, 'Style', 'text', 'String', ['Select a file ' ...
    'containing data of all samples (xlsx or xls)'], ...
    'Pos', [30 56 320 25], 'Fontsize', 10);
h_group1_but1=uicontrol(h_group1, 'Style', 'pushbutton', 'String', 'Help', ...
    'Pos', [385 58 50 23], 'Fontsize', 10);
h_group1_but1.Callback=@helpDataformat;
h_group1_text2=uicontrol(h_group1, 'Style', 'edit', 'String', '', ...
    'Pos', [30 35 330 23], 'Fontsize', 10);
h_group1_but2=uicontrol(h_group1, 'Style', 'pushbutton', 'String', 'Load', 'Fontsize', 10, ...
    'pos', [385 34 50 25]);
h_group1_but2.Callback=@loadSeparateFile;
h_group1_but3=uicontrol(h_group1, 'Style', 'pushbutton', 'String', 'Run', 'Fontsize', 10, ...
    'pos', [210 5 50 25]);
h_group1_but3.Callback=@runSeparate;

% UI Button Group 2: Data Size Reduction
h_group2=uibuttongroup(h_main, 'Title', '2. Data Size Reduction', ...
    'Fontsize', 11, 'pos', [0.05 0.02 0.9 0.75], ...
    'ForegroundColor', 'b');
hAxes1=axes(h_group2, 'pos', [0.08 0.68 0.88 0.22]);
h_group2_text1=uicontrol(h_group2, 'Style', 'text', 'String', 'Select a CSV file for preview', 'Fontsize', 10, ...
    'pos', [10 293 210 20]);
h_group2_but1=uicontrol(h_group2, 'Style', 'pushbutton', 'String', 'Load', 'Fontsize', 10, ...
    'pos', [210 293 50 23]);
h_group2_but1.Callback=@loadGraph;
h_group2_text2=uicontrol(h_group2, 'Style', 'text', 'String', 'No. of data points: ', 'Fontsize', 10, ...
    'pos', [10 175 210 25]);
h_group2_text3=uicontrol(h_group2, 'Style', 'text', 'String', '', 'Fontsize', 10, ...
    'pos', [190 175 100 25]);

h_group2_sub1=uibuttongroup(h_group2, 'Title', 'Size reduction & Baseline subtraction parameters', ...
    'Fontsize', 11, 'pos', [0.02 0.02 0.96 0.53]);
h_group2_sub1_text1=uicontrol(h_group2_sub1, 'Style', 'text', 'String', 'Original x range:', 'Fontsize', 10, ...
    'pos', [10 120 60 25]);
h_group2_sub1_text2=uicontrol(h_group2_sub1, 'Style', 'text', 'String', '0', 'Fontsize', 10, ...
    'pos', [140 120 60 25]);
h_group2_sub1_text3=uicontrol(h_group2_sub1, 'Style', 'text', 'String', '-', 'Fontsize', 10, ...
    'pos', [200 120 25 25]);
h_group2_sub1_text4=uicontrol(h_group2_sub1, 'Style', 'text', 'String', '0', 'Fontsize', 10, ...
    'pos', [235 120 60 25]);
h_group2_sub1_text5=uicontrol(h_group2_sub1, 'Style', 'text', 'String', 'New x range:', 'Fontsize', 10, ...
    'pos', [5 100 86 25]);
h_group2_sub1_text6=uicontrol(h_group2_sub1, 'Style', 'edit', 'String', '0', 'Fontsize', 10, ...
    'pos', [140 100 60 25]);
h_group2_sub1_text7=uicontrol(h_group2_sub1, 'Style', 'text', 'String', '-', 'Fontsize', 10, ...
    'pos', [200 100 25 25]);
h_group2_sub1_text8=uicontrol(h_group2_sub1, 'Style', 'edit', 'String', '0', 'Fontsize', 10, ...
    'pos', [235 100 60 25]);
h_group2_sub1_text9=uicontrol(h_group2_sub1, 'Style', 'pushbutton', 'String', 'Estimate data No. ', 'Fontsize', 10, ...
    'pos', [300 100 120 25]);
h_group2_sub1_text9.Callback=@estimateNo;
h_group2_sub1_text10=uicontrol(h_group2_sub1, 'Style', 'text', 'String', 'x-fold for reduction: ', 'Fontsize', 10, ...
	'pos', [5 75 130 25]);
h_group2_sub1_text11=uicontrol(h_group2_sub1, 'Style', 'popupmenu', 'String', ' ', 'Fontsize', 10, ...
    'pos', [140 75 60 25]);
h_group2_sub1_text12=uicontrol(h_group2_sub1, 'Style', 'text', 'String', 'Baseline subtraction using:', 'Fontsize', 10, ...
    'pos', [5 47 170 25]);
h_group2_sub1_text13=uicontrol(h_group2_sub1, 'Style', 'radiobutton', 'String', 'averaged y values btw x', 'Fontsize', 10, ...
    'pos', [30 30 160 25]);
h_group2_sub1_text14=uicontrol(h_group2_sub1, 'Style', 'edit', 'String', '0', 'Fontsize', 10, ...
    'pos', [200 30 60 25]);
h_group2_sub1_text15=uicontrol(h_group2_sub1, 'Style', 'text', 'String', '-', 'Fontsize', 10, ...
    'pos', [275 30 10 25]);
h_group2_sub1_text16=uicontrol(h_group2_sub1, 'Style', 'edit', 'String', '0', 'Fontsize', 10, ...
    'pos', [300 30 60 25]);
h_group2_sub1_text17=uicontrol(h_group2_sub1, 'Style', 'radiobutton', 'String', 'manual input', 'Fontsize', 10, ...
    'pos', [30 5 130 25], 'value', 1);
h_group2_sub1_text18=uicontrol(h_group2_sub1, 'Style', 'edit', 'String', '0', 'Fontsize', 10, ...
    'pos', [200 5 60 25]);

h_group2_sub1_but2=uicontrol(h_group2_sub1, 'Style', 'pushbutton', 'String', 'Run', 'Fontsize', 12, ...
    'pos', [370 5 50 50], 'fontweight', 'bold');
h_group2_sub1_but2.Callback=@runReduction;





%% Functions
    function loadSeparateFile(h_main, event)
        [fileName, fileLoc]=uigetfile({'*.xlsx; *.xls'},'Select file containing data of multiple samples');
        set(h_group1_text2, 'String', [fileLoc fileName])
    end

    function helpDataformat(h_main, event)
        h_FigHelp=figure;
        set(gcf, 'pos', [posCenter+10 wHorSize wVerSize], ...
            'Numbertitle', 'off', 'name', 'Help for data structure');
        h_help1=uibuttongroup(h_FigHelp, 'Title', 'Help for data structure', ...
            'Fontsize', 11, 'pos', [0.05 0.05 0.9 0.9], ...
            'ForegroundColor', 'b');
        h_help1_text1=uicontrol(h_help1, 'Style', 'Text', 'String', ...
            ['1. Data File Separation" section is to separate an excel file ' ...
            '(xlsx or xls) into CSV files that correspond to each sample. ', ...
            'For example, if you have n samples with m data points per sample, ', ...
            'this separation process results in n CSV files that has a matrix ', ...
            'm x 2. Noted that the first column is an indenpendent variable ', ...
            'and the secnd column contains signals from a detector. Excel file should ', ...
            'have the first row as a header that is going to be used for CSV filenames. '...
            'An example is a case in which XLSX data has a 100 x 4 matrix.'], ...
            'pos', [10 220 400 160], 'fontsize', 10);
        h_help1_text2=uicontrol(h_help1, 'Style', 'Text', 'String', 'Example', ...
            'pos', [10 230 65 20], 'fontsize', 10, 'fontweight', 'bold');
        h_help1_text3=uicontrol(h_help1, 'Style', 'Text', 'String', ['XLSX with 100 x 4' ...
            newline '(col1: time, col2-4: samples)'], ...
            'pos', [10 185 200 50], 'fontsize', 10);
        h_help1_text4=uicontrol(h_help1, 'Style', 'Text', 'String', ['-------->' ...
            newline 'Separation'], ...
            'pos', [210 185 65 50], 'fontsize', 10);
        h_help1_text5=uicontrol(h_help1, 'Style', 'Text', 'String', ['3 CSV files' ...
            newline '3 samples x [100 x 2]'], ...
            'pos', [275 185 200 50], 'fontsize', 10);
        Time=[1:100]';
        Sample1=normpdf(Time,50,3)+0.5*normpdf(Time,40,3); ...
        Sample2=0.3*normpdf(Time,50,3)+0.6*normpdf(Time,40,3);
        Sample3=0.4*normpdf(Time,50,3)+0.2*normpdf(Time,40,3);
        dataTable=[Time, Sample1, Sample2, Sample3];
        dataTableCell=num2cell(dataTable);
        h_table=uitable( ...
            'pos', [30 30 wHorSize/2+10 wVerSize*0.4], ...
            'Data', dataTable, ...
            'ColumnWidth', {35,55,55,55}, ...
            'ColumnName', {'Time','Sample1','Sample2','Sample3'});
        hAxes2=axes(h_help1, 'pos', [0.68 0.12 0.26 0.35]);
        plot(hAxes2, Time, dataTable(:,2:end))
        xlabel('Time')
        ylabel('Signal')
        legend({'Sample1', 'Sample2', 'Sample3'})
    end

    function runSeparate(h_main, event)
        fullPath=get(h_group1_text2, 'String');
        indBackslash=regexpi(fullPath, '\');

        if isempty(indBackslash)
            msgbox('No file location was loaded.')
            return;
        else
            indDot=regexpi(fullPath, '\.');
            folderName=fullPath(indBackslash(end)+1:indDot(end)-1);
        end

        fprintf('Data file is being loaded....\n')
        [dataXlsx, txt]=xlsread(fullPath);
        fprintf('Excel file %s was loaded.\n', fullPath)  

        csvLoc=[fullPath(1:indBackslash(end)) folderName];

        warning off
        mkdir(csvLoc);
        fprintf('Folder %s was created.\n',csvLoc)

        warning on
        for i=2:size(dataXlsx,2)
            csvwrite([csvLoc '/' txt{i} '.csv'], [dataXlsx(:,1) dataXlsx(:,i)]);
            fprintf('%s was created.\n', [txt{i} '.csv'])
        end
        fprintf('Data file separation completed.\n')
        msgbox('Data file separation completed.')
    end

    function loadGraph(h_main, event)
        [fileNameCSV, fileLocCSV]=uigetfile({'*.csv'},'Select CSV data file');
        dataGraph=csvread([fileLocCSV fileNameCSV]);
        plot(hAxes1, dataGraph(:,1), dataGraph(:,2));
        set(h_group2_text3, 'String', num2str(size(dataGraph,1)))
        set(h_group2_sub1_text2, 'String', num2str(min(dataGraph(:,1))))
        set(h_group2_sub1_text4, 'String', num2str(max(dataGraph(:,1))))
        set(h_group2_sub1_text6, 'String', num2str(min(dataGraph(:,1))))
        set(h_group2_sub1_text8, 'String', num2str(max(dataGraph(:,1))))
        
        listFoldMat=listFoldMatPreset(size(dataGraph,1)./listFoldMatPreset>=2);
        for ind=1:numel(listFoldMat)
            listFold{ind}=num2str(listFoldMat(ind));
        end

        set(h_group2_sub1_text11, 'String', listFold);
    end

    function estimateNo(h_main, event)
        xStartNew=str2num(get(h_group2_sub1_text6, 'String'));
        xEndNew=str2num(get(h_group2_sub1_text8, 'String'));
        xNo=size(dataGraph,1);
        newData=dataGraph(and(dataGraph(:,1)>=xStartNew,dataGraph(:,1)<=xEndNew),:);
        newData=newData(1:str2num(listFold{get(h_group2_sub1_text11,'value')}):end,:);
        xNoNew=size(newData,1);
        
        warning off
        msgbox(['Estimated No. of new data after the new x range is applied is' ...
            newline newline repmat(' ',1,44+length(num2str(xNoNew))) num2str(xNoNew)], 'Estimated data size')
        warning on
        
        listFoldMat=listFoldMatPreset(xNoNew./listFoldMatPreset>=2);
        for ind=1:numel(listFoldMat)
            listFold{ind}=num2str(listFoldMat(ind));
        end

        set(h_group2_sub1_text11, 'String', listFold);
    end

    function runReduction(h_main,event)
        xStartNew=str2num(get(h_group2_sub1_text6, 'String'));
        xEndNew=str2num(get(h_group2_sub1_text8, 'String'));
        xNo=size(dataGraph,1);
        newData=dataGraph(and(dataGraph(:,1)>=xStartNew,dataGraph(:,1)<=xEndNew),:);
        
        nFold=str2num(listFold{get(h_group2_sub1_text11,'value')});
        newData=newData(1:nFold:end,:);
        xNoNew=size(newData, 1);

        if xNoNew>xNo
            msgbox('No. of the origianl data points is smaller than the final No. of data points.')
            return;

        else
            msgText=['Resulting No. of data is ' num2str(xNoNew) '.'];
            if xNoNew<100
                msgText=[msgText newline '\color{Red}Resulting No. is <100.'];
            end
            msgText=[msgText newline '\color{Black}Do you want to proceed data reduction process?'];
            opts.Interpreter = 'default';
            opts.Default = 'Yes';
            selBox = questdlg(msgText,'Data reduction?',...
                'Yes','No',opts);
            if strcmp(selBox, 'Yes')
                fileNames=ls([fileLocCSV '*.csv']);
                for ind=1:size(fileNames,1)
                    if ind==1
                        warning off
                        indBackslash=regexpi(fileLocCSV, '\');
                        newFolderName=[fileLocCSV(1:indBackslash(end-1)) ...
                            fileLocCSV(indBackslash(end-1)+1:indBackslash(end)-1) ...
                            '_' num2str(round(xStartNew,0)) '_' ...
                            num2str(round(xEndNew,0)) '_N' num2str(xNoNew)];
                        mkdir(newFolderName)
                        warning on
                        fprintf('<strong>--------------------------------------------------------</strong>\n')
                        fprintf('%s is created.\n', newFolderName)
                        fprintf('<strong>--------------------------------------------------------</strong>\n')
                    end
                    csvData=csvread([fileLocCSV deblank(fileNames(ind,:))]);
                    crtManual=get(h_group2_sub1_text17, 'value');
                    valManual=str2num(get(h_group2_sub1_text18, 'String'));
                    if get(h_group2_sub1_text13, 'value')==1
                        xStartAvg=str2num(get(h_group2_sub1_text14,'String'));
                        xEndAvg=str2num(get(h_group2_sub1_text16,'String'));
                        if xStartAvg>=xEndAvg
                            msgbox('x range is not valid.')
                            return;
                        end
                        csvData(:,2)=csvData(:,2)-mean(csvData(and(csvData(:,1)>=xStartAvg, ...
                        	csvData(:,1)<=xEndAvg),2));
                    elseif crtManual==1
                        if valManual~=0
                            if ~isempty(valManual)
                                csvData(:,2)=csvData(:,2)-valManual;
                            end
                        end

                    end
                    csvData=csvData(1:nFold:end,:);
                    csvData(or(csvData(:,1)<xStartNew, csvData(:,1)>xEndNew),:)=[];
                    csvwrite([newFolderName '\' deblank(fileNames(ind,:))], csvData)
                    fprintf('%s was created.\n', deblank(fileNames(ind,:)))

                end
                fprintf('<strong>--------------------------------------------------------</strong>\n')
                fprintf('Data size reduction process was completed.\n')
                msgbox('Data size reduction process was completed.')
            end
            
        end
    end
    
end