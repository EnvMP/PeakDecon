function guiPost()
%% ========================================================
%  GUI for Peak Deconvolution Postprocessor (ver 0.1)
%  packaged with GUI Peak Deconvolution with NMF algorithm
%  (ver. 0,1)
% ========================================================
% function guiPostprocess()  
% ========================================================
% This software was intested to postprocess data processed
% by the peak deconvolution software.
% ========================================================
% == Version history ==
% 2/21/2020: ver. 0.1
% ========================================================
% Minkyu Park

clear all; close all hidden; clc
global data Res indGraph hGraph hAxes hAxesSub

indGraph=0; 

% Get the screen size to put the software panel at the center of the screen
ss = get(0,'screensize');
wHorSize=700;
wVerSize=450;
posCenter=[ss(3)/2-wHorSize/2 ss(4)/2-wVerSize/2];

% Create main figure
h_main=figure;
set(gcf, 'pos', [posCenter wHorSize wVerSize], ...
    'Numbertitle', 'off', 'name', 'PeakDecon Postprocessor ver. 0.1');
set(h_main, 'handlevisibility', 'off');

% UI Button Group 1: Data File Separation
h_group1=uibuttongroup(h_main, 'Title', 'Load mat file', ...
    'Fontsize', 11, 'pos', [0.03 0.38 0.94 0.6], ...
    'ForegroundColor', 'b');

h_group1_text1=uicontrol(h_group1, 'Style', 'text', 'String', 'Select a mat file', ...
    'Pos', [30 230 320 25], 'Fontsize', 10);
set(h_group1_text1, 'Units', 'Normalized')
h_group1_text2=uicontrol(h_group1, 'Style', 'text', 'String', 'File:', ...
    'Pos', [5 208 40 25], 'Fontsize', 10);
set(h_group1_text2, 'Units', 'Normalized')
h_group1_text3=uicontrol(h_group1, 'Style', 'edit', 'String', '', ...
    'Pos', [45 211 270 23], 'Fontsize', 10);
set(h_group1_text3, 'Units', 'Normalized')
h_group1_but1=uicontrol(h_group1, 'Style', 'pushbutton', 'String', 'Load', 'Fontsize', 10, ...
    'pos', [320 210 50 25]);
h_group1_but1.Callback=@loadMatFile;
set(h_group1_but1, 'Units', 'Normalized')
h_group1_text4=uicontrol(h_group1, 'Style', 'text', 'String', 'Select N th component:', ...
    'Pos', [15 178 150 25], 'Fontsize', 10);
set(h_group1_text4, 'Units', 'Normalized')
h_group1_text5=uicontrol(h_group1, 'Style', 'popupmenu', 'String', ' ', ...
    'Pos', [160 180 45 25], 'Fontsize', 10);
set(h_group1_text5, 'Units', 'Normalized')
h_group1_text5.Callback=@plotW;
h_group1_text6=uicontrol(h_group1, 'Style', 'text', 'String', 'Cosine similarity index', ...
    'Pos', [30 150 150 25], 'Fontsize', 10, 'FontWeight', 'bold');
set(h_group1_text6, 'Units', 'Normalized')
h_group1_text7=uicontrol(h_group1, 'Style', 'text', 'String', 'Component', ...
    'Pos', [235 150 150 25], 'Fontsize', 10, 'FontWeight', 'bold');
set(h_group1_text7, 'Units', 'Normalized')
h_group1_text8=uicontrol(h_group1, 'Style', 'pushbutton', 'String', 'Plot Split half', ...
    'Pos', [205 180 105 25], 'Fontsize', 10);
h_group1_text8.Callback=@plotSplit;
set(h_group1_text8, 'Units', 'Normalized')
h_group1_text9=uicontrol(h_group1, 'Style', 'pushbutton', 'String', 'plot in separate window', ...
    'Pos', [20 5 150 25], 'Fontsize', 10);
h_group1_text9.Callback=@plotCosSim;
h_group1_text10=uicontrol(h_group1, 'Style', 'pushbutton', 'String', 'plot in separate window', ...
    'Pos', [200 5 150 25], 'Fontsize', 10);
h_group1_text10.Callback=@plotComp;
h_group1_text11=uicontrol(h_group1, 'Style', 'pushbutton', 'String', 'Export', ...
    'Pos', [320 177 50 25], 'Fontsize', 10);
set(h_group1_text11, 'Units', 'Normalized')
h_group1_text11.Callback=@export;

h_group1_sub1=uibuttongroup(h_group1, 'Title', 'Label setup', ...
    'Fontsize', 11, 'pos', [0.6 0.03 0.47 0.98]);
h_group1_sub1_text1=uicontrol(h_group1_sub1, 'Style', 'text', 'String', 'x axis label for signals:', ...
    'Pos', [5 195 150 25], 'Fontsize', 10);
h_group1_sub1_text2=uicontrol(h_group1_sub1, 'Style', 'edit', 'String', 'x', ...
    'Pos', [155 197 80 25], 'Fontsize', 10);
h_group1_sub1_text3=uicontrol(h_group1_sub1, 'Style', 'text', 'String', 'Setup for display order and text of groups:', ...
    'Pos', [5 155 250 25], 'Fontsize', 10);
h_table=uitable(h_group1_sub1, ...
    'pos', [5 5 250 150], ...
    'ColumnWidth', {90,90,40}, ...
    'ColumnName', {'Original label','Final label','Order'}, ...
    'ColumnFormat', {'Char', 'Char', 'Numeric'}, ...
    'ColumnEditable', [false true true]);

hAxes1=axes(h_group1, 'pos', [0.05 0.25 0.2 0.35]);
hAxes2=axes(h_group1, 'pos', [0.3 0.25 0.27 0.35]);

h_group2=uibuttongroup(h_main, 'Title', 'Plot graphs', ...
    'Fontsize', 11, 'pos', [0.03 0.03 0.94 0.33], ...
    'ForegroundColor', 'b');
h_group1_but2=uicontrol(h_group2, 'Style', 'pushbutton', 'String', 'Correlation plot for components', 'Fontsize', 10, ...
    'pos', [5 100 200 25]);
h_group1_but2.Callback=@plotCorrPlot;
h_group1_but2_text1=uicontrol(h_group2, 'Style', 'text', 'String', 'Type:', 'Fontsize', 10, ...
    'pos', [210 98 50 25]);
h_group1_but2_text2=uicontrol(h_group2, 'Style', 'popupmenu', 'String', {'Pearson', 'Kendall', 'Spearman'}, 'Fontsize', 10, ...
    'pos', [260 100 80 25]);
h_group1_but2_text3=uicontrol(h_group2, 'Style', 'text', 'String', 'Tail:', 'Fontsize', 10, ...
    'pos', [340 98 50 25]);
h_group1_but2_text4=uicontrol(h_group2, 'Style', 'popupmenu', 'String', {'both', 'right', 'left'}, 'Fontsize', 10, ...
    'pos', [390 100 50 25]);
h_group1_but2_text5=uicontrol(h_group2, 'Style', 'text', 'String', 'Significant level:', 'Fontsize', 10, ...
    'pos', [460 98 100 25]);
h_group1_but2_text6=uicontrol(h_group2, 'Style', 'edit', 'String', '0.05', 'Fontsize', 10, ...
    'pos', [560 100 50 25]);

h_group1_but3=uicontrol(h_group2, 'Style', 'pushbutton', 'String', 'Matrix plot by group', 'Fontsize', 10, ...
    'pos', [5 75 200 25]);
h_group1_but3.Callback=@plotMatrix;
h_group1_but3_text1=uicontrol(h_group2, 'Style', 'text', 'String', 'Display options for diagonal plots :', 'Fontsize', 10, ...
    'pos', [213 73 210 25]);
h_group1_but3_text2=uicontrol(h_group2, 'Style', 'popupmenu', 'String', {'grpbars', 'stairs', 'hist', 'none', 'variable'}, 'Fontsize', 10, ...
    'pos', [460 75 80 25]);

h_group1_but4=uicontrol(h_group2, 'Style', 'pushbutton', 'String', 'Box plot by group', 'Fontsize', 10, ...
    'pos', [5 50 200 25]);
h_group1_but4.Callback=@plotBox;
h_group1_but4_text1=uicontrol(h_group2, 'Style', 'text', 'String', 'Style:', 'Fontsize', 10, ...
    'pos', [210 48 50 25]);
h_group1_but4_text2=uicontrol(h_group2, 'Style', 'popupmenu', 'String', {'traditional', 'compact'}, 'Fontsize', 10, ...
    'pos', [260 50 80 25]);
h_group1_but4_text3=uicontrol(h_group2, 'Style', 'text', 'String', 'Color:', 'Fontsize', 10, ...
    'pos', [340 48 50 25]);
h_group1_but4_text4=uicontrol(h_group2, 'Style', 'popupmenu', 'String', {'No', 'Yes'}, 'Fontsize', 10, ...
    'pos', [390 50 50 25]);
h_group1_but4_text5=uicontrol(h_group2, 'Style', 'text', 'String', 'x Label:', 'Fontsize', 10, ...
    'pos', [460 48 50 25]);
h_group1_but4_text6=uicontrol(h_group2, 'Style', 'popupmenu', 'String', {'inline', 'horizontal'}, 'Fontsize', 10, ...
    'pos', [520 50 90 25]);


%% Functions
    function loadMatFile(h_main, event)
        [fileName, fileLoc]=uigetfile({'*.mat'},'Select mat file containing NMF results');
        set(h_group1_text3, 'String', [fileLoc fileName])
        load([fileLoc fileName], 'data', 'Res')
        listComp=num2cell(Res.components);
        set(h_group1_text5, 'String', listComp);
        if isempty(Res.splitInd{1}) || ~isnan(Res.splitInd{1})
            bar(hAxes1, Res.cosSimilarity(Res.components));
            set(hAxes1, 'xticklabel', Res.components);
            xlabel(hAxes1, 'No. of components')
            ylabel(hAxes1, 'Cosine similarity index')
        end
        for i=1:size(data.groupLabel,1)
            dataTable{i,1}=deblank(data.groupLabel(i,:));
            dataTable{i,2}=deblank(data.groupLabel(i,:));
            dataTable{i,3}=i;
        end
        h_table.Data=dataTable;
        plotW;
    end
    
    function plotW(h_main, event)
        if ~strcmp(get(h_group1_text5, 'string'),' ')
            nComp=Res.components(get(h_group1_text5, 'value'));
            plot(hAxes2, data.x, Res.W{nComp});
            xlabel(hAxes2, getXlabel)
            for i=1:nComp
                compLabel{i}=['C' num2str(i)];
            end
            legend(hAxes2,compLabel, 'location', 'eastoutside')
        end
    end

    function plotCosSim(h_main, event)
        if isempty(Res.splitInd{1}) || ~isnan(Res.splitInd{1})
            plotBase;
            listComp=num2cell(Res.components);
            set(h_group1_text5, 'String', listComp);
            bar(Res.cosSimilarity(Res.components));
            set(gca, 'xticklabel', Res.components);
            xlabel('No. of components')
            ylabel('Cosine similarity index')
            for i=1:length(Res.components)
                fprintf('C%i, Cosine similarity index= %1.3f\n', Res.components(i),Res.cosSimilarity(Res.components(i)))
            end
            fprintf('\n')
        else
            msgbox('Split half analysis was not conducted for the loaded data')
        end
    end

    function plotComp(h_main, event)
        plotBase;
        nComp=Res.components(get(h_group1_text5, 'value'));
        plot(data.x, Res.W{nComp});
        xlabel(getXlabel)
        for i=1:nComp
            compLabel{i}=['C' num2str(i)];
        end
        legend(gca, compLabel, 'location', 'northeast')
        
        % Mw, Mn, and dispersity
        for i=1:nComp

            Mw(i)=sum(Res.W{nComp}(:,i).*data.x.^2)./sum(Res.W{nComp}(:,i).*data.x);
            Mn(i)=sum(Res.W{nComp}(:,i).*data.x)./sum(Res.W{nComp}(:,i));

            polyD(i)=Mw(i)./Mn(i);
            printTxt{i}=sprintf('C%i: Mw= %1.5f, Mn= %1.5f, Dispersity= %1.2f\n', i, Mw(i), Mn(i), polyD(i));
            fprintf('C%i: Mw= %1.5f, Mn= %1.5f, Dispersity= %1.2f\n', i, Mw(i), Mn(i), polyD(i));
        end
        fprintf('Average dispersity = %1.2f\n', mean(polyD))
        fprintf('Keep in mind that the Mw and Mn values were based on your x-axis values.\n')
        fprintf('If time is used as x-axis, the obtained Mn and Mw should be further convered to your desired unit (e.g., Da).\n\n')
        msgbox(printTxt)
    end

    function plotSplit(h_main, event)
        if ~strcmp(get(h_group1_text5, 'string'),' ')
            if isempty(Res.splitInd{1}) || ~isnan(Res.splitInd{1})
                plotBase;
                nComp=Res.components(get(h_group1_text5, 'value'));

                nThreshold=4;
                [subplotRow,subplotCol]=calcDim(nComp, nThreshold);

                for i=1:nComp
                    compLabel{i}=['C' num2str(i)];
                end

                for i=1:nComp
                    subplot(subplotRow,subplotCol,i)
                    plot(data.x, Res.W{nComp}(:,i))
                    hold on
                    plot(data.x, Res.WSplit{nComp}{1}(:,i))
                    plot(data.x, Res.WSplit{nComp}{2}(:,i))
                    xlabel(getXlabel)
                    if i==subplotCol
                        legend({'All', 'Split1', 'Split2'})
                    end

                end
                hAxesSub = findobj(hGraph{indGraph}, 'Type', 'Axes');
                linkaxes(hAxesSub, 'x')
                textDataX='';
                textDataY='Score';
                locTitle(nComp, subplotRow,subplotCol, compLabel, textDataX, textDataY)
            else
                msgbox('Split half analysis was not conducted for the loaded data')
            end
        else
            msgbox('Data was not loaded.')
        end
    end

    function plotCorrPlot(h_main, event)
        tableData=get(h_table, 'Data');
        [orderInd, returnVal]=tableCheck(tableData);
        if returnVal==1
            msgbox(['Order in the Group name setup table is not correct.' newline ...
                'Order should start from 1 and ends with N'])
            return;
        end
        
        if ~strcmp(get(h_group1_text5, 'string'),' ')
            plotBase;
            nComp=Res.components(get(h_group1_text5, 'value'));
            for i=1:nComp
                compLabel{i}=['C' num2str(i)];
            end
            typeCell=get(h_group1_but2_text2, 'string');
            typeVal=typeCell(get(h_group1_but2_text2, 'value'));
            tailCell=get(h_group1_but2_text4, 'string');
            tailVal=tailCell(get(h_group1_but2_text4, 'value'));
            alphaVal=str2num(get(h_group1_but2_text6, 'string'));
            corrplot(Res.H{Res.components(get(h_group1_text5, 'value'))}', ...
                'varNames', compLabel, ...
                'type', typeVal{1}, ...
                'tail', tailVal{1}, ...
                'alpha', alphaVal, ...
                'testR','on')
        end
        
    end

    function plotMatrix(h_main, event)
        if numel(data.groupLabel)==1
            msgbox('Data do not have groups.')
            return;
        else
            tableData=get(h_table, 'Data');
            [orderInd, returnVal]=tableCheck(tableData);
            if returnVal==2
                msgbox('Order in the Group name setup table is not correct.')
                return;
            end
            newGroupLabel=tableData(orderInd, 2);
            newGroupInd=changeGroupInd(data.groupInd, orderInd);
            
            plotBase;
            nComp=Res.components(get(h_group1_text5, 'value'));
            for i=1:nComp
                compLabel{i}=['C' num2str(i)];
            end

            dispoptCell=get(h_group1_but3_text2, 'string');
            dispoptVal=dispoptCell(get(h_group1_but3_text2, 'value'));
            gplotmatrix(Res.H{nComp}',[],newGroupInd,lines,[],[],[],...
                dispoptVal{1}, compLabel, compLabel)
            hLegend = findobj('Tag','legend');
            set(hLegend, 'String', newGroupLabel)
            
        end
    end

    function plotBox(h_main, event)
        if numel(data.groupLabel)==1
            msgbox('Data do not have groups.')
            return;
        else
            tableData=get(h_table, 'Data');
            [orderInd, returnVal]=tableCheck(tableData);
            if returnVal==1
                msgbox('Order in the Group name setup table is not correct.')
                return;
            end
            newGroupLabel=tableData(orderInd, 2);
            newGroupInd=changeGroupInd(data.groupInd, orderInd);
            
            plotBase;
            nComp=Res.components(get(h_group1_text5, 'value'));
            for i=1:nComp
                compLabel{i}=['C' num2str(i)];
            end
            
            % Calculatin of m x n for subplot
            nThreshold=4;
            [subplotRow,subplotCol]=calcDim(nComp, nThreshold);

            % Boxplot
            if get(h_group1_but4_text4, 'value')==1
                colorVal=[];
            else
                tempLines=lines;
                colorVal=tempLines(1:length(orderInd), :);
            end
                labelOriCell=get(h_group1_but4_text6, 'string');
                labelOriVal=labelOriCell{get(h_group1_but4_text6, 'value')};
                plotStyleCell=get(h_group1_but4_text2, 'string');
                plotStyleVal=plotStyleCell{get(h_group1_but4_text2, 'value')};
            for i=1:nComp
                subplot(subplotRow,subplotCol,i);
                boxplot(Res.H{nComp}(i,:),newGroupInd ,...
                    'colorgroup', colorVal, ...
                    'PlotStyle', plotStyleVal, ...
                    'LabelOrientation', labelOriVal, ...
                    'labels', newGroupLabel);
            end
            hAxesSub = findobj(hGraph{indGraph}, 'Type', 'Axes');
            linkaxes(hAxesSub, 'x')
            
            textDataX='Group';
            textDataY='Score';
            locTitle(nComp, subplotRow,subplotCol, compLabel, textDataX, textDataY)
            
        end
    end
   
    function [subplotRow,subplotCol]=calcDim(nComp, nThreshold)
        fac=factor(nComp);
        if nComp<4
            subplotRow=1;
            subplotCol=nComp;
        else
            halfInd=length(fac)/2;
            if mod(length(fac),2)==0
                subplotRow=prod(fac(1:halfInd));
                subplotCol=prod(fac(halfInd+1:end));
            else
                    subplotRow=min(prod(fac(1:ceil(halfInd))), ...
                        prod(fac(ceil(halfInd)+1:end)));
                    subplotCol=max(prod(fac(1:ceil(halfInd))), ...
                        prod(fac(ceil(halfInd)+1:end)));
            end
            facPlus=factor(nComp);

            while or(subplotRow>nThreshold, subplotCol>nThreshold)
                facPlus=factor(prod(facPlus)+1);
                halfInd=length(facPlus)/2;
                if mod(length(facPlus),2)==0
                    subplotRow=min(prod(facPlus(1:halfInd)),prod(facPlus(halfInd+1:end)));
                    subplotCol=max(prod(facPlus(1:halfInd)),prod(facPlus(halfInd+1:end)));
                else
                    subplotRow=min(prod(facPlus(1:ceil(halfInd))),prod(facPlus(ceil(halfInd)+1:end)));
                    subplotCol=max(prod(facPlus(1:ceil(halfInd))),prod(facPlus(ceil(halfInd)+1:end)));
                end
            end
        end
    end

    function [orderInd, returnVal]=tableCheck(tableData)
        [~,orderInd]=sort(cell2mat(tableData(:,3)));
        if sum(sort(orderInd)==[1:length(orderInd)]')~=length(orderInd)
            returnVal=1;
        else
            returnVal=0;
        end
    end

    function newGroupInd=changeGroupInd(groupInd, orderInd)
        tableData=get(h_table, 'Data');
        
        [~,newOrderInd]=sort(orderInd);
        refInd=1:length(newOrderInd);
        for i=1:length(newOrderInd)
            ind{i}=groupInd==refInd(newOrderInd==i);
        end
        for i=1:length(newOrderInd)
            newGroupInd(ind{i})=i;
        end
        
    end

    function locTitle(nComp, subplotRow,subplotCol, compLabel, textDataX, textDataY)
            
        indSubplot=ceil(subplotRow/2)*subplotCol-(subplotRow*subplotCol-nComp);
        posFig=get(hGraph{indGraph}, 'pos');
        posAxis=get(hAxesSub(indSubplot), 'pos');
        ntextDataX=length(textDataX);
        ntextDataY=length(textDataY);
        nTextToPx=4.8;          % Empirical value to convert the number of
                                % texts to pixel points

        hTextX=text(hAxesSub(1), 0,0, textDataX);
        set(hTextX,'Units', 'Normalized');



        posAxisX1=get(hAxesSub(1), 'pos');
        posAxisX2=get(hAxesSub(subplotCol), 'pos');
        wXAxis=posAxisX1(3);
        hYAxis=posAxisX1(4);
        figPos=get(hGraph{indGraph}, 'pos');
        set(hTextX,'Pos', [(0.5-posAxisX1(1)-ntextDataX*nTextToPx/figPos(3)/2)/wXAxis ...
            -posAxisX1(2)*0.7/hYAxis 0])

        hTextY=text(hAxesSub(subplotCol), 0,0, textDataY, 'rotation', 90);
        set(hTextY,'Units', 'Normalized');
        set(hTextY,'Pos', [-posAxisX2(1)*0.8/wXAxis ...
            (0.5-posAxisX2(2)-ntextDataY*nTextToPx/figPos(4)/2)/hYAxis 0])
        

        % Add components at the top left corner of subplot as texts
        % rather than using title function
        for i=1:nComp
            text(hAxesSub(nComp-i+1), 0.05, (1-0.15/hYAxis*wXAxis), compLabel{i}, 'Units', 'Normalized', 'fontweight', 'bold')
        end
    end

    function plotBase(h_main, event)
        figure; 
        indGraph=indGraph+1;
        hGraph{indGraph}=gcf;
        indGraph=indGraph+1;
        hGraph{indGraph}=gcf;
        hAxes{indGraph}=gca;
    end

    function labelVal=getXlabel()
        labelVal=get(h_group1_sub1_text2, 'string');
    end

    function export(h_main, event)
        if ~strcmp(get(h_group1_text5, 'string'),' ')
            clear HCell WCell
            nComp=Res.components(get(h_group1_text5, 'value'));
            for i=1:size(data.Filenames,1)
                if i==1
                    HCell{1,1}='Sample';
                    WCell{1,1}='x axis';
                    areaCell{1,1}='Sample';
                    for j=1:nComp
                        HCell{1,j+1}=['Component ' num2str(j)];
                        WCell{1,j+1}=['Component ' num2str(j)];
                        areaCell{1,j+1}=['Component ' num2str(j)];
                        areaVal(1,j)=trapz(data.x, Res.W{nComp}(:,j));
                    end
                end
                HCell{i+1,1}=deblank(data.Filenames(i,:));
                HCell{i+1,1}(end-3:end)=[];
                areaCell{i+1,1}=HCell{i+1,1};
            end
            HCell(2:end,2:nComp+1)= num2cell(Res.H{nComp}');
            WCell=[WCell; num2cell([data.x Res.W{nComp}])];
            areaCell(2:end,2:nComp+1)=num2cell(Res.H{nComp}' ...
                .* repmat(areaVal,size(Res.H{nComp}',1),1));
            fileLoc=get(h_group1_text3, 'String');
            tempLoc=regexp(fileLoc, '\');
            fileLoc=fileLoc(1:tempLoc(end));
            xlswrite([fileLoc 'Comp' num2str(nComp) '_data.xlsx'], WCell, 'W')
            xlswrite([fileLoc 'Comp' num2str(nComp) '_data.xlsx'], HCell, 'H')
            xlswrite([fileLoc 'Comp' num2str(nComp) '_data.xlsx'], areaCell, 'Area')
            msgbox(['Comp' num2str(nComp) '_data.xlsx was saved.'])
        else
            msgbox('Data was not loaded.')
        end
    end

end