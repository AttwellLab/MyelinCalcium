% AK 29/07/2016
% Plot calcium imaging data from GECIquant results for ROIs in sheaths,
% processes and somas
% synchronicity graphs


%% also need: 

%% Intersections
% Version: 1.12, 27 January 2010
% Author:  Douglas M. Schwarz

%% function [ycorr,yfit] = bf(y,varargin)

%% Functions by AK
% read_txt.m
% read_num.m
% remove_spikes.m

%%
% Example file: GC1 Results, tab one with data from GeciQuant, tab two with
% false positives identified for removal
% folder test_test contains results for when you run the code on the
% example file
%%

clear all;
clc;
close all;

%%
date = '_test'; % change this date to match the analysed experiment

%% define variables

yfitAVERAGE = 60; % the number of points of the trace that get average for baseline
noiseFactor = 2.25; % multipler for noise defined as std(tracesample(1:100))
minPEAKarea = 0.275; % min area for a detected peak to be considered real
protectPEAKS = 0.75; % values above this are excluded from the trace sample in order not to increase noise of traces with large peaks
maxPEAKnumber_estimate = 15; % if any of the sheaths has more than 10 peaks per timepoint change this number


%%
PEAKcount_allTimePoints = [];
mean_Durations_allTimePoints = [];
mean_Amplitude_allTimePoints = [];
mean_Area_allTimePoints = [];

% list xls files in the root path
rootpath    = '../result files/';
files       = dir([rootpath '*.xls']);

mainFOLDER = strcat('Test', date);
folderpath1 = [rootpath, mainFOLDER];
folderCHECK1 = exist(folderpath1,'dir')==~7; %check if the folder exists already, if not create a new one
if folderCHECK1 == 1;
   mkdir(rootpath, mainFOLDER);
end

rootpath_mainFOLDER = strcat(rootpath, mainFOLDER,'/');
folderpath2 = [rootpath_mainFOLDER, 'Results'];
folderCHECK2 = exist(folderpath2,'dir')==~7; %check if the folder exists already, if not create a new one
if folderCHECK2 == 1;
   mkdir(rootpath_mainFOLDER, 'Results');
end


% loop over the files
for fls = 1:1:length(files)
    
    % load the xls file
    [num,txt,raw]   = xlsread([rootpath files(fls).name]);
        
    [A,B] = xlsfinfo([rootpath files(fls).name]);
    sheetValid = any(strcmp(B, 'Remove'));
    
    if sheetValid==1
    [num2,txt2,raw2]   = xlsread([rootpath files(fls).name], 'Remove');
    end
    
    filename = files(fls).name(1:end-4);
    folderpath3 = [rootpath_mainFOLDER, filename];
    
    folderCHECK3 = exist(folderpath3,'dir')==~7; %check if the folder exists already, if not create a new one
    if folderCHECK3 == 1;
       mkdir(rootpath_mainFOLDER, filename);
    end
        
    intensity = num(1:end,5);
    z = 180;
    ROInumber = length(num)/z;
    n = zeros(1,ROInumber);
    dFF = zeros(1,z);
  
     
    maxPEAKS_allROI = zeros(ROInumber, maxPEAKnumber_estimate);
    PEAKdurations_allROI = zeros(ROInumber, maxPEAKnumber_estimate);
    TRANSIENTdurations_perCELL = zeros(1, maxPEAKnumber_estimate);
    realPEAKSarea = zeros (1,maxPEAKnumber_estimate);
    realPEAKSarea_allROI = zeros(ROInumber, maxPEAKnumber_estimate);
    timeOFspikes_allROI = zeros(1,16);
          
    clearvars threshold_allROI TRANSIENTamplitude_perCELL allDATA.yfit2_table
    
                 
    
        for i = 1:ROInumber;
                                   
            %% extract data from matrix
            n(i) = i-1;
            ROIintensity= intensity(n(i)*z+1:(n(i)+1)*z); % grab all intensity datapoints for a given ROI
            ROIlabel = txt(n(i)*z+3,4);
            AllROIlabel(i)= ROIlabel;
            
                        
            if (isnan(ROIintensity(1,1)) == 1)==1; % if the first value in ROIintensity is nan display as below
                clearvars peakSTART peakEND AREA realPEAKSonly maxPEAKS timeOFmaxPEAKs PEAKdurations realPEAKSarea
                display(strcat(ROIlabel, ' sheath not present in this timepoint'));
                PEAKcount = 0;
                maxPEAKS = 0;
                PEAKdurations = 0;
                timeOFmaxPEAKs = 0;
                realPEAKSarea = 0;
                timeOFspikes = zeros(1,15);
                
            else 
                baselineSize = 10;
                baseline = mean(ROIintensity(1:baselineSize)); % define which entries are the baseline for delta F/F
                dFF = (ROIintensity - baseline) / baseline; % show delta F/F
                time = (1:z); % define x axis
                [ycorr,yfit] = bf(dFF,[0,30,60,90,120,150,180],yfitAVERAGE,'pchip');    % fit baseline; pchip piecewise cubic Hermite interpolation method

                % create graphs for myelin sheaths only, i.e. ROI label must contain letter M to proceed

                c = ROIlabel(1);
                lettertofind = 'M';

                %%% SOMA section 
                
                if cellfun(@(s) any(s == lettertofind), c) == 0
                   % clear peak related variables
                   clearvars peakSTART peakEND AREA realPEAKSonly maxPEAKS timeOFmaxPEAKs PEAKdurations timeOFmaxPEAKs
                   PEAKcount = 0;
                   maxPEAKS = 0;
                   PEAKdurations = 0;
                   timeOFmaxPEAKs = 0;
                   realPEAKSarea = 0;
                   timeOFspikes = zeros(1,15);
                   
                   clearvars transientSTART transientEND AREA realTRANSIENTSonly maxTRANSIENTS timeOFmaxTRANSIENTS TRANSIENTdurations
                   display(strcat(ROIlabel, ' This is soma'));
                   TRANSIENTScount = 0;
                   maxTRANSIENTS = 0;
                   TRANSIENTdurations = 0;
                   timeOFmaxTRANSIENTS = 0;
                   realTRANSIENTSarea = 0;
                   timeOFtransients = zeros(1,15);
                   timeOF_SOMA_transients = timeOFtransients;
                   
                % soma baseline
                
                SOMAbaselineSize = 10;
                SOMAyfitAVERAGE = 120; % % the number of points of the trace that get average for soma baseline
                SOMAbaseline = mean(ROIintensity(1:SOMAbaselineSize)); % define which entries are the baseline for delta F/F
                dFF = (ROIintensity - SOMAbaseline) / SOMAbaseline; % show delta F/F
                time = (1:z); % define x axis
                [ycorr,yfit] = bf(dFF,[0,30,60,90,120,150,180],SOMAyfitAVERAGE,'pchip');    % fit baseline; pchip piecewise cubic Hermite interpolation method
                [ycorr2,yfit2] = bf(ycorr,0:2:180,2,'pchip'); 
                
                protectSOMA_PEAKS = 0.20;
                tracesample = ycorr(ycorr<protectSOMA_PEAKS);
                noise = std(tracesample(1:100)); %noise is sefined as the standard deviation of the first 100 values of baseline adjusted trace that are <0.5)
                SOMAnoiseFactor = 1.00; % multipler for noise defined as std(tracesample(1:100))
                MinPeakHeight = noise*SOMAnoiseFactor;
                MinPeakHeight = round(MinPeakHeight,2);
                 
                % find intersection between (0.001)and the peak curve, i.e. beginning and end of each peak
                % create a line that marks minimal peak height
                 
                y = linspace(MinPeakHeight,MinPeakHeight,z); % a line y=MinPeakHeight
                SOMAtransients = yfit2' - (y);
                SOMAtransients(SOMAtransients<0)=0;

                x1 = [0 z];
                y1 = [0.001 0.001]; %vector in which betweenMinMax value is repeated for the length of x axis (distance)
                [X0,Y0] = intersections(time,SOMAtransients,x1,y1);
                
                if (~isempty(X0) == 1)==0; % if treshold is not reached produce 2/3 graphs
                     
                        clearvars transientSTART transientEND AREA realTRANSIENTSonly maxTRANSIENTS timeOFmaxTRANSIENTS TRANSIENTdurations 
                        display([ROIlabel,'no transients above threshold in soma trace']);
                        TRANSIENTScount = 0;
                        maxTRANSIENTS = 0;
                        TRANSIENTdurations = 0;
                        timeOFmaxTRANSIENTS = 0;
                        realTRANSIENTSarea = 0;
                        timeOFtransients = zeros(1,15);
                        timeOF_SOMA_transients = timeOFtransients;
               
                        % plot figure 
                        figure(i);
                        subplot(3,1,1)       % add first plot in 2 x 1 grid
                        plot(time,dFF,'k', time, yfit, 'b');
                        grid on

                        xlabel('time'); ylabel('dF/F')
                        ROIlabel = txt(n(i)*z+3,4);
                        xmin = 0;
                        xmax = z;
                        set(gca,'XTick',xmin:10:xmax);
                        xlim([xmin xmax]);
                        ymin = -1;
                        title(ROIlabel); legend('dF/F','baseline fit', 'Location','southeast');

                        subplot(3,1,2)       % add second plot in 2 x 1 grid
                        plot(time,yfit2,'b',time, y, 'r');
                        grid on

                        xlabel('time'); ylabel('dF/F with baseline correction')
                        xmin = 0;
                        xmax = z;
                        set(gca,'XTick',xmin:10:xmax);
                        xlim([xmin xmax]);
                        ymin = -1;
                        noiseLABEL = strcat('noise< ', num2str(MinPeakHeight));
                        legend('ycorr',num2str(noiseLABEL), 'Location','southeast');

                        FigHandle = figure(i);
                        set(FigHandle, 'Position', [-1200, 200, 1080, 810]); % use [100, 100, resolution, resolution] to display fig in the middle
                        
                        % append yfit2 data to the matrix
                        
                        NEWyfit2_table = horzcat(allDATA.yfit2_table,yfit2);
                        allDATA.yfit2_table = NEWyfit2_table;
                        
           elseif (~isempty(X0) == 1)==1;  % if treshold is reached produce 3/3 graphs

                        % Check if soma transients are real (min area condition)
                        
                        minTRANSIENTarea = 0.24; % min area for a supratreshold soma transient to be considered real

                        clearvars TRANSIENTScount transientSTART transientEND AREA realTRANSIENTSonly realTRANSIENTstart realTRANSIENTend realTRANSIENTSarea

                        numberOFtransients = length(X0)/2;


                        for k = 1:numberOFtransients;
                        transientSTART(k) = round(X0(2*k-1));
                        transientEND(k) = round(X0(2*k));
                        AREA(k)=trapz(SOMAtransients(transientSTART(k):transientEND(k)));

                        end

                        for f = 1:numberOFtransients;

                                if AREA(f) > minTRANSIENTarea         
                                realTRANSIENTSonly(f) = 1;
                                realTRANSIENTSarea(f) = AREA(f);

                                else 
                                realTRANSIENTSonly(f) = 0;
                                realTRANSIENTSarea(f) = 0;
                                end

                        end

                        %% find beginnings and ends of real transients in soma trace

                        allTRANSIENTScount = length(realTRANSIENTSonly);
                        
                        for g = 1:allTRANSIENTScount;

                                if realTRANSIENTSonly(g) == 1          
                                realTRANSIENTstart(g) =  transientSTART(g);

                                elseif realTRANSIENTSonly(g) == 0          
                                realTRANSIENTstart(g) =  NaN;

                                end

                        end

                        for g = 1:length(realTRANSIENTSonly);

                                if realTRANSIENTSonly(g) == 1          
                                realTRANSIENTend(g) =  transientEND(g);

                                elseif realTRANSIENTSonly(g) == 0          
                                realTRANSIENTend(g) =  NaN;

                                end

                        end

                        realTRANSIENTstart(find(isnan(realTRANSIENTstart))) = []; % removes NaNs from an array
                        realTRANSIENTend(find(isnan(realTRANSIENTend))) = [];

                        TRANSIENTScount = length(realTRANSIENTend);

                        % check if real transients are present in soma trace
                        if (~isempty(realTRANSIENTend) == 1)==0;
                        display([ROIlabel,'no transients with min area in this soma trace']);
                        display(AREA');

                        clearvars maxTRANSIENTS timeOFmaxTRANSIENTS TRANSIENTdurations
                        TRANSIENTScount = 0;
                        maxTRANSIENTS = 0;
                        TRANSIENTdurations = 0;
                        timeOFmaxTRANSIENTS = 0;
                        timeOFtransients = zeros(1,15);
                        timeOF_SOMA_transients = timeOFtransients;

                        elseif (~isempty(realTRANSIENTend) == 1)==1;
                        display([ROIlabel,'real transients detected in this soma trace']);
                        display(AREA');    

                        clearvars maxTRANSIENTS timeOFmaxTRANSIENTS TRANSIENTdurations
                        
                        timeOFtransients = zeros(10,15);
                        
                            for d = 1:length(realTRANSIENTstart);

                                duration = realTRANSIENTstart(d):1:realTRANSIENTend(d);
                                maxTRANSIENT = max(SOMAtransients(duration));
                                timeOFmaxTRANSIENT = find(SOMAtransients==maxTRANSIENT);
                                maxTRANSIENTS(d) = maxTRANSIENT;
                                timeOFmaxTRANSIENTS(d) = timeOFmaxTRANSIENT;                           
                                TRANSIENTdurations(d) = length(duration);
                                timeOFtransients(d,1:length(duration)) = duration;

                            end
                            
                            timeOFtransients(~any(timeOFtransients,2),: ) = []; %keeps non zero rows
                            timeOF_SOMA_transients = timeOFtransients;
                            
                        end

                        %% create a plot

                        figure(i);                    
                        subplot(3,1,1)       % add first plot in 3 x 1 grid
                        plot(time,dFF,'k', time, yfit, 'b');
                        grid on

                        xlabel('time'); ylabel('dF/F')
                        ROIlabel = txt(n(i)*z+3,4);
                        xmin = 0;
                        xmax = z;
                        set(gca,'XTick',xmin:10:xmax);
                        xlim([xmin xmax]);
                        ymin = -1;
                        title(ROIlabel); legend('dF/F','baseline fit', 'Location','southeast');

                        subplot(3,1,2)       % add second plot in 3 x 1 grid
                        plot(time,yfit2,'b',time, y, 'r');
                        grid on

                        xlabel('time'); ylabel('dF/F with baseline correction')
                        xmin = 0;
                        xmax = z;
                        set(gca,'XTick',xmin:10:xmax);
                        xlim([xmin xmax]);
                        ymin = -1;
                        noiseLABEL = strcat('noise< ', num2str(MinPeakHeight));
                        legend('ycorr',num2str(noiseLABEL),'Location','southeast');

                        subplot(3,1,3)       % add third plot in 2 x 1 grid
                        plot(time,SOMAtransients,'k'); grid on
                        hold on

                        if (~isempty(realTRANSIENTend) == 1)==1;
                            plot(timeOFmaxTRANSIENTS, maxTRANSIENTS,'r^', 'MarkerFaceColor','r')

                        end

                        xlabel('time'); ylabel('dF/F with baseline correction')
                        xmin = 0;
                        xmax = z;
                        set(gca,'XTick',xmin:10:xmax);
                        xlim([xmin xmax]);
                        ymin = -1;

                        if (~isempty(realTRANSIENTend) == 1)==1; % display real transients legends only if real transients exist 
                            TRANSIENTcountLEGEND = strcat('Real transients: ', num2str(TRANSIENTScount));
                            legend('Transients only',TRANSIENTcountLEGEND, 'Location','southeast');
                        else 
                            legend('Transients only', 'Location','southeast');
                        end
                        hold off

                        FigHandle = figure(i);
                        set(FigHandle, 'Position', [-1200, 200, 1080, 810]);
                        
                        % save transients to a matrix
                        if TRANSIENTScount==0
                        TRANSIENTdurations_perCELL(1,1) = TRANSIENTdurations;
                        TRANSIENTamplitude_perCELL(1,1) = maxTRANSIENTS;
                        else     
                        TRANSIENTdurations_perCELL(1,1:TRANSIENTScount) = TRANSIENTdurations;                                              
                        TRANSIENTamplitude_perCELL(1,1:TRANSIENTScount) = maxTRANSIENTS;
                        end 
                        
                        TRANSIENTdurations_perCELL_s = TRANSIENTdurations*3.5; % transient durations in seconds
                        
                    
                end
                     
                %%% processes section

                elseif cellfun(@(s) any(s == lettertofind), c) == 1
                    %% find peaks
                    tracesample = ycorr(ycorr<protectPEAKS);
                    noise = std(tracesample(1:100)); %noise is defined as the standard deviation of the first 100 values of baseline adjusted trace that are <0.5)
                    MinPeakHeight = noise*noiseFactor;
                    MinPeakHeight = round(MinPeakHeight,2);


                    %% Find intersections
                    % find intersection between (0.001)and the peak curve, i.e. beginning and end of each peak

                    % create a line that marks minimal peak height
                    y = linspace(MinPeakHeight,MinPeakHeight,z); % a line y=MinPeakHeight
                    peaksonly = ycorr' - (y);
                    peaksonly(peaksonly<0)=0;

                    x1 = [0 z];
                    y1 = [0.001 0.001]; %vector in which betweenMinMax value is repeated for the length of x axis (distance)
                    [X0,Y0] = intersections(time,peaksonly,x1,y1);


                    %% Proceed with the code only if threshold for peaks is reached

                    if (~isempty(X0) == 1)==0;
                        clearvars peakSTART peakEND AREA realPEAKSonly maxPEAKS timeOFmaxPEAKs PEAKdurations timeOFmaxPEAKs
                        display([ROIlabel,'no peaks above threshold', MinPeakHeight]);
                        PEAKcount = 0;
                        maxPEAKS = 0;
                        PEAKdurations = 0;
                        timeOFmaxPEAKs = 0;
                        realPEAKSarea = 0;
                        timeOFspikes = zeros(1,15);

                        figure(i);
                        subplot(3,1,1)       % add first plot in 2 x 1 grid
                        plot(time,dFF,'k', time, yfit, 'b');
                        grid on

                        xlabel('time'); ylabel('dF/F')
                        ROIlabel = txt(n(i)*z+3,4);
                        xmin = 0;
                        xmax = z;
                        set(gca,'XTick',xmin:10:xmax);
                        xlim([xmin xmax]);
                        ymin = -1;
                        title(ROIlabel); legend('dF/F','baseline fit');

                        subplot(3,1,2)       % add second plot in 2 x 1 grid
                        plot(time,ycorr,'k',time, y, 'r');
                        grid on

                        xlabel('time'); ylabel('dF/F with baseline correction')
                        xmin = 0;
                        xmax = z;
                        set(gca,'XTick',xmin:10:xmax);
                        xlim([xmin xmax]);
                        ymin = -1;
                        noiseLABEL = strcat('noise< ', num2str(MinPeakHeight));
                        legend('ycorr',num2str(noiseLABEL));

                        FigHandle = figure(i);
                        set(FigHandle, 'Position', [-1200, 200, 1080, 810]); % use [100, 100, resolution, resolution] to display fig in the middle
                        
 
                    elseif(~isempty(X0) == 1)==1;

                        %% Check if peaks are real (min area condition)

                        clearvars peakSTART peakEND AREA realPEAKSonly realPEAKstart realPEAKend realPEAKSarea

                        numberOFpeaks = length(X0)/2;


                        for k = 1:numberOFpeaks;
                        peakSTART(k) = round(X0(2*k-1));
                        peakEND(k) = round(X0(2*k));
                        AREA(k)=trapz(peaksonly(peakSTART(k):peakEND(k)));

                        end

                        for f = 1:numberOFpeaks;

                                if AREA(f) > minPEAKarea         
                                realPEAKSonly(f) = 1;
                                realPEAKSarea(f) = AREA(f);

                                else 
                                realPEAKSonly(f) = 0;
                                realPEAKSarea(f) = 0;
                                end

                        end

                        %% find beginnings and ends of real peaks

                        allPEAKcount = length(realPEAKSonly);
                        
                        for g = 1:allPEAKcount;

                                if realPEAKSonly(g) == 1          
                                realPEAKstart(g) =  peakSTART(g);

                                elseif realPEAKSonly(g) == 0          
                                realPEAKstart(g) =  NaN;

                                end

                        end

                        for g = 1:length(realPEAKSonly);

                                if realPEAKSonly(g) == 1          
                                realPEAKend(g) =  peakEND(g);

                                elseif realPEAKSonly(g) == 0          
                                realPEAKend(g) =  NaN;

                                end

                        end

                        realPEAKstart(find(isnan(realPEAKstart))) = []; % removes NaNs from an array
                        realPEAKend(find(isnan(realPEAKend))) = [];

                        PEAKcount = length(realPEAKend);

                        %% find maximum value in a peak
                        if (~isempty(realPEAKend) == 1)==0;
                        display([ROIlabel,'no peaks with min area', MinPeakHeight]);
                        display(AREA');

                        clearvars maxPEAKS timeOFmaxPEAKs PEAKdurations
                        PEAKcount = 0;
                        maxPEAKS = 0;
                        PEAKdurations = 0;
                        timeOFmaxPEAKs = 0;
                        timeOFspikes = zeros(1,15);

                        elseif (~isempty(realPEAKend) == 1)==1;
                        display([ROIlabel,'real peaks detected', MinPeakHeight]);
                        display(AREA');    

                        clearvars maxPEAKS timeOFmaxPEAKs PEAKdurations
                        
                        timeOFspikes = zeros(10,15);
                        
                            for d = 1:length(realPEAKstart);

                                duration = realPEAKstart(d):1:realPEAKend(d);
                                maxPEAK = max(peaksonly(duration));
                                timeOFmaxPEAK = find(peaksonly==maxPEAK);
                                maxPEAKS(d) = maxPEAK;
                                timeOFmaxPEAKs(d) = timeOFmaxPEAK;                           
                                PEAKdurations(d) = length(duration);
                                timeOFspikes(d,1:length(duration)) = duration;

                            end
                            
                            timeOFspikes(~any(timeOFspikes,2),: ) = []; %keeps non zero rows
                            
                        end

                        %% create a plot

                        figure(i);                    
                        subplot(3,1,1)       % add first plot in 2 x 1 grid
                        plot(time,dFF,'k', time, yfit, 'b');
                        grid on

                        xlabel('time'); ylabel('dF/F')
                        ROIlabel = txt(n(i)*z+3,4);
                        xmin = 0;
                        xmax = z;
                        set(gca,'XTick',xmin:10:xmax);
                        xlim([xmin xmax]);
                        ymin = -1;
                        title(ROIlabel); legend('dF/F','baseline fit');

                        subplot(3,1,2)       % add second plot in 2 x 1 grid
                        plot(time,ycorr,'k',time, y, 'r');
                        grid on

                        xlabel('time'); ylabel('dF/F with baseline correction')
                        xmin = 0;
                        xmax = z;
                        set(gca,'XTick',xmin:10:xmax);
                        xlim([xmin xmax]);
                        ymin = -1;
                        noiseLABEL = strcat('noise< ', num2str(MinPeakHeight));
                        legend('ycorr',num2str(noiseLABEL));

                        subplot(3,1,3)       % add third plot in 2 x 1 grid
                        plot(time,peaksonly,'k'); grid on
                        hold on

                        if (~isempty(realPEAKend) == 1)==1;
                            plot(timeOFmaxPEAKs, maxPEAKS,'r^', 'MarkerFaceColor','r')

                        end

                        xlabel('time'); ylabel('dF/F with baseline correction')
                        xmin = 0;
                        xmax = z;
                        set(gca,'XTick',xmin:10:xmax);
                        xlim([xmin xmax]);
                        ymin = -1;

                        if (~isempty(realPEAKend) == 1)==1; % display real peaks legends only if real peaks exist 
                            PEAKcountLEGEND = strcat('Real peaks: ', num2str(PEAKcount));
                            legend('Peaks only',PEAKcountLEGEND);
                        else 
                            legend('Peaks only');
                        end


                        FigHandle = figure(i);
                        set(FigHandle, 'Position', [-1200, 200, 1080, 810]);
                        
                       
                                                
                    end
                       


                end
                
                % save figures

                        timePOINT = filename(1:3);
                        graphname = strcat(timePOINT,'_graph_', ROIlabel{1});

                        img = getframe(gcf);
                        rootpath2 = strcat(rootpath_mainFOLDER, filename,'/', graphname);
                        imwrite(img.cdata, [rootpath2 , '.bmp']);
           

            end
            
             %% load data into a matrix

            PEAKcount_allROI(i)= PEAKcount;        
            maxPEAKS_allROI(i,1:PEAKcount) = maxPEAKS;
            PEAKdurations_allROI(i,1:PEAKcount) = PEAKdurations;
            timeOFmaxPEAKs_allROI(i,1:PEAKcount) = timeOFmaxPEAKs;
            realPEAKSarea( :, ~any(realPEAKSarea,1) ) = []; %keeps non zero columns
            realPEAKSarea_allROI(i,1:PEAKcount) = realPEAKSarea;
            
            [b,n] = size(timeOFspikes);
            for e = 1:b
            ROI_timeOFspikes = horzcat(i,timeOFspikes(e,:));
            timeOFspikes_allROI = vertcat(timeOFspikes_allROI,ROI_timeOFspikes);
            end
            
            % threshold all ROI
            threshold_allROI(i) = MinPeakHeight;
                        
        end
        
        timeOFspikes_allROI = timeOFspikes_allROI(2:end,:);
        
        
       %% clean up data
       PEAKcount_allROI = PEAKcount_allROI';
       maxPEAKS_allROI( :, ~any(maxPEAKS_allROI,1) ) = []; %keeps non zero columns
       PEAKdurations_allROI( :, ~any(PEAKdurations_allROI,1) ) = []; %keeps non zero columns
       realPEAKSarea_allROI( :, ~any(realPEAKSarea_allROI,1) ) = []; %keeps non zero columns
       timeOFspikes_allROI( :, ~any(timeOFspikes_allROI,1) ) = []; %keeps non zero columns

       %% Remove peaks as specified in the file
        if sheetValid==1
           
       [FIGtoREMOVE, PEAKtoREMOVE] = read_num(num2);
       
       
       if FIGtoREMOVE==0
          disp('No sheaths were removed')
       else
            [Reason]=read_txt(txt2);
            for removal=1:length(FIGtoREMOVE)
       
            % mark peaks to be removed on the figure

            figure(FIGtoREMOVE(removal));
            subplot(3,1,3) 
            hold on
            message = ['\rightarrow ', Reason(removal)];
            plot(timeOFmaxPEAKs_allROI(FIGtoREMOVE(removal),PEAKtoREMOVE(removal)), maxPEAKS_allROI(FIGtoREMOVE(removal),PEAKtoREMOVE(removal)),'b^', 'MarkerFaceColor','b');
            text(timeOFmaxPEAKs_allROI(FIGtoREMOVE(removal),PEAKtoREMOVE(removal)), maxPEAKS_allROI(FIGtoREMOVE(removal),PEAKtoREMOVE(removal)) , message);
            legend('Peaks only', 'Real peaks', 'Removed peaks');

            % remove peaks

            PEAKcount_allROI(FIGtoREMOVE(removal)) = PEAKcount_allROI(FIGtoREMOVE(removal))-1; % removes the chosen peak from peakcount
            maxPEAKS_allROI(FIGtoREMOVE(removal),PEAKtoREMOVE(removal)) = 0; % replaces chosen peak amplitude with zero
            PEAKdurations_allROI(FIGtoREMOVE(removal),PEAKtoREMOVE(removal)) = 0; % replaces chosen peak duration with zero
            realPEAKSarea_allROI(FIGtoREMOVE(removal),PEAKtoREMOVE(removal)) = 0; % replaces chosen peak area with zero
            
            %remove peak values from time of peaks
            timeOFspikes_toREMOVE = timeOFspikes_allROI(timeOFspikes_allROI(:,1)==FIGtoREMOVE(removal),:);
            timeOFspikes_toREMOVE = timeOFspikes_toREMOVE(PEAKtoREMOVE(removal),:);
                                                
            timeOFspikes_allROI(ismember(timeOFspikes_allROI,timeOFspikes_toREMOVE,'rows')==1,:)=zeros;
            ind = find(timeOFspikes_allROI(:,1)==0);
            timeOFspikes_allROI(ind,1)=FIGtoREMOVE(removal);
            
            % overwrite corrected figure
            timePOINT = filename(1:3);
            graphname = strcat(timePOINT,'_graph_', AllROIlabel{FIGtoREMOVE(removal)});

            img = getframe(gcf);
            rootpath2 = strcat(rootpath_mainFOLDER, filename,'/', graphname);
            imwrite(img.cdata, [rootpath2 , '.bmp']);
            
            end
       end

        end
       
   
          
%% Amplitude, duration and area
                
        % average duration of all spikes for a given ROI in a timepoint
        [rws,cols] = size(PEAKdurations_allROI);
        for row = 1:rws
        row_withZEROS = PEAKdurations_allROI(row,:);
        row_durations = row_withZEROS(row_withZEROS~=0);
        row_durations_mean = mean(row_durations);
        durations_mean_allROI(row) = row_durations_mean;
        end
        
        durations_mean_allROI(isnan(durations_mean_allROI))=0;
        durations_mean_S_allROI = durations_mean_allROI*3.5;
        durations_mean_S_allROI = durations_mean_S_allROI';
      
        % average amplitude of all spikes for a given ROI in a timepoint
        maxPEAKS_allROI(maxPEAKS_allROI == 0) = NaN;
        maxPEAKSdFF_allROI = bsxfun(@plus,maxPEAKS_allROI,threshold_allROI');
        maxPEAKSdFF_allROI(isnan(maxPEAKSdFF_allROI))=0;
        
        [rws,cols] = size(maxPEAKSdFF_allROI);
        for row = 1:rws
        row_withZEROS = maxPEAKSdFF_allROI(row,:);
        row_amplitudes = row_withZEROS(row_withZEROS~=0);
        row_amplitudes_mean = mean(row_amplitudes);
        amplitudes_mean_allROI(row) = row_amplitudes_mean;
        end
        
        amplitudes_mean_allROI(isnan(amplitudes_mean_allROI))=0;
                
        % average area of all spikes for a given ROI in a timepoint
        realPEAKS_base_allROI = bsxfun(@times,threshold_allROI',PEAKdurations_allROI); %multiply the size of treshold by the duration to establish area base
        realPEAKSarea_allROI(realPEAKSarea_allROI == 0) = NaN; %remove zeros
        realPEAKSarea_WITH_base_allROI = bsxfun(@plus,realPEAKSarea_allROI,realPEAKS_base_allROI);
        realPEAKSarea_WITH_base_allROI(isnan(realPEAKSarea_WITH_base_allROI))=0;
        
        [rws,cols] = size(realPEAKSarea_WITH_base_allROI);
        for row = 1:rws
        row_withZEROS = realPEAKSarea_WITH_base_allROI(row,:);
        row_area = row_withZEROS(row_withZEROS~=0);
        row_area_mean = mean(row_area);
        areas_mean_allROI(row) = row_area_mean;
        end
        
        areas_mean_allROI(isnan(areas_mean_allROI))=0;
                       
      
        %% enter data for all timepoints (files) into a matrix
        
        % peak count
        NEW_PEAKcount_allTimePoints = horzcat(PEAKcount_allTimePoints,PEAKcount_allROI);
        PEAKcount_allTimePoints = NEW_PEAKcount_allTimePoints;
        
        % durations
        NEW_mean_Durations_allTimePoints = horzcat(mean_Durations_allTimePoints,durations_mean_S_allROI);
        mean_Durations_allTimePoints = NEW_mean_Durations_allTimePoints;
        
        % amplitudes
        NEW_mean_Amplitude_allTimePoints = horzcat(mean_Amplitude_allTimePoints,amplitudes_mean_allROI');
        mean_Amplitude_allTimePoints = NEW_mean_Amplitude_allTimePoints;
        
        % area                        
        NEW_mean_Area_allTimePoints = horzcat(mean_Area_allTimePoints,areas_mean_allROI');
        mean_Area_allTimePoints = NEW_mean_Area_allTimePoints;
        
                               
        % save filenames into a string
        ALL_timePOINTs(fls) = cellstr(timePOINT);
       
      close all;
      
      ROIlabels = AllROIlabel;
            
      clearvars PEAKcount_allROI maxPEAKS_allROI PEAKdurations_allROI realPEAKSarea_allROI AllROIlabel timeOFtransients 
            
 
end

%% data for all timepoints 

% sum of peak count per ROI
PEAKcount_sum = sum(PEAKcount_allTimePoints,2); % summ of all rows (2 denotes dimension)

% average amplitude (exclude zeros - only existing spikes)
        [rws,cols] = size(mean_Amplitude_allTimePoints);
        for row = 1:rws
        row_withZEROS = mean_Amplitude_allTimePoints(row,:);
        row_amplitudes = row_withZEROS(row_withZEROS~=0);
        row_amplitudes_mean = mean(row_amplitudes);
        Mean_ROI_amplitude(row) = row_amplitudes_mean;
        end

        Mean_ROI_amplitude(isnan(Mean_ROI_amplitude))=0;
      
% average amplitude (exclude zeros - only existing spikes)    
        [rws,cols] = size(mean_Durations_allTimePoints);
        for row = 1:rws
        row_withZEROS = mean_Durations_allTimePoints(row,:);
        row_durations = row_withZEROS(row_withZEROS~=0);
        row_durations_mean = mean(row_durations);
        Mean_ROI_duration(row) = row_durations_mean;
        end
        
        Mean_ROI_duration(isnan(Mean_ROI_duration))=0;
        
% average area (exclude zeros - only existing spikes)            
        [rws,cols] = size(mean_Area_allTimePoints);
        for row = 1:rws
        row_withZEROS = mean_Area_allTimePoints(row,:);
        row_area = row_withZEROS(row_withZEROS~=0);
        row_area_mean = mean(row_area);
        Mean_ROI_area(row) = row_area_mean;
        end
        
        Mean_ROI_area(isnan(Mean_ROI_area))=0;

%% save

savename = 'Results_matlab';
rootpathRESULTS = strcat(rootpath_mainFOLDER, 'Results','/');
xlsfilename = [rootpathRESULTS, savename, '.xls'];

xlswrite(xlsfilename, ROIlabels', 'count', 'A2'); xlswrite(xlsfilename, ALL_timePOINTs, 'count', 'B1'); 
xlswrite(xlsfilename, PEAKcount_allTimePoints, 'count', 'B2'); xlswrite(xlsfilename, PEAKcount_sum, 'count', 'I2'); 

xlswrite(xlsfilename, ROIlabels', 'amplitude', 'A2'); xlswrite(xlsfilename, ALL_timePOINTs, 'amplitude', 'B1'); 
xlswrite(xlsfilename, mean_Amplitude_allTimePoints, 'amplitude', 'B2'); xlswrite(xlsfilename, Mean_ROI_amplitude', 'amplitude', 'I2'); 

xlswrite(xlsfilename, ROIlabels', 'duration', 'A2'); xlswrite(xlsfilename, ALL_timePOINTs, 'duration', 'B1'); 
xlswrite(xlsfilename, mean_Durations_allTimePoints, 'duration', 'B2'); xlswrite(xlsfilename, Mean_ROI_duration', 'duration', 'I2'); 

xlswrite(xlsfilename, ROIlabels', 'area', 'A2'); xlswrite(xlsfilename, ALL_timePOINTs, 'area', 'B1'); 
xlswrite(xlsfilename, mean_Area_allTimePoints, 'area', 'B2'); xlswrite(xlsfilename, Mean_ROI_area', 'area', 'I2'); 
