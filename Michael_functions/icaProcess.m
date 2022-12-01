function [eegA,currentDataNum]=icaProcess(eegA,currentDataNum,icaResponse,varargin)
%EEG=eegA.ica;
currentData=eegA(currentDataNum);
EEG=currentData.eeg;
if nargin<4
    eventNames={  'S  1'  'S  1 MISS' 'S  2' 'S  2 MISS' 'S  3' 'S  3 MISS' 'S  4'  'S  4 MISS'  'S  5' 'S  5 MISS' 'qrs1' 'R128'  };
    eventDescriptions={  'Responded WG'  'Missed WG' 'Responded SG'  'Missed SG'  'Responded WNG' 'InhibitedWNG' 'Responded SNG' 'Inhibited SNG'  'Ct' 'Missed Ct' 'Heart R-Peak' 'MRI R' };
    eventNamesPlot={  'S  1'  'S  2'  'S  3 MISS'  'S  4'  };
end
if isempty(currentData.epochs)
            %eventNames={  'S  1'  'S  1 MISS' 'S  2' 'S  2 MISS' 'S  3' 'S  3 MISS' 'S  4'  'S  4 MISS'  'S  5' 'S  5 MISS' 'qrs1' 'R128'  };
            %eventDescriptions={  'Responded WG'  'Missed WG' 'Responded SG'  'Missed SG'  'Responded WNG' 'InhibitedWNG' 'Responded SNG' 'Inhibited SNG'  'Ct' 'Missed Ct' 'Heart R-Peak' 'MRI R' };
            %EEG = pop_epoch( EEG, eventNames, [-0.5  1], 'newname', char(string(EEG.setname)+'_'+string(EEG.srate)+'Hz_Epoched'), 'epochinfo', 'yes');
            [evNa,~]=eeg_eventtypes(EEG);
            eventDescriptions=eventDescriptions(ismember(eventNames,evNa));
            eventNames=eventNames(ismember(eventNames,evNa));
            for i=1:length(eventNames)
                regexprep(eventNames{i},'\W','')
                currentData.epochs.(regexprep(eventNames{i},'\W',''))=pop_epoch(EEG,eventNames(i), [-0.5  1]);
                currentData.epochs.(regexprep(eventNames{i},'\W','')).description=eventDescriptions{i};
                currentData.epochs.(regexprep(eventNames{i},'\W',''))=pop_rmbase(currentData.epochs.(regexprep(eventNames{i},'\W','')), [-500 0] ,[]);
                %     else
                %        eventNames(i)=[];
            end
            %eventNamesPlot={  'S  1'  'S  2'  'S  3 MISS'  'S  4'  };
            %currentData.epochs.('combined')=pop_epoch(EEG,eventNames(i), [-0.5  1]);
            currentData.epochs.('combined')=pop_epoch(EEG,eventNamesPlot, [-0.5  1]);
            currentData.epochs.('combined').description='Stimulycombined';
            currentData.epochs.('combined')=pop_rmbase(currentData.epochs.('combined'), [-500 0] ,[]);
end
    switch icaResponse
        case 'epoch'
            %eventNames={  'S  1'  'S  1 MISS' 'S  2' 'S  2 MISS' 'S  3' 'S  3 MISS' 'S  4'  'S  4 MISS'  'S  5' 'S  5 MISS' 'qrs1' 'R128'  };
            %eventDescriptions={  'Responded WG'  'Missed WG' 'Responded SG'  'Missed SG'  'Responded WNG' 'InhibitedWNG' 'Responded SNG' 'Inhibited SNG'  'Ct' 'Missed Ct' 'Heart R-Peak' 'MRI R' };
            %EEG = pop_epoch( EEG, eventNames, [-0.5  1], 'newname', char(string(EEG.setname)+'_'+string(EEG.srate)+'Hz_Epoched'), 'epochinfo', 'yes');
            [evNa,~]=eeg_eventtypes(EEG);
            eventDescriptions=eventDescriptions(ismember(eventNames,evNa));
            eventNames=eventNames(ismember(eventNames,evNa));
            for i=1:length(eventNames)
                regexprep(eventNames{i},'\W','')
                currentData.epochs.(regexprep(eventNames{i},'\W',''))=pop_epoch(EEG,eventNames(i), [-0.5  1]);
                currentData.epochs.(regexprep(eventNames{i},'\W','')).description=eventDescriptions{i};
                currentData.epochs.(regexprep(eventNames{i},'\W',''))=pop_rmbase(currentData.epochs.(regexprep(eventNames{i},'\W','')), [-500 0] ,[]);
                %     else
                %        eventNames(i)=[];
            end
            %eventNamesPlot={  'S  1'  'S  2'  'S  3 MISS'  'S  4'  };
            %currentData.epochs.('combined')=pop_epoch(EEG,eventNames(i), [-0.5  1]);
            currentData.epochs.('combined')=pop_epoch(EEG,eventNamesPlot, [-0.5  1]);
            currentData.epochs.('combined').description='Stimulycombined';
            currentData.epochs.('combined')=pop_rmbase(currentData.epochs.('combined'), [-500 0] ,[]);
        case 'viewComponent'
            pop_viewprops(EEG,0);
        case 'icaBlink'
            icablinkmetrics(EEG,'ArtifactChannel',EEG.data(find(strcmp({EEG.chanlocs.labels},'EOGR')),:),'MetricThresholds',[0.001 0.001 0.001],'VisualizeData','True');
            icablinkmetrics(EEG,'ArtifactChannel',EEG.data(find(strcmp({EEG.chanlocs.labels},'EOGL')),:),'MetricThresholds',[0.001 0.001 0.001],'VisualizeData','True');
            %  EEG.icaquant = icablinkmetrics(EEG,'ArtifactChannel',EEG.data(find(strcmp({EEG.chanlocs.labels},'EOGR')),:),'MetricThresholds',[0.001 0.001 0.001],'VisualizeData','True');
            %  EEG.icaquant = icablinkmetrics(EEG,'ArtifactChannel',EEG.data(find(strcmp({EEG.chanlocs.labels},'EOGL')),:),'MetricThresholds',[0.001 0.001 0.001],'VisualizeData','True');
            plotBlink=0;
        case 'erpComponentArtefact'
            figure(333);
            %figure('Units', 'normalized','PaperPositionMode','auto','InvertHardcopy','off');
            subplot(2,1,1);pop_envtopo(currentData.epochs.qrs1, [-500  996],'limcontrib',[-500 996],'compsplot',[7],'title', currentData.epochs.qrs1.description,'electrodes','off');
            %figure('Units', 'normalized','PaperPositionMode','auto','InvertHardcopy','off');
            subplot(2,1,2);pop_envtopo(currentData.epochs.R128, [-500  996],'limcontrib',[-500 996],'compsplot',[7],'title', currentData.epochs.R128.description,'electrodes','off');
        case 'erpComponentEvent'
            figure(444);
            % EEG=eegA.ica2;
            % EEG=currentData.epochs.S3MISS;
            % %EEG=currentData.epochs.combined;
            % EEG = eeg_checkset( EEG );
            %figure('Units', 'normalized','PaperPositionMode','auto','InvertHardcopy','off'); pop_envtopo(EEG, [-500  996],'limcontrib',[-500 996],'compsplot',[7],'title', EEG.description,'electrodes','off');
            
            %figure('Units', 'normalized','PaperPositionMode','auto','InvertHardcopy','off');
            subplot(2,1,1);pop_envtopo(currentData.epochs.S3MISS, [-500  996],'limcontrib',[0 500],'compsplot',[7],'title', currentData.epochs.S3MISS.description,'electrodes','off');
            %figure('Units', 'normalized','PaperPositionMode','auto','InvertHardcopy','off');
            subplot(2,1,2);pop_envtopo(currentData.epochs.combined, [-500  996],'limcontrib',[0 500],'compsplot',[7],'title', currentData.epochs.combined.description,'electrodes','off');
            
            % EEG=eegA.ica2;
            
            % EEG=currentData.epochs.combined;
            % EEG = eeg_checkset( EEG );
            % eeglab redraw;
        case 'erpCombined'
            figure('Name',sprintf('%d-Combined Erps',currentDataNum),'NumberTitle','off');
            hold on;
            smothErp=20;
            plot(currentData.epochs.S1.times,mean(smoothdata(currentData.epochs.S1.data(60,:,:),'gaussian',smothErp),3),'DisplayName',currentData.epochs.S1.description)
            plot(currentData.epochs.S2.times,mean(smoothdata(currentData.epochs.S2.data(60,:,:),'gaussian',smothErp),3),'DisplayName',currentData.epochs.S2.description)
            plot(currentData.epochs.S3MISS.times,mean(smoothdata(currentData.epochs.S3MISS.data(60,:,:),'gaussian',smothErp),3),'DisplayName',currentData.epochs.S3MISS.description)
            plot(currentData.epochs.S4.times,mean(smoothdata(currentData.epochs.S4.data(60,:,:),'gaussian',smothErp),3),'DisplayName',currentData.epochs.S4.description)
            xline(0)
            %legend('Responded WG','Responded SG','Responded WNG','InhibitedWNG','Responded SNG')
            legend
            set(gca,'YDir','reverse','YMinorTick','on','XMinorTick','on')
            ylabel('\muV')
            xlabel('Time (ms)')
            rl=refline([0 0])
            rl.Annotation.LegendInformation.IconDisplayStyle = 'off';
            rl.Color='black';
            rl.LineStyle='--';
            title({EEG.setname,'FCz(Gauss 20)'})
            hold off;
        case 'erpSingle'
            for i=eventNamesPlot
                try
            iterName=regexprep(i{1},'\W','')
            figure('Name',sprintf('%d-%s Erps',currentDataNum,iterName),'NumberTitle','off');
            pop_erpimage(currentData.epochs.(iterName),1, [60],[[]],currentData.epochs.(iterName).description,10,1,i,[],'latency' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [60] EEG.chanlocs EEG.chaninfo } );
            catch ME
                disp(getReport(ME));
                end
            end
            
%             figure(663); pop_erpimage(currentData.epochs.S3MISS,1, [60],[[]],currentData.epochs.S3MISS.description,10,1,{ 'S  3 MISS'},[],'latency' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [60] EEG.chanlocs EEG.chaninfo } );
%             figure(664); pop_erpimage(currentData.epochs.S4,1, [60],[[]],currentData.epochs.S4.description,10,1,{ 'S  4'},[],'latency' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [60] EEG.chanlocs EEG.chaninfo } );
%             %figure; pop_erpimage(currentData.epochs.S4MISS,1, [60],[[]],currentData.epochs.S4.description,10,1,{ 'S  4 MISS'},[],'latency' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [60] EEG.chanlocs EEG.chaninfo } );
%             figure(662); pop_erpimage(currentData.epochs.S2,1, [60],[[]],currentData.epochs.S2.description,10,1,{ 'S  2'},[],'latency' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [60] EEG.chanlocs EEG.chaninfo } );
%             figure(661); pop_erpimage(currentData.epochs.S1,1, [60],[[]],currentData.epochs.S1.description,10,1,{ 'S  1'},[],'latency' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [60] EEG.chanlocs EEG.chaninfo } );
        case 'deleteComponent'
            EEG.deletedComp=str2num(icaChooseDialog('Select num',num2str([find(EEG.reject.gcompreject == 1)'])));
            %pop_viewprops( EEG);
            %find(EEG.reject.gcompreject == 1)
            %EEG = pop_subcomp( EEG, find(EEG.reject.gcompreject == 1), 0);
           
            eegA(end+1).eeg=pop_subcomp( EEG,EEG.deletedComp, 0);
            eegA(end).epochs=[];
            eegA(currentDataNum)=currentData;
            currentDataNum=length(eegA);
            return;
        case 'saveSet'
             pop_saveset( EEG);
        case 'changeSet'
            %EEG =eegA(menu({'Select dataset'},arrayfun(@(x){num2str(x)},1:length(eegA))));
            currentDataNum=menu({'Select dataset'},arrayfun(@(x){num2str(x)},1:length(eegA)));
            return;
        case 'deleteSet'
%             EEG=eegA(1);
            try
                eegA(icaChooseDialog('Datasets to be deleted',num2str(1:length(eegA))))=[];
                currentDataNum=1;
                return;
            catch
                warning('Selection Error');
            end
        case 'loadNew'
            [icaChoice,dataFolder]=uigetfile('*ICA.set');
    currentDataNum=1;
    eegA=struct('eeg',[pop_loadset( icaChoice, dataFolder)],'epochs',[]);
    currentData=eegA;
    %%
%     EEG = pop_loadset( icaChoice, dataFolder);
    return;
%             eegA(1)=EEG;
        otherwise
    end
%  icaResponse=icaSelector(currentDataNum)
% currentData.eeg=EEG;
eegA(currentDataNum)=currentData;
function choice = icaChooseDialog(tiTle,defaultText)
% if isnumeric(defaultText)
%     defaultText=num2str(defaultText);
% end
d = dialog('Position',[300 300 250 150],'Name','Select');
txt = uicontrol('Parent',d,...
    'Style','text',...
    'Position',[20 80 210 40],...
    'String',tiTle...
    );

edt = uicontrol('Parent',d,...
    'String',defaultText,...
    'Style','edit',...
    'Position',[75 70 100 25],...
     'Callback',@qrs_popup_callback);

btn = uicontrol('Parent',d,...
    'Position',[89 20 70 25],...
    'String','OK',...
    'Callback',@qrs_popup_callback);

choice = defaultText;

% Wait for d to close before running to completion
uiwait(d);

    function qrs_popup_callback(popup,event)
        %idx = edt.Value;
        %popup_items = popup.String;
        %choice = char(popup_items(idx,:));
        choice = edt.String;
        %noMriFiles(idx).name;%char(popup_items(idx,:));
        delete(gcf);
    end
end
end
