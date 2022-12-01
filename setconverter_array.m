% i have to remove the first and last 2 secondsfor comparability
% i have to remove the first and last 2 secondsfor comparability
clc
clear
close all
%%
addpath('~/ownCloud/EEG/eeglab2019_1/');
addpath('~/ownCloud/EEG/Michael_functions/');
eeglab;
eventNames={  'S  1' 'S  2' 'S  3' 'S  4' 'S  5' 'R'  };
%%
[fileNames,filefolders]=uigetfile('*.vhdr','MultiSelect','on');
currentData.files=struct('name',fileNames,'folder',filefolders);
dataFolder=currentData.files(1).folder;
cd(dataFolder);
if ~isfolder('Processed')
    mkdir([dataFolder 'Processed']);
end
saveFolder=[dataFolder 'Processed'];
outputSampleRate=160;
% if isfolder('Processed')
%     saveFolder=[dataFolder 'Processed'];
% else
%     saveFolder=uigetdir(dataFolder);
% end
%%
%raSampleCell=cell(length(currentData.files),2);
for iData=1:length(currentData.files)
    disp(iData);
    currentFile=currentData.files(iData);
    EEG = pop_loadbv(currentFile.folder,currentFile.name);
    EEG.setname=currentFile.name(1:end-5);
    EEG=pop_select(EEG,'nochannel',{'ECG1','ECG2','ECG3'});
    eventNames={  'S  1'  'S  2'  'S  3'  'S  4'  'S  5' };
    EEG=eventeResponseFilter(EEG,eventNames,'S129',EEG.srate);
    
    %%EEG = eeg_checkset( EEG );
    % keyboard
    if EEG.srate~=outputSampleRate
        EEG=pop_resample(EEG, outputSampleRate)
    end
    
    %ra=[EEG.event(strcmp({EEG.event.type},'R')).latency];
    EEG.event(strcmp({EEG.event.type},'R'))=[];
    
    ignoreMargin=10;
    qrsSampleNum=100;
    %EEG.urevent=EEG.event;
    raSample=raSampleCell{iData,1};
    % raSample=ra(randperm(length(ra)-2*ignoreMargin,qrsSampleNum)+ignoreMargin);
    EEG=eeg_addnewevents(EEG,{raSample},{'R-Peak'})
    % EEG = eeg_checkset( EEG );
    % raSampleCell{iData,1}=raSample;
    % raSampleCell{iData,2}=EEG.srate;
    % save([currentFile.folder 'raSample_' char(num2str(iData)) '.txt'],'raSampleCell','-ascii')
    
    %     r128=[EEG.event(strcmp({EEG.event.type},'R128')).latency];
    %     eventI=find(strcmp({EEG.event.type},'R')&...
    %         ([EEG.event.latency]<r128(2)|...
    %         [EEG.event.latency]>r128(end)));
    %     EEG=pop_selectevent(EEG,'omitevent',eventI);
    %
    
    
    % EEG=pop_select(EEG,'nopoint',[1 r128(2)]);
    % %to avoid first and last heartbeast
    % EEG.event(strcmp({EEG.event.type},'R')&...
    %     ([EEG.event.latency]<r128(2)|[EEG.event.latency]>r128(end)))...
    %     =[];
    
    % EEG.urevent(strcmp({EEG.urevent.type},'R')&...
    %     ([EEG.urevent.latency]<r128(2)|[EEG.urevent.latency]>r128(end)))...
    %     =[];
    % EEG = eeg_checkset( EEG );
    %keyboard
    
    %iArr.R=[EEG.event(strcmp({EEG.event.type},'R')).latency];
    
    
    %     iArr.event=EEG.event;
    %     iArr.urevent=EEG.urevent;
    %     iArr.eventdescription=EEG.eventdescription;
    %     iArr.chanlocs=EEG.chanlocs;
    %     iArr.urchanlocs=EEG.urchanlocs;
    %     arr.(currentFile.name(1:5)).((currentFile.name(end-8:end-5)))...
    %         =iArr;
    %
    
    
    %    =struct('event',EEG.event,'eventdescription',EEG.eventdescription,'urevent',EEG.urevent);
    %arr.(currentFile.name(1:5))(str2double(currentFile.name(end-5)))...
    EEG=setFilter(EEG);
    %     EEG=pop_epoch(EEG,eventNames, [-1  1]);
    %     EEG=pop_rmbase(EEG, [-1000 -200] ,[]);
    pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath',saveFolder);%EEG = pop_saveset( EEG, 'filename',[EEG.setname 'noMRI.set'],'filepath','/run/media/ne550/TOSHBIG/EEG/');
end
% sRate=EEG.srate;
%save('raCell','raSampleCell')
%%
load('raCell')
%%
outputSampleRate=160;
while 1
    [fileNames,filefolders]=uigetfile('*.set','MultiSelect','on');
    currentData.files=struct('name',fileNames,'folder',filefolders);
    dataFolder=currentData.files(1).folder;
    cd(dataFolder);
    if ~isfolder('Processed')
        mkdir([dataFolder 'Processed']);
    end
    saveFolder=[dataFolder 'Processed'];
    %%
    for iData=1:length(currentData.files)
        disp(iData);
        currentFile=currentData.files(iData);
        EEG = pop_loadset(currentFile.name,currentFile.folder);
        %  r128=[EEG.event(strcmp({EEG.event.type},'R128')).latency];
        % EEG=pop_select(EEG,'nopoint',[1 r128(2)]);
        removed_channels=[];
        for i=1:size(EEG.data,1)
            naTimes=EEG.times(find(isinf(EEG.data(i,:))|isnan(EEG.data(i,:))));
            if ~isempty(naTimes)
                disp(currentFile.name)
                disp(EEG.chanlocs(i).labels)
                %disp(naTimes)
               removed_channels(end+1)=i;
            end
        end
        if ~isempty(removed_channels)
            disp(removed_channels)
        EEG = pop_interp(EEG,removed_channels, 'spherical');
        end
        
        if EEG.srate~=outputSampleRate
            EEG=pop_resample(EEG, outputSampleRate)
        end
        if isempty(EEG.urevent)
            EEG.urevent=EEG.event;
        end
        raSample=raSampleCell{iData,1};
        EEG=eeg_addnewevents(EEG,{raSample},{'R-Peak'});
        
        
        % %     iArr=arr.(currentFile.name(1:5)).(currentFile.name(end-12:end-9));
        %     iArr=arr.(currentFile.name(1:5)).(currentFile.name(strfind(currentFile.name,'seg'):strfind(currentFile.name,'seg')+3));
        %     %iArr=arr.(currentFile.name(1:5)).((currentFile.name(end-8:end-5)))
        %     for i=fieldnames(iArr)'
        %         disp(i)
        %         EEG.(i{1})=iArr.(i{1});
        %     end
        %
        %     %EEG=eeg_addnewevents(EEG,{iArr.R},{'R'});
        
        %EEG = eeg_checkset( EEG );
        EEG=pop_eegfiltnew(EEG,0,70);
        EEG=setFilter(EEG);
        %EEG = pop_firws(EEG, 'fcutoff', [1 70], 'ftype', 'bandpass', 'wtype', 'blackman', 'forder', 1376, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
        %     EEG=pop_epoch(EEG,eventNames, [-1  1]);
        %     EEG=pop_rmbase(EEG, [-1000 -200] ,[]);
        pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath',saveFolder);%EEG = pop_saveset( EEG, 'filename',[EEG.setname 'noMRI.set'],'filepath','/run/media/ne550/TOSHBIG/EEG/');
    end
end
function EEGclean=setFilter(EEG)
notch50=true;
artefactLengthSamples=289.419354838710;%289;%MB1
%artefactLengthSamples=300;%MB2
artefactLength=artefactLengthSamples/5000;
lowPassF=68;%chosen to avoid the resonance at 69Hz
NoiseFreq=[1/artefactLength:1/artefactLength:min(EEG.srate/2,lowPassF) min(EEG.srate/2,lowPassF) 50 ];
% plotChan=60;
% skipIter=false;
% while true

if notch50
    EEGclean = pop_firws(EEG, 'fcutoff', [49 51], 'ftype', 'bandstop', 'wtype', 'blackman', 'forder', 1376, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
else
    %EEGclean =pop_eegfiltnew(EEG, 'locutoff',49,'hicutoff',51,'revfilt',1);
    lineNoiseIn = struct('lineNoiseChannels', [1:EEG.nbchan],...
        'Fs', EEG.srate, ...
        'lineFrequencies', NoiseFreq,...
        'p',0.05, ...
        'fScanBandWidth', 4, ...
        'taperBandWidth',3, ...
        'taperWindowSize', 4, ...
        'taperWindowStep', 4, ...
        'tau', 100, ...
        'pad', 0, ...
        'fPassBand', [0 EEG.srate/2], ...
        'maximumIterations', 10);
    EEGclean = cleanLineNoise(EEG, lineNoiseIn);
end
EEGclean = pop_firws(EEGclean, 'fcutoff', [1 lowPassF], 'ftype', 'bandpass', 'wtype', 'blackman', 'forder', 1376, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
end

