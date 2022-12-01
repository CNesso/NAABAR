%function bcgFun(EGG,secondAlign,plotVis,menuselect)
clc
clear
close all
%%
addpath('~/ownCloud/EEG/eeglab2019_1/');
addpath('~/ownCloud/EEG/Michael_functions/');

[qrsChoice,dataFolder]=uigetfile('*.vhdr','MultiSelect','on');
if ~iscell(qrsChoice),qrsChoice={qrsChoice};end
cd(dataFolder);
meanMethodArr={'gauss','mean','smooth'};
for meanm=1:2%length(meanMethodArr)
    meanMethod=meanMethodArr{meanm}
    if ~isfolder(meanMethod)
        mkdir([dataFolder meanMethod]);
    end
    saveFolder=[dataFolder meanMethod];
    for ji=1:length(qrsChoice)
        try
        bcgFunc(qrsChoice{ji},dataFolder,meanMethod,saveFolder)
        catch ME
            er.chan=qrsChoice{ji};
            er.message=ME.message;
            warning(er.message);
            continue;
        end
    end
end
function bcgFunc(qrsChoice,dataFolder,meanMethod,saveFolder)
eeglab;

EEG = pop_loadbv(dataFolder,qrsChoice);
%EEG = eeg_checkset( EEG );
EEG.setname=qrsChoice(1:end-5)
%EEG=pop_chanedit(EEG, 'lookup','/home/ne550/ownCloud/EEG/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
EEG = eeg_checkset( EEG );
eeglab redraw;
if EEG.srate<5000
    [EEG] = pop_resample( EEG, 5000);
end

secondAlign=1;
noMenu=1;
preFilter=false;
plotVis=false;
menuSelect=false;
saveQrsToFile=true;
debaseStart=EEG.srate/2;
maxD=EEG.srate*0.05;
cutTollerance=EEG.srate*0.02;
trimPart=30;% percentage of trimmean
fadeRange=maxD*2;
minAvg=0.8;
refW=30;
eventType='qrs1';%backcompatibility overwritten in hep
plotSample=false;%'FCz';%false;%'FCz';%true;+del2
qrsLpf=25;%30
qrsHpf=0;
onlySample=0;

a=[EEG.event(strcmp({EEG.event.type},'R128')).latency];
subMa=round(mean(diff(a(1:10))));
subBreak=find(diff(a)>2*subMa);
if subBreak
    subEnd=[a(subBreak)+2*subMa];
    subBegin=[a(1)-subMa a(subBreak+1)];
    subEnd=[subEnd a(end)+subMa];in-own
    EEG.subsections=[subBegin' subEnd'];
else
    EEG.subsections=[1 EEG.pnts];
end
eventNames={  'S  1'  'S  2'  'S  3'  'S  4'  'S  5' };
EEG=eventeResponseFilter(EEG,eventNames,'S129',EEG.srate);
eeg_eventtypes(EEG)


%% ECG creation
qrsFile=dir([dataFolder qrsChoice(1:end-find(fliplr(qrsChoice=='.'))) '*.txt']);
if logical(length(qrsFile))
    qrsFile(1).num=str2num(qrsFile(1).name(end-5));
    qrsFile(1).enable=true;
else
    clear qrsFile
end
if ~isfield(EEG,'ECG')
    EEG.ECG=pop_select(EEG,'channel',{'ECG1','ECG2'});
    EEG.ECG.data(3,:)=(EEG.ECG.data(1,:)+mean(EEG.ECG.data(1,1:EEG.ECG.srate*10)))+(EEG.ECG.data(2,:)+mean(EEG.ECG.data(2,1:EEG.ECG.srate*10)));
end
if(exist('qrsFile','var'))
    if ~qrsFile(1).enable
        EEG.ECG=pop_fmrib_qrsdetect(EEG.ECG,1,'qrs1','no');
        EEG.ECG=pop_fmrib_qrsdetect(EEG.ECG,2,'qrs2','no');
        EEG.ECG=pop_fmrib_qrsdetect(EEG.ECG,3,'qrs3','no');
        % EEG.event=eegT.event;
        % EEG.eventdescription=eegT.eventdescription;
        % EEG.urevent=eegT.urevent;
        %%
    end
end
EEG=pop_select(EEG,'nochannel',{'ECG1','ECG2'});


if noMenu
    qrsA=load([dataFolder qrsFile(1).name]);
else
    completedQrsCorrection=2;
    while (completedQrsCorrection-1)
        
        if plotSample
            plotString='Plot ON';
            michael_hep;
        else
            plotString='Plot OFF';
            michael_hep;
        end
        completedQrsCorrection=menu('QRS correction','Finished','Resume','To ECG1','Input Delay',plotString,'Keyboard','Clear');
        if completedQrsCorrection==7
            clear HEP;
            if (~exist('qrsFile','var'))
                qrsFile(1).enable=false;
            end
            completedQrsCorrection=2;
        end
        if completedQrsCorrection==3
            ecg = double(EEG.ECG.data(1,:));
            ecg = ecg(:); % make sure ecg is one-column vector
            
            HEP.ecg = (ecg-min(ecg))/(max(ecg)-min(ecg));
            completedQrsCorrection=2;
        end
        if completedQrsCorrection==4
            delayQrs=input('Delay(ms): ')/1000;
            HEP.qrs=HEP.qrs+delayQrs*EEG.srate
            completedQrsCorrection=2;
        end
        if completedQrsCorrection==5
            plotSample=~plotSample;
        end
        if completedQrsCorrection==6
            keyboard;
        end
    end
    
    if isgraphics(333);close(333);end
    if isgraphics(hepgui.main);close(hepgui.main);end
    qrsA=HEP.qrs;
    clear HEP
    if saveQrsToFile
        save([dataFolder qrsChoice(1:end-find(fliplr(qrsChoice=='.'))) '_' char(num2str(chanChoice)) '.txt'],'qrsA','-ascii')
    end
end
 eegT=eeg_addnewevents(EEG,{qrsA(1:end)'},{char(eventType)});%first and last events can be excluded to increse stability

%         tic
%         eegT= fmrib_pas(EEG,qrsA','mean')
%         toc

tic
if (~exist('eegF','var')) &&  preFilter
    eegF=pop_eegfiltnew(eegT,0,qrsLpf);
end
fmri_hr_rem2

if plotVis
    menu('accept deletion',"OK");
end
close all
pop_saveset( pop_resample(eegT, 500), 'filename',[qrsChoice(1:end-find(fliplr(qrsChoice=='.'))) '_' meanMethod '.set'],'filepath',saveFolder);%EEG = pop_saveset( EEG, 'filename',[EEG.setname 'noMRI.set'],'filepath','/run/media/ne550/TOSHBIG/EEG/');
end
