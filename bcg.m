% while true
clc
clear
close all

loopProcess=false;
preFilter=false;

%%
addpath('~/ownCloud/EEG/eeglab2019_1/');
addpath('~/ownCloud/EEG/Michael_functions/');
%dataFolder='/run/media/ne550/TOSHBIG/EEG/Mb1extend/processed/';
% if isfolder(dataFolder)
% %     cd(dataFolder);
% end
newD=true;

eeglab;
%%
meanMethodArr={'gauss','mean','smooth'};
% meanMethod='gauss';
secondAlign=1;
plotVis=true;
menuSelect=false;
% dataFolder='/run/media/ne550/xchange/fin/EEG_Data_1/';
[qrsChoice,dataFolder]=uigetfile('*.vhdr','MultiSelect','on');
cd(dataFolder);
% if isfolder('Processed')
% saveFolder=[dataFolder 'Processed'];
% else
%     saveFolder=uigetdir(dataFolder);
% end
noMenu=1;
for meanm=1:length(meanMethodArr)
    meanMethod=meanMethodArr{meanm}
    if ~isfolder(meanMethod)
        mkdir([dataFolder meanMethod]);
    end
    saveFolder=[dataFolder meanMethod];
    % end
    for ji=1:length(qrsChoice)
        if (exist('HEP','var')),clear HEP;end
        if (exist('eegF','var')), clear eggF;end;
        EEG = pop_loadbv(dataFolder,qrsChoice{ji});
        EEG = eeg_checkset( EEG );
        EEG.setname=qrsChoice{ji}(1:end-5)
        EEG=pop_chanedit(EEG, 'lookup','/home/ne550/ownCloud/EEG/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
        EEG = eeg_checkset( EEG );
        eeglab redraw;
        if EEG.srate<5000
            [EEG] = pop_resample( EEG, 5000);
        end
        a=[EEG.event(strcmp({EEG.event.type},'R128')).latency];
        subMa=round(mean(diff(a(1:10))));
        subBreak=find(diff(a)>2*subMa);
        if subBreak
            subEnd=[a(subBreak)+2*subMa];
            subBegin=[a(1)-subMa a(subBreak+1)];
            subEnd=[subEnd a(end)+subMa];
            EEG.subsections=[subBegin' subEnd'];
        else
            EEG.subsections=[1 EEG.pnts];
        end
        eventNames={  'S  1'  'S  2'  'S  3'  'S  4'  'S  5' };
        EEG=eventeResponseFilter(EEG,eventNames,'S129',EEG.srate);
        eeg_eventtypes(EEG)
        qrsFile=dir([dataFolder qrsChoice{ji}(1:end-find(fliplr(qrsChoice{ji}=='.'))) '*.txt']);
        if logical(length(qrsFile))
            qrsFile(1).num=str2num(qrsFile(1).name(end-5));
            qrsFile(1).enable=true;
        else
            clear qrsFile
        end
        EEG.ECG=pop_select(EEG,'channel',{'ECG1','ECG2'});
        EEG.ECG.data(3,:)=(EEG.ECG.data(1,:)+mean(EEG.ECG.data(1,1:EEG.ECG.srate*10)))+(EEG.ECG.data(2,:)+mean(EEG.ECG.data(2,1:EEG.ECG.srate*10)));
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
        EEG=pop_select(EEG,'nochannel',{'ECG1','ECG2'});%%
        % if strcmp(qrsChoice(end-2:end),'set')
        % qrsFile=dir([dataFolder qrsChoice(1:end-4) '*.txt']);
        % else
        %  qrsFile=dir([dataFolder qrsChoice(1:end-find(fliplr(qrsChoice=='.'))) '*.txt']);
        % end
        % qrsFile=dir([dataFolder qrsChoice{ji}(1:end-find(fliplr(qrsChoice{ji}=='.'))) '*.txt']);
        % if logical(length(qrsFile))
        %     qrsFile(1).num=str2num(qrsFile(1).name(end-5));
        %     qrsFile(1).enable=true;
        % else
        %     clear qrsFile
        % end
        %% QRS deletion
        %EEG=eeg_addnewevents(EEG,{HEP.qrs'},{'qrs1'});
        debaseStart=EEG.srate/2;
        maxD=EEG.srate*0.05;
        cutTollerance=EEG.srate*0.02;
        trimPart=30;% percentage of trimmean
        fadeRange=maxD*2;
        minAvg=0.8;
        refW=30;
        eventType='qrs1';%backcompatibility overwritten in hep
        plotSample=false;%false;%'FCz';%true;+del2
        qrsLpf=30;
        qrsHpf=0;
        onlySample=0;
        
        if noMenu
            michael_hep;
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
        end
        
        if isgraphics(333);close(333);end
        if isgraphics(hepgui.main);close(hepgui.main);end
        %%
        saveQrsToFile=true;
        qrsA=HEP.qrs;
        if saveQrsToFile
            %     save([dataFolder EEG.setname 'noMRI' char(num2str(chanChoice)) '.txt'],'qrsA','-ascii')
            save([dataFolder qrsChoice{ji}(1:end-find(fliplr(qrsChoice{ji}=='.'))) '_' char(num2str(chanChoice)) '.txt'],'qrsA','-ascii')
        end
        % continue;
        eegT=eeg_addnewevents(EEG,{qrsA'},{char(eventType)});
        %eegT=eeg_addnewevents(EEG,{HEP.qrs'},{char(eventType)});
        %fmri_hr_rem(eegT,EEG.subsections,eventType,debaseStart,maxD,cutTollerance,trimPart,fadeRange,refW,plotSample,eegF,qrsHpf,onlySample,minAvg);
        %%
        tic
        if (~exist('eegF','var')) &&  preFilter
            eegF=pop_eegfiltnew(eegT,0,qrsLpf);
            %     EEG.ECG.event=EEG.event;EEG.ECG=pop_fmrib_qrsdetect(EEG.ECG,3,'qrs1','no');EEG.event=EEG.ECG.event
        end
        %currentData.noECG
        %=fmri_hr_rem(currentData.noMri,subsections,eventType,debaseStart,maxD,cutTollerance,trimPart,fadeRange,refW,plotSample,eegF,qrsLpf,onlySample,minAvg);
        %%
        % plotVis=false;
        % meanMethod='smooth';
        % secondAlign=1;
        fmri_hr_rem2
        %eegT =EEGc
        %%
        tic
        %qrsS.c=rmfield(qrsS.c,'data');
        %save([dataFolder qrsChoice(1:end-find(fliplr(qrsChoice=='.'))) 'noBCG.mat'],'qrsS')
        %movefile([dataFolder qrsChoice{ji}(1:end-find(fliplr(qrsChoice{ji}=='.'))) '*'],saveFolder)
        toc
        if plotVis
            menu('accept deletion',"OK");
        end
        pop_saveset( pop_resample(eegT, 500), 'filename',[qrsChoice{ji}(1:end-find(fliplr(qrsChoice{ji}=='.'))) '_' meanMethod '.set'],'filepath',saveFolder);%EEG = pop_saveset( EEG, 'filename',[EEG.setname 'noMRI.set'],'filepath','/run/media/ne550/TOSHBIG/EEG/');

        continue;
        %eegT =fmri_hr_rem(eegT,EEG.subsections,eventType,debaseStart,maxD,cutTollerance,trimPart,fadeRange,refW,plotSample,eegF,qrsLpf,onlySample,minAvg);
        %%
        % toc
        % EEG=pop_select(eegT,'point',EEG.subsections);
        % clear eegT;
        % EEG = eeg_checkset( EEG );
        % EEG= pop_resample(EEG, 250);
        % EEG = eeg_checkset( EEG );
        % %pop_saveset( EEG, 'filename',[EEG.setname
        % %'noQrsResamp.set'],'filepath','/run/media/ne550/TOSHBIG/EEG/');
        % EEG.noQRS=EEG;
        % %vis_artifacts(currentData.noMri,currentData.noECG);
        %vis_artifacts(EEG,eegT);
        if plotSample &  preFilter
            qrsAccept=menu('Accept Qrs deletion?','No','Yes');
            if qrsAccept==1
                keyboard;
            end
        end
        %%
        ch=60;
        gr=10;
        preLen=200;
        erpLen=800;
        chanQrs=arrayfun(@(y){cell2mat(arrayfun(@(x){[x.qrs(1:end-1)]'},qrsS.c(y,:)))},1:size(qrsS.c,1));
        chanDb=arrayfun(@(y){cell2mat(arrayfun(@(x){[x.dB]'},qrsS.c(y,:)))},1:size(qrsS.c,1));
        noQrsA=cell2mat(arrayfun(@(x){eegT.data(ch,x:(x-1)+erpLen)-mean(eegT.data(ch,x-preLen:x))},chanQrs{ch}'));
        figure;
        nexttile;hold on;
        %plot(noQrsA')
        arrayfun(@(y)plot((cell2mat(arrayfun(@(x){eegT.data(y,x:(x-1)+erpLen)-mean(eegT.data(y,x-preLen:x))},chanQrs{y}')))'),1:size(EEG.data,1));
        nexttile;hold on;
        arrayfun(@(y)plot((cell2mat(arrayfun(@(x,z){eegT.data(y,x:(x-1)+erpLen)-z},chanQrs{y}',chanDb{y}')))'),1:size(EEG.data,1));
        title('new')
        nexttile;hold on;
        arrayfun(@(y)plot((cell2mat(arrayfun(@(x){EEG.data(y,x:(x-1)+erpLen)-mean(EEG.data(y,x-preLen:x))},chanQrs{y}')))'),1:size(EEG.data,1));
        %plot((cell2mat(arrayfun(@(x){EEG.data(ch,x:(x-1)+erpLen)-mean(EEG.data(ch,x-preLen:x))},chanQrs{ch}')))')
        title('old')
        %%
        figure;
        nexttile;hold on;
        plot(cell2mat(arrayfun(@(x){eegT.data(ch,x:(x-1)+erpLen)-mean(eegT.data(ch,x-preLen:x))},qrsA)));
        title('new')
        nexttile;hold on;
        plot(cell2mat(arrayfun(@(x){EEG.data(ch,x:(x-1)+erpLen)-mean(EEG.data(ch,x-preLen:x))},qrsA)));
        title('old')
        nexttile;hold on;
        plot(mean(noQrsA))
        plot(std(noQrsA))
        nexttile;hold on;
        plot((noQrsA./std(noQrsA))')
        
        %%
        %EEG=eegT;
        %clear eegT;
        if isgraphics(111);close(111);end
        if isgraphics(112);close(112);end
        if isgraphics(hepgui.main);close(hepgui.main);end
        if loopProcess
            movefile([dataFolder qrsChoice(1:end-find(fliplr(qrsChoice=='.'))) '*'],saveFolder)
        else
            break;
        end
    end
end
% end
%%