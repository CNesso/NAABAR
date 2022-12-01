function eegO=fastMriM(eeg,eventArr,volOrSlice,exc_chan,options)
%{eeg,lpf,L,Win,etype,strig,anc_chk,trig_correct,Volumes,Slices,pre_frac,exc_chan,NPC}
arguments
    eeg
    eventArr
    volOrSlice
    exc_chan
    options.dummyExklude=6
    options.chanXworker=1
    options.oversample=20000/eeg.srate
    options.window=32
    options.anc=0
    options.trig_correct=0%search for missing trigger
    options.Volumes=410
    options.Slices=32
    options.lpf=69
    options.NPC='auto'%number of factors for OBS pca
    options.pre_frac=0 %position of artefact as fraction 0 artefact at beginning 1 at end
end
if iscell(eventArr)
    len=length(eventArr);
else
    len=1;
    eventArr={eventArr};
end
if options.dummyExklude>0
    for i=1:len
        eventArr{i}(1:options.dummyExklude)=[];
    end
end

eventArr=cellfun(@round,eventArr,'UniformOutput',false);

win=options.window;
chanXworker=options.chanXworker;
lpf=options.lpf;
anc=options.anc;
ChanGroups=ceil(size(eeg.data,1)/chanXworker);
oversample=options.oversample;
trig_correct=options.trig_correct;
Volumes=options.Volumes;
Slices=options.Slices;
NPC=options.NPC;
pre_frac=options.pre_frac;
if chanXworker==1
    pData=eeg.data;
    parfor p=1:ChanGroups
        [Lia,Locb]=ismember(exc_chan,p);
        Locb(Lia);
        size(Locb(Lia));
        evArr=eventArr;
        eegT=eegChanSel(eeg,p);
        for k=1:len
%             fprintf('\nProcessing Chan: %d - Subsection %d\n',p,k);
%             eegT=arrToEvent(eegT,evArr{k},'event');
%             eegT=pop_fmrib_fastr(eegT,lpf,oversample,win,'event',volOrSlice,anc,trig_correct,Slices,size(eegT.event,2),pre_frac,Locb(Lia),NPC);
%             eegT=fmrib_fastr(eegT,lpf,oversample,win,evArr{k},volOrSlice,anc,trig_correct,Slices,size(eegT.event,2),pre_frac,Locb(Lia),NPC);
            eegT=michael_fmrib_fastr(eegT,lpf,oversample,win,evArr{k},volOrSlice,anc,trig_correct,Slices,size(eegT.event,2),pre_frac,Locb(Lia),NPC,[k,p]);
        end
        pData(p,:)=eegT.data;
    end
    eegO=eeg;
    eegO.data=pData;
else
    pData=cell(1,ChanGroups);
    parfor p=1:ChanGroups
        i=(p-1)*chanXworker+1;
        if p==ChanGroups && rem(size(eeg.data,1),chanXworker)~=0
            WorkerRange=i:size(eeg.data,1)
        else
            WorkerRange=i:i+chanXworker-1
        end
        fprintf('\nProcessing Chan: %d\n',WorkerRange(1));
        [Lia,Locb]=ismember(exc_chan,WorkerRange);
        Locb(Lia)
        size(Locb(Lia));
        %     if len==1
        %         eegT=arrToEegEvent(eeg,WorkerRange,eventArr{1},etype);
        %         %eegT=pop_fmrib_fastr(eegT,lpf,4,win,eType,volOrSlice,anc,0,32,size(eegT.event,2),0,Loch(Lia),'auto');
        %         eegT=pop_fmrib_fastr(eegT,lpf,oversample,win,eType,volOrSlice,anc,trig_correct,Slices,size(eegT.event,2),pre_frac,Locb(Lia),'auto');
        %     else
        %         for k=1:len
        %         eegT=arrToEegEvent(eeg,WorkerRange,eventArr{i},etype);
        %         eegT=pop_fmrib_fastr(eegT,lpf,oversample,win,char(string(k)+eType),volOrSlice,anc,trig_correct,Slices,size(eegT.event,2),pre_frac,Locb(Lia),'auto');
        %         end
        %     end
        evArr=eventArr;
        eegT=eegChanSel(eeg,WorkerRange);
        for k=1:len
            eegT=arrToEvent(eegT,evArr{k},'event');
            eegT=pop_fmrib_fastr(eegT,lpf,oversample,win,'event',volOrSlice,anc,trig_correct,Slices,size(eegT.event,2),pre_frac,Locb(Lia),NPC);
        end
        pData(p)={eegT.data};
    end
    eegO=eeg;
    currChan=1;
    for i=1:size(pData,2)
        eegO.data(currChan:currChan+size(pData{i},1)-1,:)=pData{i};
        currChan=currChan+size(pData{i},1);
    end
end
end
