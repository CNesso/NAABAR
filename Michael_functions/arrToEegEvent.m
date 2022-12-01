function eegOut=arrToEegEvent(eeg,chan,eventArr,EventName)
if iscell(eventArr)
    len=length(eventArr);
else
    len=1;
    eventArr={eventArr};
end
for i=1:len
    evArr=eventArr{i};
    T = table(evArr',ones(length(evArr),1),...
        zeros(length(evArr),1),ones(length(evArr),1),...
        [1:length(evArr)]',repmat(EventName,length(evArr),1),...
        repmat('Response',length(evArr),1),...
        [1:length(evArr)]','VariableNames',fieldnames(eeg.event)');
    eegOut(i)=eeg;
    eegOut(i).data=eeg.data(chan,:);
    eegOut(i).chanlocs=eeg.chanlocs(:,chan);
    eegOut(i).nbchan=length(chan);
    eegOut(i).event=table2struct(T);
end
end
