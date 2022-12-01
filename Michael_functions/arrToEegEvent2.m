function eegOut=arrToEegEvent2(eeg,chan,eventArr,EventName)
if iscell(eventArr)
    len=length(eventArr);
    if iscell(EventName)
        if length(eventArr)>length(EventName)
            EventName=arrayfun(@(x){sprintf('%s %d',x,EventName{1})},1:length(eventArr));
        end
    else
        EventName=arrayfun(@(x){sprintf('%s %d',x,EventName)},1:length(eventArr));
        %EventName=repmat({EventName},1,length(eventArr));
    end
else
    len=1;
    eventArr={eventArr};
    if ~iscell(EventName)
        EventName={EventName}
    end
end
for i=1:len
    evArr=eventArr{i};
    evNam=EventName(i);
    if exist('T','var')
        Told=T;
    end
    T= table(evArr',ones(length(evArr),1),...
        zeros(length(evArr),1),ones(length(evArr),1),...
        [1:length(evArr)]',repmat(evNam,length(evArr),1),...
        repmat({'Response'},length(evArr),1),...
        [1:length(evArr)]','VariableNames',fieldnames(eeg.event)');
    if exist('Told','var')
        T=[T;Told];
    end
end
    eegOut=eeg;
    eegOut.data=eeg.data(chan,:);
    eegOut.chanlocs=eeg.chanlocs(:,chan);
    eegOut.nbchan=length(chan);
    eegOut.event=table2struct(T);
end
