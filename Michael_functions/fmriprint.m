function [subsections,generatedSlices,correctedVolumes]=fmriprint(eeg,subsessionMinDistanceSamples,artfactRelPosSamples,artefactLengthSamples,eType,tr,dummyVolumes,slices,subsectionsPadding,addevent)
trSamples=tr*eeg.srate;
eegPart=[];
subsections=[];
correctedVolumes={};
generatedSlices={};
lastTrigger=0;
eegPartArtefact=[];
subsessionStart=0;
for i=1:length(eeg.event)
    if strcmp(eeg.event(:,i).type,eType)
        if ((eeg.event(:,i).latency-lastTrigger)>subsessionMinDistanceSamples)
            if subsessionStart~=0
                subsections(end+1,:)=[subsessionStart,lastTrigger+trSamples+subsectionsPadding];
                correctedVolumes(end+1)={eegPart};
                generatedSlices(end+1)={eegPartArtefact};
                eegPart=[];
                eegPartArtefact=[];
            end
            subsessionStart=eeg.event(:,i).latency-subsectionsPadding;
            eegPart(end+1:end+dummyVolumes)=(eeg.event(:,i).latency+artfactRelPosSamples)-dummyVolumes*trSamples:trSamples:(eeg.event(:,i).latency+artfactRelPosSamples-1);
        end
        for j=1:slices
            eegPartArtefact(end+1)=eeg.event(:,i).latency+(artefactLengthSamples*(j-1))+artfactRelPosSamples;
        end
        eegPart(end+1)=eeg.event(:,i).latency+artfactRelPosSamples;
        lastTrigger=eeg.event(:,i).latency;
    end
end
subsections(end+1,:)=[subsessionStart,lastTrigger+trSamples+subsectionsPadding];
correctedVolumes(end+1)={eegPart};
generatedSlices(end+1)={eegPartArtefact};
if addevent
    for i=1:size(correctedVolumes,2)
        EEG = eeg_addnewevents(EEG, {correctedVolumes{i}}, {char(string(i)+'correctedVol')});
        EEG = eeg_addnewevents(EEG, {generatedSlices{i}}, {char(string(i)+'slice')});
    end
end
end