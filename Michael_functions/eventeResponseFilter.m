function varargout=eventeResponseFilter(eeg,stimuliTypes,responseType,timeWindow)
% eventNames={  'S  1'  'S  2'  'S  3'  'S  4'  'S  5' };% 'S129'  'qrs1'  };
% eegt=eventeResponseFilter(EEG,eventNames,'S129',EEG.srate);
% [eegt,eve,meve]=eventeResponseFilter(EEG,eventNames,'S129',EEG.srate);
if nargout>1
for i=[stimuliTypes,responseType]
   eventStruct.(regexprep(i{1},'\W',''))=[];
   if nargout>2
       missStruct.(regexprep(i{1},'\W',''))=[];
   end
end
end
currentStimulus={'',0,0};%type,eventnumber,latency
for i=1:length(eeg.event)
    if strcmp(eeg.event(:,i).type,responseType)
        if nargout>1
            tArr=eventStruct.(responseType);
            tArr(end+1)=eeg.event(:,i).latency;
            eventStruct.(responseType)=tArr;
        end
        if (eeg.event(:,i).latency-currentStimulus{3})<timeWindow
            if nargout>1
                tArr=eventStruct.(currentStimulus{1});
                tArr(end+1)=eeg.event(:,i).latency;
                eventStruct.(currentStimulus{1})=tArr;
            end
            currentStimulus={'',0,0};
        else
            if currentStimulus{2}==0
                continue
            end
            eeg.event(:,currentStimulus{2}).type=strjoin({eeg.event(:,currentStimulus{2}).type, 'MISS'});
            if nargout>2
                tArr=missStruct.(currentStimulus{1});
                tArr(end+1)=eeg.event(:,i).latency;
                missStruct.(currentStimulus{1})=tArr;
            end
            currentStimulus={'',0,0};
        end
    elseif ismember(eeg.event(:,i).type,stimuliTypes)
        if currentStimulus{2}~=0
            eeg.event(:,currentStimulus{2}).type=strjoin({eeg.event(:,currentStimulus{2}).type, 'MISS'});
            if nargout>2
                tArr=missStruct.(currentStimulus{1});
                tArr(end+1)=eeg.event(:,i).latency;
                missStruct.(currentStimulus{1})=tArr;
            end
        end
        currentStimulus={regexprep(eeg.event(:,i).type,'\W',''),i,eeg.event(:,i).latency};
    else
        if nargout>1
            fName=regexprep(eeg.event(:,i).type,'\W','');
            if ~isfield(eventStruct,fName);
                eventStruct.(fName)=[];
            end
            tArr=eventStruct.(fName);
            tArr(end+1)=eeg.event(:,i).latency;
            eventStruct.(fName)=tArr;
        end
    end
end
        varargout{1}=eeg;
switch nargout
    case 2
        varargout{2}=eventStruct;
    case 3
        varargout{2}=eventStruct;
        varargout{3}=missStruct;
end
end