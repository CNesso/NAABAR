function eegS=eegSubdivide(eeg,subsections)
rootSetName=eeg.setname;
for k=1:size(subsections,1)
    eegT=pop_select( eeg,'point',subsections(k,:));
    eegT=eeg_checkset(eegT);
    eegT.setname=[rootSetName '_subsection_' num2str(k)];
    fprintf('\nSubsession %d - Start: %3.3f s | End: %3.3f\n',k,subsections(k,1)/eeg.srate,subsections(k,2)/eeg.srate);
    eegS(k)=eegT;
end
end