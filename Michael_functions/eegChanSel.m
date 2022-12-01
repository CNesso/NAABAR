function eegOut=eegChanSel(eeg,chan)
    eegOut=eeg;
    eegOut.data=eeg.data(chan,:);
    eegOut.chanlocs=eeg.chanlocs(:,chan);
    eegOut.nbchan=length(chan);
end