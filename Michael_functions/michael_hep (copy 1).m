%clear HEP
if (~exist('HEP','var'))
    qrsOn=false;
    chanChoice=1;
    if(exist('qrsFile','var'))
        if qrsFile(1).enable
            qrsOn=true;
        end
    end
    if qrsOn
        HEP.qrs =load([dataFolder qrsFile(1).name]);
        if ~isempty(str2num(qrsFile(1).name(end-5)))
            chanChoice=str2num(qrsFile(1).name(end-5));
        end
    else
        %EEGBrowser( EEG.ECG);
        qrsAA=[];
        figure(333);
        for ii=1:3
            qrsA=[];
            for i=1:length(EEG.ECG.event)
                if strcmp(EEG.ECG.event(:,i).type,char("qrs"+num2str(ii))) ...
                        qrsA(end+1)=EEG.ECG.event(:,i).latency;
                end
            end
            qrsA=unique(qrsA);
            qrsAA{ii}=qrsA;
            nexttile;plot(diff(qrsA));title(num2str(ii));
        end
        chanChoice = menu('Choose a color','1:ECG1','2:ECG2','3:Difference');
        eventType=char("qrs"+num2str(chanChoice));
        qrsA=qrsAA{chanChoice};
        HEP.qrs =qrsA';
    end
    subsectionsI=~logical(EEG.times);
    % subsectionsI(1,:)=false;
    subsectionsI(cell2mat(arrayfun(@(x){[EEG.subsections(x,1):EEG.subsections(x,2)]},1:size(EEG.subsections,1))))=true;
    
    %     qrsA=[];
    %     for i=1:length(EEG.ECG.event)
    %         if strcmp(EEG.ECG.event(:,i).type,eventType) ...
    %             qrsA(end+1)=EEG.ECG.event(:,i).latency;
    %         end
    %     end
    %     qrsA=unique(qrsA);
    %clear HEP
    
    ecg = double(EEG.ECG.data(chanChoice,:));
    ecg = ecg(:); % make sure ecg is one-column vector
    
    HEP.ecg = (ecg-min(ecg))/(max(ecg)-min(ecg));
    
    %HEP.ecg = ecg/max(abs(ecg)); % normalise ecg
    
    HEP.srate = EEG.srate;
    HEP.winsec = 10; %janela de tempo mostrada no grafico
    HEP.sec_ini = 0; %posicao inicial do grafico
    HEP.ignore=[];
    %     HEP.qrs =qrsA';
    
    clear ecg qrsOn
    
end
heplab
set(hepgui.save_btn,'CallBack','delete(gcf)');
set(hepgui.main, 'WindowKeyPressFcn',[...
    'subsectionsI(HEP.ignore)=false;',...
    'qrsDiffAvg=abs(smoothdata(diff(HEP.qrs),"gaussian",length(HEP.qrs)/20)-diff(HEP.qrs));',...
    'qrsDiffAvg=[0;abs(mean(qrsDiffAvg)-qrsDiffAvg)];',...
    'mapI=1:length(qrsDiffAvg);',...
    'mapI=mapI(subsectionsI(HEP.qrs));',...
    '[~,maxQrsDiff]=max(qrsDiffAvg(subsectionsI(HEP.qrs)));',...
    'maxQrsDiffI=mapI(maxQrsDiff);',...
    'maxQrsDiff=HEP.qrs(maxQrsDiffI)/EEG.srate;',...
    'A=maxQrsDiff-HEP.winsec*0.5;',...
    'if abs(str2num(get(hepgui.win_label,''String''))-A)<0.1,',...
    'HEP.ignore(end+1)=HEP.qrs(maxQrsDiffI);',...
    'qrsDiffAvg(maxQrsDiffI)=false;',...
    '[~,maxQrsDiff]=max(qrsDiffAvg(subsectionsI(HEP.qrs)));',...
    'maxQrsDiff=HEP.qrs(mapI(maxQrsDiff))/EEG.srate;',...
    'A=maxQrsDiff-HEP.winsec*0.5;',...
    'end,',...
    'if length(A)==1,',...
    'if A(1) <= get(hepgui.slider_arrow,''Max'') & A(1)>=0,',...
    'set(hepgui.slider_arrow,''Value'',A(1));',...
    'elseif A(1)<=get(hepgui.slider_arrow,''Max'')+HEP.winsec & A(1)>=0,',...
    'set(hepgui.slider_arrow,''Value'',get(hepgui.slider_arrow,''Max''));',...
    'end,',...
    'HEP.sec_ini=round(1000*get(hepgui.slider_arrow,''Value''))/1000;',...
    'HEP.ecg_handle=heplab_ecgplot(HEP.ecg,HEP.srate,HEP.qrs,',...
    'HEP.sec_ini,HEP.ecg_handle,HEP.winsec);',...
    'end,',...
    'set(hepgui.win_label, ''String'',num2str(HEP.sec_ini));',...
    'axes(HEP.ecg_handle);',...
    'hold on;xline(maxQrsDiff);',...
    ]);
% qrsDiffAvg=abs(smoothdata(diff(HEP.qrs),"gaussian",length(HEP.qrs)/20)-diff(HEP.qrs));
% qrsDiffAvg=[0;abs(mean(qrsDiffAvg)-qrsDiffAvg)];
% mapI=1:length(qrsDiffAvg);
% mapI=mapI(subsectionsI(HEP.qrs));
% qrsDiffAvg=qrsDiffAvg(subsectionsI(HEP.qrs));
% [~,maxQrsDiff]=max(qrsDiffAvg);
% maxQrsDiff=HEP.qrs(mapI(maxQrsDiff))/EEG.srate
%
% stdI=qrsDiffAvg>std(qrsDiffAvg)*3;
% mapI=mapI(stdI);
% [~,sortQrsDiff]=sort(qrsDiffAvg(stdI),'descend');
% sortQrsDiff=HEP.qrs(mapI(sortQrsDiff))/EEG.srate
