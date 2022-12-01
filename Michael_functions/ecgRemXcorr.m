function varargout=ecgRemXcorr(eeg,qrsType,avgWindow,options)
arguments
    eeg
    qrsType
    avgWindow
    options.initialSkip=eeg.srate;
    % %%
    options.artStart=600
    options.artEnd=3550
    options.artFull=4000
    options.artPeakStart=1700%artStart
    options.artPeakEnd=2500%1200
    options.artPeakFactor=1
    options.maxDelay=200
    options.lpf=30
    options.plotFigures=flase
    options.eventTypeFigures='S 3'
    options.correctLine=true
end
initialSkip=options.initialSkip;
artStart=options.artStart;
artEnd=options.artEnd ;
artFull=options.artFull ;
artPeakStart=options.artPeakStart ;
artPeakEnd=options.artPeakEnd ;
artPeakFactor=options.artPeakFactor ;
maxDelay=options.maxDelay ;
lpf=options.lpf ;
plotFigures=options.plotFigures ;

refW=avgWindow;
eegM=eeg.data;

qrsA=[];
if (~exist('eegF','var'))
    eegF=pop_eegfiltnew(eeg,0,lpf);
end
refA=cell(size(eeg.data,1));% zeros(size(eeg.data,1),zeros(ceil(size(qrsA,2)/refW),artFull+1));
qrsCA=cell(size(eeg.data,1));
for i=1:length(eeg.event)
    if strcmp(eeg.event(:,i).type,qrsType)
        for j=1:size(subsections,1)
            if eeg.event(:,i).latency>subsections(j,1)+initialSkip & eeg.event(:,i).latency<subsections(j,2)
                qrsA(end+1)=eeg.event(:,i).latency;
                break;
            end
        end
        if eeg.event(:,i).latency+artFull>size(eeg.data,2)
            warning('Last QRS ignored. End of data preceedes end of artefact');
            break;
        end
    end
end
qrsA=unique(qrsA);
tStart = tic;
for par=1:size(eeg.data,1)
    tic
    eegChan=par
    parData=eeg.data(eegChan,:);
    lineMax=std(parData(subsections(1,1):subsections(1,2)));
    rebaseA=zeros(size(qrsA));
    qrsB=zeros(size(qrsA));
    delQ=zeros(size(qrsA));
    qrsD=zeros(size(qrsA,2),artFull+1);
    iter=1;
    refD=0;
    %
    
    for i=1:length(qrsA)
        
        sampleRange=qrsA(i)+artStart:qrsA(i)+artEnd;
        fullMean=mean(eegF.data(eegChan,qrsA(i):qrsA(i)+artFull));
        iterData=eegF.data(eegChan,sampleRange)-fullMean;
        if iter==1
            refD=iterData;
        end
        covQ=xcorr(iterData,refD)+(xcorr(iterData,refD(artPeakStart:artPeakEnd))*artPeakFactor);
        [~,delCo]=max(covQ(length(iterData)-maxDelay:length(iterData)+maxDelay));
        delQ(i)=delCo-maxDelay;
        qrsD(i,:)=eegF.data(eegChan,qrsA(i)+delQ(i):qrsA(i)+delQ(i)+artFull);
        qrsD(i,:)=qrsD(i,:)-mean(qrsD(i,:));
        rebaseA(i)=trimmean(eegF.data(eegChan,qrsA(i)-trSamples:qrsA(i)),33);
        iter=iter+1;
        if iter>refW
            iter=1;
        end
        qrsB(i)=qrsA(i)+delQ(i);
    end
    qrsC=zeros(size(qrsA));
    delQ2=zeros(size(qrsA));
    qrsD2=zeros(size(qrsA,2),artEnd+1);%artFull+1);
    iter=1;
    reff=zeros(ceil(size(qrsA,2)/refW),artEnd+1);%artFull+1);
    for i=1:length(qrsB)
        sampleRange=qrsB(i)+artStart:qrsB(i)+artEnd;
        fullMean=mean(eegF.data(eegChan,qrsB(i):qrsB(i)+artFull));
        iterData=eegF.data(eegChan,sampleRange)-fullMean;
        if iter>refW | i==1
            if i+refW-1>length(qrsB)
                refD=mean(smoothdata(qrsD(i-refW+1:i,:),'gaussian',refW));
            else
                refD=mean(smoothdata(qrsD(i:i+refW-1,:),'gaussian',refW));
            end
            iter=1;
            if i~=1
                reff(floor(i/refW),:)=mean(smoothdata(qrsD2(i-refW:i-1,:),'gaussian',refW));
            end
        end
        covQ=xcorr(iterData,refD(artStart:artEnd))+(xcorr(iterData,refD(artPeakStart:artPeakEnd))*artPeakFactor);
        [~,delCo]=max(covQ(length(iterData)-maxDelay:length(iterData)+maxDelay));
        delQ2(i)=delCo-maxDelay;
        qrsD2(i,:)=eegF.data(eegChan,qrsB(i)+delQ2(i):qrsB(i)+delQ2(i)+artEnd);%artFull);
        qrsD2(i,:)=qrsD2(i,:)-rebaseA(i);
        iter=iter+1;
        qrsC(i)=qrsB(i)+delQ2(i);
    end
    reff(end,:)=mean(smoothdata(qrsD2(end-refW+1:end,:),'gaussian',refW));
    refA{par}=reff;
    qrsCA{par}=qrsC;
    
    lineDiffA=zeros(size(qrsA));
    paraArr=zeros(size(qrsA,2),2);
    pData=cell(size(qrsA,2),1);
    iter=1;
    reffF=cell(ceil(size(qrsA,2)/refW),1);%artFull+1);
    reffLen=zeros(ceil(size(qrsA,2)/refW));
    for k=1:length(qrsC)
        if iter>refW
            reffLen(floor(k/refW))=min((qrsC(k+1-refW:k)-1)-qrsA(k-refW:k-1));
            ppData=cell2mat(pData(k-refW:k-1));
            reffF{floor(k/refW)}=mean(smoothdata(ppData(1:reffLen(floor(k/refW))),'gaussian',refW));
            iter=1;
        end
        pData{k}=parData(qrsA(k):qrsA(k)+trSamples*2)-rebaseA(k);
        if options.correctLine
            para=parData(qrsC(k):qrsC(k)+artEnd)-reff(ceil(k/refW),:);
            paraA=[parData(qrsC(k))-para(1),parData(qrsC(k)+artEnd)-para(end)];
            lineDiff=abs(paraA(1)-paraA(2));
            lineDiffA(k)=paraA(1)-paraA(2);
            paraArr(k,:)=paraA;
            if lineDiff<lineMax
                if (k>1) && (((qrsC(k)-1)-(qrsC(k-1)+artEnd+1))>0)
                    parData(qrsC(k):qrsC(k)+artEnd)=parData(qrsC(k):qrsC(k)+artEnd)-reff(ceil(k/refW),:);
                    lineIni=parData(qrsC(k-1)+artEnd)+(eeg.data(eegChan,(qrsC(k-1)+artEnd+1))-eeg.data(eegChan,(qrsC(k-1)+artEnd)));
                    lineEnd=parData(qrsC(k))-(eeg.data(eegChan,(qrsC(k)+1))-eeg.data(eegChan,(qrsC(k))));
                    paraA=[parData(qrsC(k-1)+artEnd+1)-lineIni,parData(qrsC(k)-1)-lineEnd];
                    lineDiff=abs(paraA(1)-paraA(2));
                    
                    [~,lineI]=min(paraA);
                    if lineI==2
                        lineDiff=-lineDiff;
                    end
                    lineAdd=paraA(1):lineDiff/((qrsC(k)-1)-(qrsC(k-1)+artEnd+1)):paraA(2);
                    parData(qrsC(k-1)+artEnd+1:qrsC(k)-1)=parData(qrsC(k-1)+artEnd+1:qrsC(k)-1)-lineAdd;
                else
                    %qrsC(k)*eeg.srate
                    paraAA=[qrsC(k),qrsC(k)+artEnd];
                    [~,paraI]=min(abs(paraA));
                    parData(qrsC(k):qrsC(k)+artEnd)=parData(qrsC(k):qrsC(k)+artEnd)-reff(ceil(k/refW),:)+paraA(paraI);
                end
            else
                %             qrsC(k)*eeg.srate
                paraAA=[qrsC(k),qrsC(k)+artEnd];
                [~,paraI]=min(abs(paraA));
                if paraI==1
                    paraM=2;
                else
                    paraM=1;
                end
                parData(qrsC(k):qrsC(k)+artEnd)=parData(qrsC(k):qrsC(k)+artEnd)-reff(ceil(k/refW),:)+paraA(paraI);
                parData(paraAA(paraM)-60:paraAA(paraM)+60)=smoothdata(parData(paraAA(paraM)-60:paraAA(paraM)+60),'rloess');
            end
        else
            parData(qrsC(k):qrsC(k)+artEnd)=parData(qrsC(k):qrsC(k)+artEnd)-reff(ceil(k/refW),:);
        end
        iter=iter+1;
    end
    eegM(par,:)= parData;
    toc
end
toc(tStart)
eegTT=eeg
eegTT.data=eegM
if plotFigures
    EEGBrowser(eegTT)
    eegChan=randi(64)
    reff=refA{eegChan};
    qrsC=qrsCA{eegChan};
    figure;
    nexttile;
    plot(delQ);hold on;plot(delQ2)
    numFig=refW;
    figStart=randi(length(qrsC)-numFig)
    nexttile;
    plot(qrsD2(figStart:figStart+numFig,:)');
    disp('plot from '+ string(round(qrsC(figStart)/eeg.srate))+ 's')
    for i=figStart:figStart+numFig
        nexttile;
        hold on;
        plot(eegTT.data(eegChan,qrsC(i):qrsC(i)+artFull)-mean(eegTT.data(eegChan,qrsC(i):qrsC(i)+artFull)));
        hold on;
        plot(reff(ceil(i/20),:),'--');
        plot(eeg.data(eegChan,qrsC(i):qrsC(i)+artFull)-mean(eeg.data(eegChan,qrsC(i):qrsC(i)+artFull)),'--');
    end
    %%
    %eegFF=eegTT;
    eegFF=pop_eegfiltnew(eegTT,0,lpf);
    menn=zeros(length(qrsC),9001);
    for i=1:length(qrsC)
        menn(i,:)=eegFF.data(62,qrsA(i)-1000:qrsA(i)+artFull*2)-trimmean(eegFF.data(62,qrsA(i)-5000:qrsA(i)),33);
    end
    figure;
    nexttile;
    hold on;
    plot(mean(menn(1:20,:)))
    plot(mean(menn))
    
    menn=zeros(length(qrsC),8001);
    for i=1:length(qrsC)
        menn(i,:)=eegFF.data(62,qrsC(i):qrsC(i)+artFull*2)-trimmean(eegFF.data(62,qrsC(i)-5000:qrsC(i)),33);
    end
    plot(mean(menn(1:20,:)))
    plot(mean(menn))
    % eeg = eeg_addnewevents(eeg, {qrsCA{62}}, {'qrsC64'});
    s3A=[];
    for i=1:length(eeg.event)
        if strcmp(eeg.event(:,i).type,options.eventTypeFigures)
            s3A(end+1)=eeg.event(:,i).latency;
        end
    end
    s3B=zeros(length(s3A),4001);
    for i=1:length(s3A)
        menn(i,:)=eegTT.data(62,s3A(i):s3A(i)+10000)-mean(eegTT.data(62,s3A(i)-5000:s3A(i)));
    end
    nexttile;
    plot(mean(menn))
    nexttile;
    plot(menn(1:10,:)')
end

 varargout{1}=eegTT;
end
%%