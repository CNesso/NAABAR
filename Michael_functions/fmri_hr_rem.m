%% inportant do only one section at a time!!!!
% subsections=[611345,2662745;3107835,5159235];
% debaseStart=EEG.srate/2;
% maxD=EEG.srate*0.05;
% cutTollerance=EEG.srate*0.02;
% trimPart=30;% percentage of trimmean
% fadeRange=maxD*2;
% %end of func arguments
% refW=30;
% eventType='qrs1';
% plotSample='FCz';
% qrsLpf=30;
% qrsHpf=0;
% onlySample=false;
%fmri_hr_rem(EEG,subsections,eventType,debaseStart,maxD,cutTollerance,trimPart,fadeRange,refW,plotSample,qrsLpf,qrsHpf,onlySample)
function eegO=fmri_hr_rem(EEG,subsections,eventType,debaseStart,maxD,cutTollerance,trimPart,fadeRange,refW,plotSample,qrsLpf,qrsHpf,onlySample,minAvg)
%important do only one section at a time!!!!
try
if ~isnumeric(qrsLpf)
    eegF=qrsLpf;
else
    if qrsHpf
        EEG=pop_eegfiltnew(EEG,1);%suggested not necessary
    end
    if qrsLpf
        eegF=pop_eegfiltnew(EEG,0,qrsLpf);
    else
        eegF=EEG;
    end
end
%%
wholeChan=false;
if plotSample
        wholeChan=true;
        if ischar(plotSample)
            chanNum=find(strcmp(plotSample,{EEG.chanlocs.labels}));
        elseif isnumeric(plotSample)
            chanNum=plotSample;
        else
            wholeChan=false;
            chanNum=randi(size(EEG.data,1))
        end
    chandata=eegF.data(chanNum,:);
    chandataN=EEG.data(chanNum,:);
    if isgraphics(111);clf(111);end;figure(111);
    for i=1:size(subsections,1)
        sectionStart=subsections(i,1);
        sectionEnd=subsections(i,2);
        [diffCell,modNotZero,winRem]=qrsEvent2array(EEG,eventType,refW,sectionStart,sectionEnd,true);
    end
    sectionStart=subsections(1,1);
    sectionEnd=subsections(1,2);
    [diffCell,modNotZero,winRem]=qrsEvent2array(EEG,eventType,refW,sectionStart,sectionEnd,false);
    
    cellNum=randi(length(diffCell))%9%
    
    
    qrs=diffCell{cellNum};
    if isgraphics(112);clf(112);end;figure(112);
    artRange=windowAvg(chandata,debaseStart,qrs,maxD*2,trimPart,EEG.chanlocs(chanNum).labels);
    winSub(chandataN,artRange,fadeRange,minAvg,sectionStart,EEG.chanlocs(chanNum).labels);
    winData=cellfun(@(x)winSub(chandataN,windowAvg(chandata,debaseStart,x,maxD*2,trimPart,0),fadeRange,minAvg,sectionStart,0),diffCell);
    drawnow;
    %figure;plot(cell2mat(winData.data'));
end
%%
if ~onlySample
    eegO=EEG;
    completeData=cell(1,size(subsections,1));
    for i=1:size(subsections,1)
        fprintf('\n\n Section %d\n',i);

        sectionEEG=pop_select(EEG,'point',subsections(i,:));
        sectionEegF=pop_select(eegF,'point',subsections(i,:));
        sectionStart=1;%subsections(i,1);
        sectionEnd=sectionEEG.pnts;%subsections(i,2);
        
        [diffCell,modNotZero,winRem]=qrsEvent2array(sectionEEG,eventType,refW,sectionStart,sectionEnd,false);
        
        %completeData{i}=arrayfun(@(x) qrsChannel(x,EEG,eegF,diffCell,sectionStart,sectionEnd,modNotZero,winRem,debaseStart,maxD,trimPart,fadeRange,cutTollerance,minAvg),(1:size(EEG.data,1)),'UniformOutput',false);
        sectionData=arrayfun(@(x) qrsChannel(x,sectionEEG,sectionEegF,diffCell,sectionStart,sectionEnd,modNotZero,winRem,debaseStart,maxD,trimPart,fadeRange,cutTollerance,minAvg),(1:size(EEG.data,1)),'UniformOutput',false);
        
        
        sectionStart=subsections(i,1);
        sectionEnd=subsections(i,2);
        %eegO.data(:,sectionStart:sectionEnd)=cell2mat(completeData{i}');
        eegO.data(:,subsections(i,1):subsections(i,2))=cell2mat(sectionData');
 
        fprintf(' Finisched\n');
    end
    
    if plotSample
%         EEGBrowser(eegO);
        vis_artifacts(eegO,EEG);
    end
else
    if wholeChan
    eegO=winData;
    else
       eegO=true;
    end
end
catch ME
    disp(getReport(ME));
    keyboard;
end
end
%%

function completeChan=qrsChannel(chanNum,EEG,eegF,diffCell,sectionStart,sectionEnd,modNotZero,winRem,debaseStart,maxD,trimPart,fadeLength,cutTollerance,minAvg)
try
    fadeF=arrayfun(@(x) 1/((1+exp(-(15/fadeLength*(-0.5*fadeLength+x))))),1:fadeLength);
    %setCells=cell(size(EEG.data,1),1);
    %completeData=cell(size(EEG.data,1),1);
   msg = sprintf(['Channel: %03d'],chanNum);
   reverseStr = repmat(sprintf('\b'), 1, length(msg));
   if chanNum>1
   fprintf([reverseStr msg]);
   else
     fprintf(msg);
   end
   
   
    
%     fprintf('%d.',chanNum);
    chandata=eegF.data(chanNum,:);
    chandataN=EEG.data(chanNum,:);
    chanCells=cellfun(@(x) winSub(chandataN,windowAvg(chandata,debaseStart,x,maxD,trimPart,false),fadeLength,minAvg,sectionStart,false),diffCell,'UniformOutput',false);
    if modNotZero
        lastCell=chanCells{end};
        %chanCells{end}=lastCell(end-winRem+1:end,:);
        chanCells{end}.data=lastCell.data(end-winRem+1:end,:);
        chanCells{end}.qrs=lastCell.qrs(end-winRem+1:end,:);
    end
    chanBegin=chandataN(sectionStart:chanCells{1}.qrs(1));%chanCells{1}{1}(1));
    %chanTcell=[chanCells(2:end);{{sectionEnd+1}}];
    %continue here
secondCell=cellfun(@(x){x.qrs(1),x.before},chanCells(2:end),'UniformOutput',false);

%diffLastCell=cellfun(@(x,y)(y{1}-x.qrs(end))-1,chanCells(1:end-1),secondCell,'UniformOutput',false);
diffLastCell=cellfun(@(x,y)(y{1}-x.qrs(end)),chanCells(1:end-1),secondCell,'UniformOutput',false);

beforeA=cellfun(@(x,y,z) x{2}(end-(min([max([y-length(z.data{end}),0]), length(x{2})-fadeLength])+(fadeLength-1)):end),secondCell,diffLastCell,chanCells(1:end-1),'UniformOutput',false);
beforeA=cellfun(@(x) [x(1:fadeLength).*fadeF x(fadeLength+1:end)],beforeA,'UniformOutput',false);
artA=cellfun(@(x,y,z) x.data{end}(1:min([length(x.data{end}),y])),chanCells(1:end-1),diffLastCell,'UniformOutput',false);
artA=cellfun(@(x) [x(1:end-fadeLength)  x(end-fadeLength+1:end).*(1-fadeF)],artA,'UniformOutput',false);


last=chanCells{end}.data{end}(1:min([min(diff(chanCells{end}.qrs)),length(chandataN)-chanCells{end}.qrs(end),length(chanCells{end}.data{end})])-1);
last=[last(1:(end-fadeLength)-1),last(end-(fadeLength-1):end).*fadeF];

%lasting=[chandataN(chanCells{end}.qrs(end):chanCells{end}.qrs(end)+(length(last)-1))-last,chandataN(chanCells{end}.qrs(end)+length(last)+1:end)];
lasting=chandataN(chanCells{end}.qrs(end):chanCells{end}.qrs(end)+(length(last)-1))-last;

%ending=chandataN(chanCells{end}.qrs(end)+length(last)+1:sectionEnd);
% ending=chandataN(chanCells{end}.qrs(end)+length(last)+1:sectionEnd+1);
ending=chandataN(chanCells{end}.qrs(end)+length(last):sectionEnd);


fadeLength=min([fadeLength chanCells{1}.qrs(1)-sectionStart]);
fadeF=arrayfun(@(x) 1/((1+exp(-(15/fadeLength*(-0.5*fadeLength+x))))),1:fadeLength);
beforeFirst=chanCells{1}.before(end-(fadeLength-1):end).*fadeF;
beginning=[chandataN(sectionStart:(chanCells{1}.qrs (1)-fadeLength)-1), chandataN(chanCells{1}. qrs (1)-fadeLength:chanCells{1}.qrs (1)-1)-beforeFirst];

compose=cellfun(@(x,y,z) [x zeros(1, z-length(x))]+[zeros(1, z-length(y)) y],artA,beforeA,diffLastCell,'UniformOutput',false);
compose(end+1)={lasting};

completeArt=cellfun(@(x,y) [x.data(1:end-1);{y}], chanCells, compose,'UniformOutput',false);

winData=cellfun(@(a,b)arrayfun(@(x,y)chandataN(x:x+(length(y{1})-1))-y{1},a.qrs,b,'UniformOutput',false),chanCells,completeArt,'UniformOutput',false);

%completeChan=[beginning,cell2mat([cellfun(@(x,y) cell2mat([x.data(1:end-1);{y}]'), chanCells, compose,'UniformOutput',false)]'),ending];
completeChan=[beginning,cell2mat([cellfun(@(x) cell2mat(x'),winData,'UniformOutput',false)]'),ending];
% %     chanCells{2}.qrs(1)-(chanCells{1}.qrs(end)+length(chanCells{1}.data{end}))
% %     chanCells{2}.before((end-(max([chanCells{2}.qrs(1)-(chanCells{1}.qrs(end)+length(chanCells{1}.data{end})) 0])+fadeLength)+1):end)
%     keyboard
%     intervals=arrayfun(@(x,y)chandataN(x{end}{end,1}(2):y{1}{1}(1)-1),chanCells,chanTcell,'UniformOutput',false);
%     for j=1:length(intervals)
%         if length(intervals{j})<1
%             chanCells{j}(end,2)={chanCells{j}{end,2}(1:chanCells{j+1}{1}(1)-chanCells{j}{end,1}(1))};
%         end
%     end
%     completeChan=cellfun(@(x,y)[x(:,2);{y}],chanCells,intervals,'UniformOutput',false);
%     completeChan=[chanBegin(1:end-1) cell2mat(cellfun(@(x)cell2mat(x'),completeChan,'UniformOutput',false)')];
%     
%     %         completeData{i}=completeChan;
%     %setCells{i}=chanCells;
%     
    
catch ME
    disp(getReport(ME));
    keyboard;
end

end

function outData=windowAvg(chandata,debaseStart,qrs,maxD,trimPart,plotData)
try
    srate=5000;
    %qrsStd=std(qrs)[1:beforeArtLenght
    artEnd=min(diff(qrs));
    qBase=arrayfun(@(x) trimmean(chandata(x-debaseStart:x),trimPart),qrs);
    winData=arrayfun(@(x)chandata(x: x+artEnd),qrs,'UniformOutput',false);
    winDataO=winData;
    
    [~,del1]=cellfun(@(x) max(xcov(x,winData{1}, maxD)), winData(2:end));
    del1=[0 ; del1-maxD];
    
    winData=cell2mat(arrayfun(@(x,y)chandata(x+y: x+y+artEnd),qrs,del1,'UniformOutput',false));
    winData1=winData-qBase;
    %art1=mean(smoothdata(winData1,'gaussian',length(qrs)));
    art1=gaussM(winData1);
    var1=var(winData1);
        qBaseC=arrayfun(@(x){sum((winData1(x,:)-art1).*(1./var1))/sum(1./var1)},1:size(winData1,1));
    %[~,del2]=cellfun(@(x) max(xcov(x-mean(x-art1,2),art1,maxD)), winDataO);
    [~,del2]=cellfun(@(x,y) max(xcov(x-y,art1,maxD)), winDataO,qBaseC');
    del2=del2-maxD;
        qrsCorrected=qrs+del2;
    outData.artDiff=diff(qrsCorrected);
    artEnd=min(outData.artDiff);

    winDataN=arrayfun(@(x)chandata(x: x+artEnd),qrsCorrected,'UniformOutput',false);
    winData=cell2mat(winDataN)-qBase;

    %art2=mean(smoothdata(winData,'gaussian',length(qrs)));
    art2=gaussM(winData);
    var2=var(winData);
%     [~,del3]=cellfun(@(x) max(xcov(x-mean(x-art2,2),art2,maxD)), winDataN);
%     del3=del3-maxD
%         qrsCorrected=qrsCorrected+del3;
%     outData.artDiff=diff(qrsCorrected);
%     artEnd=min(outData.artDiff);
% 
%     winDataN3=arrayfun(@(x)chandata(x: x+artEnd),qrsCorrected,'UniformOutput',false);
%     winData3=cell2mat(winDataN3)-qBase;
% 
%     art3=mean(smoothdata(winData3,'gaussian',length(qrs)));
    %art3=mean(smoothdata(winData-cellfun(@(x) trimmean(x-art2,trimPart,2),winDataN),'gaussian',length(qrs)));
    if plotData
        fig=gcf;
        fig.Name=['ECG Plots sample:'  char(string(plotData))  ' [' char(string(round(qrs(1)))) '-' char(string(round(qrs(end)))) '(samples)]'];
%         nexttile;hold on;cellfun(@plot, winDataO);title("Raw");
        nexttile;hold on;plot((1:length(winDataO{1}))/srate, (cell2mat(winDataO)-qBase)');title("Raw");ylabel('\muV');
        nexttile;hold on;plot((1:length(winData1))/srate,winData1');title("Delay(1)");
        nexttile;plot((1:length(winData1))/srate,(winData1-mean(winData1-art1,2))');title("Delay(1)&Art.debase(1)");
        
        nexttile;hold on;plot((1:length(winData))/srate,(winData'));title("Delay(2,LFP)");ylabel('\muV');
        nexttile;hold on;plot((1:length(winData1))/srate,art1);plot((1:length(winData))/srate,art2);yyaxis right;plot((1:min([length(art2) length(art1)]))/srate,art2(1:min([length(art2) length(art1)]))-art1(1:min([length(art2) length(art1)])),'Color',[0.47,0.67,0.19]);legend("1st aproximation","2nd aproximation","Difference("+string(round(mean(abs(del2))))+"samp.del.)");title("Aproximations(1,2)");set(gca,'YColor',[0.47,0.67,0.19]);%plot(art3(1:artEnd)-art2(1:artEnd))
        nexttile;hold on;plot((1:length(winData))/srate,(winData-art2)');title("Art.Subtracted(LPF)");
        
    end
        art2=art2-mean([art2(1) art2(end)]);
        %qBaseC=cellfun(@(x)sum((x-art2).*(1./var2))/sum(1./var2),winDataN);
    %outData.meanDiff=cellfun(@(x) trimmean(x-art2,trimPart,2),winDataN);
    outData.meanDiff=cellfun(@(x)sum((x-art2).*(1./var2))/sum(1./var2),winDataN);
    outData.qrs=qrsCorrected;
catch ME
    disp(getReport(ME));
    keyboard;
end
end
function [diffCell,modNotZero,winRem]=qrsEvent2array(EEG,eventType,refW,sectionStart,sectionEnd,plotGraph)
try
    qrsA=[];
    for i=1:length(EEG.event)
        if strcmp(EEG.event(:,i).type,eventType) ...
                && EEG.event(:,i).latency>sectionStart+EEG.srate ...
                && EEG.event(:,i).latency<sectionEnd-EEG.srate
            qrsA(end+1)=EEG.event(:,i).latency;
        end
    end
    qrsA=unique(qrsA);
    %qrsA(end)=[];qrsA(1)=[];%to avoid cutted artefacts
    
    diffQ=(diff(qrsA));
    diffQ=[mean(diffQ(1:refW)) diffQ];
    hRate=60./(diffQ/EEG.srate);
    if plotGraph
        fig=gcf;
        fig.Name='Section ECG Plots';
        %         nexttile;
        %         plot(diffQ/(EEG.srate/1000));
        %         title("R-R interval (ms) mean: "+string(mean(diffQ/(EEG.srate/1000))))
        nexttile;
        plot(qrsA/EEG.srate,hRate);
%         title("Heart-RatefadeF=arrayfun(@(x) 1/((1+exp(-(15/fadeLength*(-0.5*fadeLength+x))))),1:fadeLength); (bpm) mean: "+string(mean(hRate)))
        title("Heart-Rate(mean: "+string(mean(hRate))+" bpm)")
        xlabel('BPM')
        ylabel('Seconds')
    end
    winRem=mod(length(qrsA),refW);
    if winRem~=0
        modNotZero=true;
        qrsAt=qrsA(1:end-winRem);
        diffCell=mat2cell(qrsAt',[repmat(refW,1,length(qrsAt)/refW)]);
        diffCell{end+1}=qrsA(end-refW+1:end)';
    else
        modNotZero=false;
        diffCell=mat2cell(qrsA',[repmat(refW,1,length(qrsA)/refW)]);
    end
    
catch ME
    disp(getReport(ME));
    keyboard;
end
end

function outData=winSub(chandata,artRange,fadeLength,minAvg,sectionStart,plotData)
try
    srate=5000;
    qrs=artRange.qrs;
    artDiff=artRange.artDiff;
    dBase=artRange.meanDiff;
    fadeF=arrayfun(@(x) 1/((1+exp(-(15/fadeLength*(-0.5*fadeLength+x))))),1:fadeLength);
    %artDiff=artDiff-1;%adapt for length
    [artDiffS,artDiffI]=sort(artDiff,'descend');
    artEnd=artDiffS(end);
    artDiff(end+1)=artEnd;%adding length of last
    if minAvg<1
        minAvg=floor(length(artDiff)*minAvg);
    end
    beforeArtLenght=(artDiffS(1)-artEnd)*2+fadeLength;
    try
    winDataNa=cell2mat(arrayfun(@(q,b)chandata(q-beforeArtLenght:q+artDiffS(minAvg))-b,qrs,dBase,'UniformOutput',false));
    catch ME
    warning('Index exceeds the number of array elements.Using NaN buffer.');
%     disp(getReport(ME));
    chandataBuff=[NaN(1,beforeArtLenght) chandata NaN(1,artDiffS(minAvg))];
    winDataNa=cell2mat(arrayfun(@(x,q)chandataBuff((x+beforeArtLenght)-beforeArtLenght:(x+beforeArtLenght)+artDiffS(minAvg))-q,qrs,dBase,'UniformOutput',false));
    end
    winDataNa(artDiffI>minAvg,[1:beforeArtLenght,beforeArtLenght+artEnd:size(winDataNa,1)])=NaN;
%     if qrs(1)-beforeArtLenght<sectionStart
%     winDataNa(1,1:beforeArtLenght)=NaN;%relewant if dirst qrs is NOT omitted
%     end
    %art3=mean(smoothdata(winDataNa,'gaussian',length(qrs),'omitnan'),'omitnan');
    art3=gaussM(winDataNa);
    before=art3(1:beforeArtLenght)-(art3(beforeArtLenght)-(art3(beforeArtLenght+1)+(art3(beforeArtLenght+1)-art3(beforeArtLenght+2))));

    %     art4=art3(beforeArtLenght+1:end)+zeros(1,artDiffS(1)-length(art3(beforeArtLenght+1:end)));
art4=art3(beforeArtLenght+1:end);
art3=[before art4];
%     winDataNb=winDataNa(:,beforeArtLenght+1:end);
%     winDataNb=arrayfun(@(a,w) w-a(1:length(w)),art4,winDataNb);
%     winDataNb=winDataNb(:,1:length(art3(beforeArtLenght+1:end)));
artDiff(end)=[];
beforeA=arrayfun(@(x) before(end-(max([x-artDiffS(minAvg),0])+fadeLength):end),artDiff,'UniformOutput',false);
beforeA=cellfun(@(x) [x(1:fadeLength).*fadeF x(fadeLength+1:end)],beforeA,'UniformOutput',false);
artA=arrayfun(@(x) art4(1:min([x,artDiffS(minAvg)])),artDiff,'UniformOutput',false);
artA=cellfun(@(x) [x(1:end-fadeLength)  x(end-fadeLength+1:end).*(1-fadeF)],artA,'UniformOutput',false);

compose=arrayfun(@(x,y,z) [x{1} zeros(1, z-length(x{1}))]+[zeros(1, z-length(y{1})) y{1}],artA,beforeA,artDiff,'UniformOutput',false);
    if plotData
        fig=gcf;
        winData=arrayfun(@(x,y,z,k)(chandata(x:x+y-1)-z{1})-k,qrs(1:end-1),artDiff,compose,dBase(2:end),'UniformOutput',false);
        winDataX=arrayfun(@(x,y,z)chandata(x:x+y-1)-z,qrs(1:end-1),artDiff,dBase(2:end),'UniformOutput',false);
        nexttile;hold on;cellfun(@(y)plot((1:length(y))/srate,y),winDataX);title("Delay(2)&Art.debase(2)");ylabel('\muV');xlabel('Seconds');
%         nexttile;hold on;plot(art4);
%         nexttile;hold on;cellfun(@(x)plot(x),compose);beforePlot=plot([artDiffS(1)-max(cellfun(@(x)length(x),beforeA)):artDiffS(1)],before(end-max(cellfun(@(x)length(x),beforeA)):end)+(max(compose{1})*0.5));title("Aproximation+Xfade");
        nexttile;hold on;
        compPlot=cellfun(@(y)plot((1:length(y))/srate,y),compose);
        yyaxis right;beforePlot=plot([artDiffS(1)-max(cellfun(@(x)length(x),beforeA)):artDiffS(1)]/srate,before(end-max(cellfun(@(x)length(x),beforeA)):end));%,'DisplayName','Aprox.preQrs');
        %compPlot.Annotation.LegendInformation.IconDisplayStyle = 'off';
        for pi=2:length(compPlot)
        compPlot(pi).Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        %beforePlot.Color='black';
        ax = gca;
        axLim1=ax.YAxis(1).Limits;
        ax.YAxis(2).Limits = [axLim1(1)*0.5,axLim1(2)+abs(axLim1(1)*0.5)];
        beforePlot.LineStyle='--';
        xlabel('Seconds');
        legend('Aprox.postQrs','Aprox.preQrs');
        title("Aproximation&Xfade");
        nexttile;hold on;cellfun(@(y)plot((1:length(y))/srate,y),winData);title("Art.Subtracted&Art.debase(2)");xlabel('Seconds');
%         nexttile;hold on;plot(winData{2});plot(winDataX{2});plot(art4);
%         nexttile;hold on;;plot(cell2mat(winDataX'));plot(cell2mat(winData'));

        set(gcf,'Renderer',"painters");
    end
compose{end+1}=art4;
%compose{end+1}=[art4(1:end-fadeLength) art4(end-fadeLength+1:end).*(1-fadeF)];
outData.data=compose;
outData.before=before;
outData.qrs=qrs;
outData.artEnd=artEnd;
    %beforearrayfun(@(x)x(end-beforeLen:end)+fadeLength x(end-(beforeLen+fadeLength):end),before
%     if plotData
%         fig=gcf;
%         fig.Name=['ECG Plots sample:' char(string(plotData)) ' [' char(string(artRange{1}(1))) '-' char(string(artRange{end}(2))) '(samples)]'];
%         bot=mean(winDataNa)-std(winDataNa);
%         to=mean(winDataNa)+std(winDataNa);
%         nexttile;hold on;plot(winDataNa');legend(string(art3(1)-art3(end)))
%         nexttile;hold on;yyaxis left;plot(art3);plot(art3-mean([art3(1) art3(end)]));plot(to);plot(bot);rl=refline([0 0]);rl.Annotation.LegendInformation.IconDisplayStyle = 'off';rl.Color='black';rl.LineStyle='--';yyaxis right;plot(mean(winDataN-art3));legend("aproximation","mean after","filtered aproximation")
%         nexttile;hold on;cellfun(@plot,outData)
%     end
catch ME
    disp(getReport(ME));
    keyboard;
end
end
function out=gaussM(d)
mWeights=normpdf(normalize(d,'zscore','robust'));%also possible 'std'
%mWeights=1./(abs(normalize(B,'zscore','robust'))+1)
out=sum(d.*mWeights,'omitnan')./sum(mWeights,'omitnan');
end