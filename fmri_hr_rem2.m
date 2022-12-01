%%
fprintf('\n Starting BCG removal \n');

tic
morePlot=true;
a=[eegT.event(strcmp({eegT.event.type},eventType)).latency]';%substitute with qrs1
%a=qrsA;%substitute later
%diffA=diff(a)
padA=round(mean(diff(a))/2);
dBstart=0.2;%in seconds
dBstart=round(dBstart*EEG.srate);
if dBstart>padA
    padA=dBstart;
end
if a(1)<padA
    a(1)=[];
end
if a(end)+padA>size(EEG.data,2)
    a(end)=[];
end

winLen=20;
minWin=ceil(winLen*(2/3));% 20;%can also be  calculated like winLen*.75 
winI=1:winLen: length(a);

% if mod(length(a),winLen)~=0 %maybe length(a)+1?
% winI(end)=length(a)-winLen;
% end
winI(end)=length(a)-winLen;
if (~exist('plotVis','var'))
    plotVis=true;
end
if (~exist('eegF','var'))
    fprintf('\n Filtering:\n');
    %eegF=pop_eegfiltnew(EEG,1,[]);
    %eegF=pop_eegfiltnew(eegF,0,qrsLpf);
    eegF=pop_eegfiltnew(eegT,0,qrsLpf);

    %     EEG.ECG.event=EEG.event;EEG.ECG=pop_fmrib_qrsdetect(EEG.ECG,3,'qrs1','no');EEG.event=EEG.ECG.event
end

cellA=arrayfun(@(x){a(x:x+winLen)},winI);
cellAdatalim=cellfun(@(x){[min(x)-padA max(x)+padA]},cellA);
%cellAdata=cellfun(@(x){chanData(x(1):x(2))},cellAdatalim);
cellAevents=cellfun(@(x,y){x-y(1)},cellA,cellAdatalim); 
%cellAlen=cellfun(@(x){diff(x)},cellA);
cellAevents=cellfun(@(x){x(1:end-1)},cellAevents);
cellAdB=cellfun(@(x){zeros(size(x))},cellAevents);
cellAerp=cellfun(@(x){[]},cellAevents);
cellAvar=cellfun(@(x){[]},cellAevents);
cellAbefore=[];
cellAafter=[];
cellAfade=cellfun(@(x){[]},cellAevents);
%arrayfun(@(y){cellfun(@(x){eegF.data(y,x(1):x(2))},cellAdatalim)},1:size(eegF.data,1))

qrsS.c=cell2mat(arrayfun(@(y){struct('qrs',cellA,'qrsR',cellAevents,...'data',cellAdata,'len',cellAlen,
    'dataLim',cellAdatalim,'dB',cellAdB,'delay',cellAdB,...
    'erp',cellAerp,'var',cellAvar,'delayA',cellAerp,'dBa',cellAerp,...
    'before',cellAbefore,'after',cellAafter,'fade',cellAfade,...
    'data',cellfun(@(x){eegF.data(y,x(1):x(2))},cellAdatalim))},1:size(eegF.data,1))');

qrsS.o=cell2mat(arrayfun(@(y){struct(...
        'data',cellfun(@(x){eegT.data(y,x(1):x(2))},cellAdatalim))},1:size(eegT.data,1))');
qrsS.h={rmfield(qrsS.c,'data')};%qrsS.history{end+1}=rmfield(qrsS.chan,'data');

qrsS.p.padA=padA;
qrsS.p.winI=winI;
qrsS.p.winLen=winLen;
qrsS.p.minWin=minWin;
qrsS.p.dbStart=dBstart;
qrsS.p.trimPart=50;
%qrsS.p.maxDiff=repmat(max(diff(a)),size(eegF.data,1),1);
%melgio un valore singolo
qrsS.p.minDiff=min(diff(a));

qrsS.p.lastQ=repmat(a(end),size(eegF.data,1),1);
qrsS.p.maxD=eegF.srate*0.05;%check this number!!
qrsS.p.fadeLen=qrsS.p.maxD*2;

qrsS.f.randCol=@(x)cell2mat(arrayfun(@(i){x(randperm(size(x,1)),i)},1:size(x,2)));
%qrsS.f.meanFunRoot=@(x)mean(smoothdata(x,'gaussian',size(x,2)));
qrsS.f.meanFunRoot=@(x)mean(smoothdata(x,'gaussian',size(x,2)),'omitnan');
qrsS.f.meanFunRand=@(x)qrsS.f.meanFunRoot(qrsS.f.randCol(x));
%qrsS.f.meanFun=qrsS.f.meanFunRand;
if (exist('meanMethod','var'))
switch meanMethod
    case 'gauss'
        disp('gauss')
        qrsS.f.meanFun=@gaussM;
    case 'mean'
        disp('mean')
        qrsS.f.meanFun=@mean;
    case 'smooth'
        disp('smooth')
        %qrsS.f.meanFun=qrsS.f.meanFunRand;
        qrsS.f.meanFun=qrsS.f.meanFunRoot;
end
else
  qrsS.f.meanFun=@gaussM;  
end
% qrsS.f.meanFun=@mean
qrsS.f.bBase=@(a,b,trimPart)arrayfun(@(x) trimmean(a,trimPart),b);
qrsS.f.fadeF=@(fadeLength)arrayfun(@(x) 1/((1+exp(-(15/fadeLength*(-0.5*fadeLength+x))))),1:fadeLength);
qrsS.p.fadeArr=qrsS.f.fadeF(qrsS.p.fadeLen);

ran=randi(length(cellA));
chan=60;

if morePlot
figure;
nexttile;hold on;
arrayfun(@(x,y)plot(eegF.data(chan,x:x+y)),qrsS.c(chan,ran).qrs(1:end-1),diff(qrsS.c(chan,ran).qrs)-1);
end

fprintf('\n Preparation finisched\n Starting aligment (%.2f s)\n',toc);


dToc=toc;
qrsS.c=arrayfun(@(x)align1(x,qrsS.p,qrsS.f),qrsS.c);
qrsS.h{end+1}=rmfield(qrsS.c,'data');
fprintf('\n First aligment (%.2f s)\n',toc-dToc);

if morePlot
nexttile;hold on;
arrayfun(@(x,y,z)plot(eegF.data(chan,x:x+y)-z),qrsS.c(chan,ran).qrs(1:end-1),diff(qrsS.c(chan,ran).qrs)-1,qrsS.c(chan,ran).dB);
end

dToc=toc;
if (exist('secondAlign','var'))
    
    if secondAlign
        qrsS.c=arrayfun(@(x)align2(x,qrsS.p,qrsS.f),qrsS.c);
    end
    secondAlign=secondAlign
else
qrsS.c=arrayfun(@(x)align2(x,qrsS.p,qrsS.f),qrsS.c);
end
qrsS.h{end+1}=rmfield(qrsS.c,'data');
fprintf('\n Second aligment (%.2f s)\n',toc-dToc);

if morePlot
nexttile;hold on;
arrayfun(@(x,y,z)plot(eegF.data(chan,x:x+y)-z),qrsS.c(chan,ran).qrs(1:end-1),diff(qrsS.c(chan,ran).qrs)-1,qrsS.c(chan,ran).dB);
nexttile;hold on;%yyaxis 'right';
plot(qrsS.c(chan,ran).var(1:min(diff(qrsS.c(chan,ran).qrs))))
nexttile;hold on;
arrayfun(@(x,y)plot(eegF.data(chan,x:x+y)),qrsS.c(chan,ran).qrs(1:end-1),diff(qrsS.c(chan,ran).qrs)-1);

end
%% group align
maxI=0;
minI=0;
chanErpI=[];
chanErpCell=cell(1,size(qrsS.c,1));
%dBaseCell=cell(1,size(qrsS.c,1));
dBaseCell=cell(size(qrsS.c));
delayCell=cell(size(qrsS.c));
delayI=0;
delayIa=0;
chanI=nan;%randi(62)
%chan 38 on subject 22 problematic but does it justify recalculate all erps
% CONTROL big delays!!!!!!!
for i=1:size(qrsS.c,1)
%for i=chanI
   maxI=max(arrayfun(@(x)length(x.erp),qrsS.c(i,:)),[],'all');
    minI=min(arrayfun(@(x)length(x.erp),qrsS.c(i,:)),[],'all');

%     chanErpI=cell2mat(cellfun(@(x){[x nan(1,maxI-length(x))]},{qrsS.c(i,:).erp}'));
%     chanErpCell{i}=qrsS.f.meanFun(chanErpI(:,1:minI));
%     dBaseCell{i}=arrayfun(@(x) trimmean(chanErpI(x,1:1:size(chanErpCell{i},2))-chanErpCell{i},qrsS.p.trimPart),1:size(chanErpI,1));
    %dBaseCell{i}=arrayfun(@(x) trimmean(chanErpI(x,1:1:size(chanErpCell{i},2))-chanErpCell{i},qrsS.p.trimPart),1:size(chanErpI,1));
   
    if (i==chanI)
     figure;hold on;
     nexttile;hold on;
     cellfun(@plot,{qrsS.c(i,:).erp}')
    end
    
    chanErpI=cellfun(@(x){x(1:minI)},{qrsS.c(i,:).erp}');
    
    if (i==chanI) 
     nexttile;hold on;
     plot(cell2mat(chanErpI)')
    end
    
    [~,delayI]=cellfun(@(x) max(xcorr(x,chanErpI{1}, qrsS.p.maxD*2)), chanErpI(2:end));
    delayI=[0 ; delayI-qrsS.p.maxD*2];
    delayI=delayI+abs(min(delayI));
    delayIa=delayI;
    
    chanErpI=arrayfun(@(x,y){x{1}(1+y:minI)},{qrsS.c(i,:).erp}',delayI);
    
    if (i==chanI) 
     nexttile;hold on;
     cellfun(@plot,chanErpI)
    end
    
    minI=min(cellfun(@(x)length(x),chanErpI));
    minIo=minI;
    chanErpI=cellfun(@(x){x(1:minI)},chanErpI);
    chanErpCell{i}=qrsS.f.meanFun(cell2mat(chanErpI));
    
    dBaseCell(i,:)=arrayfun(@(x) {trimmean(x{1}-chanErpCell{i},qrsS.p.trimPart)},chanErpI);
    %dBaseCell{i};
    [~,delayI]=cellfun(@(x,y) max(xcorr(x-y,chanErpCell{i}, qrsS.p.maxD*2)), chanErpI,dBaseCell(i,:)');
    delayI=[delayI-qrsS.p.maxD*2];
    delayI=delayI+abs(min(delayI));
    delayI=delayI+delayIa;
    
    chanErpI=arrayfun(@(x,y){x{1}(1+y:minI)},{qrsS.c(i,:).erp}',delayI);
    
    if (i==chanI) 
     nexttile;hold on;
     cellfun(@plot,chanErpI)
    end
    

    minI=min(cellfun(@(x)length(x),chanErpI));
    chanErpI=cellfun(@(x){x(1:minI)},chanErpI);
    chanErpCell{i}=qrsS.f.meanFun(cell2mat(chanErpI));
    dBaseCell(i,:)=arrayfun(@(x) {trimmean(x{1}-chanErpCell{i},qrsS.p.trimPart)},chanErpI);%dangerous but keep
%     dBaseCell{i}
    %possible last step is to add mean to avoid a further drift
    %see:
    %[delayI round(delayI-mean(delayI))]
    %use
    %delayI=round(delayI-mean(delayI))
    %delayI=round(delayI-median(delayI));
    delayCell(i,:)=arrayfun(@(x){x},delayI);
    if (i==chanI) 
     nexttile;hold on;
     cellfun(@(x,y)plot(x-y),chanErpI,dBaseCell(i,:)')
    end
%     nexttile;plot(chanErpCell{i})

    %if (i==chanI);break;end
end
if ~isnan(chanI) 
figure;nexttile;hold on;arrayfun(@(x)plot(x.erp),qrsS.c(chan,:))
end

dToc=toc;
qrsS.c=arrayfun(@(x,y,z,w)finalCorrect(x,qrsS.p,qrsS.f,y{1},z{1},w),qrsS.c,dBaseCell,delayCell,qrsS.o);
qrsS.h{end+1}=rmfield(qrsS.c,'data');
fprintf('\n Group aligment justification (%.2f s)\n',toc-dToc);

if ~isnan(chanI) 
nexttile;hold on;arrayfun(@(x)plot(x.erp),qrsS.c(chan,:))
figure;
hold on;
cellfun(@plot,chanErpI);
end
%% mod compensation
bqrs=qrsS.c;
qrsMod=mod(length(a),winLen);
if qrsMod==0%to avoid problem when nuber of qrs are a multiplacaor of the number of windows 
    qrsMod=1;
else
    qrsMod=(winLen-(qrsMod-1));%empiricaly tested to fit with lsat window
end
%arrayfun(@(x){qrsS.c(i,end).qrs(qrsMod+1)},1:size(qrsS.c,1))
%arrayfun(@(x){x.qrs(1:end-1)},qrsS.c)
for i=1:size(qrsS.c,1)%awfull please correct
    qrsS.c(i,end).qrs(1:qrsMod)=[];
end
qrsCompl=arrayfun(@(x){[x.qrs(1:end-1)]},qrsS.c);
qrsCompl=cell2mat(arrayfun(@(x){cell2mat(qrsCompl(x,:)')},1:size(qrsCompl,1)));
%qrsS.c=arrayfun(@(x)align1(x,qrsS.p,qrsS.f),qrsS.c);
qrsS.c=bqrs;
qrsS.p.minDiff=min(diff(qrsCompl),[],'all');
qrsS.c=arrayfun(@(x)winFader(x,qrsS.p,qrsS.f),qrsS.c);
for i=1:size(qrsS.c,1)
    qrsS.c(i,end).qrs(1:qrsMod)=[];
    qrsS.c(i,end).qrsR(1:qrsMod)=[];
    qrsS.c(i,end).dB(1:qrsMod)=[];
    qrsS.c(i,end).fade(1:qrsMod)=[];
    qrsS.c(i,end).delay(1:qrsMod)=[];
    qrsS.p.lastQ(i)=qrsS.c(i,end).qrs(end);
end
qrsS.h{end+1}=rmfield(qrsS.c,'data');

%% fade this needs to be rewritten for neutral
diffI=0;
erpI=[];
beforeI=[];
afterI=[];
lastLenI=[];
completeAprox=cell(size(qrsS.c));
startFade=cell(size(qrsS.c,1),1);
completeChan=cell(size(qrsS.c,1),1);
completeChanstartstop=cell(size(qrsS.c,1),1);
for i=1:size(qrsS.c,1)
    startFade{i}=qrsS.c(i,1).before((end-qrsS.p.fadeLen)+1:end).*qrsS.p.fadeArr;
    for j=2:size(qrsS.c,2)    
        diffI=qrsS.c(i,j).qrs(1)-qrsS.c(i,j-1).qrs(end-1);
        erpI=qrsS.c(i,j-1).fade{end};
        diffI=(diffI-length(erpI));
        beforeI=qrsS.c(i,j).before((length(qrsS.c(i,j).before)-max([(diffI),0]))-qrsS.p.fadeLen+1:end);
        afterI={erpI(1:(end+min([(diffI),0]))-qrsS.p.fadeLen) ...
                erpI(((end+min([(diffI),0]))-qrsS.p.fadeLen)+1:end+min([(diffI),0])).*(1-qrsS.p.fadeArr)};
       
            beforeI={beforeI(1:qrsS.p.fadeLen).*(qrsS.p.fadeArr) ...
                beforeI(qrsS.p.fadeLen+1:end)};
            qrsS.c(i,j-1).fade{end}=[afterI{1} afterI{2}+beforeI{1} ...
                beforeI{2}];
        completeAprox{i,j-1}=cell2mat(qrsS.c(i,j-1).fade');
    end
    lastLenI=qrsS.c(i,j).qrs(end)-qrsS.c(i,j).qrs(end-1);
    qrsS.c(i,j).fade{end}=[
        qrsS.c(i,j).fade{end}(1:min([lastLenI,length(qrsS.c(i,j).fade{end})])),...
        nan(1,lastLenI-length(qrsS.c(i,j).fade{end}))];
        %completeAprox{i,j}=single([]);
        completeAprox{i,j}=cell2mat(qrsS.c(i,j).fade');
        %to control length use:
        %completeChan{i}=[startFade{i} cell2mat(completeAprox(i,:))];
        %(length(completeChan{i})-length(startFade{i}))-(qrsS.p.lastQ(i)-(qrsS.c(i,1).qrs(1)))
        completeChan{i}=[zeros(1,(qrsS.c(i,1).qrs(1)-length(startFade{i}))),...
            startFade{i},cell2mat(completeAprox(i,:)),...
            repmat(qrsS.c(i,j).fade{end}(end),1,qrsS.p.fadeLen).*(1-qrsS.p.fadeArr),...
            zeros(1,(size(eegT.data,2)-(qrsS.p.lastQ(i)+qrsS.p.fadeLen)))];
end
if morePlot
figure;plot(completeChan{60})
end

% EEGc=EEG;
% EEGc.data=EEG.data-cell2mat(completeChan);
% EEGc=eeg_addnewevents(EEGc,{HEP.qrs'},{char(eventType)})
eegTemp=eegT.data-cell2mat(completeChan);
for c=1:size(eegT.data,1)
    eegTemp(c,isnan(eegTemp(c,:)))=eegT.data(c,isnan(eegTemp(c,:)));
end
eegT.data=eegTemp;
%eegT.data=eegT.data-cell2mat(completeChan);
fprintf('\n Total BCG correction time %.2f s\n\n',toc);
if plotVis
vis_artifacts(eegT,EEG);
end
%%
function out=finalCorrect(d,p,f,db,delay,o)
    d.dB=d.dB-db;
    d.delay=d.delay+delay;
    d.qrs=d.qrs+delay;
    d.qrsR=d.qrsR+delay;
    artEnd=min(diff(d.qrs))-1;

    winData=cell2mat(arrayfun(@(x){d.data(x: x+artEnd)},d.qrsR))-d.dB;
    d.data=o.data;
    d.erp=f.meanFun(winData);
    d.var=mean((winData-d.erp).^2);
    %d.var=var(winData);
    out=d;
end
function out=gaussM(d)
mWeights=normpdf(normalize(d,'zscore','robust'));%also possible 'std'
%mWeights=1./(abs(normalize(d,'zscore','robust'))+1);
out=sum(d.*mWeights,'omitnan')./sum(mWeights,'omitnan');
%out=smoothdata(out,'gaussian',100);%(1/50)*EEG.srate)%~50Hz filter
end
function out=align1(d,p,f)
    artEnd=min(diff(d.qrs))-1;
    qBase=arrayfun(@(x) trimmean(d.data(x-p.dbStart:x),p.trimPart),d.qrsR);
    winData=arrayfun(@(x,y){d.data(x: x+artEnd)-y},d.qrsR,qBase);
    [~,delay]=cellfun(@(x) max(xcorr(x,winData{1}, p.maxD)), winData(2:end));
    delay=[0 ; delay-p.maxD];
    
    d.dB=qBase;
    d.delay=delay;
    d.qrs=d.qrs+[delay;delay(end)];
    d.qrsR=d.qrsR+delay;
    artEnd=min(diff(d.qrs))-1;

    winData=cell2mat(arrayfun(@(x,y){d.data(x: x+artEnd)},d.qrsR))-qBase;
    d.erp=f.meanFun(winData);
    d.var=mean((winData-d.erp).^2);
    %d.var=var(winData);

    out=d;
end
function out=align2(d,p,f)
    artEnd=length(d.erp)-1;
    
    qBase=arrayfun(@(x)sum((d.data(x:x+artEnd)-d.erp).*(1./d.var))/sum(1./d.var),d.qrsR);
    winData=arrayfun(@(x,y){d.data(x: x+artEnd)-y},d.qrsR,qBase);
    [~,delay]=cellfun(@(x) max(xcorr(x,d.erp,p.maxD)), winData);

    delay=delay-max(delay);
    d.delay=delay;
    d.qrs=d.qrs+[delay;delay(end)];
    d.qrsR=d.qrsR+delay;
    artEnd=min(diff(d.qrs))-1;
    
    winData=cell2mat(arrayfun(@(x,y){d.data(x: x+artEnd)},d.qrsR))-qBase;
    d.erp=f.meanFun(winData);
    d.var=mean((winData-d.erp).^2);
    %d.var=var(winData);
    d.erp=d.erp-mean(d.erp);
    qBase=arrayfun(@(x,y)sum((d.data(x+y:x+y+artEnd)-d.erp).*(1./d.var))/sum(1./d.var),d.qrsR,d.delay);
    
    d.dB=qBase;
    out=d;
end
function out=alignfinal(erp,p,f)
    [~,delay]=cellfun(@(x) max(xcorr(x,erp{1}, p.maxD)), erp(2:end));
    delay=[0 ; delay-p.maxD];
    out=delay;
end
function out=meanNA(a,len,data,f)
ma=max(len);
%pad data with nan
win=cell2mat(arrayfun(@(x,y){[data(x:x+y-1) nan(1,ma-y)]},a(1:end-1),len));%exclude last to eqyuiliza array length
%out=mean(win,'omitnan')
out=f.meanFun(win)

end

function [d,artOut,beforeOut]=winFader(d,p,f)
diffQrs=diff(d.qrs);
%minDiff=min(diffQrs);%provissiory use per channel later
minDiff=p.minDiff;
minWin=p.minWin;
artEnd=minDiff;%provissiory

sor=sort(diffQrs,'descend');
diffQrsMin=sor(min([minWin,length(sor)]))<=diffQrs;
%schift forvard to accomodate that its the following has its before
diffQrsMinShift=circshift(diffQrsMin,1);
diffQrsMinShift(1)=0;%to improve to avoid carryover
beforeLen=(max(diffQrs)-minDiff)+p.fadeLen;

beforeQrsR=[d.qrsR(diffQrsMinShift)-beforeLen d.qrsR(diffQrsMinShift)+p.fadeLen d.dB(diffQrsMinShift)];%+2 for a startpoint and an endpoint
afterQrsR=[(d.qrsR(diffQrsMin)+artEnd)-p.fadeLen d.qrsR(diffQrsMin)+diffQrs(diffQrsMin) d.dB(diffQrsMin)];%+2 for a startpoint and an endpoint

before=cell2mat(arrayfun(@(x){d.data(beforeQrsR(x,1):beforeQrsR(x,2))-beforeQrsR(x,3)},1:size(beforeQrsR,1))');
beforeMean=f.meanFun(before);

after=cell2mat(arrayfun(@(x){d.data(afterQrsR(x,1):afterQrsR(x,1)+min(afterQrsR(:,2)-afterQrsR(:,1)))-afterQrsR(x,3)},1:size(afterQrsR,1))');
afterMean=f.meanFun(after);

endings=diffQrs;
endings(diffQrsMin)=sor(min([minWin,length(sor)]));
%endings(end)=sor(minWin);%to ensure no cutting on the last-> better
%correction at the end
diffEndings=diffQrs-endings;

beforeMeanCell={beforeMean(1:end-p.fadeLen) beforeMean(end-(p.fadeLen-1):end).*(1-p.fadeArr)};
afterMeanCell={afterMean(1:p.fadeLen).*p.fadeArr afterMean(p.fadeLen+1:end)};

artMeanCell={d.erp(1:p.fadeLen).*p.fadeArr,...
    d.erp(p.fadeLen+1:(end-(p.fadeLen))),...
    d.erp(end-(p.fadeLen-1):end).*(1-p.fadeArr)};
artMeanCellArr=[beforeMeanCell{2}+artMeanCell{1},...
    artMeanCell{2},artMeanCell{3}+afterMeanCell{1},afterMeanCell{2}];
artArr=repmat({artMeanCellArr},p.winLen,1);

beforeArr=repmat({beforeMean(1:end-p.fadeLen)},p.winLen,1);%subtract the already faded part
beforeArr=arrayfun(@(x,y){{x{1}((end-y)-(p.fadeLen-1):end-y).*p.fadeArr ...
    x{1}((end-y)+1:end)}},beforeArr,diffEndings);

artFadeArr=arrayfun(@(x,y,z){[x{1}(1:y-p.fadeLen) ...
    x{1}(y-(p.fadeLen-1):y).*(1-p.fadeArr)+z{1}{1} z{1}{2}]},artArr,endings,beforeArr);
artFadeArr{end}=artMeanCellArr;% removed after from last.
beforeOut=beforeMeanCell{1};
artOut=artFadeArr;
d.fade=artFadeArr;
d.before=beforeMeanCell{1};
d.after=afterMeanCell{2};
end