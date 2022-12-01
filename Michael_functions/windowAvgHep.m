function outData=windowAvgHep(chandata,debaseStart,qrs,maxD,trimPart,plotData,subsec)
try
    qrsO=qrs;
    mapI=1:length(qrs);
    mapI=mapI(subsec(qrs));
    qrs=qrs(subsec(qrs));
    %(randi(10)+1)
    %qrsStd=std(qrs)[1:beforeArtLenght
    artEnd=round(min(diff(qrs))/2);
    qBase=arrayfun(@(x) trimmean(chandata(x-debaseStart:x),trimPart),qrs);
    winData=arrayfun(@(x)chandata(x-artEnd: x+artEnd),qrs,'UniformOutput',false);
    winDataO=winData;
    
    %[~,del1]=cellfun(@(x) max(xcov(x,winData{1}, maxD)), winData(2:end));
    [~,del1]=cellfun(@(x) max(xcov(x,winData{1}, maxD)), winData(2:end));
    del1=[0  del1-maxD];
    
    winData=cell2mat(arrayfun(@(x,y)chandata((x+y)-artEnd: x+y+artEnd),qrs,del1,'UniformOutput',false));
    winData1=winData-qBase;
    art1=mean(smoothdata(winData1,'gaussian',length(qrs)));
    
    
    
    [~,del2]=cellfun(@(x) max(xcov(x-mean(x-art1,2),art1,maxD)), winDataO);
    del2=del2-maxD;
        qrsCorrected=qrs+del2;
    outData.artDiff=diff(qrsCorrected);
    artEnd=round(min(outData.artDiff)/2);

    winDataN=arrayfun(@(x)chandata(x-artEnd: x+artEnd),qrsCorrected,'UniformOutput',false);
    winData=cell2mat(winDataN)-qBase;

    art2=mean(smoothdata(winData,'gaussian',length(qrs)));
    
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
        nexttile;hold on;cellfun(@plot, winDataO);
        nexttile;hold on;plot(winData1');
        nexttile;plot((winData-mean(winData1-art1,2))');
        
        nexttile;hold on;plot(winData');
        nexttile;hold on;plot(art1);plot(art2);yyaxis right;plot(art2(1:min([length(art2) length(art1)]))-art1(1:min([length(art2) length(art1)])));legend("1st aproximation","2nd aproximation","Difference("+string(round(mean(abs(del2))))+"samp.del.)");%plot(art3(1:artEnd)-art2(1:artEnd))
        nexttile;hold on;plot((winData-art2)');
        
    end
%         art2=art2-mean([art2(1) art2(end)]);
%     outData.meanDiff=cellfun(@(x) trimmean(x-art2,trimPart,2),winDataN);
%     outData.qrs=qrsCorrected;
qrsO(mapI)=qrsCorrected;
    %outData=qrsCorrected;
    outData=qrsO;
catch ME
    disp(getReport(ME));
    keyboard;
end
end