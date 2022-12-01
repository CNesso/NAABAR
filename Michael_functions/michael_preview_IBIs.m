% michael_preview_IBIs(HEP,subsectionsI)
function michael_preview_IBIs(HEP,subsectionsI,refW,hepgui)
% qrs = HEP.qrs(subsectionsI(HEP.qrs))./HEP.srate;
% %qrs(1)=[];
% x=qrs(2:end);
% y=diff(qrs);

x=HEP.qrs(subsectionsI(HEP.qrs))./HEP.srate;
y= diff(x);
[~,maxQ]=max(y);
y(maxQ)=[];
x([maxQ length(x)])=[];
if isgraphics(333);clf(333);end;figure(333);
if exist('refW','var')
    %stdPlot=movstd(diff(y)*HEP.srate,refW);
    subplot(2,1,1);
%     aH1=axes;
%     winAvg=movavg(diff(y)*HEP.srate,'simple',refW);
    % winAvg=movmean(diff(y)*HEP.srate,refW,'Endpoints','fill');
    %plot(x(2:end),y(1:end-1)-winAvg);
    qrsDiffAvg=abs(smoothdata(diff(x),"gaussian",length(x)/refW)-diff(x));
    qrsDiffAvg=[0;abs(mean(qrsDiffAvg)-qrsDiffAvg)];
    plot(x,qrsDiffAvg);
    ylim([0,mean(qrsDiffAvg)+2*std(qrsDiffAvg)]);
%title(sprintf('Interbeat Intervals window(meanStd: %.3g)',mean(stdPlot)));
title(sprintf('Difference to Gauss(%d) smothed curve',refW));
xlabel('Seconds')
ylabel('Sample')
%             dcm = datacursormode;
%             dcm.Enable = 'on';(5:end)
%             dcm.DisplayStyle = 'window';
    subplot(2,1,2);
end
% aH(2)=axes;
stdPlot=std(diff(y));
plot(x,y);
title(sprintf('Interbeat Intervals(std: %.3g seconds)',stdPlot));
xlabel('Seconds')
ylabel('Seconds')
aH=findobj(333,'Type','axes');
lH = findobj(aH,'Type','line');
set(lH,'hittest','off'); % so you can click on the Markers
set(aH,'ButtonDownFcn',{@michael_getCoord,HEP,hepgui});
%             dcm = 
% datacursormode on;
%             dcm.Enable = 'on';
%             dcm.DisplayStyle = 'window';
end