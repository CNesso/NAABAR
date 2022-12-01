function final=michael_getCoord(aH,evnt,HEP,hepgui)
drawnow
f = ancestor(aH,'figure');
click_type = get(f,'SelectionType');
ptH = getappdata(aH,'CurrentPoint');
delete(ptH)
if strcmp(click_type,'normal')
    %Finding the closest point and highlighting it
    lH = findobj(aH,'Type','line');
    minDist = realmax;
    finalIdx = NaN;
    finalH = NaN;
    pt = get(aH,'CurrentPoint'); %Getting click position
    for ii = lH'
        xp=get(ii,'Xdata'); %Getting coordinates of line object
        yp=get(ii,'Ydata');
        dx=daspect(aH);      %Aspect ratio is needed to compensate for uneven axis when calculating the distance
        [newDist idx] = min( ((pt(1,1)-xp).*dx(2)).^2 + ((pt(1,2)-yp).*dx(1)).^2 );
        if (newDist < minDist)
            finalH = ii;
            finalIdx = idx;
            minDist = newDist;
        end
    end
    xp=get(finalH,'Xdata'); %Getting coordinates of line object
    yp=get(finalH,'Ydata');
    final=[xp(finalIdx),yp(finalIdx)];
    fprintf('%s\n',num2str(final));
    A=xp(finalIdx);
     %figure(hepgui.main.CurrentObject.Parent)
     
hepfig(1)=findobj('Name','HEPLAB','Type','figure');
     figure(hepfig(1));
     hepedit=findobj(hepfig(1),'Style','edit');
%     if A(1) <= get(hepgui.slider_arrow,'Max') & A(1)>=0
%     set(hepgui.slider_arrow,'Value',A(1));
%     elseif A(1)<=get(hepgui.slider_arrow,'Max')+HEP.winsec & A(1)>=0
%     set(hepgui.slider_arrow,'Value',get(hepgui.slider_arrow,'Max'));
%     end
%     delete(findobj(hepgui.main,'Type','axes'));
%     HEP.sec_ini=round(1000*get(hepgui.slider_arrow,'Value'))/1000;
%     HEP.ecg_handle=heplab_ecgplot(HEP.ecg,HEP.srate,HEP.qrs,HEP.sec_ini,HEP.ecg_handle,HEP.winsec);
%     end
%     set(hepgui.win_label, 'String',num2str(HEP.sec_ini));
%     axes(HEP.ecg_handle);
    %figure(hepgui.main.CurrentObject.Parent),
    hepwin=str2num(get(hepedit(1),'String'))/2;
    set(hepedit(end),'String',num2str(A-hepwin));
    notify(hepgui.main,'SizeChanged');
    
%     ptt=pt(:,3)>0
%     out=(pt(pt(:,3)>0,[1 2]))
%     [xp',yp']  
    
    
%     ptH = plot(aH,xp(finalIdx),yp(finalIdx),'k*','MarkerSize',20);
%     setappdata(aH,'CurrentPoint',ptH);
% elseif strcmp(click_type,'alt')
%     %do your stuff once your point is selected   
%     disp('Done clicking!');
%     % HERE IS WHERE YOU CAN PUT YOUR STUFF
%     uiresume(f);
end