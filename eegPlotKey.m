function eegPlotKey(pKey)
%helper to function eegplot_readkey.
%insert at line 8:
%     else
%         eegPlotKey(evnt.Key)


%eegplot( EEG.data,'winrej',TMPREJ,'color',colPlot,'winlength',20,'command','fprintf(''Rejections saved in variable TMPREJ\n'');','events', EEG.event,'eloc_file', EEG.chanlocs );
            global TMPREJ eegSrateG
ax1 = findobj('tag','eegaxis','parent',gcf);         % axes handle
g = get(gcf,'UserData');
data = get(ax1, 'userdata');
ESpacing = findobj('tag','ESpacing','parent',gcf);   % ui handle
EPosition = findobj('tag','EPosition','parent',gcf);

if strcmp(pKey, 'uparrow')==1
afterPos=sort(TMPREJ(find(TMPREJ(:,1)/eegSrateG>(str2num(get(EPosition,'string')))+round(g.winlength/2)+1)));
if ~isempty(afterPos)
newPos=ceil(afterPos(1,1)/eegSrateG);
newWin=newPos-g.winlength/2;
set(EPosition,'string',num2str(newWin));
eegplot('drawp',0);
end
elseif strcmp(pKey, 'downarrow')==1
prePos=sort(TMPREJ(find(TMPREJ(:,1)/eegSrateG<str2num(get(EPosition,'string')))));
if ~isempty(prePos)
newPos=round(prePos(end,1)/eegSrateG);
newWin=newPos-g.winlength/2;
set(EPosition,'string',num2str(newWin));
eegplot('drawp',0);
end
end
end