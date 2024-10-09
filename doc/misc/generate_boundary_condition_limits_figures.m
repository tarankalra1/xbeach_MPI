function generate_boundary_condition_limits_figures
% Generates illustrative plots for acceptable XBeach boundary condition
% depths


%% set up x and y meshes

Hrange = [2:0.1:8];
Drange = [0:1:40];
Trange = [6:0.25:18];

[HD,DH] = meshgrid(Hrange,Drange);
[TD,DT] = meshgrid(Trange,Drange);

%% Factors
HiD = HD./DH;
[c,cg] = wavevelocity(TD,DT);
CgiC = cg./c;CgiC(isnan(CgiC)) = 1;
k = disper(2*pi./TD,DT,9.81);
KH = k.*DT;

%% plot colors
gre = rgb('greenyellow');
yel = rgb('gold');
ora = rgb('darkorange');
red = rgb('darkred');
sil = rgb('dimgray');

%% surfbeat plot
figure;
s(1) = subplot(1,2,1);hold on;
levels = [1/1;1/2;1/3;1/4;1/5;0];
cmap = [red;red;ora;gre;gre;gre];
make_subplot_local(s(1),HD,DH,HiD,[1/1;1/2;1/3;1/4;1/5;0],[red;red;ora;gre;gre;gre],gre,{'1/1';'1/2';'1/3';'1/4';'1/5'});
title('Wave height to water depth ratio (H_s/h)');
xlabel('Significant wave height H_s (m)');
ylabel('Water depth (m)');
xt = 2:8;
yt = 0:5:40;
xl = xlim;yl=ylim;
pp =[];
for i=2:length(xt)-1;
    pp(end+1) = plot([xt(i) xt(i)],[yl(1) yl(2)],'k:');
end
for i=2:length(yt)-1;
    pp(end+1) = plot([xl(1) xl(2)],[yt(i) yt(i)],'k:');
end
set(pp,'color',sil);
set(s(1),'xtick',xt,'ytick',yt);

s(2) = subplot(1,2,2);hold on;
make_subplot_local(s(2),TD,DT,CgiC,[0.95:-0.05:0.6 0],[red;red;ora;yel;gre;gre;gre;gre;gre],gre);
title('Wave group velocity to wave celerity ratio (n = c_g/c)');
xlabel('Wave period T_{m-1,0} \approx T_p/1.1 (s)');
ylabel('Water depth (m)');
xt = 6:2:18;
yt = 0:5:40;
xl = xlim;yl=ylim;
pp =[];
for i=2:length(xt)-1;
    pp(end+1) = plot([xt(i) xt(i)],[yl(1) yl(2)],'k:');
end
for i=2:length(yt)-1;
    pp(end+1) = plot([xl(1) xl(2)],[yt(i) yt(i)],'k:');
end
set(pp,'color',sil);
set(s(2),'xtick',xt,'ytick',yt);

fname = 'boundaries_SB';
saveas(gcf,[fname '.fig']);
savefigure(gcf,[fname '.png'],'size',4*[2 1],'dpi',900,'margins',1.1,'printoptions','-dpng');
savefigure(gcf,[fname '.eps'],'size',4*[2 1],'dpi',300,'margins',1.1,'printoptions','-depsc2');
close;

%% nonh plot
figure;
s(1) = subplot(1,3,1);hold on;
levels = [1/1;1/2;1/3;1/4;1/5;0];
cmap = [red;red;ora;gre;gre;gre];
make_subplot_local(s(1),HD,DH,HiD,[1/1;1/2;1/3;1/4;1/5;0],[red;red;ora;gre;gre;gre],gre,{'1/1';'1/2';'1/3';'1/4';'1/5'});
title('Wave height to water depth ratio (H_s/h)');
xlabel('Significant wave height H_s (m)');
ylabel('Water depth (m)');
xt = 2:8;
yt = 0:5:40;
xl = xlim;yl=ylim;
pp =[];
for i=2:length(xt)-1;
    pp(end+1) = plot([xt(i) xt(i)],[yl(1) yl(2)],'k:');
end
for i=2:length(yt)-1;
    pp(end+1) = plot([xl(1) xl(2)],[yt(i) yt(i)],'k:');
end
set(pp,'color',sil);
set(s(1),'xtick',xt,'ytick',yt);

s(2) = subplot(1,3,2);hold on;
make_subplot_local(s(2),TD,DT,CgiC,[0.95:-0.05:0.6 0],[red;red;ora;yel;gre;gre;gre;gre;gre],gre);
title('Wave group velocity to wave celerity ratio (n = c_g/c)');
xlabel('Wave period T_{m-1,0} \approx T_p/1.1 (s)');
ylabel('Water depth (m)');
xt = 6:2:18;
yt = 0:5:40;
xl = xlim;yl=ylim;
pp =[];
for i=2:length(xt)-1;
    pp(end+1) = plot([xt(i) xt(i)],[yl(1) yl(2)],'k:');
end
for i=2:length(yt)-1;
    pp(end+1) = plot([xl(1) xl(2)],[yt(i) yt(i)],'k:');
end
set(pp,'color',sil);
set(s(2),'xtick',xt,'ytick',yt);

s(3) = subplot(1,3,3);hold on;
make_subplot_local(s(3),TD,DT,KH,[0.5 1 1.1 1.25 2],[gre;yel;ora;red;red],gre);
title('Relative water depth (kh)');
xlabel('Wave period T_{m-1,0} \approx T_p/1.1 (s)');
ylabel('Water depth (m)');
xt = 6:2:18;
yt = 0:5:40;
xl = xlim;yl=ylim;
pp =[];
for i=2:length(xt)-1;
    pp(end+1) = plot([xt(i) xt(i)],[yl(1) yl(2)],'k:');
end
for i=2:length(yt)-1;
    pp(end+1) = plot([xl(1) xl(2)],[yt(i) yt(i)],'k:');
end
set(pp,'color',sil);
set(s(3),'xtick',xt,'ytick',yt);

fname = 'boundaries_NH';
saveas(gcf,[fname '.fig']);
savefigure(gcf,[fname '.png'],'size',4*[3 1],'dpi',900,'margins',1.1,'printoptions','-dpng');
savefigure(gcf,[fname '.eps'],'size',4*[2 1],'dpi',300,'margins',1.1,'printoptions','-depsc2');

%% nonh plot combined wave period effect
ch = get(s(3),'children');
ind =[];
for i=1:length(ch)
    if strcmpi(get(ch(i),'type'),'line')
        if all(get(ch(i),'color')==0)
            copyobj(ch(i),s(2));
        end
    elseif strcmpi(get(ch(i),'type'),'patch')
        ind(end+1) = i;
    end
end

for i=4:-1:1
    copyobj(ch(ind(i)),s(2));
end

ch = get(s(2),'children');
for i=1:length(ch)
    if strcmpi(get(ch(i),'type'),'hggroup')
        delete(ch(i));
    end
end
axes(s(2));
title('Wave period and water depth conditions');
delete(s(3));

fname = 'boundaries_NH2';
saveas(gcf,[fname '.fig']);
close
open([fname '.fig']);
savefigure(gcf,[fname '.png'],'size',4*[2 1],'dpi',900,'margins',1.1,'printoptions','-dpng');
savefigure(gcf,[fname '.eps'],'size',4*[2 1],'dpi',300,'margins',1.1,'printoptions','-depsc2');
close

%% nonhplus plot
figure;
s(1) = subplot(1,3,1);hold on;
levels = [1/1;1/2;1/3;1/4;1/5;0];
cmap = [red;red;ora;gre;gre;gre];
make_subplot_local(s(1),HD,DH,HiD,[1/1;1/2;1/3;1/4;1/5;0],[red;red;ora;gre;gre;gre],gre,{'1/1';'1/2';'1/3';'1/4';'1/5'});
title('Wave height to water depth ratio (H_s/h)');
xlabel('Significant wave height H_s (m)');
ylabel('Water depth (m)');
xt = 2:8;
yt = 0:5:40;
xl = xlim;yl=ylim;
pp =[];
for i=2:length(xt)-1;
    pp(end+1) = plot([xt(i) xt(i)],[yl(1) yl(2)],'k:');
end
for i=2:length(yt)-1;
    pp(end+1) = plot([xl(1) xl(2)],[yt(i) yt(i)],'k:');
end
set(pp,'color',sil);
set(s(1),'xtick',xt,'ytick',yt);

s(2) = subplot(1,3,2);hold on;
make_subplot_local(s(2),TD,DT,CgiC,[0.95:-0.05:0.6 0],[red;red;ora;yel;gre;gre;gre;gre;gre],gre);
title('Wave group velocity to wave celerity ratio (n = c_g/c)');
xlabel('Wave period T_{m-1,0} \approx T_p/1.1 (s)');
ylabel('Water depth (m)');
xt = 6:2:18;
yt = 0:5:40;
xl = xlim;yl=ylim;
pp =[];
for i=2:length(xt)-1;
    pp(end+1) = plot([xt(i) xt(i)],[yl(1) yl(2)],'k:');
end
for i=2:length(yt)-1;
    pp(end+1) = plot([xl(1) xl(2)],[yt(i) yt(i)],'k:');
end
set(pp,'color',sil);
set(s(2),'xtick',xt,'ytick',yt);

s(3) = subplot(1,3,3);hold on;
make_subplot_local(s(3),TD,DT,KH,[0.5 1 2 3 3.5 4 5],[gre;gre;gre;yel;ora;red;red],gre);
title('Relative water depth (kh)');
xlabel('Wave period T_{m-1,0} \approx T_p/1.1 (s)');
ylabel('Water depth (m)');
xt = 6:2:18;
yt = 0:5:40;
xl = xlim;yl=ylim;
pp =[];
for i=2:length(xt)-1;
    pp(end+1) = plot([xt(i) xt(i)],[yl(1) yl(2)],'k:');
end
for i=2:length(yt)-1;
    pp(end+1) = plot([xl(1) xl(2)],[yt(i) yt(i)],'k:');
end
set(pp,'color',sil);
set(s(3),'xtick',xt,'ytick',yt);

fname = 'boundaries_NHplus';
saveas(gcf,[fname '.fig']);
savefigure(gcf,[fname '.png'],'size',4*[3 1],'dpi',900,'margins',1.1,'printoptions','-dpng');
savefigure(gcf,[fname '.eps'],'size',4*[2 1],'dpi',300,'margins',1.1,'printoptions','-depsc2');

%% nonhplus plot combined wave period effect
ch = get(s(3),'children');
ind =[];
for i=1:length(ch)
    if strcmpi(get(ch(i),'type'),'line')
        if all(get(ch(i),'color')==0)
            copyobj(ch(i),s(2));
        end
    elseif strcmpi(get(ch(i),'type'),'patch')
        ind(end+1) = i;
    end
end

for i=4:-1:1
    copyobj(ch(ind(i)),s(2));
end

ch = get(s(2),'children');
for i=1:length(ch)
    if strcmpi(get(ch(i),'type'),'hggroup')
        delete(ch(i));
    end
end
axes(s(2));
title('Wave period and water depth conditions');
delete(s(3));

fname = 'boundaries_NHplus2';
saveas(gcf,[fname '.fig']);
close
open([fname '.fig']);
savefigure(gcf,[fname '.png'],'size',4*[2 1],'dpi',900,'margins',1.1,'printoptions','-dpng');
savefigure(gcf,[fname '.eps'],'size',4*[2 1],'dpi',300,'margins',1.1,'printoptions','-depsc2');
close

function make_subplot_local(axh,X,Y,Z,levels,colors,bgcolor,labels)

if~exist('labels','var')
    for i=1:length(levels)
        labels{i} = num2str(levels(i),'%0.2f');
    end
end

axes(axh);
holdstatus = get(axh,'NextPlot');
set(axh,'NextPlot','add');

% % [c,h] = contour(X,Y,Z,levels,'k');
% c = contourc(X(1,:),Y(:,1),Z,levels);

% [xp,yp,lp] = contour_split(c);
[c,h] = contourf(X,Y,Z,levels,'k');
ch = get(h,'children');
xp=[];
yp =[];
lp = [];
for i=1:length(ch)
    if strcmpi(get(ch(i),'type'),'patch')
        xp{end+1} = get(ch(i),'xdata');
        yp{end+1} = get(ch(i),'ydata');
        lp(end+1) = get(ch(i),'FaceVertexCData');
    end
end
% ordering is not correct for plotting patches, so correct
[~,~,lptemp] = contour_split(c,'nansplit',0);
[~,ind] = ismember(lptemp,lp);
% reorder
xp = xp(ind);
yp = yp(ind);
lp = lp(ind);

delete(h);



p = patch([min(X(:)) min(X(:)) max(X(:)) max(X(:))],[min(Y(:)) max(Y(:)) max(Y(:)) min(Y(:))],bgcolor);
set(p,'edgecolor','none');
for i=1:length(xp)
    xt = xp{i};
    yt = yp{i};
    rem = isnan(xt) | isnan(yt);
    xt(rem) = [];
    yt(rem) = [];
    
%     xt = [xt(1);xt(:);xt(end)];
%     yt = [min(Y(:));yt(:);min(Y(:))];
    p = patch(xt,yt,'b');
    ind = find(levels==lp(i));
    if ~isempty(ind)
        set(p,'facecolor',colors(ind,:),'edgecolor','none','facealpha',1);
    else
        set(p,'facecolor','none','edgecolor','none','facealpha',1);
    end
end
[c,h] = contour(X,Y,Z,levels,'k');
cl = clabel(c,h,'labelspacing',1000000);
for i=1:length(cl)
    str = get(cl(i),'String');
    strn = str2double(str);
    ind = find(abs(levels-strn)==min(abs(levels-strn)));
    set(cl(i),'String',labels{ind});
end
box on;
grid on;
axis square;
set(axh,'NextPlot',holdstatus);




