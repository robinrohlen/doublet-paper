function [x,y]=plot_broken_axis(x,y,mrksize,vals_broken,upper_vals,tick_size,set_xlim,col,scat,broken_axis,symboltype)

if nargin<9
    scat=1;
end

if nargin<10
    broken_axis=1;
end

if nargin<11
    symboltype=0;
end

find_ind=find(y>vals_broken(1) & y<vals_broken(2));

x(find_ind)=[];
y(find_ind)=[];

find_ind=find(y>=vals_broken(2));

divide_factor=upper_vals(1)/(vals_broken(1)+tick_size);

y(find_ind)=y(find_ind)./divide_factor;

hold on;
if scat==1
    scatter(x,y,mrksize,'MarkerFaceColor',col,'MarkerEdgeColor','k','MarkerFaceAlpha',1.0,'MarkerEdgeAlpha',1.0);
    plot(x,y,':','Color',col);
else
    if symboltype==0
        scatter(x,y,mrksize,'MarkerFaceColor',col,'MarkerEdgeColor','k','MarkerFaceAlpha',1.0,'MarkerEdgeAlpha',1.0);
    else
        plot(x,y,'-','Color',col,'LineWidth',3,'MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',mrksize);
    end
end
hold off;

xlim(set_xlim);

num_ticks=0:tick_size:(upper_vals(2)/divide_factor);
yticks(num_ticks);
ylim([0 upper_vals(2)/divide_factor]);
get_yticklabels=get(gca,'YTickLabels');

for i=find(num_ticks>=vals_broken(1)+tick_size)
    get_yticklabels{i}=num2str(round(num_ticks(i)*divide_factor,-1));%num2str(round((num_ticks(i)-add_diff-tick_size)*divide_factor,-1));
end

set(gca,'YTickLabels',get_yticklabels);

if broken_axis==1
    xtick=get(gca,'XTick');
    t1=text(set_xlim(1),vals_broken(1)+tick_size/2+0.5,'//','FontSize',16);
    set(t1,'rotation',270);
end

hold on;
plot(set_xlim,[vals_broken(1)+tick_size/2 vals_broken(1)+tick_size/2],'k:');
hold off;

g=gcf;
g.Renderer='painters';
