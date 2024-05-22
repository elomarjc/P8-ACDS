xlabel('Longitude [\circ]')
ylabel('Latitude [\circ]')
%set(gca,'DataAspectRatio',[1 1 1])
%set(gca,'xtick',[-180 0 180])
set(gca,'xtick',[1  72 144 216 288],'XTickLabel',{'-180','-90','0','90','180'})
set(gca,'ytick',[1 45 90 135 180],'YTickLabel',{'-90','-45','0','45','90'})

