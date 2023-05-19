function makesimplecontourplot(F,em,ex,num)
 [C,h]=contour(em,ex,F,num,'k'); set(h,'linewidth',2); 
 set(gca,'linewidth',2,'fontsize',12)
 axis([min(em) max(em) min(ex) max(ex)])
 h=xlabel('Emission (nm)'); set(h,'fontsize',12)
 h=ylabel('Excitation (nm)');  set(h,'fontsize',12)
end