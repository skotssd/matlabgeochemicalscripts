function makecontourplot(F,em,ex)
 colormap('jet')
 h=surf(em,ex,F);
 shading interp
 hold on; [C,h]=contour3(em,ex,F,3,'k'); set(h,'linewidth',2); 
 set(gca,'linewidth',2,'fontsize',12)
 axis([min(em) max(em) min(ex) max(ex) 0 max(max(F))])
 view([0 90])
 hold on; plot3([250 600],[450 450],[0 0],'k','linewidth',2)
 hold on; plot3([600 600],[220 450],[0 0],'k','linewidth',2)
 h=xlabel('Emission (nm)'); set(h,'fontsize',12)
 h=ylabel('Excitation (nm)');  set(h,'fontsize',12)
 %add scale bar
 colorbar 
end