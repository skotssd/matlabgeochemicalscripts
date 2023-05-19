function makesurfaceplot(F,em,ex)
 colormap('jet'); shading interp
 h=surf(em,ex,F);
 set(gca,'linewidth',2)
 view([-26 48])
 axis([min(em) max(em) min(ex) max(ex) 0 max(max(F))*1.1])
 h=xlabel('Emission (nm)'); set(h,'fontsize',12)
 h=ylabel('Excitation (nm)');  set(h,'fontsize',12)
 h=zlabel('Intensity (arb.)'); set(h,'fontsize',12)
end