function [] = save_fig_and_crop(ps)


%%%% save the figure in tightly bound pdf

ratio =  ps(4) / ps(3)
paperWidth = 15;
paperHeight = paperWidth*ratio;


set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth+0.5 paperHeight+0.5]);


print(gcf, '-dpdf', '../../figures/Tails_Sch_1.pdf');
!pdfcrop ../../figures/Tails_Sch_1.pdf ../../figures/Tails_Sch_1.pdf


return