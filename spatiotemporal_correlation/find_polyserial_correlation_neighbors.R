library(polycor)

corfun = function(d) { polyserial(d$near,d$cr, std.err=TRUE) }
bootfun = function(d,i) { d2 = d[i,]; return(polyserial(d2$near,d2$cr, std.err=FALSE)) }

PCdf = NULL
thisPC = NULL

for (geno in c('twist','wt')) {
  for (bin in 6:10) {
    for (i in 1:2) {
      
      I = pulses$behavior == i & pulses$genotype == geno & pulses$amplitude == bin;
      
      if (length(I[I]) > 6) {

        rho = corfun(pulses[I,]);
        #       rho = polyserial(pulses$cr[I], pulses$near[I], std.err = TRUE)
#         bootstat = boot(pulses[I,], bootfun, 1000);
#         ci = boot.ci(bootstat);
        
        r = rho$rho;
        std = sqrt(diag( rho$var ));
        
        thisPC = data.frame(rho = r,
                            std = std,
                            wald_p = 1-pchisq((r/std)^2,df=1),
#                            cih = ci$basic[5], cil = ci$basic[4],
                            genotype = factor(geno),
                            amplitude = factor(bin),
                            behavior = factor(cluster_names[i]) );
        PCdf = rbind(thisPC,PCdf);
      }
      else {
        thisP = data.frame(rho = NA, std = NA,
                           genotype = factor(geno),
                           amplitude = factor(bin),
                           behavior = factor(cluster_names[i]) );
        PCdf = rbind(thisPC,PCdf);
      }
    }
  }
}

PCdf$amplitude = factor(PCdf$amplitude, levels = rev(levels(PCdf$amplitude)) )
PCdf$behavior = factor(PCdf$behavior, levels = rev(levels(PCdf$behavior)) )

## ggplot R value

library(ggplot2)

p.rho = ggplot(data = PCdf, aes(colour=behavior, x=amplitude, y=rho))
limits = aes( ymax = rho + std, ymin = rho - std)
p.rho = p.rho + geom_point(size=5) 
p.rho = p.rho + geom_errorbar(limits, width=0.2)
p.rho = p.rho + geom_hline( yintercept = 0 )
p.rho = p.rho + theme(axis.text = element_text(size=20),
                        title = element_text(size=20),
                        strip.text.y = element_text(size = 24, colour = "blue"),
                        legend.text = element_text(size = 20),
                        legend.title = element_text(size = 24)
)

p.rho = p.rho + guides(fill = guide_legend(reverse=TRUE))
p.rho = p.rho + facet_grid(genotype ~ behavior)
p.rho

## ggplot p.value (plot -log(P))

p.pval = ggplot(data = PCdf, aes(colour=behavior, x=amplitude, y=-log10(wald_p)))
# limits = aes( ymax = cih, ymin = cil)
p.pval = p.pval + geom_point(size=5) 
# p.pval = p.pval + geom_errorbar(limits, width=0.2)
p.pval = p.pval + geom_hline( yintercept = -log10(0.05) )
p.pval = p.pval + theme(axis.text = element_text(size=20),
              title = element_text(size=20),
              strip.text.y = element_text(size = 24, colour = "blue"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 24)
)

p.pval = p.pval + guides(fill = guide_legend(reverse=TRUE))
p.pval = p.pval + facet_grid(genotype ~ behavior)
p.pval

### Wald test for polyserial
