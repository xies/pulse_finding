####

library(ppcor)

CORRdf = NULL
thisCORR = NULL

for (geno in c('wt','twist')) {
  for (beha in 1:2) {
    
    I = pulses$behavior == beha & pulses$genotype == geno & pulses$amplitude > 6;
    df = pulses[I,]; df$genotype = NULL;
    
    prho = pcor.test(df$near,df$cr,df$amplitude)
    phi = qnorm(0.0001360578/2);
    error = sqrt( phi^2 / (prho$n-3 + phi^2) );

    thisCORR = data.frame(
      prho = prho$estimate,
      pvalue = prho$p.value,
      error = error,
      genotype = geno, behavior = behavior_names[beha]
    )
    
    CORRdf = rbind(thisCORR,CORRdf)
  }
}

## ggplot rho value

p = ggplot(data = CORRdf,
           aes(colour=genotype, x = behavior, y=prho))
# limits = aes( ymax = prho + error, ymin = prho - error)
p = p + geom_point(size=5) 
# p = p + geom_errorbar(limits, width=0.2)
p = p + theme(axis.text = element_text(size=20),
              title = element_text(size=20),
              strip.text.y = element_text(size = 24, colour = "blue"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 24)
)

p = p + guides(fill = guide_legend(reverse=TRUE))
p


## ggplot rho value

p = ggplot(data = CORRdf,
           aes(colour=genotype, x = behavior, y=prho))
p = p + geom_point(size=5)
p = p + geom_hline( yintercept = 0 )
p = p + theme(axis.text = element_text(size=20),
              title = element_text(size=20),
              strip.text.y = element_text(size = 24, colour = "blue"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 24)
)
p = p + facet_grid(genotype ~ .)
p = p + guides(fill = guide_legend(reverse=TRUE))
p

## ggplot p.value (plot -log(P))

p.pval = ggplot(data = CORRdf, aes(colour=behavior, x=behavior, y=-log10(pvalue)))
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
p.pval = p.pval + facet_grid(genotype ~ .)
p.pval

