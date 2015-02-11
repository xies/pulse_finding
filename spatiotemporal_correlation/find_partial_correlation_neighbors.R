####

CORRdf = NULL
thisCORR = NULL

for (geno in c('wt','twist')) {
  for (beha in 1:2) {
    
    I = pulses$behavior == beha & pulses$genotype == geno & pulses$amplitude > 6;
    df = pulses[I,]; df$genotype = NULL;
    
    partial_rho = pcor(c('near','cr', 'amplitude'), var(df))
    prho_pvalue = pcor.test(partial_rho, 1, length(I[I]))
    
    thisCORR = data.frame(
      prho = partial_rho,
      pvalue = prho_pvalue$pvalue,
      genotype = geno, behavior = beha
    )
    
    CORRdf = rbind(thisCORR,CORRdf)
  }
}

p = ggplot(data = CORRdf,
           aes(colour=genotype, x = behavior, y=prho))
p = p + geom_point(size=5)
p = p + theme(axis.text = element_text(size=20),
              title = element_text(size=20),
              strip.text.y = element_text(size = 24, colour = "blue"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 24)
)

p = p + guides(fill = guide_legend(reverse=TRUE))
p