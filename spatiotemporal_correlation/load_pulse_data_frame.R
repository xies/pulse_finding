wt = read.csv('~/Dropbox (MIT)/wt.txt',header=FALSE)
control = read.csv('~/Dropbox (MIT)/Pulse exports/control.txt',header=FALSE)
twist = read.csv('~/Dropbox (MIT)/twist.txt',header=FALSE)

colnames(wt)[1] = 'near'; colnames(twist)[1] = 'near'; colnames(control)[1] = 'near';
colnames(wt)[2] = 'cr'; colnames(twist)[2] = 'cr'; colnames(control)[2] = 'cr';
colnames(wt)[3] = 'range'; colnames(twist)[3] = 'range'; colnames(control)[3] = 'range';
colnames(wt)[4] = 'behavior'; colnames(twist)[4] = 'behavior'; colnames(control)[4] = 'behavior';
colnames(wt)[5] = 'amplitude'; colnames(twist)[5] = 'amplitude'; colnames(control)[5] = 'amplitude';
# colnames(wt)[6] = 'persistence'; colnames(twist)[6] = 'persistence';

pulses = data.frame(near=wt$near,cr=wt$cr/7.27,behavior=wt$behavior,range=wt$range,
                genotype='wt',amplitude=wt$amplitude)
pulses = rbind(pulses, data.frame(near=twist$near,cr=twist$cr/7.27,behavior=twist$behavior,range=twist$range,
                          genotype='twist',amplitude=twist$amplitude) )
pulses = rbind(pulses, data.frame(near=control$near,cr=control$cr/7.27,behavior=control$behavior,range=control$range,
                                  genotype='control',amplitude=control$amplitude) )

###

P = NULL
thisP = NULL

for (bin in 10) {
  for (i in 1:2) {
    
    I = pulses$behavior == i & pulses$genotype == 'wt' & pulses$amplitude == bin;
  
    if (length(I[I]) > 6) {
#             rho =cor(pulses$near[I], pulses$range[I], alternative = 't', conf.level=0.95)
      rho = cor(pulses$cr[I],pulses$near[I])
#       rho = polyserial(pulses$cr[I], pulses$near[I], std.err = TRUE)
      bootstat_wt = boot(pulses[I,],f,1000)
      ci = boot.ci(bootstat_wt)
      thisP = data.frame(rho = rho,
#                          std = sqrt(diag(rho$var))[1],
                         cih = ci$basic[5], cil = ci$basic[4],
                         genotype = factor('wt'),
                         amplitude = factor(bin), behavior = factor(cluster_names[i]) )
      
      P = rbind(thisP,P)
    }
    else {
      thisP = data.frame(rho = NA, cil=NA, cih=NA,
                         genotype = factor('wt'),
                         amplitude = factor(bin), behavior = factor(cluster_names[i]) )
      P = rbind(thisP,P)
    }
    
    I = pulses$behavior == i & pulses$genotype == 'twist' & pulses$amplitude == bin;
    if (length(I[I]) > 6) {
      #         rho = cor.test(pulses$near[I], pulses$cr[I], alternative = 't', conf.level=0.95)
      rho = cor(pulses$range[I],pulses$near[I])
#       rho = polyserial(pulses$cr[I], pulses$near[I], std.err = TRUE)
      bootstat_twist = boot(pulses[I,],f,1000)
      ci = boot.ci(bootstat_twist)
      thisP = data.frame(rho = rho,
#                          std = sqrt(diag(rho$var))[1],
                         cih = ci$basic[5], cil = ci$basic[4],
                         genotype = factor('twist'),
                         amplitude = factor(bin), behavior = factor(cluster_names[i]) )
      P = rbind(thisP,P)
    }
    else {
      thisP = data.frame(rho = NA, cil = NA, cih = NA,
                         genotype = factor('twist'),
                         amplitude = factor(bin), behavior = factor(cluster_names[i]) )
      P = rbind(thisP,P)
    }
  }
}

P$amplitude = factor(P$amplitude, levels = rev(levels(P$amplitude)))
P$behavior = factor(P$behavior, levels = rev(levels(P$behavior)))

p = ggplot(data = subset(P, behavior != 'Unconstricting'),
           aes(colour=behavior, x=amplitude, y=rho))
limits = aes( ymax = cih, ymin = cil)
p = p + geom_point(size=5) 
p = p + geom_errorbar(limits, width=0.2)
p = p + geom_hline( yintercept = 0 )
p = p + theme(axis.text = element_text(size=20),
              title = element_text(size=20),
              strip.text.y = element_text(size = 24, colour = "blue"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 24)
)

p = p + guides(fill = guide_legend(reverse=TRUE))
p = p + facet_grid(genotype ~ behavior)
p

### Wald test for polyserial

wald_stat = 
  
  qchisq(.95, df=7)


###



df = pulses;
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, colour="red", geom=geom, width=0.2, ...)
}

p = ggplot(data = subset(df, behavior < 3 & amplitude < 6 & genotype == 'twist'),
           aes(x=near, y = cr) )
p = p + geom_point(size = 5, color = 'blue')
p = p + stat_sum_df( mean_cl_boot , geom = 'smooth',conf.int = 0.95)
p = p + facet_grid( behavior  ~ amplitude )
p = p + coord_cartesian(ylim=c(0, .8))
p = p + scale_x_continuous(breaks = seq(0 , 10, 1))
p = p + theme(axis.text = element_text(size=20),
              title = element_text(size=20),
              strip.text.y = element_text(size = 24, colour = "blue"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 24)
              )
p


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