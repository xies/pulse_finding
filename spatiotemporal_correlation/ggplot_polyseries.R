wt = read.csv('~/Dropbox/wt.txt',header=FALSE)
twist = read.csv('~/Dropbox/twist.txt',header=FALSE)

colnames(wt)[1] = 'near'; colnames(twist)[1] = 'near';
colnames(wt)[2] = 'cr'; colnames(twist)[2] = 'cr';
colnames(wt)[3] = 'range'; colnames(twist)[3] = 'range';
colnames(wt)[4] = 'behavior'; colnames(twist)[4] = 'behavior';
colnames(wt)[5] = 'amplitude'; colnames(twist)[5] = 'amplitude';
# colnames(wt)[6] = 'persistence'; colnames(twist)[6] = 'persistence';

pulses = data.frame(near=wt$near,cr=wt$cr,behavior=wt$behavior,range=wt$range,
                genotype='wt',amplitude=wt$amplitude)
pulses = rbind(pulses, data.frame(near=twist$near,cr=twist$cr,behavior=twist$behavior,range=twist$range,
                          genotype='twist',amplitude=twist$amplitude) )

###

P = NULL
thisP = NULL

for (bin in 1:10) {
  for (i in 1:3) {

    I = pulses$behavior == i & pulses$genotype == 'wt' & pulses$amplitude == bin;
    
    if (length(I[I]) > 4) {
#       rho = rcorr(pulses$near[I], pulses$cr[I],type = 'spearman')
      rho = polyserial(pulses$cr[I], pulses$near[I],std.err = TRUE)
      thisP = data.frame(rho = rho$rho, std = rho$var,
                         genotype = factor('wt'),
                         amplitude = factor(bin), behavior = factor(cluster_names[i]) )
      
      P = rbind(thisP,P)
    }
    else {
      thisP = data.frame(rho = NA, std = NA,
                         genotype = factor('wt'),
                         amplitude = factor(bin), behavior = factor(cluster_names[i]) )
      P = rbind(thisP,P)
    }
    
    I = pulses$behavior == i & pulses$genotype == 'twist' & pulses$amplitude == bin;
    if (length(I[I]) > 4) {
#             rho = rcorr(pulses$near[I], pulses$cr[I],type = 'spearman')
      rho = polyserial(pulses$cr[I], pulses$near[I], std.err = TRUE)
      thisP = data.frame(rho = rho$rho, std = rho$var,
                         genotype = factor('twist'),
                         amplitude = factor(bin), behavior = factor(cluster_names[i]) )
      P = rbind(thisP,P)
    }
    else {
      thisP = data.frame(rho = NA, std = NA,
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
limits = aes( ymax = rho + std, ymin = rho - std)
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

###

df = pulses;
df$near = df$near;

stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, colour="red", geom=geom, width=0.2, ...)
}

p = ggplot(data = subset(df, behavior < 3 & amplitude == 10 & genotype == 'wt'),
           aes(x=near, y = cr) )
           
p = p + geom_point(size = 5, color = 'blue')
p = p + stat_sum_df( mean_cl_boot , geom = 'smooth',conf.int = 0.95)

p = p + facet_grid( behavior  ~ . )
p = p + coord_cartesian(ylim=c(0, 7.5))
p = p + scale_x_continuous(breaks = seq(0 , 10, 1))
p = p + theme(axis.text = element_text(size=20),
              title = element_text(size=20),
              strip.text.y = element_text(size = 24, colour = "blue"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 24)
              )
p

