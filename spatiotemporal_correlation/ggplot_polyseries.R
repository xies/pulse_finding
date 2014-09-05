wt = read.csv('~/Dropbox/wt.txt',header=FALSE)
twist = read.csv('~/Dropbox/twist.txt',header=FALSE)

colnames(wt)[1] = 'near'; colnames(twist)[1] = 'near';
colnames(wt)[2] = 'cr'; colnames(twist)[2] = 'cr';
colnames(wt)[3] = 'behavior'; colnames(twist)[3] = 'behavior';
colnames(wt)[4] = 'amplitude'; colnames(twist)[4] = 'amplitude';
colnames(wt)[5] = 'persistence'; colnames(twist)[5] = 'persistence';

pulses = data.frame(near=wt$near,cr=wt$cr,behavior=wt$behavior,
                genotype='wt',amplitude=wt$amplitude)
pulses = rbind(pulses, data.frame(near=twist$near,cr=twist$cr,behavior=twist$behavior,
                          genotype='twist',amplitude=twist$amplitude) )

###

P = NULL
thisP = NULL

for (bin in 1:10) {
  for (i in 1:3) {

    I = pulses$behavior == i & pulses$genotype == 'wt' & pulses$amplitude == bin;
    
    if (length(I[I]) > 4) {
      rho = polyserial(pulses$cr[I], pulses$near[I],std.err = TRUE)
      thisP = data.frame(rho = rho$rho, var= rho$var, genotype = factor('wt'),
                       amplitude = factor(bin), behavior = factor(cluster_names[i]) )
      P = rbind(thisP,P)
    }
    else {
      thisP = data.frame(rho = NA, var= NA, genotype = factor('wt'),
                         amplitude = factor(bin), behavior = factor(cluster_names[i]) )
      P = rbind(thisP,P)
    }
    
    I = pulses$behavior == i & pulses$genotype == 'twist' & pulses$amplitude == bin;
    if (length(I[I]) > 4) {
      rho = polyserial(pulses$cr[I], pulses$near[I], std.err = TRUE)
      thisP = data.frame(rho = rho$rho, var= rho$var, genotype = factor('twist'),
                         amplitude = factor(bin), behavior = factor(cluster_names[i]) )
      P = rbind(thisP,P)
    }
    else {
      thisP = data.frame(rho = NA, var= NA, genotype = factor('twist'),
                         amplitude = factor(bin), behavior = factor(cluster_names[i]) )
      P = rbind(thisP,P)
    }
  }
}

P$amplitude = factor(P$amplitude, levels = rev(levels(P$amplitude)))
P$behavior = factor(P$behavior, levels = rev(levels(P$behavior)))

p = ggplot(data = subset(P, behavior != 'Unconstricting'),
           aes(colour=behavior, x=amplitude, y=rho))
limits = aes( ymax = rho + var, ymin = rho - var)
p = p + geom_point(size=5) + geom_errorbar(limits, width=0.2)
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
df$near = factor(df$near);

p = ggplot(data = subset(df, behavior < 3 & amplitude > 6),
           aes(x=near,y = cr,color = amplitude) )
p = p + geom_point(size = 5)
p = p + facet_grid( genotype  ~ . )
p = p + theme(axis.text = element_text(size=20),
              title = element_text(size=20),
              strip.text.y = element_text(size = 24, colour = "blue"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 24)
              )
p



