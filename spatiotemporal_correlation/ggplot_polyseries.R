wt = read.csv('~/Dropbox/wt.txt',header=FALSE)
twist = read.csv('~/Dropbox/twist.txt',header=FALSE)

colnames(wt)[1] = 'near'; colnames(twist)[1] = 'near';
colnames(wt)[2] = 'cr'; colnames(twist)[2] = 'cr';
colnames(wt)[3] = 'behavior'; colnames(twist)[3] = 'behavior';
colnames(wt)[4] = 'amplitude'; colnames(twist)[4] = 'amplitude';

df = data.frame(near=wt$near,cr=wt$cr,behavior=wt$behavior,
                genotype='wt',amplitude=wt$amplitude)
df = rbind(df, data.frame(near=twist$near,cr=twist$cr,behavior=twist$behavior,
                          genotype='twist',amplitude=twist$amplitude) )

###

P = NULL
thisP = NULL

for (bin in 1:10) {
  for (i in 1:3) {

    I = df$behavior == i & df$genotype == 'wt' & df$amplitude == bin;
    if (length(I[I]) > 4) {
      rho = polyserial(df$cr[I], df$near[I],std.err = TRUE)
      print('foo')
      thisP = data.frame(rho = rho$rho, var= rho$var, genotype = factor('wt'),
                       amplitude = factor(bin), behavior = factor(cluster_names[i]) )
      P = rbind(thisP,P)
    }
    
    I = df$behavior == i & df$genotype == 'twist' & df$amplitude == bin;
    if (length(I[I]) > 4) {
      rho = polyserial(df$cr[I], df$near[I], std.err = TRUE)
      thisP = data.frame(rho = rho$rho, var= rho$var, genotype = factor('twist'),
                         amplitude = factor(bin), behavior = factor(cluster_names[i]) )
      P = rbind(thisP,P)
    }
  }
}

p = ggplot(data = P, aes(colour=behavior, y=rho, x=amplitude))
limits = aes( ymax = rho + var, ymin = rho - var)
p = p + geom_point() + geom_errorbar(limits, width=0.2)
p = p + theme(axis.text = element_text(size=20),
                    title = element_text(size=20),
                    strip.text.y = element_text(size = 24, colour = "blue"),
                    legend.text = element_text(size = 20),
                    legend.title = element_text(size = 24)
                    axis.sides = "tr")
p = p + guides(fill = guide_legend(reverse=TRUE))
p = p + facet_grid(genotype ~ .)
# p + scale_x_discrete(limits = rev(levels(factor(c('twist','wt')))) )
p

##

p = ggplot(data = subset(df,amplitude == 10 & behavior == 1),
           aes(x=near,y = cr,color=amplitude) )
p = p + geom_point()
p = p + facet_grid( genotype + behavior ~ amplitude )
p


