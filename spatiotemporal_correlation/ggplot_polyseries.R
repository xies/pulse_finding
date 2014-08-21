P = NULL
thisP = NULL

for (i in 1:3) {
  rho = polyserial(df$cr[df$behavior == i & df$gentype == 'wt'],
             df$near[df$behavior == i & df$gentype == 'wt'], std.err = TRUE)
  thisP = data.frame(rho = rho$rho, var= rho$var, gentype = factor('wt'),
                     behavior = factor(cluster_names[i]))
  P = rbind(thisP,P)
  
  rho = polyserial(df$cr[df$behavior == i & df$gentype == 'twist'],
             df$near[df$behavior == i & df$gentype == 'twist'], std.err = TRUE)
  thisP = data.frame(rho = rho$rho, var= rho$var, gentype = factor('twist'),
                     behavior = factor(cluster_names[i]))
  P = rbind(thisP,P)
}

p = ggplot(data = P, aes(colour=behavior, y=rho, x=gentype))
limits = aes( ymax = rho + var, ymin = rho - var)
p = p + geom_point(size = 3) + geom_errorbar(limits, width=0.2)
p = p + theme(axis.text = element_text(size=20),
                    title = element_text(size=20),
                    strip.text.y = element_text(size = 24, colour = "blue"),
                    legend.text = element_text(size = 20),
                    legend.title = element_text(size = 24)
                    axis.sides = "tr")
p = p + guides(fill = guide_legend(reverse=TRUE))
p + scale_x_discrete(limits = rev(levels(factor(c('twist','wt')))) )
p

##

df = data.frame(near=wt$near,cr=wt$cr,behavior=wt$behavior,gentype='wt')
df = rbind(df, data.frame(near=twist$near,cr=twist$cr,behavior=twist$behavior,gentype='twist') )
p = ggplot(data = df,aes(x=near,y = cr, group = behavior, color = behavior))
p = p + geom_point()
p = p + facet_grid(gentype ~ .)
p

