### Load pulse data (csv) into dataframe

wt = read.csv('~/Dropbox (MIT)/Pulse export/wt_full.csv',header=TRUE)
control = read.csv('~/Dropbox (MIT)/Pulse export/control.txt',header=TRUE)
twist = read.csv('~/Dropbox (MIT)/Pulse export/twist_full.csv',header=TRUE)

behavior_names = c('Ratcheted','Unratcheted','Unconstricting');

colnames(wt)[1] = 'near'; colnames(twist)[1] = 'near'; colnames(control)[1] = 'near';
colnames(wt)[2] = 'cr'; colnames(twist)[2] = 'cr'; colnames(control)[2] = 'cr';
colnames(wt)[3] = 'range'; colnames(twist)[3] = 'range'; colnames(control)[3] = 'range';
colnames(wt)[4] = 'behavior'; colnames(twist)[4] = 'behavior'; colnames(control)[4] = 'behavior';
colnames(wt)[5] = 'amplitude'; colnames(twist)[5] = 'amplitude'; colnames(control)[5] = 'amplitude';
colnames(wt)[6] = 'center';

pulses = data.frame(near=wt$near,cr=wt$cr,behavior=wt$behavior,range=wt$range,
                genotype='wt',amplitude=wt$amplitude,center=wt$center)
pulses = rbind(pulses, data.frame(near=twist$near,cr=twist$cr,behavior=twist$behavior,range=twist$range,
                          genotype='twist',amplitude=twist$amplitude,center=twist$center) )
# pulses = rbind(pulses, data.frame(near=control$near,cr=control$cr,behavior=control$behavior,range=control$range,
#                                   genotype='control',amplitude=control$amplitude) )

### Plot raw relationships

df = pulses;
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, colour="red", geom=geom, width=0.2, ...)
}

p = ggplot(data = subset(df, behavior < 3 & amplitude < 6 & genotype != 'control'),
           aes(x=near, y = cr) )
p = p + geom_point(size = 5, color = 'blue')
p = p + stat_sum_df( mean_cl_boot , geom = 'smooth',conf.int = 0.95)
p = p + facet_grid( genotype  ~ amplitude )
p = p + coord_cartesian(ylim=c(0, 0.82))
p = p + scale_x_continuous(breaks = seq(0 , 10, 1))
p = p + theme(axis.text = element_text(size=20),
              title = element_text(size=20),
              strip.text.y = element_text(size = 24, colour = "blue"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 24)
)
p
