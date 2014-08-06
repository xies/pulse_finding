fits2plot = f[f$behavior != "N/A",]

p_xy = ggplot(data=fits2plot,aes(x = x, y = y))
p_xy = p_xy + stat_density2d( geom="tile", aes(fill = ..density..), contour = FALSE)
p_xy = p_xy + geom_point(color="magenta")
p_xy = p_xy + facet_grid(behavior ~ .)
p_xy = p_xy + labs( x = expression(
  paste("Centroid-x (",mu,"m)") ) )
p_xy = p_xy + labs( y = expression(
  paste("Centroid-y (",mu,"m)") ) )
p_xy = p_xy + theme(axis.text = element_text(size=18),
      title = element_text(size=18),
      strip.text.y = element_text(size = 20, colour = "blue"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 20) )
p_xy = p_xy + ggtitle('Pulse spatial locations')
p_xy

##
p_y = ggplot(data=fits2plot,aes(x=y,
                                color=behavior))
p_y = p_y + geom_density()
p_y = p_y + facet_grid(embryoID ~ .)
p_y = p_y + labs( x = expression(
  paste("Centroid-y (",mu,"m)") ) ) + ggtitle('Pulse location')
p_y = p_y + theme(axis.text = element_text(size=14),
                    title = element_text(size=18),
                    strip.text.y = element_text(size = 18, colour = "blue"),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 18) )
p_y

##
p_x = ggplot(data=fits2plot,aes(x=x,
                                color=behavior))
p_x = p_x + geom_density()
p_x = p_x + facet_grid(embryoID ~ .)
p_x = p_x + labs( x = expression(
  paste("Centroid-x (",mu,"m)") ) ) + ggtitle('Pulse location')
p_x = p_x + theme(axis.text = element_text(size=14),
                  title = element_text(size=18),
                  strip.text.y = element_text(size = 18, colour = "blue"),
                  legend.text = element_text(size = 18),
                  legend.title = element_text(size = 18) )
p_x

###

# Px = matrix(seq(1:9), 3)
# Py = matrix(seq(1:9), 3)

Pks = NULL
this_f = NULL

x = fits2plot$x;
y = fits2plot$y;
l = fits2plot$behavior;
for (i in 1:3) {
  for (j in 1:3) {
    x_ks_stat = ks.test(x[l == cluster_names[i]], x[l == cluster_names[j]])
    y_ks_stat = ks.test(y[l == cluster_names[i]], y[l == cluster_names[j]])
    thisf = data.frame( fitID = raw[,1],
                        x = raw[,2] - xref[embryoID], y = raw[,3] - yref[embryoID],
                        t = raw[,4])
    thisf = data.frame( behavior1 = cluster_names[i], behavior2 = cluster_names[j],
                        p_value = x_ks_stat$p.value,
                        direction = 'x')
    Pks = rbind(Pks,thisf)
    thisf = data.frame( behavior1 = cluster_names[i], behavior2 = cluster_names[j],
                        p_value = y_ks_stat$p.value,
                        direction = 'y')
    Pks = rbind(Pks,thisf)
    
  }
}

p_ks = ggplot(data = Pks, aes(x=behavior1,y=behavior2,fill=p_value))
p_ks = p_ks + geom_tile() + labs( x = "Pulse behavior", y = "Pulse behavior" )
p_ks = p_ks + theme(axis.text = element_text(size=14),
                  title = element_text(size=18),
                  strip.text.y = element_text(size = 18, colour = "blue"),
                  legend.text = element_text(size = 18),
                  legend.title = element_text(size = 18) )
p_ks = p_ks + facet_grid( direction ~ .)
p_ks
