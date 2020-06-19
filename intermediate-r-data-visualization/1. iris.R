#Read data in.
dat <- iris

#Plot using 'base' graphics.
plot(x = dat$Sepal.Length, y = dat$Petal.Length,
     xlab = "Sepal Length", ylab = "Petal Length")

#Install ggplot2.
library(ggplot2)

#Get details of ggplot2 package.
library(help = "ggplot2")  

#Plot using qplot. 
qplot(x = Sepal.Length, y = Petal.Length, data = dat, 
      xlab = "Sepal Length", ylab = "Petal Length")

qplot(x = Sepal.Length, y = Petal.Length, data = dat, 
      xlab = "Sepal Length", ylab = "Petal Length",
      color = Species)

#But the journals charge extra for color figures.
#Can we use shapes to distinguish species?
qplot(x = Sepal.Length, y = Petal.Length, data = dat, 
      xlab = "Sepal Length", ylab = "Petal Length", 
      shape = Species)

#Check the boxplots for all species simultaneously.
qplot(x = Species, y = Sepal.Length, data = dat, geom = "boxplot",
      ylab = "Sepal Length")


#----------------------------------------
#----------------------------------------
#Underlying grammar is not clear. We'll use ggplot2.
#Specify what goes on which axis. Specify the data.
ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length))

#Add geometrical representation of data.
ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length)) +
  geom_point()

ggplot(data = dat, aes(x = log10(Sepal.Length), y = Petal.Length)) +
  geom_point()

#Add axis labels.
ggplot(data = dat, aes(x = log10(Sepal.Length), y = Petal.Length)) +
  geom_point() +
  xlab("Sepal Length (log10)") +
  ylab("Petal Length")

#Explore options with geom_point.
ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length)) +
  geom_point(shape =1)

#Explore options with geom_point.
ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length)) +
  geom_point(aes(shape = Species))

ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length, shape = Species)) +
  geom_point()

#Explore other geometrical mapping options.
ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length)) +
  geom_line()

#Explore other geometrical mapping options.
ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_line()

#Add a trendline to data.
ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point() + geom_smooth(method = lm)

#Limit the axis.

ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point() + geom_smooth(method = lm) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10))

#Sepcify where the breaks should be.
ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point() + geom_smooth(method = lm) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10))+
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10))

ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point() + geom_smooth(method = lm) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10))+
  scale_x_continuous(breaks = 0:5*2)


#Move legend to top.
ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point() + geom_smooth(method = lm) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10))+
  scale_x_continuous(breaks = 0:5*2)+
  theme(legend.direction = "horizontal", legend.position = "top",
        legend.title = element_blank())+
  xlab("Sepal Length") +
  ylab("Petal Length")


#Add border to plot.
ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point() + geom_smooth(method = lm) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10))+
  scale_x_continuous(breaks = 0:5*2)+
  theme(legend.direction = "horizontal", legend.position = "top",
        legend.title = element_blank(),
        panel.border = element_rect(size = c(1,1,1,1), color = "black", 
                                    fill = NA))+
  xlab("Sepal Length") +
  ylab("Petal Length")


#Specifying limits with xlim instead of coord_cartesian. 
#Note the different behavior.
ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point() + geom_smooth(method = lm) +
  xlim(5, 10) + ylim(0, 10)+
  theme(legend.direction = "horizontal", legend.position = "top",
        legend.title = element_blank(),
        panel.border = element_rect(size = c(1,1,1,1), color = "black", 
                                    fill = NA))+
  xlab("Sepal Length") +
  ylab("Petal Length")

ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point() + geom_smooth(method = lm) +
  coord_cartesian(xlim = c(5, 10), ylim = c(0, 10)) +
  theme(legend.direction = "horizontal", legend.position = "top",
        legend.title = element_blank(),
        panel.border = element_rect(size = c(1,1,1,1), color = "black", 
                                    fill = NA))+
  xlab("Sepal Length") +
  ylab("Petal Length")
  
#----------------------------------------
#----------------------------------------
#Adjust more theme elements.
#Change axis text, axis label, label direction.

#Make figures for publication.
p1 <- ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point(size = 0.5) + geom_smooth(method = lm) +
  coord_cartesian(xlim = c(4, 8.5), ylim = c(0, 8))+
  scale_x_continuous(breaks = 0:5*2)+
  theme(legend.direction = "horizontal", legend.position = "top",
        legend.title = element_blank(),
        panel.border = element_rect(size = c(1,1,1,1), color = "black", 
                                    fill = NA))+
  xlab("Sepal Length") +
  ylab("Petal Length")

p2 <- ggplot(data= dat, aes(x = Species, y = Sepal.Length)) +
  geom_boxplot() + ylab("Sepal Length")  +
  theme(panel.border = element_rect(size = c(1,1,1,1), color = "black", 
                                    fill = NA))

#Package to arrange figures.
library(cowplot)
theme_set(theme_grey())

#Arrange plots.
p <- plot_grid(p1, p2, labels = c("a", "b"))

#Save plot to file.
ggsave("Iris.pdf", 
       plot = p,
       width = 17, height = 8, units = "cm")

#----------------------
#Plots with multiple facet panels.
#Facets row wise.
ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length)) +
  geom_point() + geom_smooth(method = lm) +
  facet_grid(Species~.)+
  coord_cartesian(xlim = c(1, 8), ylim = c(1, 8)) +
  theme(legend.direction = "horizontal", legend.position = "top",
        legend.title = element_blank(),
        panel.border = element_rect(size = c(1,1,1,1), color = "black", 
                                    fill = NA))+
  xlab("Sepal Length") +
  ylab("Petal Length")


#Facets column wise.
ggplot(data = dat, aes(x = Sepal.Length, y = Petal.Length)) +
  geom_point() + geom_smooth(method = lm) +
  facet_grid(.~Species)+
  coord_cartesian(xlim = c(1, 8), ylim = c(1, 8)) +
  theme(legend.direction = "horizontal", legend.position = "top",
        legend.title = element_blank(),
        panel.border = element_rect(size = c(1,1,1,1), color = "black", 
                                    fill = NA))+
  xlab("Sepal Length") +
  ylab("Petal Length")
