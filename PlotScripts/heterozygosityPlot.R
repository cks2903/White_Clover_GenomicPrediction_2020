library(ggplot2)
library("lattice")
library(viridis)

setwd("/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/data/")
data <- read.table("heterozygosity_per_individual.csv", header = TRUE, sep = ",")


p <- ggplot(data, aes(heterozygisty)) +
  geom_histogram(aes(fill=..count..), bins = 50) +
  scale_x_continuous(limits = c(0, 0.4), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 28), expand = c(0, 0)) +
  scale_fill_viridis(option = "magma", direction = -1) +
  ylab("Frequency") +
  xlab("Heterozygosity") +
  theme_classic() +
  theme(aspect.ratio=1,
        legend.position="none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),)

ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/Clover_GP_paper/figures/heterozygosity.pdf', plot = p, width = 15, height = 12, unit = 'cm')

