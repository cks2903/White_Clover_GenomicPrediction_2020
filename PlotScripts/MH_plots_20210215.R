#########################################
#########################################
# Make a Manhattan plot of any .pvals file
# inputs are:
# 1) .pvals file
# 2) name of trait
#########################################
#########################################

library(ggplot2)
library(dplyr)
args=commandArgs(trailingOnly = TRUE)

# load file

data <- read.table(args[1], header=TRUE, sep=",")
datafiltered1=data[-which(data$mafs<0.05),]

trait=args[2]

# Plotting
{
  #prepare for plotting
  datafiltered <- subset(datafiltered1, chromosomes >0)
  datafiltered$p.bh <- p.adjust(datafiltered$scores, "BH")
  p.bonferroni <- 1/nrow(datafiltered)*0.05
  
  
  don <- datafiltered %>% 
    
    # Compute chromosome size
    group_by(chromosomes) %>% 
    summarise(chr_len=max(positions)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(datafiltered, ., by=c("chromosomes"="chromosomes")) %>%
    
    # Add a cumulative position of each SNP
    arrange(chromosomes, positions) %>%
    mutate( BPcum=positions+tot) %>%
    
    mutate( is_highlight=ifelse(p.bh<0.05, "yes", "no")) 
  
  
  # Prepare X axis 
  axisdf = don %>% group_by(chromosomes) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  
  
  # plotting
  p <- ggplot(don, aes(x=BPcum, y=-log10(scores))) +
    
    # Show all points
    geom_point( aes(color=as.factor(chromosomes)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(c("gray33", "gray67"), 16 )) +
    
    # custom X axis:
    scale_x_continuous(label = axisdf$chromosomes, breaks= axisdf$center, expand = c(0, 0)) +
    #scale_y_continuous(limits = c(0,7.5), expand = expand_scale(mult = c(0, .1))) +     # remove space between plot area and x axis
    ylim(0,7.5) +
    # Add thresholds
    geom_hline(yintercept=-log10(p.bonferroni), linetype="dashed", color = "#f04546") + 
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
     panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 22),
      title = element_text(size = 20)
    ) +
    labs(y = ("-log10(p-value)"), x = ("Chromosome"), title = (trait))
  
  p
}
ggsave(paste('ManhattanPlot_',trait,Sys.Date(),'.pdf',sep="_"), plot = p, width = 75, height = 15, unit = 'cm')
ggsave(paste('ManhattanPlot_',trait,Sys.Date(),'.png',sep="_"), plot = p, width = 75, height = 15, unit = 'cm')



