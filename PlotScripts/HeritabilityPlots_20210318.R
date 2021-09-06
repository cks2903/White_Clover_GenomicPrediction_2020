library(ggplot2)
library(wesanderson)
library(ggpubr)

# full model
df <- data.frame(x=c(28.3,22.4,30.9,23.3,13.3,1.5,2.9,1.2,1.8,3.7,0.9,1.1,0.8,1.6,0.9), y=c(5,4,3,2,1,5,4,3,2,1,5,4,3,2,1))
df$leftI <- c(21.9,16.8,24.1,17.8,9.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
df$RightI <- c(34.6,28.4,37.6,29.0,18.0,3.1,5.3,2.7,3.7,6.8,2.0,2.9,2.2,4.1,2.3)
df$group=c("gpd","gpdCor","gpi","gpiCor","iSize","gpd","gpdCor","gpi","gpiCor","iSize","gpd","gpdCor","gpi","gpiCor","iSize")
df$group=as.factor(df$group)
df$type=c(rep("Clover",5),rep("Rhizobium",5),rep("Interaction",5))

head(df)

df$group=factor(df$group,levels=c("iSize","gpiCor","gpi","gpdCor","gpd"))
df$type=factor(df$type,levels=c("Clover","Rhizobium","Interaction"))


p <- ggplot(df, aes(x=x,y=group,colour= group)) +
        geom_point(size=3) +   
        geom_errorbarh(aes(xmax = RightI, xmin = leftI),height=0.001,size=1) +
  scale_color_manual(values=c(wes_palette("Rushmore1")[4],"cadetblue","brown4",wes_palette("Rushmore1")[2],wes_palette("Rushmore1")[3])) +
  
  labs( x="% variance explained") +
  facet_wrap(~type, ncol = 1) +
  theme(axis.title.y=element_blank(),axis.ticks.y=element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

p


# heritability, n=145
h2_df=data.frame(x=c(0.68,0.23,0.41,0.31,0.83),y=c(5,4,3,2,1))
h2_df$leftI <- c(0.32,0.00,0.00,0.25,0.57)
h2_df$RightI <- c(1.00,0.60,0.82,1.00,1.00)
h2_df$group=c("gpd","gpdCor","gpi","gpiCor","iSize")
h2_df$group=as.factor(h2_df$group)


h2_df$group=factor(h2_df$group,levels=c("iSize","gpiCor","gpi","gpdCor","gpd"))

p2 <- ggplot(h2_df, aes(x=x,y=group,colour= group)) +
  geom_point(size=3) +   
  geom_errorbarh(aes(xmax = RightI, xmin = leftI),height=0.001,size=1) +
  scale_color_manual(values=c(wes_palette("Rushmore1")[4],"cadetblue","brown4",wes_palette("Rushmore1")[2],wes_palette("Rushmore1")[3])) +
  theme(axis.title.y=element_blank(),
                axis.ticks.y=element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  labs( x="h2")  
p2

combined=ggarrange(p, p2,common.legend = TRUE,
          ncol = 2, nrow = 1,align="v")

combined

ggsave("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/variance_her_20210903.png", plot=combined,width =20, height =9,units = "cm"
       )
ggsave("/Volumes/NAT_MBG-PMg/Cathrine/Nchain/Genomic_prediction_yield_July2020/V2_LessHarsh_SaraQualityFilter/Article_GP_gpd/variance_her_20210903.pdf", plot=combined,width =20, height =9,units = "cm"
)
