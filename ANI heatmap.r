```{r ANI heatmap}

#load ANI data
Saureus_ANI <- read.csv("fastani_phage.csv",header = TRUE, sep = ",", encoding = "UTF-8") #check.names = FALSE)

Saureus_ANI$Query = factor(Saureus_ANI$Query,
                           levels = unique(Saureus_ANI$Query))
Saureus_ANI$ANI

Isolates<-unique(Saureus_ANI$Query)


library(tidyr)
Saureus_ANI_RESULTS<-spread(Saureus_ANI, Query, ANI)
#Encoding(rownames(Saureus_ANI_RESULTS)) <- "UTF-8"
rownames(Saureus_ANI_RESULTS) <- Saureus_ANI_RESULTS[,1]
Saureus_ANI_matrix <- as.matrix(Saureus_ANI_RESULTS[,-1])
print(rownames(Saureus_ANI_matrix))

library(pheatmap)
pheatmap(Saureus_ANI_matrix)

library(dplyr)
library(reshape)

#Extra Colours = "#88CCEE", "#4897DC", "#004488", "#01539C", "#537DA2", 

library(viridis)
pheatmap(Saureus_ANI_matrix,
         color= colorRampPalette(c( "#DADADA", "#A4C3DE","#66A2D8", "#1F68A9", "#0B4171", "#222255"))(6),
         breaks=c(90, 92, 94, 96, 98, 100),
         legend_breaks = c(90,92, 94, 96, 98, 100),
         legend=TRUE,
         fontsize = 14, 
         fontsize_row = 14, 
         fontsize_col = 14,
         filename ="Phage_heatmap.jpg")
         dev.off()

```