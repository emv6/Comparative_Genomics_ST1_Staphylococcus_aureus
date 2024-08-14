---
title: "Virulence Heatmaps"
author: "Emma_Voss"
date: "2024-04-14"
---
library(reshape)
library(dplyr)
library(pheatmap)
library(viridis)

```{r VFDB Heatmap Adherence}

#Next step is using python code as heatmap is too large so collapsing heat map by matches of lines.Input the formatted combined matrix 
adherence_Heatmap <- read.csv("Virulence/adherencegroups.csv", fill = TRUE, header = TRUE, sep = ",",check.names=FALSE)
format_adherence_Heatmap_metadata <- adherence_Heatmap %>% select(Group, Host) #Saving the metadata into a separate data frame 
rownames(adherence_Heatmap) <- adherence_Heatmap$Group #Saying Group is the row names
rownames(format_adherence_Heatmap_metadata) <- format_adherence_Heatmap_metadata$Group #Saying Group is the row names
format_adherence_Heatmap_metadata <- format_adherence_Heatmap_metadata[-c(1)] #Deleting the column 1 which is group 
pformat_adherence_Heatmap_matrix <- as.matrix(adherence_Heatmap[,-c(1:2)]) #Deleting the Group, Source and Human from the format matrix 
rownames(pformat_adherence_Heatmap_matrix) <- adherence_Heatmap[,1] #Deleting the column 1 which is group 
adherence_Heatmap <- mapply(pformat_adherence_Heatmap_matrix, FUN=as.numeric)
adherence_Heatmap <- matrix(data=adherence_Heatmap, ncol=length(colnames(pformat_adherence_Heatmap_matrix)), nrow=length(row.names(pformat_adherence_Heatmap_matrix)))
row.names(adherence_Heatmap) <- row.names(pformat_adherence_Heatmap_matrix)
colnames(adherence_Heatmap) <- colnames(pformat_adherence_Heatmap_matrix)


pheatmap(adherence_Heatmap)


#test heatmap
pheatmap(adherence_Heatmap, annotation_row = format_adherence_Heatmap_metadata, labels_row = format_adherence_Heatmap_metadata$Group, cluster_rows = F)



#set color palette
Host_colour = list(
  Host = c("Human"="#EB1A17", "Bovine"="#38AD23"))
#Legend
Legend_labels <- c("Absence", "Presence")


Adherence <- pheatmap(adherence_Heatmap,
         color=c("red", "blue"),
         legend_breaks = c(0, 1),
         annotation_colors = Host_colour,
         legend_labels = Legend_labels,
         annotation_row = format_adherence_Heatmap_metadata, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         show_rownames = TRUE,
         fontsize = 12, 
         fontsize_row = 10, 
         fontsize_col = 10,
         cellwidth = 30,
         angle_col = 45,
         cellheight = 30,
         width=13,
         height=14,
         main = "Virulence Genes: Adherence",
         filename ="FigS5.jpg")

```

```{r VFDB Heatmap Enterotoxin}
#Next step is using python code as heatmap is too large so collapsing heat map by matches of lines.Input the formatted combined matrix 
Enterotoxin_Heatmap <- read.csv("Virulence/enterotoxingroups.csv", fill = TRUE, header = TRUE, sep = ",",check.names=FALSE)
format_Enterotoxin_Heatmap_metadata <- Enterotoxin_Heatmap %>% select(Group, Host) #Saving the metadata into a separate data frame 
rownames(Enterotoxin_Heatmap) <- Enterotoxin_Heatmap$Group #Saying Group is the row names
rownames(format_Enterotoxin_Heatmap_metadata) <- format_Enterotoxin_Heatmap_metadata$Group #Saying Group is the row names
format_Enterotoxin_Heatmap_metadata <- format_Enterotoxin_Heatmap_metadata[-c(1)] #Deleting the column 1 which is group 
pformat_Enterotoxin_Heatmap_matrix <- as.matrix(Enterotoxin_Heatmap[,-c(1:2)]) #Deleting the Group, Source and Human from the format matrix 
rownames(pformat_Enterotoxin_Heatmap_matrix) <- Enterotoxin_Heatmap[,1] #Deleting the column 1 which is group 
Enterotoxin_Heatmap <- mapply(pformat_Enterotoxin_Heatmap_matrix, FUN=as.numeric)
Enterotoxin_Heatmap <- matrix(data=Enterotoxin_Heatmap, ncol=length(colnames(pformat_Enterotoxin_Heatmap_matrix)), nrow=length(row.names(pformat_Enterotoxin_Heatmap_matrix)))
row.names(Enterotoxin_Heatmap) <- row.names(pformat_Enterotoxin_Heatmap_matrix)
colnames(Enterotoxin_Heatmap) <- colnames(pformat_Enterotoxin_Heatmap_matrix)


pheatmap(Enterotoxin_Heatmap)


#test heatmap
pheatmap(Enterotoxin_Heatmap, annotation_row = format_Enterotoxin_Heatmap_metadata, labels_row = format_Enterotoxin_Heatmap_metadata$Group, cluster_rows = F)



#set color palette
Host_colour = list(
  Host = c("Human"="#EB1A17", "Bovine"="#38AD23"))
#Legend
Legend_labels <- c("Absence", "Presence")

Enterotoxin <- pheatmap(Enterotoxin_Heatmap,
         color=c("red", "blue"),
         legend_breaks = c(0, 1),
         annotation_colors = Host_colour,
         legend_labels = Legend_labels,
         annotation_row = format_Enterotoxin_Heatmap_metadata, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         show_rownames = TRUE,
         fontsize = 10, 
         fontsize_row = 10, 
         fontsize_col = 8,
         cellwidth = 30,
         angle_col = 45,
         cellheight = 30,
         width=8,
         height=6,
         main = "Virulence Genes: Enterotoxins", 
         filename ="FigS6.jpg")

```

```{r VFDB Heatmap Exoenzyme}

#Next step is using python code as heatmap is too large so collapsing heat map by matches of lines.Input the formatted combined matrix 
Exoenzyme_Heatmap <- read.csv("Virulence/exoenzymegroups.csv", fill = TRUE, header = TRUE, sep = ",",check.names=FALSE)
format_Exoenzyme_Heatmap_metadata <- Exoenzyme_Heatmap %>% select(Group, Host) #Saving the metadata into a separate data frame 
rownames(Exoenzyme_Heatmap) <- Exoenzyme_Heatmap$Group #Saying Group is the row names
rownames(format_Exoenzyme_Heatmap_metadata) <- format_Exoenzyme_Heatmap_metadata$Group #Saying Group is the row names
format_Exoenzyme_Heatmap_metadata <- format_Exoenzyme_Heatmap_metadata[-c(1)] #Deleting the column 1 which is group 
pformat_Exoenzyme_Heatmap_matrix <- as.matrix(Exoenzyme_Heatmap[,-c(1:2)]) #Deleting the Group, Source and Human from the format matrix 
rownames(pformat_Exoenzyme_Heatmap_matrix) <- Exoenzyme_Heatmap[,1] #Deleting the column 1 which is group 
Exoenzyme_Heatmap <- mapply(pformat_Exoenzyme_Heatmap_matrix, FUN=as.numeric)
Exoenzyme_Heatmap <- matrix(data=Exoenzyme_Heatmap, ncol=length(colnames(pformat_Exoenzyme_Heatmap_matrix)), nrow=length(row.names(pformat_Exoenzyme_Heatmap_matrix)))
row.names(Exoenzyme_Heatmap) <- row.names(pformat_Exoenzyme_Heatmap_matrix)
colnames(Exoenzyme_Heatmap) <- colnames(pformat_Exoenzyme_Heatmap_matrix)


pheatmap(Exoenzyme_Heatmap)


#test heatmap
pheatmap(Exoenzyme_Heatmap, annotation_row = format_Exoenzyme_Heatmap_metadata, labels_row = format_Exoenzyme_Heatmap_metadata$Group, cluster_rows = F)



#set color palette
Host_colour = list(
  Host = c("Human"="#EB1A17", "Bovine"="#38AD23"))
#Legend
Legend_labels <- c("Absence", "Presence")

Exoenzyme <- pheatmap(Exoenzyme_Heatmap,
         color=c("red", "blue"),
         legend_breaks = c(0, 1),
         annotation_colors = Host_colour,
         legend_labels = Legend_labels,
         annotation_row = format_Exoenzyme_Heatmap_metadata, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         show_rownames = TRUE,
         fontsize = 10, 
         fontsize_row = 10, 
         fontsize_col = 8,
         cellwidth = 30,
         angle_col = 45,
         cellheight = 30,
         width=12,
         height=8,
         main = "Virulence Genes: Exoenzymes",
         filename ="FigS7.jpg")

```

```{r VFDB Heatmap Exotoxin}

#Next step is using python code as heatmap is too large so collapsing heat map by matches of lines.Input the formatted combined matrix 
Exotoxin_Heatmap <- read.csv("Virulence/exotoxingroups.csv", fill = TRUE, header = TRUE, sep = ",",check.names=FALSE)
format_Exotoxin_Heatmap_metadata <- Exotoxin_Heatmap %>% select(Group, Host) #Saving the metadata into a separate data frame 
rownames(Exotoxin_Heatmap) <- Exotoxin_Heatmap$Group #Saying Group is the row names
rownames(format_Exotoxin_Heatmap_metadata) <- format_Exotoxin_Heatmap_metadata$Group #Saying Group is the row names
format_Exotoxin_Heatmap_metadata <- format_Exotoxin_Heatmap_metadata[-c(1)] #Deleting the column 1 which is group 
pformat_Exotoxin_Heatmap_matrix <- as.matrix(Exotoxin_Heatmap[,-c(1:2)]) #Deleting the Group, Source and Human from the format matrix 
rownames(pformat_Exotoxin_Heatmap_matrix) <- Exotoxin_Heatmap[,1] #Deleting the column 1 which is group 
Exotoxin_Heatmap <- mapply(pformat_Exotoxin_Heatmap_matrix, FUN=as.numeric)
Exotoxin_Heatmap <- matrix(data=Exotoxin_Heatmap, ncol=length(colnames(pformat_Exotoxin_Heatmap_matrix)), nrow=length(row.names(pformat_Exotoxin_Heatmap_matrix)))
row.names(Exotoxin_Heatmap) <- row.names(pformat_Exotoxin_Heatmap_matrix)
colnames(Exotoxin_Heatmap) <- colnames(pformat_Exotoxin_Heatmap_matrix)


pheatmap(Exotoxin_Heatmap)


#test heatmap
pheatmap(Exotoxin_Heatmap, annotation_row = format_Exotoxin_Heatmap_metadata, labels_row = format_Exotoxin_Heatmap_metadata$Group, cluster_rows = F)



#set color palette
Host_colour = list(
  Host = c("Human"="#EB1A17", "Bovine"="#38AD23"))
#Legend
Legend_labels <- c("Absence", "Presence")

Exotoxin <- pheatmap(Exotoxin_Heatmap,
         color=c("red", "blue"),
         legend_breaks = c(0, 1),
         annotation_colors = Host_colour,
         legend_labels = Legend_labels,
         annotation_row = format_Exotoxin_Heatmap_metadata, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         show_rownames = TRUE,
         fontsize = 10, 
         fontsize_row = 10, 
         fontsize_col = 8,
         cellwidth = 30,
         angle_col = 45,
         cellheight = 30,
         width=14,
         height=10,
         main = "Virulence Genes: Exotoxins",
         filename ="FigS8.jpg")

```

```{r VFDB Heatmap Haemolysin}
#Next step is using python code as heatmap is too large so collapsing heat map by matches of lines.Input the formatted combined matrix 
Haemolysin_Heatmap <- read.csv("Virulence/haemolysingroups.csv", fill = TRUE, header = TRUE, sep = ",",check.names=FALSE)
format_Haemolysin_Heatmap_metadata <- Haemolysin_Heatmap %>% select(Group, Host) #Saving the metadata into a separate data frame 
rownames(Haemolysin_Heatmap) <- Haemolysin_Heatmap$Group #Saying Group is the row names
rownames(format_Haemolysin_Heatmap_metadata) <- format_Haemolysin_Heatmap_metadata$Group #Saying Group is the row names
format_Haemolysin_Heatmap_metadata <- format_Haemolysin_Heatmap_metadata[-c(1)] #Deleting the column 1 which is group 
pformat_Haemolysin_Heatmap_matrix <- as.matrix(Haemolysin_Heatmap[,-c(1:2)]) #Deleting the Group, Source and Human from the format matrix 
rownames(pformat_Haemolysin_Heatmap_matrix) <- Haemolysin_Heatmap[,1] #Deleting the column 1 which is group 
Haemolysin_Heatmap <- mapply(pformat_Haemolysin_Heatmap_matrix, FUN=as.numeric)
Haemolysin_Heatmap <- matrix(data=Haemolysin_Heatmap, ncol=length(colnames(pformat_Haemolysin_Heatmap_matrix)), nrow=length(row.names(pformat_Haemolysin_Heatmap_matrix)))
row.names(Haemolysin_Heatmap) <- row.names(pformat_Haemolysin_Heatmap_matrix)
colnames(Haemolysin_Heatmap) <- colnames(pformat_Haemolysin_Heatmap_matrix)


pheatmap(Haemolysin_Heatmap)


#test heatmap
pheatmap(Haemolysin_Heatmap, annotation_row = format_Haemolysin_Heatmap_metadata, labels_row = format_Haemolysin_Heatmap_metadata$Group, cluster_rows = F)



#set color palette
Host_colour = list(
  Host = c("Human"="#EB1A17", "Bovine"="#38AD23"))
#Legend
Legend_labels <- c("Absence", "Presence")

Haemolysin <- pheatmap(Haemolysin_Heatmap,
         color=c("red", "blue"),
         legend_breaks = c(0, 1),
         annotation_colors = Host_colour,
         legend_labels = Legend_labels,
         annotation_row = format_Haemolysin_Heatmap_metadata, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         show_rownames = TRUE,
         fontsize = 10, 
         fontsize_row = 10, 
         fontsize_col = 8,
         cellwidth = 30,
         angle_col = 45,
         cellheight = 30,
         width=8,
         height=5,
         main = "Virulence Genes: Haemolysin",
         filename ="FigS9.jpg")

```

```{r VFDB Heatmap Immune}

#Next step is using python code as heatmap is too large so collapsing heat map by matches of lines.Input the formatted combined matrix 
Immunemod_Heatmap <- read.csv("Virulence/Immunemodulationgroups.csv", fill = TRUE, header = TRUE, sep = ",",check.names=FALSE)
format_Immunemod_Heatmap_metadata <- Immunemod_Heatmap %>% select(Group, Host) #Saving the metadata into a separate data frame 
rownames(Immunemod_Heatmap) <- Immunemod_Heatmap$Group #Saying Group is the row names
rownames(format_Immunemod_Heatmap_metadata) <- format_Immunemod_Heatmap_metadata$Group #Saying Group is the row names
format_Immunemod_Heatmap_metadata <- format_Immunemod_Heatmap_metadata[-c(1)] #Deleting the column 1 which is group 
pformat_Immunemod_Heatmap_matrix <- as.matrix(Immunemod_Heatmap[,-c(1:2)]) #Deleting the Group, Source and Human from the format matrix 
rownames(pformat_Immunemod_Heatmap_matrix) <- Immunemod_Heatmap[,1] #Deleting the column 1 which is group 
Immunemod_Heatmap <- mapply(pformat_Immunemod_Heatmap_matrix, FUN=as.numeric)
Immunemod_Heatmap <- matrix(data=Immunemod_Heatmap, ncol=length(colnames(pformat_Immunemod_Heatmap_matrix)), nrow=length(row.names(pformat_Immunemod_Heatmap_matrix)))
row.names(Immunemod_Heatmap) <- row.names(pformat_Immunemod_Heatmap_matrix)
colnames(Immunemod_Heatmap) <- colnames(pformat_Immunemod_Heatmap_matrix)


pheatmap(Immunemod_Heatmap)


#test heatmap
pheatmap(Immunemod_Heatmap, annotation_row = format_Immunemod_Heatmap_metadata, labels_row = format_Immunemod_Heatmap_metadata$Group, cluster_rows = F)



#set color palette
Host_colour = list(
  Host = c("Human"="#EB1A17", "Bovine"="#38AD23"))
#Legend
Legend_labels <- c("Absence", "Presence")


Immune <- pheatmap(Immunemod_Heatmap,
         color=c("red", "blue"),
         legend_breaks = c(0, 1),
         annotation_colors = Host_colour,
         legend_labels = Legend_labels,
         annotation_row = format_Immunemod_Heatmap_metadata, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         show_rownames = TRUE,
         fontsize = 10, 
         fontsize_row = 10, 
         fontsize_col = 8,
         cellwidth = 30,
         angle_col = 45,
         cellheight = 30,
         width=18,
         height=8,
         main = "Virulence Genes: Immune Modulation",
         filename ="FigS10.jpg")

```


```{r VFDB Heatmap Intracellular}
#Next step is using python code as heatmap is too large so collapsing heat map by matches of lines.Input the formatted combined matrix 
IntracellularAdhesion_Heatmap <- read.csv("Virulence/Intracellular_adhesiongroups.csv", fill = TRUE, header = TRUE, sep = ",",check.names=FALSE)
format_IntracellularAdhesion_Heatmap_metadata <- IntracellularAdhesion_Heatmap %>% select(Group, Host) #Saving the metadata into a separate data frame 
rownames(IntracellularAdhesion_Heatmap) <- IntracellularAdhesion_Heatmap$Group #Saying Group is the row names
rownames(format_IntracellularAdhesion_Heatmap_metadata) <- format_IntracellularAdhesion_Heatmap_metadata$Group #Saying Group is the row names
format_IntracellularAdhesion_Heatmap_metadata <- format_IntracellularAdhesion_Heatmap_metadata[-c(1)] #Deleting the column 1 which is group 
pformat_IntracellularAdhesion_Heatmap_matrix <- as.matrix(IntracellularAdhesion_Heatmap[,-c(1:2)]) #Deleting the Group, Source and Human from the format matrix 
rownames(pformat_IntracellularAdhesion_Heatmap_matrix) <- IntracellularAdhesion_Heatmap[,1] #Deleting the column 1 which is group 
IntracellularAdhesion_Heatmap <- mapply(pformat_IntracellularAdhesion_Heatmap_matrix, FUN=as.numeric)
IntracellularAdhesion_Heatmap <- matrix(data=IntracellularAdhesion_Heatmap, ncol=length(colnames(pformat_IntracellularAdhesion_Heatmap_matrix)), nrow=length(row.names(pformat_IntracellularAdhesion_Heatmap_matrix)))
row.names(IntracellularAdhesion_Heatmap) <- row.names(pformat_IntracellularAdhesion_Heatmap_matrix)
colnames(IntracellularAdhesion_Heatmap) <- colnames(pformat_IntracellularAdhesion_Heatmap_matrix)


pheatmap(IntracellularAdhesion_Heatmap)


#test heatmap
pheatmap(IntracellularAdhesion_Heatmap, annotation_row = format_IntracellularAdhesion_Heatmap_metadata, labels_row = format_IntracellularAdhesion_Heatmap_metadata$Group, cluster_rows = F)



#set color palette
Host_colour = list(
  Host = c("Human"="#EB1A17", "Bovine"="#38AD23"))
#Legend
Legend_labels <- c("Absence", "Presence")





Intracellular <- pheatmap(IntracellularAdhesion_Heatmap,
         color=c("red", "blue"),
         legend_breaks = c(0, 1),
         annotation_colors = Host_colour,
         legend_labels = Legend_labels,
         annotation_row = format_IntracellularAdhesion_Heatmap_metadata, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         show_rownames = TRUE,
         fontsize = 10, 
         fontsize_row = 10, 
         fontsize_col = 8,
         cellwidth = 30,
         angle_col = 45,
         cellheight = 30,
         main = "Virulence Genes: Intracellular Adhesion",
         width=7,
         height=3, 
         filename ="FigS11.jpg")

```

```{r VFDB Heatmap VII Secretion System}

#Next step is using python code as heatmap is too large so collapsing heat map by matches of lines.Input the formatted combined matrix 
VIISecretion_Heatmap <- read.csv("Virulence/TypeVIISecretionsystemgroups.csv", fill = TRUE, header = TRUE, sep = ",",check.names=FALSE)
format_VIISecretion_Heatmap_metadata <- VIISecretion_Heatmap %>% select(Group, Host) #Saving the metadata into a separate data frame 
rownames(VIISecretion_Heatmap) <- VIISecretion_Heatmap$Group #Saying Group is the row names
rownames(format_VIISecretion_Heatmap_metadata) <- format_VIISecretion_Heatmap_metadata$Group #Saying Group is the row names
format_VIISecretion_Heatmap_metadata <- format_VIISecretion_Heatmap_metadata[-c(1)] #Deleting the column 1 which is group 
pformat_VIISecretion_Heatmap_matrix <- as.matrix(VIISecretion_Heatmap[,-c(1:2)]) #Deleting the Group, Source and Human from the format matrix 
rownames(pformat_VIISecretion_Heatmap_matrix) <- VIISecretion_Heatmap[,1] #Deleting the column 1 which is group 
VIISecretion_Heatmap <- mapply(pformat_VIISecretion_Heatmap_matrix, FUN=as.numeric)
VIISecretion_Heatmap <- matrix(data=VIISecretion_Heatmap, ncol=length(colnames(pformat_VIISecretion_Heatmap_matrix)), nrow=length(row.names(pformat_VIISecretion_Heatmap_matrix)))
row.names(VIISecretion_Heatmap) <- row.names(pformat_VIISecretion_Heatmap_matrix)
colnames(VIISecretion_Heatmap) <- colnames(pformat_VIISecretion_Heatmap_matrix)


pheatmap(VIISecretion_Heatmap)


#test heatmap
pheatmap(VIISecretion_Heatmap, annotation_row = format_VIISecretion_Heatmap_metadata, labels_row = format_VIISecretion_Heatmap_metadata$Group, cluster_rows = F)

#set color palette
Host_colour = list(
  Host = c("Human"="#EB1A17", "Bovine"="#38AD23"))
#Legend
Legend_labels <- c("Absence", "Presence")

VII_Secretion <- pheatmap(VIISecretion_Heatmap,
         color=c("red", "blue"),
         legend_breaks = c(0, 1),
         annotation_colors = Host_colour,
         legend_labels = Legend_labels,
         annotation_row = format_VIISecretion_Heatmap_metadata, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         show_rownames = TRUE,
         fontsize = 10, 
         fontsize_row = 10, 
         fontsize_col = 8,
         cellwidth = 30,
         angle_col = 45,
         cellheight = 30,
         main = "Virulence Genes: VII Secretion System",
         width=10,
         height=3, 
         filename = "Fig12.jpg")



```

```{r VFDB Heatmap Others}

#Next step is using python code as heatmap is too large so collapsing heat map by matches of lines.Input the formatted combined matrix 
Others_Heatmap <- read.csv("Virulence/Othersgroups.csv", fill = TRUE, header = TRUE, sep = ",",check.names=FALSE)
format_Others_Heatmap_metadata <- Others_Heatmap %>% select(Group, Host) #Saving the metadata into a separate data frame 
rownames(Others_Heatmap) <- Others_Heatmap$Group #Saying Group is the row names
rownames(format_Others_Heatmap_metadata) <- format_Others_Heatmap_metadata$Group #Saying Group is the row names
format_Others_Heatmap_metadata <- format_Others_Heatmap_metadata[-c(1)] #Deleting the column 1 which is group 
pformat_Others_Heatmap_matrix <- as.matrix(Others_Heatmap[,-c(1:2)]) #Deleting the Group, Source and Human from the format matrix 
rownames(pformat_Others_Heatmap_matrix) <- Others_Heatmap[,1] #Deleting the column 1 which is group 
Others_Heatmap <- mapply(pformat_Others_Heatmap_matrix, FUN=as.numeric)
Others_Heatmap <- matrix(data=Others_Heatmap, ncol=length(colnames(pformat_Others_Heatmap_matrix)), nrow=length(row.names(pformat_Others_Heatmap_matrix)))
row.names(Others_Heatmap) <- row.names(pformat_Others_Heatmap_matrix)
colnames(Others_Heatmap) <- colnames(pformat_Others_Heatmap_matrix)


pheatmap(Others_Heatmap)


#test heatmap
pheatmap(Others_Heatmap, annotation_row = format_Others_Heatmap_metadata, labels_row = format_Others_Heatmap_metadata$Group, cluster_rows = F)



#set color palette
Host_colour = list(
  Host = c("Human"="#EB1A17", "Bovine"="#38AD23"))
#Legend
Legend_labels <- c("Absence", "Presence")

Others <- pheatmap(Others_Heatmap,
         color=c("red", "blue"),
         legend_breaks = c(0, 1),
         annotation_colors = Host_colour,
         legend_labels = Legend_labels,
         annotation_row = format_Others_Heatmap_metadata, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_legend= TRUE,
         annotation_names_row = TRUE,
         show_rownames = TRUE,
         fontsize = 10, 
         fontsize_row = 10, 
         fontsize_col = 8,
         cellwidth = 30,
         angle_col = 45,
         cellheight = 30,
         main = "Virulence Gene Category: Others",
         width=18,
         height=5,
         filename ="FigS13.jpg")

```