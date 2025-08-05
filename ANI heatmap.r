```{r ANI heatmap}

## Libraries
library(tidyverse)
library(pheatmap)

#load ANI data
Saureus_ANI <- read.csv("fastani_phage.csv",header = TRUE, sep = ",", encoding = "UTF-8") #check.names = FALSE)

# Checking isolate names
Saureus_ANI$Query <- trimws(as.character(Saureus_ANI$Query))
Saureus_ANI$Reference <- trimws(as.character(Saureus_ANI$Reference))

# Build list of all isolates in consistent sorted order
all_isolates <- sort(unique(c(Saureus_ANI$Query, Saureus_ANI$Reference)))

# Create square matrix with rows = columns = isolates
ANI_matrix <- matrix(NA, nrow = length(all_isolates), ncol = length(all_isolates),
                     dimnames = list(all_isolates, all_isolates))

# Fill matrix from fast ANI csv
for (i in 1:nrow(Saureus_ANI)) {
  q <- Saureus_ANI$Query[i]
  r <- Saureus_ANI$Reference[i]
  val <- Saureus_ANI$ANI[i]
  ANI_matrix[q, r] <- val
}

# Fill lower triangle to make it symmetric, if upper triangle exists
ANI_matrix[lower.tri(ANI_matrix)] <- t(ANI_matrix)[lower.tri(ANI_matrix)]

# Ensure all values are numeric
ANI_matrix <- apply(ANI_matrix, 2, as.numeric)
rownames(ANI_matrix) <- colnames(ANI_matrix) <- all_isolates

# Numeric Matrix 
ANI_matrix <- as.matrix(ANI_matrix)

# Generate heatmap without clustering
pheatmap(ANI_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("#DADADA", "#A4C3DE", "#66A2D8", "#1F68A9", "#0B4171", "#222255"))(6),
         breaks = c(90, 92, 94, 96, 98, 100),
         legend_breaks = c(90, 92, 94, 96, 98, 100),
         fontsize_row = 12,
         fontsize_col = 12,
         display_numbers = FALSE,
         filename = "Phage_heatmap.jpg")

dev.off()
```

```
