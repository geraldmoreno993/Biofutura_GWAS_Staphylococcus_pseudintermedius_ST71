
setwd("/home/ubigem/Documentos/proyectoSP/panaroo/results_panaroo")

# Instala el paquete si no lo tienes
#install.packages("data.table")
         
# Carga el archivo Rtab
library(data.table)
data <- fread("gene_presence_absence.Rtab")


# Supongamos que tu objeto se llama 'mi_matriz'
write.table(data, file = "matriz0_1_beta.csv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)




#TRASPONER LA TABLA




# compute principal components on the accessory portion of the pangenome
PC<-prcomp(df_selected[, c(8:5864)])
# you can use "Country" instead of "Host"
#PCi<-data.frame(PC$x, Host=dafr$Host)
PCi<-data.frame(PC$x, Host=df_selected$Host)
# plot PCA labeled by Host (or Phylogroup)
ggplot(PCi, aes(x=PC1,y=PC2,fill=Host)) +
  geom_point(size = 5, alpha = 0.5, shape = 21) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw()

         