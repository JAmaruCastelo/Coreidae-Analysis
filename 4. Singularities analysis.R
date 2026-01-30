# librerias -----------------------------------------------------------------
library (terra)
library(cluster)
library(vegan)
library(ggplot2)
library (ENMTools)

# abrimos los archivos necesarios para el analisis ---------------------------
predict="D:/0. Tesis Doctoral/1. Resultados/4. Singularidad/1. Material needed/Presente_world_predict_all.tiff"
areas_compared="D:/0. Tesis Doctoral/1. Resultados/2. Analisis de centros/2. Centros de diversidad y endemismos/A. Primer analisis/Richnes/spatial_RICHNESS_SET1.tif"


# richness --------------------------------------------------------
rich=rast(areas_compared)
rich_vec <- as.polygons(rich >= 18, dissolve = TRUE, values=T)
rich_vec <- disagg(rich_vec[rich_vec$spatial_RICHNESS_SET1==1])

# Cortamos las predicciones de cada una de las especies en las areas identificadas ----------------------------
predict=rast(predict)
rich_vec$name=c(1,2,3,4)
corte=crop(predict, rich_vec)
corte=mask(corte, rich_vec)
categorized=rasterize(rich_vec, corte, field="name") # raster categorizado por tipo de area

# Sampleamos la cantidad de puntos necesarios para el analisis -----------------
puntos=spatSample(rich_vec, size=100,  strata="name") # sacamos 100 puntos de cada area necesaria para el analisis
value_suitability=extract(corte,puntos, ID=F) # estraemos los valores de suitabilidad de cada especie en los puntos

grupo=puntos$name

# extraigo la distancia ------------------------------------------------------
suit_std <- decostand(value_suitability, method = "hellinger")
mat=values(predict)
rs <- rowSums(mat, na.rm=TRUE)
mat.h <- mat
mat.h[rs > 0, ] <- sqrt(mat[rs > 0, ] / rs[rs > 0])
mat.h[rs == 0, ] <- NA
spp.h <- predict
values(spp.h) <- mat.h

# calculo el PCA de las capas ambientales ----------------------------------------



# calculo los dos ejes del NMDS - - - - - - - 
nmds <- metaMDS(suit_std, distance = "bray", k = 3, trymax = 20) # consigo los dos ejes que me expresen mejor las imagenes






# graficamos los resultados del NMDS -----
resultados=as.data.frame(nmds$points)
resultados$grupo=as.factor(grupo)

p=ggplot(data=resultados, aes(x=MDS1, y=MDS2))+
  geom_point(aes(color=grupo))+
  geom_hex(fill = NA, bins = 30,aes(group=grupo,color=grupo))+
  theme_classic()
  
ggsave("D:/5. Tesis doctoral/3.- Results/Turnover/singularity between centers/nmds_richness.svg", plot = p, width = 18, height = 15, units = "cm")

plot(categorized)

# Analisis de Permanova, para ver si son estadisticamente
dist_bray <- vegdist(suit_std, method="bray")
perm <- adonis2(dist_bray ~ resultados$grupo, permutations = 999)

# chequeo de homogeneidad de varianza de los datos --------- 
disp <- betadisper(dist_bray, resultados$grupo)
anova(disp) #### significativo, capta diferencias de dispersion y no solo de centroides


# verificamos si cada punto seleccionado en base a la matriz de disimilitus puede
# ser asociado sin problema unicamente a uno de los centros
sil<- silhouette(grupo, dist_bray)
puntos$sil=sil[,3]

# reescalamos las siluetas de cada uno de los datos ----------
rasters<- geodata::worldclim_global("bio", res=2.5, path="D:/4. GIS/Climatic Bio Variables/Presente") # climatic 
amb=crop(rasters, rich_vec)
amb=mask(amb, rich_vec) 
pca_amb<-raster.pca(amb, 4)
val_test=values(pca_amb$rasters, na.rm=T)
valamb=extract(pca_amb$rasters, puntos)
valamb <- valamb[ , -1]
pred1 <- knn.reg(train = valamb, y = puntos$sil,test = val_test,k = 5)

ok <- complete.cases(values(pca_amb$rasters))
full_1 <- rep(NA, nrow(pca_amb$rasters))
full_1[ok] <- pred1$pred
r_out_1 <- pca_amb$rasters$PC1
terra::values(r_out_1) <- full_1

plot(puntos,"sil")
plot(r_out_1)

# verificiamos con una que calcula
inter_intra_pixel <- function(i){
  
  mi_grupo <- grupo[i]
  
  idx_intra <- which(grupo == mi_grupo)
  idx_inter <- which(grupo != mi_grupo)
  
  intra_mean <- mean(Dmat[i, idx_intra])
  inter_mean <- mean(Dmat[i, idx_inter])
  
  diff <- inter_mean - intra_mean
  ratio <- inter_mean / (intra_mean + 1e-12)  # opcional
  
  return(c(intra = intra_mean,
           inter = inter_mean,
           diff = diff,
           ratio = ratio))
} # esta funcion calcula en differencias o en ratio
inter_intra_pixel2 <- function(i){
  mi_grupo <- grupo[i]
  idx_intra <- which(grupo == mi_grupo) # los indices que son mi grupo
  idx_inter <- which(grupo != mi_grupo) # los indices que no son mi grupo
  intra_mean <- mean(Dmat[i, idx_intra]) # el promedio de las diatnacia de mi grupo
  # calculo mi distancia promedio a los otros grupos
  medias_inter_por_grupo <- tapply(Dmat[i, idx_inter],grupo[idx_inter],mean,na.rm = TRUE)
  
  medias_todos <- c(medias_inter_por_grupo, setNames(intra_mean, mi_grupo))
  medias_todos
  
  # Grupo con mayor media total
  grupo_max <- names(medias_todos)[which.min(medias_todos)]
  max_mean  <- min(medias_todos)
  
  diff  <- max_mean - intra_mean #### la mayor media total menos la media interna
  ratio <- max_mean / (intra_mean + 1e-12)
  
  return(list(intra = intra_mean,max_mean = max_mean,
              grupo_max = grupo_max,diff = diff,ratio = ratio))} # esta funcion calcula cual es la mayor media

# con la primera funcion
Dmat <- as.matrix(dist_bray)
res_inter_intra <- t(sapply(seq_len(nrow(Dmat)), inter_intra_pixel2))
res_inter_intra <- as.data.frame(res_inter_intra)

puntos$grupo=res_inter_intra$grupo_max

plot(puntos, "grupo")

# con la segunda función-------------------------------------
# presenta las condiciones ambientales que permiten a otras especies de otras
# areas tambien sobrevivir,,,,,,,
#### en lugar de eso tambien se podria conocer 
res_inter_intra <- t(sapply(seq_len(nrow(Dmat)), inter_intra_pixel2))
res_inter_intra <- as.data.frame(res_inter_intra)




puntos$grupo=res_inter_intra$ratio
plot(puntos, "grupo")









# pairwise comparison----------------------------------------------------

grupo=valores$name
ss=suit_std[grupo %in% c("uno", "dos"),]
name=grupo[grupo %in% c("uno", "dos")]
dist_bray <- vegdist(ss, method="bray")
perm <- adonis2(dist_bray ~ name, permutations = 999)
ss=suit_std[grupo %in% c("tres", "cuatro"),]
name=grupo[grupo %in% c("tres", "cuatro")]
dist_bray <- vegdist(ss, method="bray")
perm <- adonis2(dist_bray ~ name, permutations = 999)
perm
# perm de uno con 2 y verificamos que salga correctamente



# contribucion de espeices, opt
simper_res <-simper (suit_std, valores$name)
summary(simper_res)


