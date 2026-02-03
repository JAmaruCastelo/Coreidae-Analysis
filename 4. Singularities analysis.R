# librerias -----------------------------------------------------------------
library (terra)
library(cluster)
library(vegan)
library(ggplot2)
library (ENMTools)
library (ks)
# abrimos los archivos necesarios para el analisis ---------------------------
predict="D:/0. Tesis Doctoral/1. Resultados/4. Singularidad/1. Material needed/Presente_world_predict_all.tiff"


# abrimos la prediccion de cada una de las especies --------------------------------------------------------
predict=rast(predict)

# convierto los valores de la prediccion en datos de Hellinger ------------------------------------------------
mat=values(predict)
rs <- rowSums(mat, na.rm=TRUE)
mat.h <- mat
mat.h[rs > 0, ] <- sqrt(mat[rs > 0, ] / rs[rs > 0])
mat.h[rs == 0, ] <- NA
spp.h <- predict  # estamos guardando como spp.h los valiores convertidos de hellinger
values(spp.h) <- mat.h

# Calculo el PCA de los valores de hellinger y me quedo con 3 caras ------------------
PCA_allv=ENMTools::raster.pca(spp.h, 5)
# verifico el % de variacion explicado por cada una de las variables -----------
prop_varianza <- PCA_allv$pca.object$sdev^2 / sum(PCA_allv$pca.object$sdev^2) # make the calculation of the variance
informacion_pca<- data.frame("Standar Deviation"=PCA_allv$pca.object$sdev,"prop_variance"=prop_varianza, "variance_acumulate"=cumsum(prop_varianza)) # store the summary of the deviation
# se trabajan con 5 ejes que ejemplifican toda la variacion existente con 5 es de 71.94% de la varianza

# archivos para guardar y abrir mas adelante --------------------------------------
writeRaster(spp.h, "D:/0. Tesis Doctoral/1. Resultados/4. Singularidad/2. Hellinger transformation/Species_Hellinger_raster2.tif", overwrite=T)
writeRaster(PCA_allv$rasters, "D:/0. Tesis Doctoral/1. Resultados/4. Singularidad/2. Hellinger transformation/Species_PCA_raster.tif", overwrite=T)

######################### SEGUNDA PARTE ####################################################
# abrimos todos los archivos realizados ---------------------------------------------------------
spp.h<-rast("D:/0. Tesis Doctoral/1. Resultados/4. Singularidad/2. Hellinger transformation/Species_Hellinger_raster.tif")
PCA_allv=rast("D:/0. Tesis Doctoral/1. Resultados/4. Singularidad/2. Hellinger transformation/Species_PCA_raster.tif")

# Para no colapsar la memoria seleccionamos 50000 grillas de todo el objeto de estudio------------
samp <- spatSample(PCA_allv, size=50000, method="random", na.rm=T)


# realizamos un kde multivariado con los 5 ejes
kde <- kde(samp[, -5]) # implementa como maximo 4 ejes


Kmax <- 15

wss <- sapply(1:Kmax, function(k){
  kmeans(samp, centers=k, nstart=50, iter.max=100)$tot.withinss
})






plot(1:Kmax, wss, type="b", xlab="K", ylab="Within-cluster SS")



### ultimo intento
is_peak <- function(mat, i, j){
  neigh <- mat[
    max(1,i-1):min(nrow(mat),i+1),
    max(1,j-1):min(ncol(mat),j+1)
  ]
  mat[i,j] == max(neigh)
}

peak_idx <- which(
  sapply(seq_len(length(z)), function(k){
    i <- ((k-1) %% nrow(z)) + 1
    j <- ((k-1) %/% nrow(z)) + 1
    is_peak(z, i, j)
  }) &
    z > quantile(z, 0.85),   # filtra ruido
  arr.ind = TRUE
)

centros <- cbind(
  kde$x[peak_idx[,1]],
  kde$y[peak_idx[,2]]
)

centros

image(kde$x, kde$y, kde$z,
      xlab="PC1", ylab="PC2")

contour(kde$x, kde$y, kde$z, add=TRUE)

points(centros[,1], centros[,2],
       pch=19, cex=1.3)


points(centros[,1], centros[,2],
       pch=19, cex=1.2)

# este funciono adecuadamente
centros


dist_fun <- function(x){
  if(any(is.na(x))) return(NA)
  
  d <- sqrt( (centros[,1] - x[1])^2 + 
               (centros[,2] - x[2])^2 )
  min(d)
}


pc <- c(TT$PC1, TT$PC2)
singularity <- app(pc, dist_fun)


plot(singularity)
## funcion que proyecta el Kde a el mapoa geografico

kde_fun <- function(x){
  
  if(any(is.na(x))) return(NA)
  
  ix <- findInterval(x[1], kde$x)
  iy <- findInterval(x[2], kde$y)
  
  ix <- pmax(pmin(ix, length(kde$x)), 1)
  iy <- pmax(pmin(iy, length(kde$y)), 1)
  
  kde$z[ix, iy]
}



# voltear la escala pra que te de la inversa
# de los valores

kde_geo <- app(pc, kde_fun)

plot(kde_geo)











# una vez que ya tengo los cuatro centros puedo calcular la distancia 

# encontrar los centros de cada una de las areas. 
center <- global (PCA_allv$rasters, mean, na.rm=T)

singularity <- app(pcs, function(x) {
  sqrt(sum((x - center)^2))
})### singularidad regional si en promedio la distancia al centro al que pertenece es mucho mas cercano o se parece mas a las otras areas.
RGB[RGB < 0] <- 0
RGB[RGB > 255] <- 255
plotRGB(PCA.corte,r="PC"1, g="PC2", b="PC3")





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


