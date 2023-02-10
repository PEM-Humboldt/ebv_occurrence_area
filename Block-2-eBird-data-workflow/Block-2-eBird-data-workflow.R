### Cargar librerias

# Nombrar librerias
packagesList<-list("rstudioapi", "magrittr", "dplyr", "plyr", "purrr", "raster", "terra", "auk", "sf","tools",
                   "gdalUtilities", "tibble", "reshape2", "vapour", "unmarked")

# Instalar librerias que no esten descargadas
install<- packagesList[!packagesList %in% rownames(installed.packages())]
lapply(install, function(x) install.packages(x) )

# Cargar librerias
lapply(packagesList, library, character.only = TRUE)

# Establecer directorio de trabajo (getSourceEditorContext establece el directorio donde está guardado el script)
dirfolder<- "D:/IAvH 2022/Bloque2/draft_ebird_bloque2_06_12_22"
setwd(dirfolder)

# Nombrar dirección de archivo de ocurrencias
file_data<- "ebd_CO_relFeb-2022.txt"

# Nombrar dirección archivo de listas
file_lists<- "ebd_sampling_relFeb-2022.txt"

# Creación archivo de referencia ebird 
ebd <- auk_ebd(file = file_data, file_sampling = file_lists) 




## Definir filtros
# Definir filtro área de estudio

# definir parametros espaciales
dir_basemap = "Caldas.shp"
resolution_basemap<- 1000
crs_basemap<- CRS("+init=epsg:3395")

info_layer<- vapour_layer_info(dir_basemap)
extentBase<- vapour::vapour_read_extent(dir_basemap) %>% {c(min(unlist(map(.,1))), max(unlist(map(.,2))), min(unlist(map(.,3))), max(unlist(map(.,4))))} %>%
  extent %>% st_bbox(crs= info_layer$projection$Proj4) %>% st_as_sfc() %>% st_transform(crs = crs_basemap) %>% st_bbox() %>% extent()

rasterbase = raster(extentBase,crs = crs_basemap, res= resolution_basemap )

# rasterizar área de estudio
tname2 = "grilla_base.tif"; t_file = writeStart(rasterbase, filename = tname2,  overwrite = T); writeStop(t_file);
gdalUtilities::gdal_rasterize(dir_basemap, tname2, burn =1, at=T)
raster_area = rast(t_file) %>% {terra::mask(setValues(., seq(ncell(.))), .)} 

# ajustar covariables área de estudio
setwd(dirfolder)
folder_covars<- "defaultcovars"
dir_covars<- list.files(folder_covars, pattern =  "\\.tif$") 

covars_raster<- lapply(dir_covars, function(x) {setwd(dirfolder); setwd(folder_covars); print(x)
  tname2<- tempfile(fileext = '.tif')
  t_file<- writeStart(rasterbase, filename = tname2,  overwrite=T); writeStop(t_file);
  
  gdalUtilities::gdalwarp(srcfile = x, dstfile= tname2, 
                          tr= res(rasterbase), te = st_bbox(extent(rasterbase)),
                          overwrite=TRUE, r= "near")
  rast(tname2) %>% setNames(gsub(".tif", "",x))
}) %>% setNames(sapply(., function(x) gsub(" ", "_", names(x)) ))






cells_area<- 

cells(raster_area)




# crear tabla covariables
covars_data<-  data.frame(Pixel= cells(raster_area))

for(i in names(covars_raster)){
  covars_data[,i]<-covars_raster[[i]][covars_data$Pixel]
}

# crear una copia del area en formato 4326. deben proyectarse en 4326, pues es la proyecion de los datos de ebird
area_4326<-  raster(raster_area) %>% projectRaster(crs = st_crs(4326)$proj4string, method = "ngb")
ebd_bbox<- raster(raster_area) %>% projectRaster(crs = st_crs(4326)$proj4string) %>% st_bbox()




# Filtrar por área de estudio y esfuerzo de muestreo
filter_bbox<- ebd %>% auk_bbox(ebd_bbox) %>% auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  auk_distance(distance = c(0,10), distance_units= "km") %>% auk::auk_duration(duration = c(0, 300)) %>% auk_date(date = c("2010-01-01", "2020-12-31"))  %>% auk_complete()



##  Generar awk filtrado
# Nombrar archivos de salida
dir.create("outputs_auk_ebird"); 
output_folder<- basename(file_path_sans_ext(dir_basemap))
output_studyarea_ebd<- paste(output_folder, file_data, sep= "_")
output_studyarea_lists<- paste(output_folder, file_lists, sep= "_")
  
## Generar archivo
setwd(paste0(dirfolder, "/outputs_auk_ebird"))
dir.create(output_folder); setwd(output_folder)

# filtro auk
auk_filter(filter_bbox, file = output_studyarea_ebd, file_sampling = output_studyarea_lists, overwrite = T)



#### Cargar archivo auk por especie
## Filtrar por especies
sp<- "Zonotrichia capensis"

# Matriz de especies por lista
setwd(paste(dirfolder, "outputs_auk_ebird", output_folder, sep= "/"))
ebird_data <- auk_zerofill(output_studyarea_ebd, output_studyarea_lists, species = sp , collapse = T, complete= F) %>%
  mutate(observation_count= ifelse(observation_count == "X", 1, observation_count)) %>% 
  mutate(observation_count= as.numeric(observation_count))



# Cortar puntos por mascara área de estudio
ebird_spatial<- ebird_data %>%  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)


ebird_spatial_mask<- as.data.frame(ebird_spatial) %>% mutate(Pixel= raster::extract(area_4326, ebird_spatial) ) %>% 
  dplyr::filter(!is.na(Pixel))  %>%  st_as_sf() %>% st_transform(crs_basemap)


ebird_data_mask<- st_drop_geometry(ebird_spatial_mask)


plot(raster_area)
plot(ebird_spatial_mask[, "geometry"], add= T)

# Generar matriz pixeles por puntos y covariables
summ_sp<- ebird_data_mask %>% mutate(month= format(observation_date,'%Y') ) %>% 
  group_by(Pixel, month) %>% dplyr::summarise(occurrence= ifelse(sum(observation_count)>0,1,0)  ) %>% 
  data.frame(xyFromCell(raster_area, .$Pixel)) %>% list(covars_data) %>% join_all()


matriz_input<- reshape2::acast(summ_sp, Pixel~month, value.var = "occurrence", fill=0)

matriz_input_covars<- list(rownames_to_column(as.data.frame(matriz_input), "Pixel"), dplyr::select(summ_sp, c("Pixel", names(covars_raster))  ) ) %>% 
  join_all() %>% dplyr::select(-"Pixel")


# modelo de ocupación


Data.1 <- unmarkedFrameOccu(y = dplyr::select(matriz_input_covars, -names(covars_raster)), siteCovs = dplyr::select(matriz_input_covars, names(covars_raster)), obsCovs = NULL)  
DataMod <- occu(~1 ~ 1, Data.1)


ests <- plogis(coef(DataMod))		# Escala original

# Obtener los resultados a una escala normal con el error estandard

psiSE <- backTransform(DataMod, type="state")
pSE <- backTransform(DataMod, type="det")


###Cuales valores obtenemos?

# Intervalos de confianza

(ciPsi <- confint(psiSE))
(ciP <- confint(pSE))

# Poner lso resultados en una tabla

resultsTable <- rbind(psi = c(ests[1], ciPsi), p = c(ests[2], ciP))
colnames(resultsTable) <- c("Estimate", "lowerCI", "upperCI")

resultsTable

# Hacer un grafico

plot(1:2, resultsTable[,"Estimate"], xlim=c(0.5, 2.5), ylim=0:1, 
     col=c("blue", "darkgreen"), pch=16, cex=2, cex.lab=1.5,
     xaxt="n", ann=F)
axis(1, 1:2, labels=c(expression(psi), expression(italic(p))), cex.axis=1.5)
arrows(1:2, resultsTable[,"lowerCI"], 1:2, resultsTable[,"upperCI"], 
       angle=90, length=0.1, code=3, col=c("blue", "darkgreen"), lwd=2)








########### Modelos de ocupacion
Data.1 <- unmarkedFrameOccu(y = dplyr::select(matriz_input_covars, -names(covars_raster)), siteCovs = dplyr::select(matriz_input_covars, names(covars_raster)), obsCovs = NULL)  
DataMod <- occu(~1 ~ 1, Data.1)

m1<-occu(~1 ~altitud,Data.1)

newData2 <- distinct(dplyr::select(covars_data, names(covars_raster)))
E.psi <- predict(m1, type="state", newdata=newData2, appendData=TRUE)


# Expected occupancy over range of BOSQUE
newData2 <- data.frame(altitud=  sort(unique(covars_data$altitud)))

E.psi <- predict(m1, type="state", newdata=newData2, appendData=TRUE)

plot(Predicted ~ altitud, E.psi, type="l", col="blue",
     xlab="Altitud",
     ylab="Probabilidad de ocupacion esperada")


match_pixels<- list(covars_data, E.psi) %>% join_all()
  
  
occupancy_model_raster<- raster_area
occupancy_model_raster[match_pixels$Pixel]<- match_pixels$Predicted

plot(occupancy_model_raster)


save.image("C:/Users/LENOVO/Desktop/Block-2-eBird-data-workflow.RData")
load("C:/Users/LENOVO/Desktop/Block-2-eBird-data-workflow.RData")

