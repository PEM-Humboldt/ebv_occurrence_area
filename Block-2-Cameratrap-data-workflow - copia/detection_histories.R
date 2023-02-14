packagesList<-list("dplyr", "plyr",  "pbapply", "camtrapR", "tidyverse", "readr", "lubridate", "sf", "terra", "raster", "vapour")
lapply(packagesList, library, character.only = TRUE)

# Load records data
dir_sites_putumayo<- "D:/Humboldt_provisional/PRUEBAS_CODIGO/ocupacion_B2/sites_putumayo.csv"
sp_rec = read.csv2(dir_sites_putumayo) %>% dplyr::select(-"X")

# List species name
unique(sp_rec$scientificName) 


# Change date order if needed (correct format yyyy-mm-dd)
sp_rec$eventDate = lubridate::parse_date_time(x = sp_rec$eventDate, order = c("dmy", "Ymd","dmY"))


# Load covariates data 
covars = read.csv("D:/Humboldt_provisional/PRUEBAS_CODIGO/ocupacion_B2/covars_putumayo_2.csv")

# Select detection covariates
covars = dplyr::select(covars, c("Cam.Site", "No_spp", "events", "Lat_Y", "Long_X", "Instal.Date", "Last.eventDate", "Cam.Days"))

# Change date order, if needed
covars$Instal.Date = lubridate::parse_date_time(x = covars$Instal.Date, order = c("dmy", "Ymd","dmY")) 
covars$Last.eventDate = lubridate::parse_date_time(x = covars$Last.eventDate, order = c("dmy", "Ymd","dmY"))

# Load study area as shapefile
dir_basemap = "D:/Humboldt_provisional/PRUEBAS_CODIGO/ocupacion_B2/Putumayo.shp"
crs_basemap<- CRS("+init=epsg:3395")
res=1000


info_layer<- vapour_layer_info(dir_basemap)
extentBase<- vapour::vapour_read_extent(dir_basemap) %>% {c(min(unlist(map(.,1))), max(unlist(map(.,2))), min(unlist(map(.,3))), max(unlist(map(.,4))))} %>%
  extent %>% st_bbox(crs= info_layer$projection$Proj4) %>% st_as_sfc() %>% st_transform(crs = crs_basemap) %>% st_bbox() %>% extent()

rasterbase = raster(extentBase,crs = crs_basemap, res= res )

tname2 = tempfile(fileext = '.tif'); t_file = writeStart(rasterbase, filename = tname2,  overwrite = T); writeStop(t_file);
gdalUtilities::gdal_rasterize(dir_basemap, tname2, burn =1, at=T)
raster_area = rast(t_file) %>% {terra::mask(setValues(., seq(ncell(.))), .)} 

cell_area = terra::cells(raster_area)

# Rasterize data
cameras_spatial<- covars %>% mutate(x= Long_X, y= Lat_Y) %>% dplyr::filter( (!is.na(x)) |  (!is.na(y)) ) %>% st_as_sf(coords = c("x", "y"), crs = info_layer$projection$EPSG) %>%
  st_as_sf() %>% st_transform(crs_basemap)

cameras_pixels<- st_drop_geometry(cameras_spatial) %>% mutate(pixel_ID= raster::extract(raster_area, cameras_spatial)[,2] )


# Join covers and records data
sp_rec1 = list(sp_rec, cameras_pixels) %>% join_all()

# Create survey lenght matrix for each site
survey_lenght = plyr::ddply(sp_rec, c("Cam.Site"), 
                            dplyr::summarise,  
                            Instal.Date = as.Date(Instal.Date),
                            Last.eventDate = max(as.Date(eventDate)),
                            Cam.Days = Last.eventDate-Instal.Date)


# Create camera operation name of stations
cam_op = cameraOperation(CTtable = cameras_pixels,
                         setupCol = "Instal.Date",
                         retrievalCol = "Last.eventDate",
                         stationCol = "pixel_ID",
                         cameraCol = "Cam.Site",
                         byCamera = FALSE,
                         allCamsOn = TRUE,
                         camerasIndependent = FALSE,
                         hasProblems  = FALSE)
str(cam_op)

# Add a column with date and time 
sp_rec1 = sp_rec1 %>% mutate(DateTimeOriginal = paste(eventDate, eventTime))
sp_rec1$DateTimeOriginal = as.POSIXlt(sp_rec1$DateTimeOriginal) # must be format="%Y-%m-%d:%H:%M:%S"

# Create list with all the species in the survey
sp_prior = c("Cabassous unicinctus", "Saimiri cassiquiarensis", "Cebus yaracus", 
              "Cuniculus paca", "Dasyprocta fuliginosa", "Dasypus novemcinctus", 
              "Didelphis marsupialis", "Eira barbara", "Leopardus pardalis", 
              "Leopardus wiedii", "Mazama americana ", "Myrmecophaga tridactyla", 
              "Nasua nasua", "Pecari tajacu", "Philander andersoni", 
              "Procyon cancrivorus", "Tamandua tetradactyla", "Tupinambis cuzcoensis", 
              "Lontra longicaudis ", "Puma concolor ", "Galictis vittata", 
              "Cuniculus Paca", "Leopardus pardalis ", "Mazama americana", 
              "Tamandua tetradactyla ", "Dasypus novemcinctus ", "Dasyprocta fuliginosa ", 
              "Cuniculus paca ", "Eira barbara ", "Philander andersoni ", "Nasua nasua ", 
              "Didelphis marsupialis ", "Panthera onca ", "Leopardus tigrinus", 
              "Procyon cancrivorus ", "Potos flavus ", "Metachirus nudicaudatus")

# Filter matrix by species list
sp_prior1 = sp_rec1 %>% filter(scientificName %in%sp_prior)

# Show all unique entries of species cam.sites that are not in the covars sites, and delete it
not_in_sprec = unique(sp_rec1$Cam.Site)[!unique(sp_rec1$Cam.Site) %in% covars$Cam.Site]
sp_rec2 = filter(sp_rec1, !Cam.Site %in% not_in_sprec)

# Validate all unique entries of species cam.sites are in the covers sites
not_in_sprec1 = unique(sp_rec2$Cam.Site)[!unique(sp_rec2$Cam.Site) %in% covars$Cam.Site]

# Create detection history for one specie, Cuniculus paca is used as an example
DH_paca = detectionHistory (recordTable = sp_rec2,
                            species = "Cuniculus paca",
                            camOp = cam_op,                      
                            stationCol = "pixel_ID", 
                            speciesCol = "scientificName",
                            recordDateTimeCol = "DateTimeOriginal",
                            recordDateTimeFormat  = "%Y-%m-%d%H:%M:%S",
                            occasionLength = 6,  #change to colaps diferent dates        
                            day1 = "station",  #first day of survey; if we want to specify a date put in "survey"
                            datesAsOccasionNames = F, 
                            includeEffort = F, #careful if trapping effort is thought to influence detection probability, it can be returned by setting includeEffort = TRUE. 
                            scaleEffort = F)            #maybe wise using T, explore later 

str(DH_paca)
DH_paca$detection_history
write.csv(DH_paca, "D:/Humboldt_provisional/PRUEBAS_CODIGO/ocupacion_B2/data_out/DH_paca.csv")
