
#This R file loads functions and datasets needed for other scripts in this repository. 

##Load required libraries
packages <- c("ggpubr", "rnaturalearth","corrplot", "ncdf4", "rerddap", "heatwaveR", "dplyr", "readr", "readxl", "sf", "sp", "ggridges", "ggplot2", "viridis", "ggpubr", "tidyr", "forcats",  "mapproj", "factoextra", "EFAtools", "psych", "reshape2", "devtools", "FactoMineR", "usmap", "fmsb")
invisible(lapply(packages, library, character.only= TRUE))

#1) Functions/data for GetOISST: 

##Functions: 
###Function 1: take nc file and put into df 
OISST_prep <- function(filename){
  nc <- ncdf4::nc_open(filename) ####Open the NetCDF connection
  res <- ncvar_get(nc, varid = "sst") ####Read data from the netCDF file based off SST, extract the SST values and add the lon/lat/time dimension names
  dimnames(res) <- list(lon = nc$dim$longitude$vals,
                        lat = nc$dim$latitude$vals,
                        t = nc$dim$time$vals)
  res <- as.data.frame(reshape2::melt(res, value.name = "temp"), row.names = NULL) %>% 
    mutate(t = as.Date(as.POSIXct(t, origin = "1970-01-01 00:00:00")),
           temp = round(temp, 2)) ####Melt: Convert the data into a long dataframe
  ####as.POSIXct stores date and time in seconds with the # of seconds beginning at 1 January 1970
  res <- res %>% drop_na(temp)
  nc_close(nc)
  return(res)
}

###Function 2: Get OISST files from folder and apply the OISST_prep to all of them
OISST_df <- function(mypath){
  files1 <- list.files(mypath)
  files1
  setwd(mypath)
  SST_data <- list()
  for(i in 1:6){
    SST_data[[i]] <- OISST_prep(files1[[i]])
  } 
  SST_data <- as.data.frame(do.call(rbind, SST_data))
  SST_data$lat_lon <- paste(SST_data$lat, SST_data$lon, sep = "_")
  return(SST_data)
}

###Function 3: Finds overlap between SST and community radius dataset
overlap_df <-function(df1, df2){
  df_edt <- df1 %>% distinct(lat_lon, .keep_all= TRUE) %>% select(lat_lon, lat, lon)
  df_edt_sites = st_as_sf(df_edt, coords=c("lon","lat")) ####Turn coordinates into geometry object for SST dataset 
  overlap <- st_intersection(df_edt_sites, df2) ####Find geometry intersection between SST dataset and radius around coastal communities dataset 
  keeperz <- as.data.frame(overlap)
  keeperz <- unique(overlap$lat_lon)
  df_final <- df1[df1$lat_lon %in% keeperz== TRUE,] #Retains rows only with lat_lon within radius1
  return(df_final)
} 

###Function 4: Visualization to confirm accuracy of pulling OISST data
radius_coords <-function(df1, df2, df3){
  final_coords <- df1[!duplicated(df1$lat_lon),] ####Keep unique coordinates in original SST dataset filtered for keeperz
  final_coords = st_as_sf(final_coords, coords= c("lon", "lat"))  
  plot(st_geometry(df2), axes = TRUE, col= "red")
  plot(st_geometry(df3), axes = TRUE, add= TRUE)
  plot(st_geometry(final_coords), axes = TRUE, add= TRUE, col= "blue") ####Plot keeper coordinates from large SST dataset 
}

###Function 5: needed for MHW analysis
block_only<-function(df){
  clim <- ts2clm(data = df[order(df$t),], climatologyPeriod = c("1983-01-01", "2012-12-31")) #### First calculate the climatologies by specifying a 30 year time period 
  event <- detect_event(data = clim) ####Then detect the MHW events
  block_average(event) ####Return only the event metric dataframe of results with yearly means 
}

###Function 6: Get rid of land points: we do this function just in case but in reality we do not need it (none of the areas had land points so virtually this function does nothing) 
land_pts_away <-function(df){
  df$lat_lon <-paste(df$lat, df$lon, sep = "_")
  landy<-aggregate(total_days ~ lat_lon + lat + lon, FUN = sum, df) ####A land point won't have any days for total MHW days
  keeperz <- unique(landy$lat_lon)
  df_new <-df[df$lat_lon %in% keeperz == TRUE,] ####Keep only the lat_lon points that have MHW statistics for them 
  return(df_new)
}

##Datasets: US mainland, HI, AK, and Australia communities with coordinates, form radius around communities  
###US mainland: 
####Loading data that has communities that are not connected with SST yet. Data accessed (https://www.st.nmfs.noaa.gov/data-and-tools/social-indicators/) for SE, NE, and West Coast US communities for an arbitrary year to get community names and their associated lat/lon 
US_main_communities <- read_csv("~/Documents/GitHub/Nye-research-technicians/MHW.MCS/Data/US.main.communities.csv")
only_US_main_orig <- US_main_communities %>% select("Community Name", State, Region, Latitude, Longitude) %>% rename(Community.Name= "Community Name")
only_US_main_sf = st_as_sf(only_US_main_orig,coords=c("Longitude","Latitude")) ####Converts to sf object, coords states name of columns holding coordinates, geometry column is just one point 
####Plot and draw radius/buffer around communities
plot(st_geometry(only_US_main_sf), axes = TRUE)
radius1_US_main <- st_buffer(only_US_main_sf, dist = 1) ###111 km radius around each community   

###AK:
###Loading data that has AK communities that are not connected with SST yet. Data accessed (https://www.st.nmfs.noaa.gov/data-and-tools/social-indicators/) for AK communities for year 2016 to get community names and their associated lat/lon 
US_AK_communities <- read_csv("~/Documents/GitHub/Nye-research-technicians/MHW.MCS/Data/US.AK.communities.csv")
only_US_AK_orig <- US_AK_communities %>% select("Community Name", State, Region, Latitude, Longitude) %>% rename(Community.Name= "Community Name") 
only_US_AK_sf = st_as_sf(only_US_AK_orig,coords=c("Longitude","Latitude")) ####Converts to sf object, coords states name of columns holding coordinates
radius1_AK <- st_buffer(only_US_AK_sf, dist = 1)

###HI: 
###Loading data that has HI communities that are not connected with SST yet. Data accessed (https://www.st.nmfs.noaa.gov/data-and-tools/social-indicators/) for HI communities for year 2016 to get community names and their associated lat/lon 
US_HI_communities <- read_csv("~/Documents/GitHub/Nye-research-technicians/MHW.MCS/Data/US.HI.communities.csv")
only_US_HI_orig <- US_HI_communities %>% select("Community Name", State, Region, Latitude, Longitude) %>% rename(Community.Name= "Community Name")
only_US_HI_sf = st_as_sf(only_US_HI_orig,coords=c("Longitude","Latitude")) 
radius1_HI <- st_buffer(only_US_HI_sf, dist = 1) 

###Australia:
###Loading data that has local government areas that are not connected with SST yet. Shapefiles accessed through Australian Bureau of Statistics (https://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/1270.0.55.003July%202015?OpenDocument#Data) for 2016 and manipulated within QGIS version 3.20.0.
only_Aus_cs_orig <- read_csv("/Users/sallydowd/Documents/GitHub/Nye-research-technicians/MHW.MCS/Data/Aus_LGAs.csv")
only_Aus_cs_orig <- only_Aus_cs_orig %>% select(-c(vertex_index, vertex_part, vertex_part_index, fid, HubName, HubDist, angle)) %>% rename("Longitude"= "xcoord", "Latitude"= "ycoord")
only_Aus_sf = st_as_sf(only_Aus_cs_orig,coords=c("Longitude","Latitude")) 
radius1_Aus <- st_buffer(only_Aus_sf, dist = 1) 

#2) CalculateMHW.Metrics.Rmd

##Functions: 
###Function 1: Find intersection b/w MHW metrics and communities 
metrics_df <- function(df1, df2, df3, number){
  df1$input_N<-NA
  df1$input_mean<-NA
  df1$input_var<-NA
  df1$input_max<-NA
  df1$input_min<-NA
  for(i in 1:nrow(df1)){ #each row is a US community 
    intersects <- as.data.frame(st_intersection(df2[i,], df3)) 
    df1$input_N[i]<-nrow(intersects) #Rows correspond to # of lat/lon points 
    intersects <- intersects[,number] #select column with the metric you're interested in, now when you call intersects df it'll be to calculate stats around metric of interest
    df1$input_mean[i]<-mean(intersects)
    df1$input_var[i]<-var(intersects)
    df1$input_max[i]<-max(intersects)
    df1$input_min[i]<-min(intersects)  
  }
  df1 <- df1[order(df1$input_mean, decreasing = TRUE) ,]
  return(df1)
}

###Function 2: Plots with coordinates, communities, and event metric for MHW or MCS
graph1 <- function(dataset1, dataset2, metric){
  ggplot() + geom_sf(data= dataset1, cex= 0.7) + geom_point(data= dataset2, aes(x=lon, y=lat, colour= metric)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.text.align= 0, legend.title= element_text(size = 12), legend.text = element_text(size= 10), axis.text=element_text(size=10), axis.title=element_text(size=12))   
} ####Add in radius with red point: geom_sf(data=radius3[1,], cex=.5, alpha= 0.05) + geom_sf(data= dataset2[1,], cex= 2, col="red")

###Function 3: Plots with communities, radius, and event metric for MHW or MCS
US_shp <- st_read(dsn = "/Users/sallydowd/Library/CloudStorage/GoogleDrive-sallycdowd@gmail.com/My Drive/Nye.lab/Coastal.vulnerability/US.data/Shapefiles/tl_2021_us_state", layer = "tl_2021_us_state")
US_shp_df <- fortify(US_shp) #Need to get US shapefile to use in graph 2 function
Aus_shp <- st_read(dsn = "/Users/sallydowd/Library/CloudStorage/GoogleDrive-sallycdowd@gmail.com/My Drive/Nye.lab/Coastal.vulnerability/Aus.data/ShapeFiles/Australia", layer = "LGA_2016_AUST")
#Aus_shp <- st_read(dsn="/Users/sallydowd/Google Drive/My Drive/Nye.lab/Coastal.vulnerability/Aus.data/ShapeFiles/Australia", layer = "aust_cd66states")
Aus_shp_df <- fortify(Aus_shp)
graph2 <- function(dataset1, shapefile, metric){
ggplot() + geom_polygon(data = shapefile, aes(x=long, y = lat, group = group), fill= "ghostwhite") +  geom_path(data = shapefile, aes(x = long, y = lat, group = group), 
color = "black", size = 0.2) + coord_equal() +geom_polygon() + coord_equal() + geom_point(data= dataset1, aes(x= Longitude, y= Latitude, colour= metric), size=2) + standard_theme + xlab("Longitude") + ylab("Latitude") + theme(panel.background = element_rect(fill= "gray")) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x= element_blank(), axis.title.y= element_blank())
}

aus_graph <- function(dataset1, metric){
  ggplot(ozmap_states) + geom_sf(fill= "white") + geom_point(data= dataset1, aes(x= Longitude, y= Latitude, colour= metric), size=2) + standard_theme + xlab("Longitude") + ylab("Latitude") + theme(panel.background = element_rect(fill= "gray")) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x= element_blank(), axis.title.y= element_blank()) + theme(legend.position= "bottom") + labs(colour= "Hazard")
}

###Function 3.5
graph2.5 <- function(dataset1, dataset2, metric){
  ggplot() + geom_point(data= dataset1, aes(x= Longitude, y= Latitude, colour= metric)) + geom_point(data= dataset2, aes(x=lon, y=lat), shape=".") + standard_theme + xlab("Longitude") + ylab("Latitude")
}

###Function 4: Transform lat and lon data so US map contains mainland, AK, and HI 
us_map_df <- function(dataset1){
  US_metric_mod <- dataset1 %>% rename("lon"= "Longitude", "lat"= "Latitude") %>% select(lat, everything()) %>% select(lon, everything()) #move latitude and longitude to the front of dataset 
  US_transformed <- usmap_transform(US_metric_mod) ####Discarded datum unknown error: I think this is ok and makes sense b/c dataset doesn't have a coordinate reference system (just have lat/lon  
  US_transformed_edt <- US_transformed %>% mutate(lon.1 = st_coordinates(.)[,1], lat.1 = st_coordinates(.)[,2])
  return(US_transformed_edt)
}

US_map_combo <- function(dataset1, metric, colour2){
  plot_usmap() + geom_point(data= dataset1, aes(x = lon.1, y = lat.1, size= metric, colour= metric)) + standard_theme + labs(colour= colour2) + theme(panel.background = element_rect(fill= "gray"), axis.text.x = element_blank(),
                                                                                                                                                                                                                      axis.text.y = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) + scale_colour_viridis(option= "C") + guides(size=FALSE) + theme(legend.title= element_text(size = 18), legend.text = element_text(size= 16)) + labs(colour= colour2)
}

US_map_zoom_in <- function(states, dataset1, metric, colour2){
    plot_usmap(include= states) + geom_point(data= dataset1, aes(x = lon.1, y = lat.1, size= metric, colour= metric)) + standard_theme + theme(panel.background = element_rect(fill= "lightgray"), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) + scale_colour_viridis(option= "C") + guides(size=FALSE) + theme(legend.position= c(0.20, 0.80), legend.background= element_rect(fill= "gray")) + labs(colour= colour2) +  theme(legend.title= element_text(size = 18), legend.text = element_text(size= 16))
}


US_map_combo_cat <- function(dataset1, metric, colour2){
  plot_usmap() + geom_point(data= dataset1, aes(x = lon.1, y = lat.1, colour= as.factor(metric)), size=2 ) + theme(panel.background = element_rect(fill= "gray"), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) + theme(legend.title= element_text(size = 18), legend.text = element_text(size= 16)) + labs(colour= colour2) + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position= "bottom") + scale_color_manual(values= cols)
}

###Create a standard ggplot theme 
standard_theme <- theme_bw() + theme(panel.border = element_rect(fill=NA, colour = "black", size=1)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.text.align= 0, legend.title= element_text(size = 14), legend.text = element_text(size= 12), axis.text=element_text(size=12), axis.title=element_text(size=14)) 

#3) PCA.Rmd

##Functions: 
###Function 1: use principal function to calculate PCA, this is code from NOAA- will need to cite 
principal_PCA <- function(dataset){
  print("CORRELATION MATRIX")
  print(cor(dataset))
  print("KMO and Bartletts test of Sphericity")
  print(KMO(cor(dataset)))
  print(cortest.bartlett(dataset))
  coms_pca_rotated <- principal(dataset, nfactors=1, rotate="varimax", scores=TRUE)
  print(loadings(coms_pca_rotated))
  print(summary(coms_pca_rotated))
  variables_num <- length(coms_pca_rotated$values)
  largest_eigenvalue <- max(coms_pca_rotated$values) ####Eigenvalues explain variance for PCs, if there are multiple eigenvalues than there are multiple PCs (we didn't force a single factor solution then, right?)
  armours_theta <- (variables_num/(variables_num-1)*(1-1/largest_eigenvalue))
  print("ARMOURS THETA")
  print(armours_theta)
  return(coms_pca_rotated)
}

###Function 2: Get factor scores for PC1, turn into categorical score 
Indicator_data <- read_csv("/Users/sallydowd/Documents/GitHub/Nye-research-technicians/MHW.MCS/Data/US.indicator.data.final.csv")
Aus_indicators <- read_csv("/Users/sallydowd/Documents/GitHub/Nye-research-technicians/MHW.MCS/Data/Aus.indicator.data.final.csv")
number_US <- 3255
number_Aus <- 209
Community_ID <- Indicator_data %>% select(Community.Name) %>% mutate(X= 1:number_US)
LGA_ID <- Aus_indicators %>% select(LGA) %>% mutate(X= 1:number_Aus)
scores_PCA <- function(df1, PCA_df, index, number, community_df) {
  ####Get PC1 scores in dataset with type of index
  df1 <- df1 %>% mutate(PC1_scores= PCA_df$scores, PC_index= index) %>% select(PC1_scores, PC_index) 
  ####Reverse PC1 scores for labor force and housing characteristics index based of PC_index column
  if(df1$PC_index == "labor_force"  || df1$PC_index == "housing_charcs") {
    df1$PC1_scores <- df1$PC1_scores* -1
  } else { df1$PC1_scores <- df1$PC1_scores* 1}
  scores <- as.data.frame(df1)
}

#4) CalculateRisk.Rmd 

##Datasets to load: 
US_scores <- read.csv("/Users/sallydowd/Documents/GitHub/Nye-research-technicians/MHW.MCS/Data/US.PC1.scores.csv")
US_total_cumi_MHW_ld_final <- read.csv("/Users/sallydowd/Documents/GitHub/Nye-research-technicians/MHW.MCS/Data/US.total.cumi.2016.final.V2.1.csv")
US_total_dur_MHW_ld_final <- read.csv("/Users/sallydowd/Documents/GitHub/Nye-research-technicians/MHW.MCS/Data/US.total.dur.2016.final.V2.1.csv")
Aus_scores <- read.csv("/Users/sallydowd/Documents/GitHub/Nye-research-technicians/MHW.MCS/Data/Aus.PC1.scores.csv")
Australia_cumi_MHW_ld_final <- read.csv("/Users/sallydowd/Documents/GitHub/Nye-research-technicians/MHW.MCS/Data/Australia.cumi.2016.final.V2.1.csv")
Australia_dur_MHW_ld_final <- read.csv("/Users/sallydowd/Documents/GitHub/Nye-research-technicians/MHW.MCS/Data/Australia.dur.2016.final.V2.1.csv")

##Combine MHW and indicator datasets 
###US: 
US_scores_MHW <- US_total_cumi_MHW_ld_final[,c(3:7,9)] %>% left_join(US_scores[,2:9], by= "Community.Name") %>% left_join(US_total_dur_MHW_ld_final[,c(3,9)], by= "Community.Name")
i<- c(6:14)
US_scores_MHW[, i] <- apply(US_scores_MHW[ , i], 2, function(x) as.numeric(as.character(x))) 
US_scores_MHW[933, 1] <- "Paso Robles, CA"

###Australia: 
Aus_scores_MHW <- Australia_cumi_MHW_ld_final[,c(5,7,12:13,15)] %>% left_join(Aus_scores[,c(2:8)], by= "LGA") %>% left_join(Australia_dur_MHW_ld_final[,c(4,14)], by = "LGA")
i<- c(5:12)
Aus_scores_MHW[, i] <- apply(Aus_scores_MHW[ , i], 2, function(x) as.numeric(as.character(x))) 

##Functions
###Function 1: Calculate categorical score for PC1 scores
get_categorical_score <- function(x) {
  cat_score <- case_when(x >= 1 ~ 4,
                         x >= 0.5 ~ 3, #0.5 >= x < 1 
                         x >= 0 ~ 2, #0 >= x < 0.5 
                         x >=  -99999 ~ 1, #if negative, -9999 >= x <0
                         TRUE ~ 0) 
  return(cat_score) 
}

###Function 2: Calculate risk from categorized data 
risk_precat <- function(df, group, exposure){
  newdf <- df %>% group_by({{ group }}) %>% mutate(risk_cumi= cumi_quartile*{{ exposure }}*avg_cat_vuln)
  plot <- ggplot(newdf, aes(x= risk_cumi)) + geom_histogram() + xlab("Risk") + ylab("Count") + standard_theme
  print(plot)
  return(newdf)
}

###Function 3: US, Risk categorization method 1: calculate quantiles for US and assign to categories for overall risk
risk_quantiles <- function(df1, plot_title){
  ####a) Filter out communities that have a 0 for risk 
  #df2 <- df1 %>% filter(risk_cumi!=0)
  df2 <- df1
  ####b) Calculate quantiles and assign categories 
  df2$Quartile <- cut(x= df2$risk_cumi, breaks= quantile(df2$risk_cumi),include.lowest=TRUE,labels=FALSE)
  quartile_plot <- df2 %>% ggplot(aes(x= as.factor(Quartile), y= risk_cumi)) + geom_boxplot() + scale_x_discrete(labels= c("Low", "Medium", "Medium-high", "High")) + labs(title= plot_title) + xlab("Quartile") + ylab("Risk") + standard_theme
  print(quartile_plot)
  df2_cat <- df2 %>% mutate("Category"= ifelse(Quartile== 1, "Low", ifelse(Quartile==2, "Medium", ifelse(Quartile==3, "Medium-high", ifelse(Quartile==4, "High", NA)))))
  return(df2_cat)
  }
 
###Function 4: Australia, Risk categorization method 1: calculate quantiles for Aus and assign to categories for overall risk
risk_quantiles_Aus <- function(df1, plot_title){
  df2 <- df1 
  ####b) Calculate quantiles and assign categories 
  df2$Quartile <- cut(x= df2$risk_cumi, breaks= quantile(df2$risk_cumi),include.lowest=TRUE,labels=FALSE)
  quartile_plot <- df2 %>% ggplot(aes(x= as.factor(Quartile), y= risk_cumi)) + geom_boxplot() + scale_x_discrete(labels= c("Low", "Medium", "Medium-high", "High")) + labs(title= plot_title) + xlab("Quartile") + ylab("Risk") + standard_theme
  print(quartile_plot)
  df2_cat <- df2 %>% mutate("Category"= ifelse(Quartile== 1, "Low", ifelse(Quartile==2, "Medium", ifelse(Quartile==3, "Medium-high", ifelse(Quartile==4, "High", NA)))))
}

###Function 5: Calculate quantiles and assign categories 
quantiles_risk_components <- function(df1, column, plot_title){
  ####a) Calculate quantiles and assign categories 
  df1$Quartile <- cut(x= column, breaks= quantile(column), include.lowest=TRUE,labels=FALSE)
  quartile_plot <- df1 %>% ggplot(aes(x= as.factor(Quartile), y= column)) + geom_boxplot() +
    scale_x_discrete(labels=seq(1,4,1)) + labs(title= plot_title) + xlab("Quartile") + ylab("% employed in fishing") + standard_theme 
  print(quartile_plot)
  df1_cat <- df1 %>% mutate("Category"= ifelse(Quartile== 1, "Low", ifelse(Quartile==2, "Medium", ifelse(Quartile==3, "Medium-high", ifelse(Quartile==4, "High", NA)))))
  return(df1_cat)
}
