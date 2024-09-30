


# Make a grid across the map area
grid <- rsa_country_sf %>%
  sf::st_make_grid(cellsize = c(cell_size, cell_size),
                   offset = c(sf::st_bbox(rsa_country_sf)$xmin,
                              sf::st_bbox(rsa_country_sf)$ymin)) %>%
  sf::st_sf() %>%
  dplyr::mutate(cellid = dplyr::row_number())

taxon.sf = taxa_Fabacae %>%
  dplyr::select(decimalLatitude,decimalLongitude,
                species,coordinateUncertaintyInMeters,year,speciesKey) %>% #select occurrence data
  dplyr::filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
  dplyr::filter(coordinateUncertaintyInMeters<=res*1000) %>% 
  sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
               crs = 4326) %>%
  st_join(grid) %>% 
  as.data.frame() %>% 
  dplyr::select(-geometry) %>% 
  dplyr::mutate(occurrenceStatus=1)



taxon_cube<-b3gbi::process_cube(taxon_sf,grid_type = "custom",cols_cellCode = "cellid", cols_year = "year",
                 cols_occurrences = "occurrenceStatus",cols_species = "species",cols_speciesKey = "speciesKey",
                 cols_minCoordinateUncertaintyInMeters = "coordinateUncertaintyInMeters")



ts_obs_rich <- obs_richness_ts(taxon_cube, first_year=1990)
plot(ts_obs_rich)
list_species(taxon_cube)
ts_turnover<-occ_turnover_ts(taxon_cube)

plot(ts_turnover)
#ts_density<-occ_density_ts(taxon_cube)
plot(ts_density)

ts_evenness<-pielou_evenness_ts(taxon_cube)

plot(ts_evenness)

period<-sort(unique(taxon_cube$data$year))

sbs.taxon_list<-list()
names(sbs.taxon_list)<-paste0("year","_",period[-c(1:12)])

sbs.fun<-function(y){
  sbs.taxon<-taxon_cube$data %>%
    filter(year==y) %>% 
    dplyr::select(scientificName,cellCode,obs) %>%
    group_by(scientificName,cellCode) %>%
    summarise(across(obs, sum), .groups = "drop") %>%
    pivot_wider(names_from = scientificName, values_from = obs) %>%
    arrange(cellCode) %>% 
    column_to_rownames(var = "cellCode") 
  return(sbs.taxon)
}

sbs.taxon_list<-map(period[-c(1:12)],sbs.fun)
for(y in period[-c(1:12)]){
  
  
  sbs.taxon_list[[paste0("year","_",y)]]<-sbs.taxon
  # 
  # species_list<-unique(names(sbs.taxon))
  # 
  # if(!exists("taxon_status_list")){
  #   full_species_list<-sort(unique(taxon_cube$data$scientificName))
  #   taxon_status_list<-taxon_status(species_list = full_species_list,
  #                                   source = "WCVP",
  #                                   region = "South Africa")
  # }
  # 
  # intro.sf<-taxon_cube$data %>% 
  #   filter(year==y) %>% 
  #   left_join(taxa_list_status,
  #             by = c("scientificName" = "taxon"))
  # 
  # 
  # status.sf <- intro.sf %>%
  #   group_by(cellCode) %>%
  #   summarise(
  #     total_intro_obs = sum(obs[introduction_status == "introduced"], na.rm = TRUE),
  #     total_native_obs = sum(obs[introduction_status == "native"], na.rm = TRUE),
  #     .groups = "drop"
  #   ) %>% 
  #   mutate(across(c(total_intro_obs, total_native_obs), ~ ifelse(.==0,NA,.))) %>% 
  #   mutate(intro_native=total_intro_obs/total_native_obs) %>% 
  #   arrange(cellCode)
  # if (!exists("eicat_score_list")){
  #   eicat_score_list=eicat_impact(eicat_data = eicat_data,species_list = species_list,
  #                                 fun="max")
  }
 
  eicat_score<-eicat_score_list[species_list,]
  
  siteScore<-status.sf$intro_native
  
  impact_metrics<-list(year=y,sbs.taxon=sbs.taxon,eicat_score=eicat_score,siteScore)
  impact_list[[as.character(y)]]<-impact_metrics
}
sbs.taxon<-taxon_cube$data %>%
  filter(year==y) %>% 
  dplyr::select(scientificName,cellCode,obs) %>%
  group_by(scientificName,cellCode) %>%
  summarise(across(obs, sum), .groups = "drop") %>%
  pivot_wider(names_from = scientificName, values_from = obs) %>%
  arrange(cellCode) %>% 
  column_to_rownames(var = "cellCode") 


species_list<-unique(names(sbs.taxon))

full_species_list<-sort(unique(taxon_cube$data$scientificName))
taxon_status_list<-taxon_status(species_list = full_species_list,
                                source = "WCVP",
                                region = "South Africa")


test_data<-data.frame(name=species_list[1:100],
                      status=sample(c("introduced","native"),100,replace = T))

taxon_status_list<-taxon_status(species_list = species_list,
                                source = "manual",
                                status_data = test_data)



intro.sf<-taxon_cube$data %>% 
  filter(year==y) %>% 
  left_join(taxa_list_status,
            by = c("scientificName" = "taxon"))


status.sf <- intro.sf %>%
  group_by(cellCode) %>%
  summarise(
    total_intro_obs = sum(obs[introduction_status == "introduced"], na.rm = TRUE),
    total_native_obs = sum(obs[introduction_status == "native"], na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(across(c(total_intro_obs, total_native_obs), ~ ifelse(.==0,NA,.))) %>% 
  mutate(intro_native=total_intro_obs/total_native_obs) %>% 
  arrange(cellCode)


#intro.sf<-left_join(intro.sf,status.sf, by="cellCode")
  
eicat_impact<-function(eicat_data,
                       species_list,
                       col_impact=NULL,
                       col_name=NULL,
                       fun="max"){
  
  
  if(all(c("impact_category","scientific_name")%in%names(eicat_data))){
    eicat_data<-eicat_data
  } else if((!is.null(col_impact)&!is.null(col_name))){
    eicat_data <- eicat_data %>% 
      rename(all_of(c(impact_category=col_impact, scientific_name=col_name)))
    
  } else{ stop("required column is not given")}
  
  
  if(fun=="max"){
    f<-function(x) max(x,na.rm = T)
  } else if(fun=="min"){
    f<-function(x) min(x,na.rm = T)
  } else if(fun=="mean"){
    f<-function(x) mean(x,na.rm = T)
  } else {stop("`fun` should be max, min or mean character")}
  
  category_M = eicat_data %>% 
    mutate(impact_category=substr(impact_category,1,2)) %>% 
    filter(impact_category %in% c("MC","MN","MO","MR","MV")) %>% 
    select(scientific_name,
           impact_category) %>% 
    mutate(category_value = case_when(
      impact_category == "MC" ~ 0,
      impact_category == "MN" ~1,
      impact_category == "MO" ~2,
      impact_category == "MR" ~3,
      impact_category == "MV" ~4,
      TRUE ~ 0  # Default case, if any value falls outside the specified ranges
    )) %>% 
    group_by(scientific_name,impact_category) %>%
    #choose the first trait value if there are multiples trait for a species
    summarise(across(category_value, first), .groups = "drop") %>% 
    # reshape to wide format to have specie by trait dataframe
    pivot_wider(names_from = impact_category, values_from = category_value) %>% 
    # select species that are only present in gbif data
    filter(scientific_name %in% species_list) %>% 
    #convert species names to row names
    column_to_rownames(var = "scientific_name") 
  
  
  category_M<-data.frame("fun"=apply(category_M,1,f))
  #names(category_M)<-fun
  na.df<-as.data.frame(matrix(NA,
                              nrow = length(setdiff(species_list,rownames(category_M))),
                              ncol = ncol(category_M)))
  row.names(na.df)<-setdiff(species_list,rownames(category_M))
  names(na.df)<-names(category_M) # column names
  category_M<-rbind(category_M,na.df)
  category_M <- category_M %>%
    dplyr::mutate(rowname = row.names(.)) %>%  
    dplyr::arrange(rowname) %>%               
    dplyr::select(-rowname)                    
  return(category_M)
    
    
}

eicat_score=eicat_impact(eicat_data = eicat_data,species_list = species_list,
                         fun="max")
eicat<-eicat_score[species_list,]

siteScore<-status.sf$intro_native

abdundance_impact = sweep(sbs.taxon,2,eicat,FUN = "*")
impactScore = siteScore*abdundance_impact
impact<-sum(impactScore,na.rm = TRUE)

# specieImpact<-colSums(impactScore,na.rm = TRUE)
# specieImpact<-specieImpact[specieImpact>0]
# 
# siteImpact<-rowSums(impactScore,na.rm = TRUE)
# siteImpact<-siteImpact[siteImpact>0]



  

