


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
  #filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
  #filter(coordinateUncertaintyInMeters<=res*1000) %>% 
  filter(!is.na(decimalLatitude)) %>% 
  sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
               crs = 4326) %>%
  st_join(grid) %>% 
  as.data.frame() %>% 
  select(-geometry) %>% 
  mutate(occurrenceStatus=1)



taxon_cube<-process_cube(taxon.sf,grid_type = "custom",cols_cellCode = "cellid", cols_year = "year",
                 cols_occurrences = "occurrenceStatus",cols_species = "species",cols_speciesKey = "speciesKey",
                 cols_minCoordinateUncertaintyInMeters = "coordinateUncertaintyInMeters")



ts_obs_rich <- obs_richness_ts(taxon_cube, first_year=1990)
plot(ts_obs_rich)
list_species(taxon_cube)
ts_turnover<-occ_turnover_ts(taxon_cube)

plot(ts_turnover)
ts_density<-occ_density_ts(taxon_cube)
plot(ts_density)

ts_evenness<-pielou_evenness_ts(taxon_cube)

plot(ts_evenness)

sbs.taxon<-taxon_cube$data %>%
  # select species name, trait and trait value
  dplyr::select(scientificName,cellCode,obs) %>%
  #group by Species and trait
  group_by(scientificName,cellCode) %>%
  #choose the first trait value if there are multiples trait for a species
  summarise(across(obs, sum), .groups = "drop") %>%
  # reshape to wide format to have specie by trait dataframe
  pivot_wider(names_from = scientificName, values_from = obs) %>%
  arrange(cellCode) %>% 
  # convert species names to row names
  column_to_rownames(var = "cellCode") 


species_list<-sort(unique(taxon_cube$data$scientificName))

# Download WCVP for native taxa in area of interest
native_list <- rWCVP::wcvp_checklist() %>%
  filter(area_code_l3 %in% rWCVP::get_wgsrpd3_codes("South Africa")) %>%
  filter((accepted_name %in% species_list) & occurrence_type=="native")

# create taxa list
taxa_list<-data.frame(taxon=species_list)

# create new dataframe with introduction status
taxa_list_status<-taxa_list%>%
  mutate(introduction_status = ifelse(taxon%in%native_list$accepted_name,
                                      "native","introduced"))


intro.sf<-taxon_cube$data %>% 
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


intro.sf<-left_join(intro.sf,status.sf, by="cellCode")
  
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
  names(category_M)<-fun
  na.df<-as.data.frame(matrix(NA,nrow = length(setdiff(species_list,rownames(category_M))),
                              ncol = ncol(category_M)))
  row.names(na.df)<-setdiff(species_list,rownames(category_M))
  names(na.df)<-names(category_M) # column names
  category_M<-rbind(category_M,na.df)
  category_M<-category_M[order(rownames(category_M)),]
  return(category_M)
    
    
}

eicat_score=eicat_impact(eicat_data = eicat_data,species_list = species_list,
                         fun="max")
eicat_score
site_pressure<-status.sf$intro_native

abdundance_impact = sweep(sbs.taxon,2,eicat_score,FUN = "*")
impactScore = site_pressure*abdundance_impact


specieImpact<-colSums(impactScore,na.rm = TRUE)
specieImpact<-specieImpact[specieImpact>0]

siteImpact<-rowSums(impactScore,na.rm = TRUE)
siteImpact<-siteImpact[siteImpact>0]



  

