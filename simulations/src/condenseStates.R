# All 44 locations
locations_44 <- c("Australia","Belgium","Brazil","Cambodia","Canada","ChinaAnhui","ChinaBeijing","ChinaChongqing","ChinaFujian","ChinaGuangdong","ChinaHenan","ChinaHongKong","ChinaHubei","ChinaHunan","ChinaJiangsu","ChinaJiangxi","ChinaShandong","ChinaSichuan","ChinaYunnan","ChinaZhejiang","Finland","France","Germany","India","Iran","Italy","Japan","Luxembourg","Mexico","Nepal","Netherlands","NewZealand","Nigeria","Portugal","Singapore","SouthKorea","Spain","Sweden","Switzerland","Taiwan","Thailand","UK","USA","Vietnam")

### If making new smaller regions sets, leave Hubei in as its own area, this script is not designed to remove any of the covariates

# A smaller 31 location set with only one for China
locations_31 <- c(locations_44[1:5],locations_44[13],"ChinaOther",locations_44[21:44])
map_44_to_31 <- as.list(match(locations_31,locations_44))
map_44_to_31[[7]] <- which(grepl("China",locations_44) & (!grepl("Hubei",locations_44)))

# China only
locations_15 <- locations_44[grepl("China",locations_44)]
map_44_to_15 <- match(locations_15,locations_44)

# The covariates 
covariates <- c("intraContinent_distances","airTravel","Hubei_assymetry")

### If adding more, make sure that the state space size is unique and the name matches this format exactly

# Assemble
location_list <- list(locations_31,locations_15)
map_list <- list(map_44_to_31,map_44_to_15)
xml_names <- c("real_data_31_state_analysis.xml","real_data_15_state_analysis.xml")

full_xml <- scan("xml/original/real_data_44_state_analysis.xml",sep="\n",what=character(),blank.lines.skip=FALSE)

for (i in 1:length(location_list)) {
  this_xml <- full_xml
  
  nstates <- length(location_list[[i]])
  
  # Change DataType
  dt_open <- grep('<generalDataType id="sampleLoc.dataType">',this_xml)
  dt_close <- grep('</generalDataType>',this_xml)
  before <- this_xml[1:dt_open]
  after <- this_xml[dt_close:length(this_xml)]
  new_dt <- paste0('\t\t<state code="',location_list[[i]],'"/>')
  if ( any(lengths(map_list[[i]]) > 1) ) {
    for (j in which(lengths(map_list[[i]]) > 1)) {
      new_dt <- c(new_dt,paste0('\t\t<alias state="',location_list[[i]][j],'" code="',locations_44[map_list[[i]][[j]]],'"/>'))
    }
  }
  this_xml <- c(before,new_dt,after)
  
  # Change dimensions
  dimtext <- "dimension=\"44\""
  this_line <- this_xml[grep(dimtext,this_xml)]
  this_line <- gsub(dimtext,paste0('dimension="',nstates,'"'),this_line)
  this_xml[grep(dimtext,this_xml)] <- this_line

  dimtext <- "dimension=\"1892\""
  this_line <- this_xml[grep(dimtext,this_xml)]
  this_line <- gsub(dimtext,paste0('dimension="',2*choose(nstates,2),'"'),this_line)
  this_xml[grep(dimtext,this_xml)] <- this_line

  # Fix covariates
  choose_44_2 <- choose(44,2)
  for (j in 1:3) {
    covariate_line_index <- grep(covariates[j],this_xml)
    this_covariate <- this_xml[covariate_line_index]
    this_covariate <- gsub("\"/>","",this_covariate)
    this_covariate <- strsplit(this_covariate,"\"")[[1]]
    this_covariate <- this_covariate[length(this_covariate)]
    old_covariate_string <- this_covariate
    this_covariate <- as.numeric(strsplit(this_covariate,"[[:space:]]")[[1]])
    
    covariate_matrix <- matrix(NA,44,44)
    covariate_matrix[lower.tri(covariate_matrix)] <- this_covariate[1:choose_44_2]
    covariate_matrix <- t(covariate_matrix)
    covariate_matrix[lower.tri(covariate_matrix)] <- this_covariate[choose_44_2+(1:choose_44_2)]
    
    new_matrix <- matrix(NA,nstates,nstates)
    for (k in 1:nstates) {
      for (l in 1:nstates) {
        if (l != k) {
          k_states <- map_list[[i]][[k]]
          l_states <- map_list[[i]][[l]]
          if (length(l_states) == 1 && length(k_states) == 1) {
            new_matrix[k,l] <- covariate_matrix[k_states,l_states]
          } else if (length(l_states) == 1) {
            k_mean <- mean(covariate_matrix[k_states,l_states])
            new_matrix[k,l] <- k_mean
          } else if (length(k_states) == 1) {
            l_mean <- mean(covariate_matrix[k_states,l_states])
            new_matrix[k,l] <- l_mean
          } else {
            the_mean <- 0.0
            the_count <- 0
            for (kidx in 1:length(k_states)) {
              for (lidx in 1:length(l_states)) {
                the_mean <- the_mean + covariate_matrix[k_states[kidx],l_states[lidx]]
                the_count <- the_count + 1
              }
            }
            the_mean <- the_mean/the_count
            new_matrix[k,l] <- the_mean
          }
        }
      }
    }
    lt <- lower.tri(new_matrix)
    first <- t(new_matrix)[lt]
    last <-new_matrix[lt]
    new_covariates <- paste0(c(first,last),collapse=" ")
    
    preamble <- strsplit(this_xml[covariate_line_index],"value")[[1]][1]
    new_covariate_string <- paste0(preamble,' value="',new_covariates,'"/>')
    
    this_xml[covariate_line_index] <- new_covariate_string
  }
  
  cat(this_xml,sep="\n",file=paste0("xml/condensed_original/",xml_names[i]))
}
