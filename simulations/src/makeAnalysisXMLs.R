pibuss_xml <- list.files("xml/piBUSS",full.names=TRUE)

beastjar <- "~/git_repos/beast-mcmc/build/dist/beast.jar"
beaglepath <- "/usr/local/lib"

for (i in 1:length(pibuss_xml)) {
  # Simulate
  system(paste0("java  -Djava.library.path=",beaglepath," -jar ",beastjar," -working -overwrite ",pibuss_xml[i]))
  
  # Get simulated data
  seq_xml <- list.files("xml/piBUSS",full.names=TRUE)
  seq_xml <- seq_xml[grepl("sequence",seq_xml)]
  simulated_data <- scan(seq_xml,sep="\n",what=character(),blank.lines.skip=TRUE)
  
  # Get XML template for analysis
  id <- basename(pibuss_xml[i])
  id <- strsplit(id,"_frac")[[1]][1]
  if ( id == "nstates_44" ) {
    template_file <- "xml/original/real_data_44_state_analysis.xml"
  } else {
    template_file <- list.files("xml/condensed_original/",full.names=TRUE)
    template_file <- template_file[grepl(gsub("nstates_","",id),template_file)]
  }
  template <- scan(template_file,sep="\n",what=character(),blank.lines.skip=FALSE)
  
  # Make datatype codes upper-case because they come out of piBUSS that way
  dt_open <- grep('<generalDataType id="sampleLoc.dataType">',template)
  dt_close <- grep('</generalDataType>',template)
  for (l in (dt_open+1):(dt_close-1)) {
    line <- template[l]
    has_code <- grepl("code",line)
    has_alias <- grepl("alias",line)
    line <- strsplit(line,'"')[[1]]
    if (has_code) {
      line[grep("code",line)+1] <- toupper(line[grep("code",line)+1])
    }
    if (has_alias) {
      line[grep("alias",line)+1] <- toupper(line[grep("code",line)+1])
    }
    template[l] <- paste0(line,collapse='"')
  }
  
  # Locate where we need to put simulated data and add it
  attribute_ends <- grep('</attr>',template)
  sample_loc_starts <- grep('<attr name="sampleLoc">',template)
  sample_ori_starts <- grep('<attr name="oriLoc">',template)
  
  taxa <- simulated_data[grepl('<taxon idref=',simulated_data)]
  taxa <- gsub('\"/>','',taxa)
  taxa <- do.call(rbind,strsplit(taxa,'idref=\"'))[,2]
  all_states <- c()
  for (taxon in taxa) {
    charstate <- simulated_data[grep(taxon,simulated_data,fixed=TRUE)+1]
    charstate <- trimws(charstate)
    all_states <- c(all_states,charstate)
    charstate <- paste0("\t\t\t\t",charstate)
    
    taxon_line <- grep(paste0('taxon id="',taxon),template,fixed=TRUE)
    
    taxon_loc <- min(sample_loc_starts[sample_loc_starts > taxon_line])
    template[taxon_loc+1] <- charstate
    
    taxon_ori <- min(sample_ori_starts[sample_ori_starts > taxon_line])
    template[taxon_ori+1] <- charstate
  }
  
  cat(template,sep="\n",file=paste0("xml/simulated_datasets_to_analyze/",basename(pibuss_xml[i])))
  system(paste0("rm ",seq_xml))
}