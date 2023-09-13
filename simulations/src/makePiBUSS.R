filepaths <- c("xml/original/real_data_44_state_analysis.xml",list.files("xml/condensed_original",full.names=TRUE))
xmls <- lapply(filepaths,scan,sep="\n",what=character(),blank.lines.skip=FALSE)

# piBUSS block
pibuss_text <- 
'	<beagleSequenceSimulator id="simulator" parallel="false" output="XML">
		<partition from="1" to="1">
			<treeModel idref="treeModel"/>
			<glmSubstitutionModel idref="sampleLoc.model"/>
			<siteModel idref="sampleLoc.siteModel"/>
			<strictClockBranchRates idref="sampleLoc.branchRates"/>
			<frequencyModel idref="frequencyModel"/>
		</partition>
	</beagleSequenceSimulator>

	<report fileName="sequences.xml">
		<beagleSequenceSimulator idref="simulator"/>
	</report>
'

# Random-effects stuff
weak_mean <- 0.0
weak_sd <- 0.1

strong_abs_mean <- log(2.0)
strong_sd <- 0.1

percent_strong <- c(0.01,0.5,1.0)

cols <- viridis::rocket(2,begin=0.3,end=0.8)

pdf(file="simulating_distributions.pdf",width=4,height=3)
  par(mai=c(0.8,0.8,0.01,0.01))
  curve(dlnorm(x,-strong_abs_mean,weak_sd),0,3,1001,add=FALSE,col=cols[2],lwd=2,
        xlab="Rate-scale random-effect",ylab="Density")
  curve(dlnorm(x,weak_mean,weak_sd),0,3,1001,add=TRUE,col=cols[1],lwd=2)
  curve(dlnorm(x,strong_abs_mean,weak_sd),0,3,1001,add=TRUE,col=cols[2],lwd=2)
dev.off()
  
set.seed(42)
for (i in 1:3) {
  nstates <- basename(filepaths[i])
  nstates <- gsub("real_data_","",nstates)
  nstates <- gsub("_state_analysis.xml","",nstates)
  nstates <- as.numeric(nstates)
  refxdim <- 2 * choose(nstates,2)
  for (j in 1:3) {
    sim_name <- paste0("nstates_",nstates,"_frac_",percent_strong[j])
    
    mean_vec <- rep(weak_mean,refxdim)
    sd_vec <- rep(weak_sd,refxdim)
    
    # Allocate the specified number of strong random-effects
    is_strong <- logical(refxdim)
    is_strong[sample.int(refxdim,round(percent_strong[j]*refxdim))] <- TRUE
    
    mean_vec[is_strong] <- ((-1)^rbinom(sum(is_strong),1,0.5)) * strong_abs_mean
    sd_vec[is_strong] <- strong_sd
    
    # Draw the random-effects
    refx_drawn <- rnorm(refxdim,mean_vec,sd_vec)
    cat(refx_drawn,sep="\n",file=paste0("true_refx/",sim_name,".txt"))
    
    # Put things in the XML
    this_xml <- xmls[[i]]
    
    # Replace fixed-effect values with posterior means
    from <- '<parameter id="glmCoefficients" value="0.1 0.8 2.0"/>'
    to <- '<parameter id="glmCoefficients" value="0.04 0.76 0.27"/>'
    this_xml[grep(from,this_xml)] <- to
    
    # Replace clock rate with posterior mean
    from <- '<parameter id="sampleLoc.clock.rate" value="3.5" lower="0.0"/>'
    to <- '<parameter id="sampleLoc.clock.rate" value="4.02" lower="0.0"/>'
    this_xml[grep(from,this_xml)] <- to
    
    # Replace random-effects
    from <- paste0('<parameter id="glmRandCoefficients" dimension="',refxdim,'" value="0.0"/>')
    to <- paste0('<parameter id="glmRandCoefficients" dimension="',paste0(refx_drawn,collapse=""),'" value="0.0"/>')
    this_xml[grep(from,this_xml)] <- to
    
    # Remove detritus
    end_useful <- grep('</siteModel>',this_xml)
    this_xml <- this_xml[1:end_useful]
    
    # put simulator call in
    this_pibuss_text <- gsub("sequences.xml",
                             paste0("sequences_",sim_name,".xml"),
                             pibuss_text,fixed=TRUE)
    this_xml <- c(this_xml,this_pibuss_text)
    
    cat(this_xml,file=paste0("xml/piBUSS/",sim_name,".xml"),sep="\n")
  }
}