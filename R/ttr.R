
##############################################################################################
################################### Input functions ##########################################
##############################################################################################



plotTTRs = function(TTRs, RTdat, profNum, sizeCut, RTCut, yspan, xspan, chr, plotLoess, plotPoly) { 					# Define function to plot TTRs
	# Separate lists of positive and negative slope TTRs
	hTTR = data.frame(TTRs[1])
	lTTR = data.frame(TTRs[2])
	# Add column with TTR genomic length to each list
	hTTR$size = hTTR$highEnd - hTTR$highStart
	lTTR$size = lTTR$lowEnd - lTTR$lowStart
	# Add column with RT change across TTR to each list
	hTTR$RT = hTTR$highRTend - hTTR$highRTstart
	lTTR$RT = lTTR$lowRTstart - lTTR$lowRTend
	# Subsets of each list with given region of chromosome (xspan) that meet size and RT change criteria
	hTTRb = subset(hTTR, CHR == chr & highEnd < xspan[2] & highStart > xspan[1] & size > sizeCut & RT > RTCut)
	lTTRb = subset(lTTR, CHR == chr & lowEnd < xspan[2] & lowStart > xspan[1] & size > sizeCut & RT > RTCut)

	# Subset of dataset from which TTRs were called with given region of chromosome (xspan)	
	RTc = subset(RTdat, CHR == chr & POSITION < xspan[2] & POSITION > xspan[1] )
	
	# Prepare plot area
	par(mar=c(2,4.1,2,1), mfrow=c(2,1), pch=16, xaxt="n", xaxs="i", yaxs="i", lwd=2,  lty=1)
	if(plotLoess) { 													# Will plot smoothed line; requires collectLoess=T in "defineTTRs" function
	par(mfrow=c(3,1))
	plot(RTc[,profNum]~RTc$POSITION, ylim=yspan, pch=19, cex=0.5, col="gray")
	lines(xC$RT~xC$POS,  col="gray30",lwd=4)
	points(xH$RT~xH$POS, col="yellow", pch=19,cex=0.75)
	points(xL$RT~xL$POS, col="blue4", pch=19, cex=0.75)
}
	if(plotPoly) {														# Will plot TTR regions rather than boundary lines
		plot(RTc[,profNum]~RTc$POSITION, ylim=yspan, pch=19, cex=0.1, col="gray")
		lines(xC$RT~xC$POS,  col="gray30",lwd=4)

		for (i in 1:nrow(hTTRb)) {		 	
			hStart = hTTRb$highStart[i]; hEnd = hTTRb$highEnd[i]
			polygon(x=c(hStart, hStart, hEnd, hEnd), y=c(5, -5, -5, 5), lwd=1, col=rgb(0.90,0.90,0.0,0.2))
		}
		for (i in 1:nrow(lTTRb)) { 
			lStart = lTTRb$lowStart[i]; lEnd = lTTRb$lowEnd[i]
			polygon(x=c(lStart, lStart, lEnd, lEnd), y=c(5, -5, -5, 5), lwd=1, col=rgb(0,0,1,0.2))
		}
		} else { 
			plot(RTc[,profNum]~RTc$POSITION, pch=19, cex=0.5, col="grey")
			abline(v=hTTRb$highStart, col="yellow2")	
			abline(v=hTTRb$highEnd, col="yellow4")
			abline(v=lTTRb$lowStart, col="blue2")
			abline(v=lTTRb$lowEnd, col="blue4")
		}
	}

##############################################################################################
##################### Define & plot TTRs from RT datasets ####################################
##############################################################################################

#Read in list of normalized replication timing (RT) datasets
setwd("/Directory/")
DS = NULL; DS <- read.table("Qnorm RT dataset list_wo avgs.txt", header=T, nrows=1e3, comment.char = "");head(DS);dim(DS);summary(DS)

# Set cutoffs to define Timing Transition Regions (TTRs)
cutoffs = NULL
cutoffs$slope = 2.75e-6			# Slope of early TTR border
cutoffs$slopeL = 1e-06			# Slope of late TTR border

# Read and call TTRs for each RT dataset in list
RTi = NULL; RTi <- read.table("http://bio.fsu.edu/~dvera/share/gilbert/pooledqnorm/hg19-720k_IMR90_qnormToPool.txt", header=F, nrows=1e6, comment.char = "");head(RTi); dim(RTi)																											# Read example dataset
chrs = NULL; chrs = levels(RTi$CHR)[1:23]																		# Create list of chromosomes
types = NULL; types = c("hg19_720k_array","hg19_2M_array") 														# Create list of dataset formats
for(type in types){																								# For each dataset format
	D = NULL; Dsets = NULL; Dsets = DS$Location[DS$Data_Type == type]											# Subset of RT datasets with current format
	
	# Depending on dataset format, define smoothing span "sp"
	if(type == "hg19_720k_array"){
		sp = NULL; sp = 120
		}else{
		sp = NULL; sp = 185																						# Empirically determined optimum span for 2M probe array
	}
	
	for(D in 1:length(Dsets)){																					# For each dataset
		RTx = NULL; RTx = read.table(as.character(Dsets[D]), header=F, nrows=3e6, comment.char = "")			# Read dataset
		names(RTx)[1:2] = c("CHR","POSITION")																	# Name chromosome and position columns
		cat("Current dataset: ", as.character(Dsets[D]), "\n")													# Output current dataset name
		TTRs = NULL; TTRs = defineTTRs(																			# Call TTRs in current dataset
			RTa=RTx,			# Dataset containing profile
			profNum=3,			# Profile column number
 			lspan=sp,			# Span for loess smoothing
			collectLoess=F)		# Collect loess curves for plotting?		
		write.table(TTRs[1], file=paste(type, "high slope TTRs_", D, ".txt", sep=" "), row.names=F, quote=F)	# Write list of TTRs in current dataset with positive slope
		write.table(TTRs[2], file=paste(type, "low slope TTRs_", D, ".txt", sep=" "), row.names=F, quote=F)		# Write list of TTRs in current dataset with negative slope
	}
}

# Define and plot TTRs for one profile
chrs = NULL; chrs = levels(RTi$CHR)[1]		# Set chromosome
TTRs = NULL; TTRs = defineTTRs(				# Call TTRs on set chromosome for following profile
	RTa=RTi,								# Dataset containing profile
	profNum=3,								# Profile column number
 	lspan=120,								# Data points for loess smoothing span
	collectLoess=T)							# Collect loess curves for plotting?	

plotTTRs(TTRs=TTRs,							# TTR output from defineTTRs
	RTdat = RTi,							# RT dataset to plot
 	profNum = 3,							# Profile number to plot
 	sizeCut = 0,							# Minimum TTR genomic length
 	RTCut = 0,								# Minimum RT change across TTR
	yspan = c(-2,2),						# y-axis limits (RT)
	xspan = c(10e6,25e6),					# x-axis limits (bp)
	chr = chrs,								# Chromosome to plot
	plotLoess = T,							# Plot the loess smoothed line; requires collectLoess=T in "defineTTRs" function
	plotPoly = T)							# Plot TTR regions rather than boundary lines

##############################################################################################
############################# Create combined TTR list #######################################
##############################################################################################	

#Read in list of normalized replication timing (RT) datasets
DS = NULL; DS <- read.table("Qnorm RT dataset list_wo avgs.txt", header=T, nrows=1e3, comment.char = "");head(DS);dim(DS);summary(DS)

TTRi = NULL
for(type in levels(DS$Data_Type)){																	# For each dataset format
	D = NULL; Dsets = NULL; Dsets = DS$Location[DS$Data_Type == type]								# Subset of names from dataset list with current format
	
	for(D in 1:length(Dsets)){																		# For each dataset with current format
		cat("Current dataset: ", as.character(Dsets[D]), "\n")										# Output name of current dataset
		
		dsName = NULL; dsName = strsplit(as.character(Dsets[D]),"_")[[1]][2]						# Extract dataset name
		repNum = NULL; repNum = strsplit(as.character(Dsets[D]),"_")[[1]][4]						# Extract replicate number
		
		# Read positive slope TTRs for current dataset
		Rv = NULL; Rv = read.table(paste(type, "high slope TTRs_", D, ".txt", sep=" "), header=T, nrows=1e5, comment.char = "")
		# Read negative slope TTRs for current dataset
		Fw = NULL; Fw = read.table(paste(type, "low slope TTRs_", D, ".txt", sep=" "), header=T, nrows=1e5, comment.char = "")
		
		# Add columns for dataset number (SET), format (TYPE), name (SOURCE), replicate number (REP), and orientation (FLIP = T for positive slope TTRs)
		Rv1 = NULL; Rv1 = data.frame(SET=D,TYPE=type,SOURCE=dsName,REP=repNum,CHR=Rv$CHR,POSITION=Rv$highEnd,FLIP=T,SIZE=(Rv$highEnd-Rv$highStart),RT_CHANGE=(Rv$highRTend-Rv$highRTstart),RT=Rv$highRTend)
		Fw1 = NULL; Fw1 = data.frame(SET=D,TYPE=type,SOURCE=dsName,REP=repNum,CHR=Fw$CHR,POSITION=Fw$lowStart,FLIP=F,SIZE=(Fw$lowEnd-Fw$lowStart),RT_CHANGE=(Fw$lowRTstart-Fw$lowRTend),RT=Fw$lowRTstart)

		# Merge lists of positive and negative slope TTRs and add them to combined list of previous datasets
		TTRi = rbind(TTRi, Rv1, Fw1)
	}
}

hT = TTRi[TTRi$SIZE >= 250e3,]														# Filter TTRs smaller than 250kb
H = cbind(1:nrow(hT),hT); names(H)[1]="ID"											# Add column with unique identification number (ID) for each TTR

#Write to file
write.table(H, file="Hg19 TTR early borders_250kb_unaligned.txt", row.names=F, quote=F)

##############################################################################################
######################## Average TTRs common to multiple datasets ############################
##############################################################################################	

# Reclassify TTRs by cell type rather than sample name
levels(H$SOURCE)
Ho = NULL; Ho = H[order(H$CHR, H$POSITION),]; head(Ho)
Ho$SOURCE[Ho$SOURCE == "Name1" | Ho$SOURCE == "Name2"] <- "Type1"
Ho$SOURCE[Ho$SOURCE == "NameA" | Ho$SOURCE == "NameB"] <- "Type2"
Ho = droplevels(Ho); summary(Ho); levels(Ho$SOURCE)

# Combine any boundaries with the same orientation within 105kb of each other
chrs = NULL; chrs = levels(Ho$CHR)
HL = NULL
for(chr in chrs){																			# For each chromosome
	cat("Current chromosome: ", as.character(chr), "\n")
	Ht = NULL; Ht = Ho[Ho$CHR == chr & Ho$FLIP == T,]										# Subset of positive slope TTRs on current chromosome
	Hf = NULL; Hf = Ho[Ho$CHR == chr & Ho$FLIP == F,]										# Subset of negative slope TTRs on current chromosome
	A = NULL; A = rbind(A,Ht[1,])															# Start list of positive slope TTRs to average beginning with the first
	for(i in 2:nrow(Ht)){																	# For each TTR after the first
		if(Ht$POSITION[i] - Ht$POSITION[i-1] <= 105e3){										# If the current early TTR border is within 105kb of the previous early border
			A = rbind(A,Ht[i,])																# Add current TTR to list for averaging
			}else{								
			S = NULL; S = levels(droplevels(A$SOURCE))										# Otherwise, extract cell types of TTRs to average
			HL = rbind(HL, data.frame(														# and add to data frame "HL" 
				CHR=chr,																	# chromosome
				POSITION=mean(A$POSITION),													# average early TTR border position
				FLIP=T,																		# orientation
				Type1="Type1" %in% S,														# and cell type (whether TTR is present in each cell type) information 
				Type2="Type2" %in% S														# for the averaged boundary
				))
			A = NULL; A = rbind(A,Ht[i,])													# and restart the list of TTRs to average beginning with the current TTR
		}
	}
	S = NULL; S = levels(droplevels(A$SOURCE))												# Extract cell types of TTRs remaining on list to average
	HL = rbind(HL, data.frame(																# and add to data frame "HL"
		CHR=chr,																			# chromosome
		POSITION=mean(A$POSITION),															# average early TTR border position
		FLIP=T,																				# orientation
		Type1="Type1" %in% S,																# and cell type (whether TTR is present in each cell type) information 
		Type2="Type2" %in% S																# for the averaged boundary
		))

	A = NULL; A = rbind(A,Hf[1,])															# Start list of negative slope TTRs to average beginning with the first
	for(i in 2:nrow(Hf)){																	# For each TTR after the first
		if(Hf$POSITION[i] - Hf$POSITION[i-1] <= 105e3){										# If the current early TTR border is within 105kb of the previous early border
			A = rbind(A,Hf[i,])																# Add current TTR to list for averaging
			}else{
			S = NULL; S = levels(droplevels(A$SOURCE))										# Otherwise, extract cell types of TTRs to average
			HL = rbind(HL, data.frame(														# and add to data frame "HL"
				CHR=chr,																	# chromosome
				POSITION=mean(A$POSITION),													# average early TTR border position
				FLIP=F,																		# orientation
				Type1="Type1" %in% S,														# and cell type (whether TTR is present in each cell type) information 
				Type2="Type2" %in% S														# for the averaged boundary				
				))
			A = NULL; A = rbind(A,Hf[i,])													# and restart the list of TTRs to average beginning with the current TTR
		}
	}
	S = NULL; S = levels(droplevels(A$SOURCE))												# Extract cell types of TTRs remaining on list to average
	HL = rbind(HL, data.frame(																# and add to data frame "HL"
		CHR=chr,																			# chromosome
		POSITION=mean(A$POSITION),															# average early TTR border position
		FLIP=F,																				# orientation
		Type1="Type1" %in% S,																# and cell type (whether TTR is present in each cell type) information 
		Type2="Type2" %in% S																# for the averaged boundary	
		))	
}

HLo = HL[order(HL$CHR,HL$POSITION),]
dim(HLo); head(HLo)

# Combine any of the remainging boundaries within 140kb of each other
HL1 = NULL
for(chr in chrs){																			# For each chromosome
	cat("Current chromosome: ", as.character(chr), "\n")
	Ht = NULL; Ht = HLo[HLo$CHR == chr,]													# Subset of boundaries on current chromosome
	A = NULL; A = rbind(A,Ht[1,])															# Start list of boundaries to average beginning with the first
	for(i in 2:nrow(Ht)){																	# For each boundary after the first
		if(Ht$POSITION[i] - Ht$POSITION[i-1] <= 140e3){										# If the current boundary is within 140kb of the previous boundary
			A = rbind(A,Ht[i,])																# Add current boundary to list for averaging
			}else{
			HL1 = rbind(HL1, data.frame(													# Otherwise, add to data frame "HL1"
				CHR=chr,																	# chromosome
				POSITION=mean(A$POSITION),													# average boundary position
				Type1=length(A$Type1[A$Type1==T]) > 0,										# whether a TTR is present in cell type "Type1"
				Type1_FLIP=if(length(A$Type1[A$Type1==T]) == 0){"NONE"}else if(length(A$FLIP[A$Type1==T & A$FLIP==T]) > 0 & length(A$FLIP[A$Type1==T & A$FLIP==F]) > 0){"BOTH"}else if(length(A$FLIP[A$Type1==T & A$FLIP==T]) > 0){"YES"}else{"NO"},				# orientation in cell type "Type1"
				Type2=length(A$Type2[A$Type2==T]) > 0,										# whether a TTR is present in cell type "Type2"
				Type2_FLIP=if(length(A$Type2[A$Type2==T]) == 0){"NONE"}else if(length(A$FLIP[A$Type2==T & A$FLIP==T]) > 0 & length(A$FLIP[A$Type2==T & A$FLIP==F]) > 0){"BOTH"}else if(length(A$FLIP[A$Type2==T & A$FLIP==T]) > 0){"YES"}else{"NO"}				# and orientation in cell type "Type2"
				))
			A = NULL; A = rbind(A,Ht[i,])													# and restart the list of boundaries to average beginning with the current one
		}
	}
	HL1 = rbind(HL1, data.frame(															# For the remaining boundaries to average, add to data frame "HL1"
		CHR=chr,																			# chromosome
		POSITION=mean(A$POSITION),															# average boundary position
		Type1=length(A$Type1[A$Type1==T]) > 0,												# whether a TTR is present in cell type "Type1"
		Type1_FLIP=if(length(A$Type1[A$Type1==T]) == 0){"NONE"}else if(length(A$FLIP[A$Type1==T & A$FLIP==T]) > 0 & length(A$FLIP[A$Type1==T & A$FLIP==F]) > 0){"BOTH"}else if(length(A$FLIP[A$Type1==T & A$FLIP==T]) > 0){"YES"}else{"NO"},							# orientation in cell type "Type1"
		Type2=length(A$Type2[A$Type2==T]) > 0,												# whether a TTR is present in cell type "Type2"
		Type2_FLIP=if(length(A$Type2[A$Type2==T]) == 0){"NONE"}else if(length(A$FLIP[A$Type2==T & A$FLIP==T]) > 0 & length(A$FLIP[A$Type2==T & A$FLIP==F]) > 0){"BOTH"}else if(length(A$FLIP[A$Type2==T & A$FLIP==T]) > 0){"YES"}else{"NO"}								# and orientation in cell type "Type2"
		))
}
HL1o = HL1[order(HL1$CHR,HL1$POSITION),]; dim(HL1o)
