defineTTRs = function(RTa, profNum, lspan, threads = 10 , eslope = 2.75e-6 , lslope = 0.1 , collectLoess = FALSE) {	# Define function to call TTRs
	library(parallel)
	highTTR = NULL
	lowTTR = NULL


	cat("Defining TTRs for", names(RTa)[profNum], ":", "\n")# Output name of current dataset

	chrs<-unique(RTa[,1])

	fttrl<-mclapply(chrs, function(chr){
	# For each chromosome
		cat("Current chromosome: ", chr, "\n")
		# Output name of current chromosome
		RTb   = RTa[which(RTa$CHR == chr),]
		# Subset of data in current dataset from current chromosome


		# smooth chromosome data
		loessSpan  = lspan/dim(RTb)[1]# Define smoothing span (proportion of data points)
		loessProf  = loess(RTb[,profNum] ~ RTb$POSITION, span=loessSpan)$fitted	# Smooth data from current chromosome using defined span
		cat("loess done\n")
		loessProf = loessProf * 1.59 / IQR(loessProf) 	# Normalize y-range IQR to 1.59 

		# Index of positions with > 20kb gaps
		aboveGap   = (c(0,RTb$POSITION)-c(RTb$POSITION,0)) < -20e3
		# RT slope for positions with no gaps
		loessSlope = (diff(loessProf)/diff(RTb$POSITION))[!aboveGap]

		# Create data frame with slope and RT for all probe positions
		df = data.frame(loessSlope, RT=loessProf[!aboveGap], POS=RTb$POSITION[!aboveGap])	 
		# Remove last two rows
		df = df[1:(nrow(df)-2),]

		# Add columns to indicate postions with slopes above early border cutoff (=TRUE)
		df$highSlope = df$loessSlope > eslope
		df$highSlope[1] = FALSE
		df$lowSlope  = df$loessSlope < -eslope
		df$lowSlope[1] = FALSE

		# Add columns to indicate postions with slopes above late border cutoff (=TRUE)
		df$highSlopeL = df$loessSlope > lslope
		df$highSlopeL[1] = FALSE
		df$lowSlopeL = df$loessSlope < -lslope
		df$lowSlopeL[1] = FALSE

		# Add placeholder columns to indicate postions of early and late TTR borders
		df$highStart = FALSE
		df$highEnd = FALSE
		df$lowStart  = FALSE
		df$lowEnd  = FALSE

		# For positive slopes, define as early TTR border if above cutoff but not next slope to right
		df$highEnd	<- (c(F, df$highSlope) & !c(df$highSlope, F))[1:nrow(df)]
		# and define as late TTR border if above cutoff but not next slope to left
		df$highStart <- (c(df$highSlopeL, F) & !c(F, df$highSlopeL))[1:nrow(df)]
		# For negative slopes, define as early TTR border if above cutoff but not next slope to left
		df$lowStart  <- (c(df$lowSlope, F) & !c(F, df$lowSlope))[1:nrow(df)]
		# and define as late TTR border if above cutoff but not next slope to right
		df$lowEnd	<- (c(F, df$lowSlopeL) & !c(df$lowSlopeL, F))[1:nrow(df)]
		# For positive slopes, eliminate extra early border calls within TTRs
		# Create list of early TTR border positions
		


		HE = NULL; HE = grep(TRUE, df$highEnd)
		# Start list of late TTR border positions with a value > the right end of the chromosome
		HS = NULL; HS = 1e10
		# For each early TTR border starting with the farthest to the right and working backwards
		
		print(length(HE))
		

		for(i in length(HE):1){	

			# If current early border is left of the current late border (indicates start of new TTR), move to next late TTR border by adding it to list "HS"
			# next late border will be furthest to right but still left of current early border
			if(HE[i] < HS[length(HS)]){HS = c(HS, max(grep(TRUE, df$highStart)[grep(TRUE, df$highStart) < HE[i]]))
				# Otherwise set current early border to zero
				}else{HE[i] = 0
				}
		}
		# Eliminate early TTR borders set to zero
		HE = subset(HE, HE != 0)
		# Eliminate starting (artificial) value in late TTR border list
		HS = subset(HS, HS != 1e10)
		# Invert late TTR border list to account for working backwards
		HS = HS[length(HS):1]	

		# For negative slopes, eliminate extra early border calls within TTRs
		# Create list of early TTR border positions
		LS = NULL; LS = grep(TRUE, df$lowStart)	
		# Start list of late TTR border positions with a value < the left end of the chromosome
		LE = NULL; LE = 0
		# For each early TTR border
		for(i in 1:length(LS)){	

		# If current early border is right of the current late border (indicates start of new TTR), move to next late TTR border by adding it to list "LE"
		# next late border will be furthest to left but still right of current early border
		if(LS[i] > LE[length(LE)]){LE = c(LE, min(grep(TRUE, df$lowEnd)[grep(TRUE, df$lowEnd) > LS[i]]))
			# Otherwise set current early border to zero
			}else{LS[i] = 0
			}
		}
		# Eliminate early TTR borders set to zero
		LS = subset(LS, LS != 0)
		# Eliminate starting (artificial) value in late TTR border list
		LE = subset(LE, LE != 0)

		# Extract positions of early and late TTR borders
		# For positive slopes
		aH = df$POS[HS]
		bH = df$POS[HE]
		# For negative slopes
		aL = df$POS[LS]
		bL = df$POS[LE]

		#print(length(aH))
		#print(length(bH))
		#print(length(aL))
		#print(length(bL))


		# Extract RT of early and late TTR borders
		# For positive slopes
		aRTh = df$RT[HS]
		bRTh = df$RT[HE]
		# For negative slopes
		aRTl = df$RT[LS]
		bRTl = df$RT[LE]

		# Create data frames with TTR data on current chromosome
		z1 = data.frame(CHR=chr, highStart=aH, highEnd=bH, highRTstart=aRTh, highRTend=bRTh)# For positive slopes
		z2 = data.frame(CHR=chr, lowStart=aL, lowEnd=bL, lowRTstart=aRTl, lowRTend=bRTl)	# For negative slopes

		ttrl<-list(z1,z2)
		return(ttrl)
		# Add TTR data from current chromosome to previous chromosomes
		
		# Collect smoothing data for plotting 
		if(collectLoess) {
			xC <<- df
			xH <<- df[df$highSlope,]
			xL <<- df[df$lowSlope,]
		}
	},mc.cores=threads,mc.preschedule=F)
	print(fttrl)
	highTTR<-do.call(rbind,lapply(fttrl,"[",1))
	lowTTR<-do.call(rbind,lapply(fttrl,"[",2))
}