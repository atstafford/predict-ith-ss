# Function to set up 11 new data table in a list
dataSetup <- function(rawdata) {
  
  # Dataframe identifying start and stop codon, and chromosome for each bin
  start.stop <- rawdata[,c(1:3)]
  start.stop$bin <- 1:nrow(start.stop)
  start.stop <- start.stop[, c(4, 1, 2, 3)]
  
  # Vector holding sample ID
  samples <- colnames(rawdata)[-c(1:3)] 
  
  # Vector holding patient identifiers, ie the string before the '.' in sample IDs
  patients <- unique(sub("\\..*", "", colnames(rawdata)))[-c(1:3)] 
  
  # Number of samples
  noSamples <- length(samples)
  
  # Number of patients
  noPatients <- length(patients)
  
  # Number of bins
  noBins <- length(start.stop$bin)
  
  # Visualisation may require a vector to identify of chromosome ends and chromosome midpoints 
  chr.ends <- cumsum(table(start.stop$chr))
  list <- list()
  l <- 1
  for ( i in 1:length(chr.ends) ) {  
    if ( i == 1 ) { 
      list[[l]] <- chr.ends[i]/2 
      l <- l+1
    }
    else { 
      list[[l]] <- chr.ends[i-1] + ((chr.ends[i]-chr.ends[i-1])/2)
      l <- l+1
    }
  }
  chr.mid <- unlist(list)
  chr.ends <- data.frame(start=c(0,chr.ends[-22]), end=(chr.ends), col=c('G','W')) 
  
  # A dataframe with a bin per row and a patient per column, with values indicating clonality. There must only be 2 multiregion samples
  #(0=notCNA, 1=subclonalCNA, 2=clonalCNA)
  clonal.data <- rawdata[,-c(1:3)]
  k <- 1
  i <- 1
  l <- 1
  clonal <- list()
  for ( k in 1:nrow(clonal.data) ) {
    while ( i < ncol(clonal.data) ) {
      if ( clonal.data[k,i] != 2 | clonal.data[k,i+1] != 2) { #if one of the samples has a mutation, proceed
        if ( clonal.data[k,i] == clonal.data[k,i+1] ) { #both the same = clonal
          clonal[[l]] <- 2
          i <- i + 2
          l <- l + 1
        }
        else { #different = subclonal
          clonal[[l]] <- 1
          i <- i + 2
          l <- l + 1
        }
      }
      else { #neither sample has a mutation
        clonal[[l]] <- 0 
        i <- i + 2
        l <- l + 1
      }
    }
    i <- 1
  }
  clonal.data <- data.frame(t(matrix(unlist(clonal), ncol=2694)))
  colnames(clonal.data) <- patients
  clonal.data[] <- lapply(clonal.data, factor, levels=unique(unlist(clonal.data)))
  
  # Cannot use averaged raw copy number across the MR samples in case of variable MR sample per patient.
  # Will therefore can use a PIC score for each patient. However this is a continuous, not binary, measure of clonality
  # A dataframe with a bin per row and a patient per column, with values indicating pic score
  pic.data <- list()
  for (k in 1:length(patients)) {
    # grep to return indecies in my.data which to the unique patient names
    # paste("^",patient_names[1],sep="") : ^to ensure search at start of colname to prevent mises
    multiregion_samples <- grep(paste("^",patients[k],sep=""), colnames(rawdata)) 
    
    # count number of cols in nultiregion_samples (will be n in PIC calc)
    number_samples <- length(multiregion_samples)
    
    # store PIC output in PIC list
    l <- patients[k]
    pic.data[[l]] <- PIC(rawdata, number_samples, multiregion_samples)
  }
  pic.data <- data.frame(pic.data) 
  
  # Calulate average PIC per bin as a measure of how diverse that bin tends to be across patients
  ave.pic <- start.stop[,c(1:2)]
  ave.pic$avePic <- rowSums(pic.data)/(ncol(pic.data))
  
  # Calculate patient diversity (ITH), as the summed pic scored per patient divded by the maximum possible pic sum if all were diverse
  ith <- data.frame(sumPIC = colSums(pic.data, na.rm = TRUE), 
                    max.PIC = max(pic.data), 
                    binCount = colSums(!is.na(pic.data)))
  
  ith$max.totalPIC <- ith$max.PIC * ith$binCount
  ith$PIC.frac <- ith$sumPIC / ith$max.totalPIC
  
  # A dataframe of the propotion of genome gained/lost per sample, alongside ith
  pga <- data.frame(t(rawdata[,-c(1:3)]), check.names = FALSE)
  pga <- data.frame(prop.gain=apply(pga,1,function(x) sum(x == 3)/ncol(pga)),
                    prop.loss=apply(pga,1,function(x) sum(x == 1)/ncol(pga)))
  pga$prop.aneu <- pga$prop.gain + pga$prop.loss
  pga <- cbind( data.frame(ith=ith[rep(seq_len(nrow(ith)), each = 2), 5]), 
                pga)
  
  # Dataframe detailing the frequency of gains/losses and whether they are subclonal or clonal
  data <- rawdata[,-c(1:3)]
  GaLo.Clo <- data.frame(bin = start.stop$bin, chr1 = NA, chr2 = NA,
                         clonal.aneu = NA, subclonal.aneu = NA, gain = NA, loss = NA,
                         clonal.gain = NA, clonal.loss = NA, clonal.noCNA = NA, subclonal.gain = NA, subclonal.loss = NA)
  
  GaLo.Clo$chr1 <- as.numeric(start.stop[match(GaLo.Clo$bin, start.stop$bin), 2]) 
  GaLo.Clo$chr2 <- as.numeric(start.stop[match(GaLo.Clo$bin, start.stop$bin), 2]) 
  GaLo.Clo$chr2 <- factor(GaLo.Clo$chr2,levels = rev(c('1','2','3','4','5','6','7','8','9','10','11',
                                                       '12','13','14','15','16','17','18','19','20','21','22'))) #to order xaxis facets 
  
  for ( k in 1:nrow(data) ) { # for a bin
    i <- 1
    clonal.all <- clonal.gain <- clonal.loss <- clonal.noCNA <- subclonal.all <- subclonal.gain <- subclonal.loss <- subclonal.noCNA <- 0
    
    while ( i < ncol(data) ) { # for a patient
      
      # If its clonal
      if ( data[k,i] == data[k,i+1] ) { 
        if ( data[k,i] == 3 ) { # if both are gains
          clonal.gain <- clonal.gain + 1
          i <- i + 2
        }
        else if ( data[k,i] == 1 ) { # if both are losses
          clonal.loss <- clonal.loss + 1
          i <- i + 2
        }
        else if ( data[k,i] == 2 ) { # if both are diploid
          clonal.noCNA <- clonal.noCNA + 1
          i <- i + 2
        }
      }
      
      # if its subclonal
      else {
        if ( data[k,i] == 3 | data[k,i+1] == 3 ) { # if one is a gain
          subclonal.gain <- subclonal.gain + 1
          i <- i + 2
        }
        else if ( data[k,i] == 1 | data[k,i+1] == 1 ) { # if one is a loss
          subclonal.loss <- subclonal.loss + 1
          i <- i + 2
        }
        else {
          i <- i + 2
        }
      }
    }
    
    GaLo.Clo$clonal.gain[k] <- clonal.gain
    GaLo.Clo$clonal.loss[k] <- clonal.loss
    GaLo.Clo$clonal.noCNA[k] <- clonal.noCNA
    GaLo.Clo$subclonal.gain[k] <- subclonal.gain
    GaLo.Clo$subclonal.loss[k] <- subclonal.loss
  }
  
  GaLo.Clo$clonal.aneu <- GaLo.Clo$clonal.gain + GaLo.Clo$clonal.loss
  GaLo.Clo$subclonal.aneu <- GaLo.Clo$subclonal.gain + GaLo.Clo$subclonal.loss
  GaLo.Clo$gain <- GaLo.Clo$clonal.gain + GaLo.Clo$subclonal.gain
  GaLo.Clo$loss <- GaLo.Clo$clonal.loss + GaLo.Clo$subclonal.loss
  GaLo.Clo$pcSubclonal.gain <- GaLo.Clo$subclonal.gain / GaLo.Clo$gain # subclonalCNA/totalCNA for gain
  GaLo.Clo$pcSubclonal.loss <- GaLo.Clo$subclonal.loss / GaLo.Clo$loss # subclonalCNA/totalCNA for loss
  
  # Count of noCNA, subclonalCNA, and clonalCNA by patient
  patientClo <- as.data.frame(t(sapply(clonal.data, table))) 
  patientClo <- patientClo[, c('0','1','2')]
  colnames(patientClo) <- c('noCNA','subclonal','clonal')
  patientClo$CNA <- patientClo$subclonal + patientClo$clonal
  patientClo$patient <- rownames(patientClo)
  
  newData <- list('start.stop'=start.stop, 'samples'=samples, 'patients'=patients, 'noSamples'=noSamples, 'noPatients'=noPatients,
                  'noBins'=noBins,'chr.mid'=chr.mid,  'chr.end'=chr.ends,
                  'clonal.data'=clonal.data, 'pic.data'=pic.data, 'ave.pic'=ave.pic, 'ith'=ith,
                  'pga'=pga, 'GaLo.clo'=GaLo.Clo, 'patientClo'=patientClo)
  return(newData)
}

# Function to generate matrices for dip vs aneu, dip vs gain, dip vs loss, gain vs loss, and dip vs gain vs loss
genMatrices <- function(rawdata) {
  
  # Create list to store correlation matrices
  matrices.list <- rep(list(data.frame((rawdata[,-c(1:3)]))),5)
  names(matrices.list) = c('diploid.aneu', 'diploid.gain', 'diploid.loss', 'loss.gain', 'loss.dip.gain')
  
  # Convert to numeric matrices
  matrices.list <- lapply(matrices.list, function(x) {
    y <- data.frame(apply(x, 2, as.numeric))
    rownames(y) <- rownames(x)
    y
  })
  
  # In the diploid/aueploid matrix make 0=diploid, 1=aneuploid.
  matrices.list$diploid.aneu[matrices.list$diploid.aneu == 1 | matrices.list$diploid.aneu == 3] <- 1
  matrices.list$diploid.aneu[matrices.list$diploid.aneu == 2] <- 0
  
  # In the diploid/gain matrix make 0=diploid, 1=gain.
  matrices.list$diploid.gain[matrices.list$diploid.gain == 1] <- NA
  matrices.list$diploid.gain[matrices.list$diploid.gain == 2] <- 0
  matrices.list$diploid.gain[matrices.list$diploid.gain == 3] <- 1
  
  # In the diploid/loss matrix make 0=diploid, 1=loss. 
  matrices.list$diploid.loss[matrices.list$diploid.loss == 1] <- 1
  matrices.list$diploid.loss[matrices.list$diploid.loss == 2] <- 0
  matrices.list$diploid.loss[matrices.list$diploid.loss == 3] <- NA
  
  # In the loss/gain matrix make 0=loss, 1=gain. 
  matrices.list$loss.gain[matrices.list$loss.gain == 1] <- 0
  matrices.list$loss.gain[matrices.list$loss.gain == 2] <- NA
  matrices.list$loss.gain[matrices.list$loss.gain == 3] <- 1
  
  # # In the loss/dip/gain matrix make 1=loss, 2=diploid, 3=gain. Requires no changes
  
  return(matrices.list)
}

# Function to calc pic score per patient across multiple samples
PIC <- function(data, sample.no, index) { 
  # based on table of unique patient nos
  # returns PIC per row (bin), calc across cols as given by index
  # PIC formula: 1- ((CN1/n)^2 + (CN2/n)^2 + (CN3/n)^2), where CN1 is no. of counts of copy number1
  PIC <- 1 - ((rowSums(data[,index]==1)/sample.no)^2 + 
                (rowSums(data[,index]==2)/sample.no)^2 + 
                (rowSums(data[,index]==3)/sample.no)^2)
  return(PIC)
}

# Function to collect pvalue from a regression
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Function to assess ploidy and recentre if average CN is 2.8 or more
ploidyRecentre <- function (cn.list) {
  # Requires list with dataframes per patient and cols chr | start | stop | sample1...
  for ( i in 1:length(cn.list) ) {
    wd <- cn.list[[i]][,-c(1:3)]
    if ( mean(colMeans(wd, na.rm = TRUE)) >= 2.8 ) {
      cn.list[[i]] <- cbind(cn.list[[i]][,c(1:3)] ,cn.list[[i]][,-c(1:3)] - 1)
    }
    else {
      cn.list[[i]] <- testcn.list[[i]]
    }
  }
  return(cn.list)
}

# Function to bin CN data to match training dataset (Cross et al)
alignBins <- function(bins, cn.list) {
  # bins needs to be a dataframe holding: bin | chr | start | stop, for the bins to align to
  # cn.list is the output from ploidyRecentre
  
  # Create dataframe for each patient and put into list
  cnBinned.list <- rep(list(bins), length(cn.list))
  
  # Convert to numeric
  cnBinned.list <- lapply(cnBinned.list, function (x) {
    x[] <- apply(x,2,as.numeric)
    x
  })
  
  # Add empty columns to hold output
  for ( i in 1:length(cnBinned.list) ) {
    cnBinned.list[[i]][,c(5:(1+ncol(cn.list[[i]])))] <- NA
  }
  
  # Bin incoming dataset (cn.list)
  for ( p in 1:length(cnBinned.list) ) {
    for ( b in 1:nrow(cnBinned.list[[p]]) ) {
      chr <- cnBinned.list[[p]]$chr[b]
      start <- cnBinned.list[[p]]$start[b] 
      stop <- cnBinned.list[[p]]$stop[b]
      bin <- cnBinned.list[[p]]$bin[b]
      
      wd <- cn.list[[p]][cn.list[[p]]$chr == chr,]
      
      for ( r in 1:nrow(wd) ) {
        
        ## if start is before the earliest row that that chromosome
        if ( start < wd$start[r] ) {
          
          # if stop is in row r of wd or else...
          if ( dplyr::between(stop, wd$start[r], wd$stop[r]) ) { 
            cnBinned.list[[p]][b,c(5:ncol(cnBinned.list[[p]]))] <- wd[r,c(4:ncol(wd))]
            break
          }
          
          # ...or else stop is beyond row r of wd
          else if ( stop > wd$stop[r] ) { 
            fraction <- list()
            f <- 1
            
            fraction[[f]] <- (wd$stop[r]-start)/(stop-start) # what fraction of the bin is the region in
            f <- f + 1
            
            for ( x in 1:(nrow(wd)-r) ) {
              # stop is within x extra wd rows
              if ( dplyr::between(stop, wd$start[r+x], wd$stop[r+x]) ) { 
                fraction[[f]] <- (wd$stop[r+x]-wd$start[r+x])/(stop-start)
                f <- f + 1
                break
              }
              
              # stop is within a later row of wd
              else if ( stop > wd$stop[r+x] ) { 
                fraction[[f]] <- (wd$stop[r+x]-wd$start[r+x])/(stop-start)
                f <- f + 1
              }
            }
            
            # majority of bin sits in row r+f (q) of wd
            q <- r + which.max(fraction)
            cnBinned.list[[p]][b,c(5:ncol(cnBinned.list[[p]]))] <- wd[q,c(4:ncol(wd))]
            break
          }
          
          # stop is also before start of first row of wd
          else if ( stop < wd$stop[r] ) {
            cnBinned.list[[p]][b,c(5:ncol(cnBinned.list[[p]]))] <- NA
            break
          }
        }
        
        ## if start is in row r of wd
        else if ( dplyr::between(start, wd$start[r], wd$stop[r]) ) { 
          
          # if stop is in row r of wd or else...
          if ( dplyr::between(stop, wd$start[r], wd$stop[r]) ) { 
            cnBinned.list[[p]][b,c(5:ncol(cnBinned.list[[p]]))] <- wd[r,c(4:ncol(wd))]
            break
          }
          
          # ...or else stop is beyond row r of wd
          else if ( stop > wd$stop[r] ) { 
            
            # if there are no more rows of wd
            if ( r == nrow(wd) ) {
              cnBinned.list[[p]][b,c(5:ncol(cnBinned.list[[p]]))] <- wd[r,c(4:ncol(wd))]
              break
            }
            
            # if there are more rows of wd
            else if ( r < nrow(wd) ) {
              fraction <- list()
              f <- 1
              
              fraction[[f]] <- (wd$stop[r]-start)/(stop-start) # what fraction of the bin is the region in
              f <- f + 1
              
              for ( x in 1:(nrow(wd)-r) ) {
                # stop is within x extra wd rows
                if ( dplyr::between(stop, wd$start[r+x], wd$stop[r+x]) ) { 
                  fraction[[f]] <- (wd$stop[r+x]-wd$start[r+x])/(stop-start)
                  f <- f + 1
                  break
                }
                
                # stop is within a later row of wd
                else if ( stop > wd$stop[r+x] ) { 
                  fraction[[f]] <- (wd$stop[r+x]-wd$start[r+x])/(stop-start)
                  f <- f + 1
                }
              }
              
              # majority of bin sits in row r+f (q) of wd
              q <- r + which.max(fraction)
              cnBinned.list[[p]][b,c(5:ncol(cnBinned.list[[p]]))] <- wd[q,c(4:ncol(wd))]
              break
            }
          }
        }
        
        ## if start is in a later row of wd
        else {
          next
        }
        
      }
    }
  }
  return(cnBinned.list)
}
  
  




