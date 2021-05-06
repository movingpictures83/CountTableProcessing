# Clear R workspace
rm(list=ls(all=TRUE))

require(stringr)


colsplit <- function(string, pattern, names) {
  m <- length(string) 
  n = length(names)
 
  df = vector("list", m)
  for (i in 1:m) {
     vars <- strsplit(toString(string[i]), pattern)
     df[[i]] <- vars[[1]][1:n]
     dim(df[[i]]) <- n
     names(df[[i]]) <- names
  }
  dim(df) <- m  

  df
}

input <- function(inputfile) {
   input2 <<- read.table(paste(inputfile, ".shared", sep=""), header=TRUE, sep="\t")
   tmp <- colnames(input2)[2]
   tmp2 <- input2[2]
   input2 <<- input2[4:ncol(input2)]
   rownames(input2) <<- t(tmp2)
   otu_names <<- read.table(paste(inputfile, ".taxonomy", sep=""), header=TRUE, sep="\t")
   otu_names <<- otu_names[,3]
}

mysub <- function(x) {
   gsub("\\(\\d+\\)", "", x)
}

run <- function() {
col_names <- c('Kingdom','Phylum', 'Class', 'Order', 'Family', 'Genus')
taxonomy <- colsplit(otu_names,";",col_names)

m <- dim(taxonomy)
n <- dim(taxonomy[[1]])
for (i in 1:m) {
   for (j in 1:n) {
      taxonomy[[i]][j] <- mysub(taxonomy[[i]][j]) 
   }
}


# Had to do this manually now TMC
column6 <- vector(mode="numeric", length=dim(taxonomy))

for (i in 1:dim(taxonomy)[1]) {
  taxa <- taxonomy[[i]]
  #TMC new version of Mothur uses Root in spot 1
  n <- length(taxa)
  OTU_kingdom = taxa[1]
  if (is.na(OTU_kingdom)) {
     OTU_kingdom = 'unclassified' 
  }
  OTU_phylum = taxa[2]
  if (is.na(OTU_phylum)) {
     OTU_phylum = 'unclassified'
  }
  OTU_class = taxa[3]
  if (is.na(OTU_class)) {
     OTU_class = 'unclassified'
  }
  OTU_order = taxa[4]
  if (is.na(OTU_order)) {
     OTU_order = 'unclassified'
  }
  OTU_family = taxa[5]
  if (is.na(OTU_family)) {
     OTU_family = 'unclassified'
  }
  OTU_genus = taxa[6]
  if (is.na(OTU_genus)) {
     OTU_genus = 'unclassified'
  }

  if (length(grep('unclassified', OTU_genus)) != 0) {
    if (length(grep('unclassified', OTU_family)) == 0) {
      taxonomy[[i]][6] <- paste('Family', OTU_family, sep=".")
    }
    else {
      if (length(grep('unclassified', OTU_order)) == 0) {
        taxonomy[[i]][6] <- paste('Order', OTU_order, sep=".")
        taxonomy[[i]][5] <- paste('Order', OTU_order, sep=".")
      } else {
        if (length(grep('unclassified', OTU_class)) == 0) {
          taxonomy[[i]][6] <- paste('Class', OTU_class, sep=".")
          taxonomy[[i]][5] <- paste('Class', OTU_class, sep=".")
          taxonomy[[i]][4] <- paste('Class', OTU_class, sep=".")
        }
        else {
          if (length(grep('unclassified', OTU_phylum)) == 0) {
            taxonomy[[i]][6] <- paste('Phylum', OTU_phylum, sep=".")
            taxonomy[[i]][5] <- paste('Phylum', OTU_phylum, sep=".")
            taxonomy[[i]][4] <- paste('Phylum', OTU_phylum, sep=".")
            taxonomy[[i]][3] <- paste('Phylum', OTU_phylum, sep=".")
          }
          else {
            taxonomy[[i]][6] <- paste('Kingdom', OTU_kingdom, sep=".")
            taxonomy[[i]][5] <- paste('Kingdom', OTU_kingdom, sep=".")
            taxonomy[[i]][4] <- paste('Kingdom', OTU_kingdom, sep=".")
            taxonomy[[i]][3] <- paste('Kingdom', OTU_kingdom, sep=".")
            taxonomy[[i]][2] <- paste('Kingdom', OTU_kingdom, sep=".")
          }
        }
      }
    }
  }
  column6[i] <- taxonomy[[i]][6]
}
tmp_colnames <- column6 #t(column6)
#Fix colnames  (This is pretty slow, maybe it can be done faster if we do it differently)
t_colnames <- c();
for(i in 1:length(column6)){
	current <- which(t_colnames[,1]==taxonomy[[i]][6])
	if(length(current)>0){
		tmp_colnames[i] <- sprintf("%s.%04d", tmp_colnames[i], as.integer(t_colnames[current,2]))
		t_colnames[current,2] <- as.integer(t_colnames[current,2]) + 1
	}else{
		tmp_colnames[i] <- paste(tmp_colnames[i], "0001", sep=".")
		t_colnames <- rbind(t_colnames, c(taxonomy[[i]][6], 2));
	}
}
colnames(input2) <- tmp_colnames;

sample_list <- strsplit(as.vector(rownames(input2)),'\\.')
samples <- cbind(matrix(unlist(sample_list), ncol = 1, byrow = TRUE), "")
#Now calculate the duplicates

preNAP = "";
for (i in 1:nrow(samples)) {
	if(samples[i,1] == preNAP){
		samples[i,2] = "X"
	}
	preNAP = samples[i,1]
}
#Remove the duplicates from input and add the row names
input2 <<- input2[which(samples[,2]!="X"),]
rownames(input2) <- samples[which(samples[,2]!="X"),1]
}

output <- function(outputfile) {
fileOut <- paste(outputfile, sep="")
write.table(input2, file=fileOut, sep=",", append=FALSE, col.names=NA, na="")

}
