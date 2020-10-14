#!/usr/bin/env Rscript

### Author: Ludovic V. Mallet, PhD
### 2016.03.22
### licence: GPLv3
### Version: Alpha.2
### Garanteed with misatkes. <- Including this one.
### Important notes:

# Simplified for the purpose of this pipeline by Alexis Simon

library(getopt)

spec <- matrix(c(
        'genome'         , 'i', 1, "character", "Multifasta of the genome assembly (required)",
        'host_learn'     , 'r', 2, "character", "Host training set (optional)",
        'conta_learn'    , 'c', 1, "character", "Contaminant training set (optional) if missing and sliding window parameters are given, the sliding windows composition will be compared to the whole genome composition to contrast potential HGTs (prokaryotes and simple eukaryotes only)",
        'win_step'       , 't', 2, "int", "Step of the sliding windows analysis to locate the contaminant (optional) default: 500bp or 100bp",
        'win_size'       , 'w', 2, "int", "Length of the sliding window to locate the contaminant (optional) default: 5000bp ",
        'outputdir'      , 'W', 1, "character", "path to outputdir directory, should end with '/'",
        'dist'           , 'd', 2, "character", "Divergence metric used to compare profiles: (KL), JSD or Eucl",
        'manual_threshold' , 'm', 0, "logical", "You will be asked to manually set the thresholds",
        'help'           , 'h', 0, "logical",   "This help"
),ncol=5,byrow=T)       

opt = getopt(spec);

if (!is.null(opt[["help"]]) || is.null(opt[["genome"]])) {
    cat(paste(getopt(spec, usage=T),"\n"));
    quit(save="no", status = 1)
}

if (!is.null(opt[["outputdir"]]) ) {
    working_dir = opt[["outputdir"]]
}else{
    working_dir = "./"
}

genome_fasta = opt[["genome"]]

conta_sample_fasta = opt[["conta_learn"]]
if (is.null(opt[["host_learn"]]) {host_sample_fasta <- ""}
if (is.null(opt[["dist"]]) {dist <- "KL"}
if (is.null(opt[["win_step"]])) {win_step <- 500}
if (is.null(opt[["win_size"]])) {win_size <- 5000}

### compute profiles:
data=list()

input_path=paste(working_dir, basename(genome_fasta),".mcp_hostwindows_vs_host_",basename(host_sample_fasta),"_",dist,".dist",sep="")
data[["host"]]=read.delim(file=input_path, header=F)

input_path=paste(working_dir, basename(genome_fasta),".mcp_hostwindows_vs_conta_",basename(conta_sample_fasta),"_",dist,".dist",sep="")
data[["conta"]]=read.delim(file=input_path, header=F)

if (! is.null(opt[["manual_threshold"]])) {
  ### Ask the trusty human to set the thresholds:
  threshold_conta=0
  repeat{
    x11()
    plot(density(data[["conta"]][,4],na.rm=TRUE),lwd=2)
    abline(v=threshold_conta,col="red")
    Sys.sleep(10)
    
    print("Give a different threshold value for the contaminant threshold. Give the same value to confirm it.")
    new_threshold=readLines(con="stdin", 1)

    new_threshold <- as.numeric(new_threshold)
    
    if(new_threshold == threshold_conta){
      break
    }
    threshold_conta=new_threshold
  }

  threshold_host=0
  repeat{
    x11()
    plot(density(data[["host"]][,4],na.rm=TRUE), lwd=2)
    abline(v=threshold_host,col="red")
    Sys.sleep(10)
    print("Give a different threshold value for the host threshold. Give the same value to confirm it.")
    new_threshold=readLines(con="stdin", 1)
    new_threshold <- as.numeric(new_threshold)
    
    if(new_threshold == threshold_host){
      break
    }
    threshold_host=new_threshold
  }

}else{

  ### Humans are not worthy to set the thresholds, stats will guess it:
  conta_threshold_name=paste(working_dir,basename(genome_fasta),"_vs_",basename(conta_sample_fasta),"_conta_threshold.pdf",sep="")
  pdf(file=conta_threshold_name)
  des_conta=density(data[["conta"]][which(!is.nan(data[["conta"]][,4]) ),4] )
  plot(des_conta,lwd=2,main="")
  steep=des_conta[["y"]][seq(which.max(des_conta[["y"]]),0)]
  i=1
  while(steep[i+1]<steep[i]){i=i+1}
  conta_min=(which.max(des_conta[["y"]])-i)
  abline(v=des_conta[["x"]][conta_min],col="blue",lwd=2)
  threshold_conta=des_conta[["x"]][conta_min]
  dev.off()
  print(paste("Threshold conta:", threshold_conta))
  print(paste("Please inspect that the automatic threshold for the contaminant was set-up properly : ",conta_threshold_name))

  host_threshold_name=paste(working_dir,basename(genome_fasta),"_vs_",basename(host_sample_fasta),"_host_threshold.pdf",sep="")
  pdf(paste(host_threshold_name,sep=""))
  des_host=density(data[["host"]][which(!is.nan(data[["host"]][,4]) ),4] )
  plot(des_host,lwd=2,main="")
  steep=des_host[["y"]][seq(which.max(des_host[["y"]]),length(des_host[["y"]]))]
  i=1
  while(steep[i+1]<steep[i]){i=i+1}
  host_min=(which.max(des_host[["y"]])+i)
  abline(v=des_host[["x"]][host_min],col="blue",lwd=2)
  threshold_host=des_host[["x"]][host_min]
  dev.off()
  print(paste("Threshold host:", threshold_host))
  print(paste("Please inspect that the automatic threshold for the host was set-up properly : ",host_threshold_name))
}

### Perform the split over the double threshold

data[["Select_conta"]]=which((data[["conta"]][,4] <= threshold_conta)*(data[["host"]][,4] >= threshold_host)>=1)
data[["windows_conta"]]=data[["conta"]][data[["Select_conta"]],]

### Regroups contiguous windows into islands and write a GFF file of the positions of the targeted species
write(x="##gff-version 2", file = paste(working_dir, basename(genome_fasta),"_contaminant_",basename(conta_sample_fasta),".gff",sep=""),ncolumns=1)
start_index=data[["Select_conta"]][1]
for(i in seq(length(data[["Select_conta"]]))){
  if(ifelse(is.na(data[["Select_conta"]][i+1]), test <- 0 , test <- data[["Select_conta"]][i+1]) == data[["Select_conta"]][i]+1 ){ #2 windows have a consecutive index number, so they are physically touching. data[["conta"]][i,1] describe the contig, and regions are only spanning one contig max. we assert that ssuccessive indexes are in the same contig to group them. If everything of this if fine, we will then group the windows in a region, by etending the end of it.
    end_index=data[["Select_conta"]][i+1]
  }else{ #well, the window index i+1 is not touching the window i, so i is the last of its island, and the island can be written. i+1 is the start of a new island.
    end_index=data[["Select_conta"]][i]
    line=paste(sep="\t", as.character(data[["conta"]][start_index,1]),"SignatureGohtam\tregion",data[["conta"]][start_index,2],data[["conta"]][end_index,c(3)],"\t.\t.\t.")
    write(x=line,append=TRUE, file = paste(working_dir, basename(genome_fasta),"_contaminant_",basename(conta_sample_fasta),".gff",sep=""),ncolumns=1)
    start_index=data[["Select_conta"]][i+1]
  }
}

print("Done")






