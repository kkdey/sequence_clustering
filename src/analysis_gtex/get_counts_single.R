get.counts.single = function(bamfile, region){
  region <- multiseq:::split_region(region)
  locus.length <- region$end - region$start + 1
  
  cmd <- paste0("| awk -v s='1' 'BEGIN{start=0; count=0}", 
                "{st=$4; if ($4>=s){if (start==st) count+=1; ", "else {if (start>0) print start, count; start=st; count=1} }}' ")
  
  command <- paste0("samtools view ", data, " ", 
                    cmd)
  print(command)
  con <- pipe(command, open = "r")
  v <- rep(0, locus.length)
  while (length(oneLine <- readLines(con, n = 1)) > 
           0) {
    oneLine <- as.numeric(unlist(strsplit(oneLine, 
                                          " ")))
    print(oneLine)
    v[oneLine[1] - region$start + 1] <- oneLine[2]
  }
  close(con)
  return(v)
}