library(multiseq)

setwd("results/analysis_gtex")

args = commandArgs(TRUE)
file = as.character(args[1])

region = unlist(strsplit(file, "_"))
region = region[-(1:2)]
region = paste(region, "_")


load(file.path("../../data/gtex", file))


res_multiseq = list()

for(i in 1:7){
  data_temp = rbind(reads[[i]][[2]], reads[[i+1]][[2]])
  g_temp = rep(0:1, each = 100)
  res_multiseq[[i]] = multiseq(data_temp, g_temp, reflect = TRUE, lm.approx = TRUE)
}

tissue1 = NULL
tissue2 = NULL
for(i in 1:7){
  tissue1[i] = c(tissue1, reads[[i]][[1]])
  tissue2[i] = c(tissue2, reads[[i+1]][[1]])
}


for(i in 1:7){
  pdfname = paste("multiseq_effect_100", tissue1[i], tissue2[i], region, "_")
  pdf(paste0(pdfname, ".pdf"), height = 6, width = 9)
  plot(res_multiseq[[i]]$effect.mean, ylab = "effect", main = paste("Effect size estimates between tissues", tissue1[i], "and", tissue2[i], " "))
  dev.off()
}

save(res_multiseq, tissue1, tissue2, file = paste0("multiseq_effect_100_", region, ".Robj"))