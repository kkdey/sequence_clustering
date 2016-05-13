library(multiseq)


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
  temp = reads[[i]][[1]]
  temp = gsub(" .*$", "", temp)
  tissue1 = c(tissue1, temp)
}
for(i in 2:8){
  temp = reads[[i]][[1]]
  temp = gsub(" .*$", "", temp)
  tissue2 = c(tissue2, temp)
}


for(i in 1:7){
  pdfname = paste("multiseq_effect_100", tissue1[i], tissue2[i], region, "_")
  pdf(paste0(pdfname, ".pdf"), height = 6, width = 9)
  plot(res_multiseq[[i]]$effect.mean, ylab = "effect", main = paste("Effect size estimates between tissues", tissue1[i], "and", tissue2[i], " "), type = 'l')
  abline(h = 0, lty = 2, col = 2)
  dev.off()
}

save(res_multiseq, tissue1, tissue2, file = paste0("multiseq_effect_100_", region, ".Robj"))