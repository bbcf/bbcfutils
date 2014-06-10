#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
pdfname = args[1]

data <- read.table(file('stdin'), header=FALSE, quote="")
colnames(data) = c("frequency","count")
maxcount = max(data$count)
pdata = rep(0, maxcount)
for (i in 1:nrow(data)){
    cnt = data$count[i]
    freq = data$freq[i]
    if (freq != 0){ freq = log(freq,10)}
    pdata[cnt] = freq
}
pdata = data.frame(count=1:maxcount,frequency=pdata)

pdf(pdfname)
library(ggplot2)
ggplot(pdata, aes(x=count,y=frequency)) +
geom_area(aes(y=frequency), fill="blue") +
scale_x_log10(name="Read count") +
scale_y_continuous(name="log(Frequency)")
dev.off()

