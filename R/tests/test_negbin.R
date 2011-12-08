
create_fake_dataset <- function(n){
    samples = c("g1.1","g1.2","g1.3","g2.1","g2.2","g2.3","g3.1","g3.2","g3.3")
    samples = c(paste("counts",samples,sep="."),paste("rpkm",samples,sep="."))
    means1 = sample(100:200,n,replace=T); thetas1 = sample(2:5,n,replace=T)/10
    means2 = sample(100:200,n,replace=T); thetas2 = sample(2:5,n,replace=T)/10
    means3 = sample(500:700,n,replace=T); thetas3 = sample(8:12,n,replace=T)/10
    features = paste(rep("feat",n), seq(n), sep="")
    data = data.frame(row.names=samples)
    for (i in 1:n){
        line = c(rnegbin(3,means1[i],thetas1[i]),rnegbin(3,means2[i],thetas2[i]),rnegbin(3,means3[i],thetas3[i]))
        data[features[i]] = c(line, line/3)
    }
    data = as.data.frame(signif(t(data),2))
    write.table(data,"tests/testing_files/data.txt", sep=",", row.names=T, col.names=T, quote=F)
}

test0 <- function(){
    data_file = "tests/testing_files/data.txt"
    design_file = "tests/testing_files/design.txt"
    contrast_file = "tests/testing_files/contrast.txt"
    output_file = "tests/testing_files/output_glm.txt"
    if (file.exists(output_file)) {unlink(output_file)} #deletes the file if it already exists
    comparisons = main(data_file,",",design_file=design_file, contrast_file=contrast_file, output_file=output_file)
}

test1 <- function(){
    data_file = "tests/testing_files/mult_genes.csv"
    design_file = "tests/testing_files/design_mef.txt"
    contrast_file = "tests/testing_files/contrast_mef.txt"
    comparisons = main(data_file,"\t",design_file=design_file, contrast_file=contrast_file)
}
