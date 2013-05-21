library(biomaRt)
library(topGO)
ensembl_url = "sep2011.archive.ensembl.org"
biomart = "ENSEMBL_MART_ENSEMBL"
attribute_go = "go_id"
attribute_gene = "external_gene_id"
filter_go = "with_go"

def.pdf = "TopGO_plots.pdf"
def.table = "TopGO_table.txt"
def.nterms = 10

getIdsInTerm = function(x,subset,data,genome) {
    I = match(subset,unlist(genesInTerm(data,x)),nomatch=0)>0
    paste(sort(genome[subset[I],1]),collapse=', ')
}

single_topGo = function( geneList, genome, gene2GO, allGenes, pdf, table, nterms ) {
    pdf(pdf,paper="a4",height=8,width=11)
    output = table
    tab = list()
    append = FALSE
    subset = allGenes[geneList == 1]
    for (ontol in c("BP","CC","MF")) {
        data = new("topGOdata",description=ontol,ontology=ontol,
          allGenes=geneList,annot=annFUN.gene2GO,gene2GO=gene2GO)
        result = list(classicFisher=runTest(data, statistic="fisher", algo="classic"),
          elimFisher=runTest(data, statistic="fisher", algo="elim"))
        tab = GenTable(data,
          classicFisher=result$classicFisher, elimFisher=result$elimFisher, 
          orderBy="elimFisher", ranksOf="elimFisher", topNodes=nterms)
        showSigOfNodes(data,score(result$elimFisher), first=nterms, useInfo="all")
        tab = cbind(tab,genes=sapply(tab$GO.ID, getIdsInTerm, subset, data, genome))
        write(c("",ontol,""),file=output,append=append)
        append = TRUE
        write.table(tab,file=output,quote=F,sep="\t",row.names=F,append=append)
    }
    dev.off()
}


multi_topGo = function( filename, assembly_id, pdf=def.pdf, table=def.table, nterms=def.nterms ) {
    ensembl = useMart(biomart,host=ensembl_url)
    lsds = listDatasets(ensembl)
    dataset = lsds$dataset[which(lsds$version == assembly_id)]
    ensembl = useDataset(as.character(dataset),mart=ensembl)
    attr1 = attribute_go
    attr2 = attribute_gene
    filt = filter_go
    genome = getBM(attr=c("ensembl_gene_id",attr1,attr2),
      filter=c(filt,"biotype"),
      values=list(TRUE,"protein_coding"),mart=ensembl)
    gene2GO = split(genome[,attr1],genome[,"ensembl_gene_id"])
    I = which(!duplicated(genome[,"ensembl_gene_id"]))
    genome = data.frame(gene_name=genome[I,attr2],row.names=genome[I,"ensembl_gene_id"])
    allGenes = row.names(genome)
    if (any(!nchar(genome[,1])))
      genome[which(!nchar(genome[,1])),1] = allGenes[which(!nchar(genome[,1]))]

    id_set = scan(filename,what=character())
    name0 = gsub("[.].*","",gsub(".*/","",filename))
    I = grep("#",id_set)
    id_sets = list()
    if (length(I) == 0) {
        id_sets[[name0]] = id_set
    } else {
        I = c(I,length(id_set)+1)
        n = 1
        i = 1
        for (j in I[-1]) {
            name = gsub("[[:space:]]+","_",gsub("#","",id_set[i]))
            if (nchar(name) == 0) name = paste(name0,n,sep='.')
            id_sets[[name]] = id_set[(i+1):(j-1)]
            i = j
            n = n+1
        }
    }
    pdflist = list()
    tablelist = list()
    for (listName in names(id_sets)) {
        subset = id_sets[[listName]]
        geneList = factor(as.numeric(match(allGenes,subset,nomatch=0)>0))
        if (length(levels(geneList)) == 2) {
            names(geneList) = allGenes
            pdflist[[listName]] = gsub(".pdf",paste("_",listName,".pdf",sep=""),pdf)
            tablelist[[listName]] = gsub(".txt",paste("_",listName,".txt",sep=""),table)
            single_topGo(geneList,genome,gene2GO,allGenes,
                         pdflist[[listName]],tablelist[[listName]],nterms)
	}
    }
    list(pdflist,tablelist)
}


