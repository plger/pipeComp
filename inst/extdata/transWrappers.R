sanityNormWrapper <- function(x, nthreads=1, SanityBin="Sanity", keep=TRUE, tmpdir="."){
        f <- tempfile(tmpdir=tmpdir)
        system(paste("mkdir -p",f))
        cmd <- paste(SanityBin,ifelse(keep,"-e",""),"-d",f,"-n",nthreads,"-f",file.path(f,"input.tab"))
        if(is(x,"SingleCellExperiment")){
                e <- data.frame(GeneID=row.names(x), as.matrix(counts(x)))
        }else{
                e <- data.frame(GeneID=row.names(x), as.matrix(GetAssayData(x, assay = "RNA", slot = "counts")))
        }
        write.table(e, file.path(f,"input.tab"), row.names=FALSE, col.names=TRUE, sep="\t")
        system(cmd)
        e <- read.delim(file.path(f,"log_transcription_quotients.txt"),header=TRUE,row.names=1)
        e <- as.matrix(e)
        if(all(dim(e)==dim(x))){
                dimnames(e) <- dimnames(x)
        }else{
                x <- x[row.names(e),colnames(e)]
        }
        if(is(x,"SingleCellExperiment")){
                logcounts(x) <- as.matrix(e)
        }else{
            x <- SetAssayData(x, slot="data", new.data=as.matrix(e))
            x <- SetAssayData(x, slot="scale.data", as.matrix(GetAssayData(x)))
        }
        x
}


