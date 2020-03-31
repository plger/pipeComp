#' checkPipelinePackages
#' 
#' Checks whether the packages required by a pipeline and its alternative 
#' methods are available.
#'
#' @param alternatives A named list of alternative parameter values
#' @param pipDef An object of class `PipelineDefinition`.
#'
#' @return Logical.
#' @export
#' 
#' @importFrom utils installed.packages
#' @examples
#' checkPipelinePackages(list(argument1="mean"), scrna_seurat_pipeline())
checkPipelinePackages <- function(alternatives, pipDef=NULL){
  fns <- unlist(alternatives[sapply(alternatives, class)=="character"])
  fns <- lapply(fns, FUN=function(x){
    if(exists(x) && is.function(get(x))) return(get(x))
    ""
  })
  fns <- paste(unlist(fns),collapse="\n")
  if(!is.null(pipDef)) fns <- paste(fns, paste(pd@functions, collapse="\n"), 
                                    paste(pd@evaluation, collapse="\n"))
  pkg <- gregexpr("library\\(([[:alnum:]])+\\)", fns)
  pkg <- unique(regmatches(fns, pkg)[[1]])
  pkg <- gsub("\\)","",gsub("^library\\(","",pkg))
  pkg <- gsub('"',"",pkg)
  misspkg <- setdiff(pkg, row.names(installed.packages()))
  if(length(misspkg)>0) message("The following packages appear to be missing:",
       paste(misspkg, collapse=", "))
  return(length(misspkg)==0)
}



#' parsePipNames
#' 
#' Parses the names of analyses performed through `runPipeline` to extract a 
#' data.frame of parameter values (with decent classes).
#'
#' @param x The names to parse, or a data.frame with the names to parse as 
#' row.names. All names are expected to contain the same parameters.
#' @param setRowNames Logical; whether to set original names as row.names of
#' the output data.frame (default FALSE)
#' @param addcolumns An optional data.frame of `length(x)` rows to cbind to the
#' output.
#'
#' @return A data.frame
#' 
#' @importFrom utils type.convert
#' @export
#'
#' @examples
#' my_names <- c("param1=A;param2=5","param1=B;param2=0")
#' parsePipNames(my_names)
parsePipNames <- function(x, setRowNames=FALSE, addcolumns=NULL){
  if(is.data.frame(x) || is.matrix(x)){
    if(!is.null(addcolumns)){
      addcolumns <- cbind(x,addcolumns)
    }else{
      addcolumns <- x
    }
    x <- row.names(x)
  }
  x2 <- lapply(strsplit(x,";"),FUN=function(x) x)
  if(length(unique(sapply(x2,length)))>1) 
    stop("The different names do not have the same number of components.")
  n <- sapply(strsplit(x2[[1]],"="),FUN=function(x) x[1])
  y <- sapply(strsplit(unlist(x2),"="),FUN=function(x) x[2])
  y <- as.data.frame(matrix(y, ncol=length(n), byrow=TRUE))
  colnames(y) <- n
  for(i in seq_len(ncol(y))) y[[i]] <- type.convert(y[[i]])
  if(setRowNames) row.names(y) <- x
  if(!is.null(addcolumns)){
    row.names(addcolumns) <- NULL
    y <- cbind(y,addcolumns)
  }
  y
}

# run function `x` on object `o`; if there is no function `x`, run `alt` passing
# `x` as second argument
.runf <- function(x, o, alt=NULL, ...){
  if(exists(x) && is.function(get(x))){
    return(get(x)(o, ...))
  }else{
    if(is.null(alt)) stop("Function '",x,"' not found in environment!")
    return(alt(o, x, ...))
  }
}


#' buildCombMatrix
#' 
#' Builds a matrix of parameter combinations from a list of alternative values.
#'
#' @param alt A named list of alternative parameter values
#' @param returnIndexMatrix Logical; whether to return a matrix of indices, 
#' rather than a data.frame of factors.
#'
#' @return a matrix or data.frame
#' @export
#'
#' @importFrom data.table data.table setorder
#' @examples
#' buildCombMatrix(list(param1=LETTERS[1:3], param2=1:2))
buildCombMatrix <- function(alt, returnIndexMatrix=FALSE){
  eg <- data.table(expand.grid(lapply(alt, FUN=seq_along)))
  eg <- setorder(eg)
  if(returnIndexMatrix) return(as.matrix(eg))
  eg <- as.data.frame(eg)
  for(f in names(alt)){
    eg[,f] <- factor(alt[[f]][eg[,f]], levels=alt[[f]])
  }
  eg
}

#' @importFrom data.table data.table setorder
.checkCombMatrix <- function(eg, alt){
  if(is.null(dim(eg))) 
    stop("`eg` should be a matrix or data.frame of indices or factors")
  if(!all(names(alt) %in% colnames(eg))) 
    stop("The columns of `eg` do not correspond to the arguments.")
  eg <- eg[,names(alt)]
  if(!is.matrix(eg) || !is.numeric(eg)){
    for(f in colnames(eg)){
      if(is.character(eg[,f])) eg[,f] <- factor(eg[,f])
      if(is.factor(eg[,f])){
        if(!all(levels(eg[,f])==alt[[f]])) 
          stop("If `eg` contains factors, the levels should be identical to 
                 the values of the corresponding element of `alternatives`")
        eg[,f] <- as.numeric(eg[,f])
      }
    }
  }
  eg <- as.matrix(as.data.frame(setorder(data.table(eg))))
  if(any(is.na(eg))) stop("Final `eg` contains missing values!")
  eg
}


#' getQualitativePalette
#'
#' Returns a qualitative color palette of the given size. If less than 23 colors
#'  are required, the colors are based on Paul Tol's palettes. If more, the 
#'  `randomcoloR` package is used.
#'
#' @param nbcolors A positive integer indicating the number of colors
#'
#' @return A vector of colors
#'
#' @export
#' @importFrom randomcoloR distinctColorPalette
#' @examples
#' getQualitativePalette(5)
getQualitativePalette <- function(nbcolors){
  nbcolors <- round(nbcolors)
  switch(as.character(nbcolors),
         "1"=c("#4477AA"),
         "2"=c("#4477AA", "#CC6677"),
         "3"=c("#4477AA", "#DDCC77", "#CC6677"),
         "4"=c("#4477AA", "#117733", "#DDCC77", "#CC6677"),
         "5"=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677"),
         "6"=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499"),
         "7"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677",
               "#AA4499"),
         "8"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
               "#CC6677","#AA4499"),
         "9"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
               "#CC6677", "#882255", "#AA4499"),
         "10"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", 
                "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"),
         "11"=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", 
                "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", 
                "#AA4499"),
         "12"=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", 
                "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", 
                "#882255", "#AA4499"),
         "13"=c("#882E72", "#B178A6", "#1965B0", "#5289C7", "#7BAFDE", 
                "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", 
                "#F1932D", "#E8601C", "#DC050C"),
         "14"=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", 
                "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", 
                "#F6C141", "#F1932D", "#E8601C", "#DC050C"),
         "15"=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", 
                "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", 
                "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA"),
         "16"=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", 
                "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", 
                "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA", "black"),
         distinctColorPalette(nbcolors)
  )
}


.getTrueLabelsFromNames <- function(x){
  if(is.null(names(x))) return(NULL)
  tl <- sapply(strsplit(names(x),".",fixed=TRUE), FUN=function(x) x[[1]])
  names(tl) <- names(x)
  tl
}

#' farthestPoint
#'
#' Identifies the point farthest from a line passing through by the first and 
#' last points. Used for automatization of the elbow method.
#'
#' @param y Monotonically inscreasing or decreasing values
#' @param x Optional x coordinates corresponding to `y` (defaults to seq)
#'
#' @return The value of `x` farthest from the diagonal.
#' @export
#'
#' @examples
#' y <- 2^(10:1)
#' plot(y)
#' x <- farthestPoint(y)
#' points(x,y[x],pch=16)
farthestPoint <- function(y, x=NULL){
  if(is.null(x)) x <- seq_len(length(y))
  d <- apply( cbind(x,y), 1, 
              a=c(1,y[1]), b=c(length(y),rev(y)[1]), 
              FUN=function(y, a, b){
                v1 <- a-b
                v2 <- y-a
                abs(det(cbind(v1,v2)))/sqrt(sum(v1*v1))
              })
  order(d,decreasing=TRUE)[1]
}
