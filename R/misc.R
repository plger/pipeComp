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
  for(i in 1:ncol(y)) y[[i]] <- type.convert(y[[i]])
  if(setRowNames) row.names(y) <- x
  if(!is.null(addcolumns)){
    row.names(addcolumns) <- NULL
    y <- cbind(y,addcolumns)
  }
  y
}

# run function `x` on object `o`; if there is no function `x`, run `alt` passing `x` as second argument
.runf <- function(x, o, alt=NULL, ...){
  if(exists(x) && is.function(get(x))){
    return(get(x)(o, ...))
  }else{
    if(is.null(alt)) stop("Function '",x,"' not found in environment!")
    return(alt(o, x, ...))
  }
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
getQualitativePalette <- function(nbcolors){
  nbcolors <- round(nbcolors)
  if(nbcolors>22){
    distinctColorPalette(nbcolors)
  }
  switch(as.character(nbcolors),
         "1"=c("#4477AA"),
         "2"=c("#4477AA", "#CC6677"),
         "3"=c("#4477AA", "#DDCC77", "#CC6677"),
         "4"=c("#4477AA", "#117733", "#DDCC77", "#CC6677"),
         "5"=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677"),
         "6"=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499"),
         "7"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499"),
         "8"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499"),
         "9"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"),
         "10"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"),
         "11"=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"),
         "12"=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499"),
         "13"=c("#882E72", "#B178A6", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C"),
         "14"=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C"),
         "15"=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA"),
         "16"=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA", "black"),
         "17"=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
         "18"=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
         "19"=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788", "black"),
         "20"= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
         "21"= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
         "22"= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788", "black"),
         stop("Unknown nbcolors")
  )
}
