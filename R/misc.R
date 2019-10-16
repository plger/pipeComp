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
  if(is.data.frame(x)){
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
