#'extrait
#' #(changer le nom)
#'
#' #function additional to prepare file for Gum.Output package
#' #extract value between I1 and I2 in column n#rang
#' #equivalent to limit data using a uniform prior U{I1,I2} on variable rang
#'
#' @param data [matrix],[mcmc] (**required**) a matrix of MCMC output
#' @param rang [numeric] (**with default value**) the column number
#' @param I1 [numeric] (**with default value**) the inferior limit
#' @param I2 [numeric] (**with default value**) the superior limit
#'
#' @return a matrix
#'
#' @export
#'

extrait<-function(data,rang=5,I1=-200,I2=100){
  Names<-coda::varnames(data)
  data_pl<-data[data[,rang]>I1&data[,rang]<I2]
  dim(data_pl)<-redim(data_pl,ncol(data))
  colnames(data_pl)<-Names
  data_pl
}

###########################################
#' redim
#'
#' calculate the dimension values to transform a data vector in matrix with n.col column
#'
#' @param data [vector]  (**required**) data
#' @param file [n.col] (**required**) column number
#'
#' @return a couple of value n.row, n.col
#'
#' @noRd
#'
redim<-function(data,n.col){
  c(length(data)/n.col,n.col)
}
