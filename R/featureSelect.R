#' Feature Selection
#' 
#' A sample feature selection framework (described in "Examples" section,
#' used as default in iterClust framework). Customized frameworks can be defined
#' following rules specified in "Usage", "Arguments" and "Value" sections.
#' 
#' @param dset (numeric matrix) features in rows and observations in columns
#' @param iteration (positive integer) specifies current iteration
#' @param feature (character array) specifies user defined features,
#' facilitating feature selection
#' 
#' @return a character array, contains features selected
#' 
#' @keywords featureSelect
#' @examples
#' featureSelect <- function(dset, iteration, feature) return(rownames(dset))
#' 
#' @author DING, HONGXU (hd2326@columbia.edu)
#' 
#' @export

featureSelect <- function(dset, iteration, feature) return(rownames(dset))
