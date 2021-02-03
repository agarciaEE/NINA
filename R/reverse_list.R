#' @title REVERSE LIST  FUNCTION
#'
#' @param ll list of lists
#'
#' @description Reverts the subset of lists
#'
#' @return List of lists
#'
#'@details Returns an error if \code{filename} does not exist.
#'
#' @examples
#' \dontrun{
#' accident_2015 <- fars_read("Project/data/accident_2015.csv.bz2")
#' }
#'
#' @importFrom stats setNames
#'
#' @export
reverse_list <-  function(ll) {
  nms <- unique(unlist(lapply(ll, function(X) names(X))))
  ll <- lapply(ll, function(X) setNames(X[nms], nms))
  ll <- apply(do.call(rbind, ll), 2, as.list)
  lapply(ll, function(X) X[!sapply(X, is.null)])
}
