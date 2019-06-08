#' Encounters of zebra mussels in Lake Burgen, MN
#'
#' Number of zebra mussel encounters made by two observers in a double-observer removal design.
#'
#' @format A data frame with 53940 rows and 10 variables:
#' @format
#' \describe{
#'     \item{Observer name}{The name of the observer that made the detection}
#'     \item{Transect number}{The transect the detection was made on}
#'     \item{size}{The number of individual mussels in the detection event.}
#' }
#' @source \href{https://github.com/troutinthemilk/Zeebs-Season2/blob/master/Tutorial/EncountersDoubleTutorial.xlsx}{github}, based on data from \href{https://conservancy.umn.edu/handle/11299/201572}{Data Repository for University of Minnesota}
#'
#' @examples
#' data(DoubleEncounter)
"DoubleEncounter"



#' Encounters of zebra mussels in Lake Burgen, MN
#'
#' Number of zebra mussel encounters made by two observers in a double-observer distance design.
#'
#' @format
#' \describe{
#'     \item{Observer name}{The name of the observer that made the detection}
#'     \item{Transect number}{The transect the detection was made on}
#'     \item{Distance}{The distance of the detection event from the transect line.}
#'     \item{size}{The number of individual mussels in the detection event.}
#' }
#'
#' @source \href{https://github.com/troutinthemilk/Zeebs-Season2/blob/master/Tutorial/EncountersDistanceTutorial.xlsx}{github}, based on data from \href{https://conservancy.umn.edu/handle/11299/201572}{Data Repository for University of Minnesota}
#'
#' @examples
#' data(DistanceEncounter)
"DistanceEncounter"


#' Transect information for zebra mussels in Lake Burgen, MN
#'
#' Transects surveyed by two observers in a double-observer removal design.
#'
#' @format
#' \describe{
#'     \item{primary}{The name of the primary observer on this transect}
#'     \item{secondary}{The name of the secondary observer on this transect}
#'     \item{length}{The length of the current transect}
#'     \item{Transect number}{The current transect}
#' }
#'
#' @source \href{https://github.com/troutinthemilk/Zeebs-Season2/blob/master/Tutorial/EncountersDoubleTransect.xlsx}{github}, based on data from \href{https://conservancy.umn.edu/handle/11299/201572}{Data Repository for University of Minnesota}
#'
#' @examples
#' data(DoubleTransect)
"DoubleTransect"



#' Transect information for zebra mussels in Lake Burgen, MN
#'
#' Transects surveyed by two observers in a double-observer distance survey with removal.
#'
#' @format
#' \describe{
#'     \item{primary}{The name of the primary observer on this transect}
#'     \item{secondary}{The name of the secondary observer on this transect}
#'     \item{length}{The length of the current transect}
#'     \item{Transect number}{The current transect}
#' }
#'
#' @source \href{https://github.com/troutinthemilk/Zeebs-Season2/blob/master/Tutorial/EncountersDoubleTransect.xlsx}{github}, based on data from \href{https://conservancy.umn.edu/handle/11299/201572}{Data Repository for University of Minnesota}
#'
#' @examples
#' data(DistanceTransect)
"DistanceTransect"
