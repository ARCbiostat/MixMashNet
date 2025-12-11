#' Bacteremia dataset
#'
#' @description
#' Example dataset used in MixMashNet examples. This dataset contains 7240
#' patients with clinical suspicion of bacteremia who underwent blood culture testing
#' at the Vienna General Hospital
#'
#' @docType data
#' @name bacteremia
#' @usage data(bacteremia)
#' @format A data frame with 7420 rows and 16 variables:
#' \describe{
#'   \item{AGE}{Age (numeric).}
#'   \item{WBC}{White blood cell (numeric).}
#'   \item{NEU}{Neutrophil counts (numeric).}
#'   \item{HGB}{Hemoglobin (numeric).}
#'   \item{PLT}{Platelet count (numeric).}
#'   \item{CRP}{C-reactive protein (numeric).}
#'   \item{APTT}{Activated partial thromboplastin time (numeric).}
#'   \item{FIB}{Fibrinogen (numeric).}
#'   \item{CREA}{Creatinine (numeric).}
#'   \item{BUN}{Blood urea nitrogen (numeric).}
#'   \item{GLU}{Glucose (numeric).}
#'   \item{ALAT}{High-sensitivity C-reactive protein (numeric).}
#'   \item{GBIL}{Total bilirubin (numeric).}
#'   \item{ALB}{Albumin (numeric).}
#'   \item{SEX}{Sex with 0=male and 1=female.}
#'   \item{BloodCulture}{Positive blood culture results with 0=no and 1=yes.}
#' }
#' @examples
#' data(bacteremia)
#' str(bacteremia)
NULL
