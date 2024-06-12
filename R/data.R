
#' A dataset from the HF-ACTION trial
#'
#' @description The Heart Failure: A Controlled Trial Investigating Outcomes of Exercise Training
#' (HF-ACTION) study was conducted between 2003--2007 to investigate whether
#' adding exercise training to the usual care of heart failure patients improves their cardiovascular outcomes (O'Conner et al., 2009).
#' This dataset is for a high-risk subgroup consisting of 741 nonischemic patients with
#' baseline cardiopulmonary test duration less than or equal to 12 minutes
#' and analyze recurrent hospitalizations as well as overall survival.
#'
#'
#'
#' @format A data frame with 2,132 rows and 4 variables:
#' \describe{
#'   \item{id}{Unique patient ID.}
#'   \item{time}{Event time (years).}
#'   \item{status}{Event type; 1 = hospitalization, 2 = death, 0 = censoring.}
#'   \item{trt}{1 = exercise training, 0 = usual care.}
#'   }
#' @references O'CONNOR, C. M., WHELLAN, D. J., LEE, K. L., KETEYIAN, S. J.,
#'  COOPER, L. S., ELLIS, S. J., LEIFER, E. S.,
#' KRAUS, W. E., KITZMAN, D. W., BLUMENTHAL, J. A. et al. (2009).
#' Efficacy and safety of exercise training in
#'  patients with chronic heart failure: Hf-action randomized controlled trial.
#'  J. Am. Med. Assoc. 301, 1439--1450.
"liver"
