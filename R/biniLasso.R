
#' Convert Numeric to Factor
#'
#' @param data a data frame
#' @param cols a vector of numeric column names (characters).
#' @param method either "quantile" (default) for using quantiles of a given column as its cut-points, or "fixed" to used a user provided vector of cut-points. See cuts_list input.
#' @param n_bins an integer value specifying number of bins to convert the numeric columns. Only needed for quantile method.
#' @param cuts_list a list of cut-points corresponding to each numeric column in the same order as entries of cols input.
#'
#' @returns data_cat: a data frame including original numeric columns as well as the converted categorical columns.
#' @returns x: a matrix including dummy variables corresponding to converted categorical variables.
#' @returns x_cuts: a list of cut-points corresponding to each entry of cols and in the same order.
#' @export
#'
#' @examples
#' simData <- data.frame(t = rexp(1000),
#'                       event = rbinom(1000, 1, 0.6),
#'                       x1 = rnorm(1000),
#'                       x2 = rnorm(1000),
#'                       x3 = rnorm(1000))
#'
#' simData_converted_1 <-
#'   num_to_cat(data = simData,
#'              cols = c("x1", "x2", "x3"),
#'              method = "quantile",
#'              n_bins = 30)
#'
#' simData_converted_2 <-
#'   num_to_cat(data = simData,
#'              cols = c("x1", "x2", "x3"),
#'              method = "fixed",
#'              cuts_list = list(c(0, 0.5),
#'                               c(-0.5, 0, 0.5),
#'                               c(-0.5, -0.1, 0.05, 0.1, 0.5)))
#'
num_to_cat <-
  function(data,
           cols,
           method = "quantile",
           n_bins = NULL,
           cuts_list = NULL) {
           # penalty.factor = FALSE) {

    len_flag <- FALSE

    # if (penalty.factor) penFac <- numeric(0)

    data_bins <- data[ , cols]
    for (nf in 1 : length(cols)) {
      if (method == "quantile") {
        x_tmp <- sort(unique(data[ , cols[nf]]))
        x_tmp <- x_tmp[- c(1, length(x_tmp))]
        if (length(x_tmp) < n_bins) {
          len_flag <- TRUE
          n_bins_tmp <- length(x_tmp)
        }
        else n_bins_tmp <- n_bins
        x_cuts_tmp <- unique(stats::quantile(x_tmp,
                                             probs = c(1 : (n_bins_tmp - 1)) / n_bins_tmp))
        x_bounds_tmp <- c(-Inf, x_cuts_tmp, Inf)
        cut_names_tmp <- paste0("_", c(1 : n_bins_tmp))
      }
      if (method == "fixed") {
        x_cuts_tmp <- cuts_list[nf][[1]]
        x_bounds_tmp <- c(-Inf, x_cuts_tmp, Inf)
        cut_names_tmp <- paste0("_", c(1 : (length(cuts_list[nf][[1]]) + 1)))
      }
      data_bins[ , paste0(cols[nf], "_bin")] <-
        cut(data[ , cols[nf]],
            breaks = x_bounds_tmp,
            labels = cut_names_tmp)
      data_bins[ , paste0(cols[nf], "_bin")] <-
        ifelse(data_bins[ , cols[nf]] == min(data_bins[ , cols[nf]]), "min",
               ifelse(data_bins[ , cols[nf]] == max(data_bins[ , cols[nf]]), "max",
                      data_bins[ , paste0(cols[nf], "_bin")]))
      data_bins[ , paste0(cols[nf], "_bin")] <- factor(data_bins[ , paste0(cols[nf], "_bin")])
      data_bins[ , paste0(cols[nf], "_bin")] <- relevel(data_bins[ , paste0(cols[nf], "_bin")],
                                                        ref = "min")
      # if (penalty.factor) {
      #   if (length(cut_names_tmp) > 2) penFac <- c(penFac, c(0, rep(1, length(cut_names_tmp) - 2), 0))
      #   else penFac <- c(penFac, rep(0, length(cut_names_tmp)))
      # }

      if (nf == 1) x_cuts <- list(x_cuts_tmp)
      if (nf == 2) x_cuts <- list(c(x_cuts, list(x_cuts_tmp)))
      if (nf > 2) x_cuts <- list(c(x_cuts[[1]], list(x_cuts_tmp)))
    }

    x <- stats::model.matrix(stats::as.formula(paste0(" ~ ", paste(paste0(cols, "_bin"), collapse = " + "))),
                      data = data_bins)[ , -1]
    x <- x[ , -which(grepl("max$", colnames(x)))]
    colnames(data_bins)[grepl("bin$", colnames(data_bins))] <-
      stringr::str_replace_all(colnames(data_bins)[grepl("bin$", colnames(data_bins))], "bin", "cat")

    if (len_flag) warning("Not enough unique values in listed numeric columns to convert all of them to exactly n_bins dummy variables. The dummy variables were adjusted according to the available unique values in each numeric column.")

    return(list(data_cat = data_bins,
                x = x,
                x_cuts = x_cuts))

    # if (penalty.factor) {
    #   return(list(data_cat = data_bins,
    #               x = x,
    #               x_cuts = x_cuts,
    #               penalty.factor = penFac))
    # }
    # else {
    #   return(list(data_cat = data_bins,
    #               x = x,
    #               x_cuts = x_cuts))
    # }

  }



#' Extract optimal cut-points from the fitted glmnet model.
#'
#' @param glm_fit the fitted glmnet object
#' @param cols a vector of numeric column names (characters)
#' @param x_cuts a list of cut-points corresponding to each entry of cols and in the same order.
#'
#' @returns a list of optimal cut-points corresponding to the columns in cols and in the same order.
#'
cuts_extractor <-
  function(glm_fit,
           cols,
           x_cuts) {

    beta_nonZero <- names(glm_fit$beta[ , 1][glm_fit$beta[, 1] != 0])
    for (nf in 1 : length(cols)) {
      X_cuts_ind_tmp <- as.numeric(unlist(lapply(beta_nonZero[grepl(cols[nf], beta_nonZero)],
                                                 function(x) {
                                                   locs <- stringr::str_locate_all(x, "_bin")[[1]]
                                                   locs <- locs[nrow(locs) , "end"]
                                                   substr(x, start = locs + 1, stop = nchar(x))
                                                 })))
      X_cuts_ind_tmp <- X_cuts_ind_tmp[c(1, which(diff(X_cuts_ind_tmp) > 1) + 1)]
      x_bounds_tmp <- x_cuts[[1]][nf][[1]]
      if (nf == 1) x_cuts_opt <- list(unique(x_bounds_tmp[X_cuts_ind_tmp]))
      if (nf == 2) x_cuts_opt <- list(c(x_cuts_opt, list(unique(x_bounds_tmp[X_cuts_ind_tmp]))))
      if (nf > 2) x_cuts_opt <- list(c(x_cuts_opt[[1]], list(unique(x_bounds_tmp[X_cuts_ind_tmp]))))
    }
  return(x_cuts_opt)
}


#' Optimal cut-points finder based on biniLasso method
#'
#' @param x input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format (inherit from class "sparseMatrix" as in package Matrix). Requirement: nvars >1; in other words, x should have 2 or more columns.
#' @param y response variable. Quantitative for family="gaussian", or family="poisson" (non-negative counts). For family="binomial" should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; for a factor, the last level in alphabetical order is the target class). For family="multinomial", can be a nc>=2 level factor, or a matrix with nc columns of counts or proportions. For either "binomial" or "multinomial", if y is presented as a vector, it will be coerced into a factor. For family="cox", preferably a Surv object from the survival package: see Details section for more information. For family="mgaussian", y is a matrix of quantitative responses.
#' @param method either "biniLasso" or its sparse version "Sparse biniLasso" based on unilasso method. See Safari et al. (2025) for more details.
#' @param family any glmnet family option is available.
#' @param lasso_rule Lasso rule for finding optimal lambda; either "1se" (default) or "min".
#' @param lasso_nfolds The number of lambda values - default is 10.
#' @param penalty.factor Separate penalty factors can be applied to each coefficient. This is a number that multiplies lasso lambda to allow differential shrinkage. Can be 0 for some variables, which implies no shrinkage, and that variable is always included in the model. Default is 1 for all variables . Also, any penalty.factor that is set to inf is converted to an exclude, and then internally reset to 1. Note: the penalty factors are internally rescaled to sum to nvars, and the lasso lambda sequence will reflect this change.
#' @param cols a vector of numeric column names (characters).
#' @param x_cuts a list of cut-points corresponding to each entry of cols and in the same order.
#'
#' @returns a list of optimal cut-points corresponding to the columns in cols and in the same order.
#' @export
#'
#' @examples
#'
#' simData <- data.frame(t = rexp(1000),
#'                       event = rbinom(1000, 1, 0.6),
#'                       x1 = rnorm(1000),
#'                       x2 = rnorm(1000),
#'                       x3 = rnorm(1000))
#'
#' simData_converted_1 <-
#'   num_to_cat(data = simData,
#'              cols = c("x1", "x2", "x3"),
#'              method = "quantile",
#'              n_bins = 30)
#'
#' opt_cuts_finder(x = simData_converted_1$x,
#'                 y = survival::Surv(simData$t, simData$event),
#'                 method = "biniLasso",
#'                 family = "cox",
#'                 lasso_rule = "1se",
#'                 lasso_nfolds = 10,
#'                 penalty.factor = NULL,
#'                 cols = c("x1", "x2", "x3"),
#'                 x_cuts = simData_converted_1$x_cuts)
#'
#'
opt_cuts_finder <-
  function(x,
           y,
           method = "biniLasso",
           family = "cox",
           lasso_rule = "1se",
           lasso_nfolds = 10,
           penalty.factor = NULL,
           cols,
           x_cuts) {

    if (method == "both") method <- c("biniLasso", "Sparse biniLasso")
    if (is.null(penalty.factor)) penalty.factor <- rep(1, ncol(x))

    if ("biniLasso" %in% method) {
      bini_cv <- glmnet::cv.glmnet(x = x, y = y,
                          family = family, nfolds = lasso_nfolds,
                          penalty.factor = penalty.factor)
      bini_fit <- glmnet::glmnet(x = x, y = y,
                        family = family,
                        lambda = ifelse(lasso_rule == "min",
                                        bini_cv$lambda.min,
                                        bini_cv$lambda.1se),
                        penalty.factor = penalty.factor)
      x_cuts_bini_opt <- cuts_extractor(glm_fit = bini_fit,
                                        cols = cols,
                                        x_cuts = x_cuts)
      x_cuts_bini_opt <- x_cuts_bini_opt[[1]]
      names(x_cuts_bini_opt) <- cols
      x_cuts_bini_opt <- dplyr::tibble(method = "biniLasso",
                                opt_cuts = list(x_cuts_bini_opt))
    }
    if ("Sparse biniLasso" %in% method) {
      ubini_cv <- uniLasso::cv.uniLasso(x = x, y = y,
                            family = family, nfolds = lasso_nfolds,
                            penalty.factor = penalty.factor)
      ubini_fit <- uniLasso::uniLasso(x = x, y = y,
                          family = family,
                          lambda = ifelse(lasso_rule == "min",
                                          ubini_cv$lambda.min,
                                          ubini_cv$lambda.1se),
                          penalty.factor = penalty.factor)
      x_cuts_ubini_opt <- cuts_extractor(glm_fit = ubini_fit,
                                         cols = cols,
                                         x_cuts = x_cuts)
      x_cuts_ubini_opt <- x_cuts_ubini_opt[[1]]
      names(x_cuts_ubini_opt) <- cols
      x_cuts_ubini_opt <- dplyr::tibble(method = "Sparse biniLasso",
                                 opt_cuts = list(x_cuts_ubini_opt))
    }

    if (all(method == "biniLasso")) return(x_cuts_bini_opt)
    if (all(method == "Sparse biniLasso")) return(x_cuts_ubini_opt)
    if (all(method == c("biniLasso", "Sparse biniLasso"))) return(x_cuts_bini_opt %>%
                                                                    dplyr::bind_rows(x_cuts_ubini_opt))

  return(x_cuts_opt)
}


