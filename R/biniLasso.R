
#' Create cumulative binarized features
#'
#' @param x input numeric covariate for cumulative binarization
#' @param breaks breaking values to be used for binarization
#' @param labels labels to be used for the created cumulative binarized columns
#'
#' @returns a matrix containing cumulative binarized features corresponds to x
#'
c_binarization <-
  function (x,
            breaks,
            labels) {
    x_length <- length(x)
    n_breaks <- length(breaks)
    matrix_binarization <- base::matrix(0, nrow = x_length,
                                        ncol = n_breaks,
                                        dimnames = list(character(),
                                                        labels))
    for (i in 1:n_breaks) {
      matrix_binarization[, i][which(x >= breaks[i])] <- 1
    }
    return(matrix_binarization)
  }


#' Convert Numeric to Factor
#'
#' @param data a data frame
#' @param cols a vector of numeric column names (characters).
#' @param method either "quantile" (default) for using quantiles of a given column as its cut-points, or "fixed" to used a user provided vector of cut-points. See cuts_list input.
#' @param n_bins an integer value specifying number of bins to convert the numeric columns. Only needed for quantile method.
#' @param cuts_list a list of cut-points corresponding to each numeric column in the same order as entries of cols input.
#'
#' @returns data_cat a data frame including original numeric columns as well as the converted categorical columns.
#' @returns x a matrix including dummy variables corresponding to converted categorical variables.
#' @returns x_cuts a list of cut-points corresponding to each entry of cols and in the same order.
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
#'   cumBinarizer(data = simData,
#'              cols = c("x1", "x2", "x3"),
#'              method = "quantile",
#'              n_bins = 30)
#'
#' simData_converted_2 <-
#'   cumBinarizer(data = simData,
#'              cols = c("x1", "x2", "x3"),
#'              method = "fixed",
#'              cuts_list = list(c(0, 0.5),
#'                               c(-0.5, 0, 0.5),
#'                               c(-0.5, -0.1, 0.05, 0.1, 0.5)))
#'
cumBinarizer <-
  function(data,
           cols,
           method = "quantile",
           n_bins = NULL,
           cuts_list = NULL) {

    len_flag <- FALSE
    x <- matrix(1, ncol = 1, nrow = nrow(data))

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
                                             probs = c(1 : (n_bins_tmp)) / (n_bins_tmp + 1)))
        cut_names_tmp <- paste0("_bin", c(1 : n_bins_tmp))
      }
      if (method == "fixed") {
        x_cuts_tmp <- cuts_list[nf][[1]]
        cut_names_tmp <- paste0("_bin", c(1 : (length(cuts_list[nf][[1]]))))
      }

      x <- cbind(x,
                 c_binarization(x = data[ , cols[nf]],
                                breaks = x_cuts_tmp,
                                labels = paste0(cols[nf], cut_names_tmp)))

      if (nf == 1) x_cuts <- list(x_cuts_tmp)
      if (nf == 2) x_cuts <- list(c(x_cuts, list(x_cuts_tmp)))
      if (nf > 2) x_cuts <- list(c(x_cuts[[1]], list(x_cuts_tmp)))
    }

    if (len_flag) warning("Not enough unique values in listed numeric columns to convert all of them to exactly n_bins dummy variables. The dummy variables were adjusted according to the available unique values in each numeric column.")

    return(list(x = x[ , -1],
                x_cuts = x_cuts))
  }


#' Extract optimal cut-points from the fitted glmnet model.
#'
#' @param glm_fit the fitted glmnet object
#' @param cols a vector of numeric column names (characters)
#' @param x_cuts a list of cut-points corresponding to each entry of cols and in the same order.
#' @param lambda_opt the value to be used as Lasso optimal lambda value for coefficients extracted. This can be left unspecified (NULL) if lambda value already has been specified in the glm_fit (regular LAsso fit), or alternatively, the glm_fit can contain model fit for a range of lambda values (mainly for uniLasso fit), but for coefficient extraction, the fit corresponding to this lambda value will be used.
#'
#' @returns a list of optimal cut-points corresponding to the columns in cols and in the same order.
#'
cuts_extractor <-
  function(glm_fit,
           cols,
           x_cuts,
           lambda_opt = NULL) {

    if (is.null(lambda_opt)) beta_nonZero <- names(glm_fit$beta[ , 1][glm_fit$beta[, 1] != 0])
    else {
      lambda_inx <- which(glm_fit$lambda == lambda_opt)
      beta_nonZero <- names(glm_fit$beta[ , lambda_inx][glm_fit$beta[, lambda_inx] != 0])
    }
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
#' @param method either "biniLasso" or its sparse version "miniLasso" based on unilasso method. See Safari et al. (2025) for more details.
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
#'   cumBinarizer(data = simData,
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

    if (method == "both") method <- c("biniLasso", "miniLasso")
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
    if ("miniLasso" %in% method) {
      ubini_cv <- uniLasso::cv.uniLasso(x = x, y = y,
                            family = family, nfolds = lasso_nfolds,
                            penalty.factor = penalty.factor)
      ubini_fit <- uniLasso::uniLasso(x = x, y = y,
                          family = family,
                          penalty.factor = penalty.factor)
      x_cuts_ubini_opt <- cuts_extractor(glm_fit = ubini_fit,
                                         cols = cols,
                                         x_cuts = x_cuts,
                                         lambda_opt = ifelse(lasso_rule == "min",
                                                             ubini_cv$lambda.min,
                                                             ubini_cv$lambda.1se))
      x_cuts_ubini_opt <- x_cuts_ubini_opt[[1]]
      names(x_cuts_ubini_opt) <- cols
      x_cuts_ubini_opt <- dplyr::tibble(method = "miniLasso",
                                 opt_cuts = list(x_cuts_ubini_opt))
    }

    if (all(method == "biniLasso")) return(x_cuts_bini_opt)
    if (all(method == "miniLasso")) return(x_cuts_ubini_opt)
    if (all(method == c("biniLasso", "miniLasso"))) return(x_cuts_bini_opt %>%
                                                                    dplyr::bind_rows(x_cuts_ubini_opt))

  return(x_cuts_opt)
}


#' Fit a GLM using detected optimal cut-points
#'
#' @param data a data frame
#' @param optCuts a matrix containing a vector of optimal cut-points (as a list) corresponding to each covariate (may be null for some X), and a column of covariate names.
#' @param y response variable. A matrix with two columns for a Cox family model, and a vector otherwise
#' @param family any glmnet family option is available.
#' @param col_cuts column name of optimal cut-points in optCuts
#' @param col_x column name of covariate names in optCuts
#'
#' @returns data frame including created categorical covariates
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
#' simData_converted <-
#'   cumBinarizer(data = simData,
#'                cols = c("x1", "x2", "x3"),
#'                method = "quantile",
#'                n_bins = 30)
#'
#' simData_sim_optCuts <-
#'  opt_cuts_finder(x = simData_converted$x,
#'                  y = survival::Surv(simData$t, simData$event),
#'                  method = "biniLasso",
#'                  family = "cox",
#'                  lasso_rule = "1se",
#'                  lasso_nfolds = 10,
#'                  penalty.factor = NULL,
#'                  cols = c("x1", "x2", "x3"),
#'                  x_cuts = simData_converted$x_cuts)
#'
#' biniFit(data = simData_converted,
#'         optCuts = .,
#'         y = Surv(simData_converted$tte, simData_converted$vital_status),
#'         family = "cox",
#'         col_cuts = "opt_cuts",
#'         col_x = "opt_cuts_id")
#'
#'
biniFit <- function(data,
                    optCuts,
                    y,
                    family,
                    col_cuts = "opt_cuts",
                    col_x = "opt_cuts_id") {
  optCuts %<>%
    rowwise %>%
    mutate(na_flag = all(is.na(as.numeric(unlist(!!sym(col_cuts)))))) %>%
    ungroup %>%
    filter(! na_flag) %>%
    select(! na_flag)
  if (nrow(optCuts) > 0) {
    cols <- optCuts[[col_x]]
    data_converted <-
      cumBinarizer(data = data,
                   cols = cols,
                   method = "fixed",
                   cuts_list = optCuts[[col_cuts]])
  } else data_converted <- NULL
  if (family == "cox") {
    if (nrow(optCuts) > 0) {
      dataFit <-
        as.data.frame(data_converted$x) %>%
        mutate(time = y[ , 1],
                     event = y[ , 2])
    } else {
      dataFit <- data.frame(time = y[ , 1],
                            event = y[ , 2])
    }
    bini_fit <- survival::coxph(formula = survival::Surv(time, event) ~ .,
                                data = dataFit, x = TRUE)
  }
  else {
    if (nrow(optCuts) > 0) {
      dataFit <-
        as.data.frame(data_converted$x) %>%
        mutate(y = y)
    } else {
      dataFit <- data.frame(y = y)
    }
    bini_fit <- stats::glm(formula = y ~ .,
                           data = dataFit,
                           family = family)
  }

  return(list(data = data_converted$data_cat,
              fit = bini_fit,
              dataFit = dataFit))
}



#' Find a fixed number of cut-points for each predictor
#'
#' @param x input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format (inherit from class "sparseMatrix" as in package Matrix). Requirement: nvars >1; in other words, x should have 2 or more columns.
#' @param y response variable. A matrix with two columns for a Cox family model, and a vector otherwise
#' @param lambda_grid a vector of sequence of lambda values. Can be left NULL.
#' @param nCuts_per_group integer, number of cut-points per predictor
#' @param miniLasso logical, whether to use miniLasso approach
#'
#' @returns a vector of cut-points (character)
#'
cuts_finder <- function(x = X_mat_tmp,
                        y = Surv(simData_mod_sc1$Y, simData_mod_sc1$delta),
                        lambda_grid = NULL,
                        nCuts_per_group = 2,
                        miniLasso = FALSE) {
  if (is.null(lambda_grid)) lambda_grid <- exp(seq(-10, 1, length.out = 100))
  if (miniLasso) {
    fit_init <- uniLasso(x, y, family = "cox", lambda = lambda_grid)
  }
  else {
    fit_init <- glmnet(x, y, family = "cox", lambda = lambda_grid)
  }
  coef_init <- as.matrix(coef(fit_init))
  nz_counts <- apply(coef_init, 2, function(col) sum(col != 0))

  if (! any(nz_counts >= nCuts_per_group)) warning("The default lambda range returned less than nCuts_per_group non-zero coefficients. Consider passing a wider lambda range.")

  lambda_idx <- names(which.min(nz_counts[nz_counts >= nCuts_per_group]))

  if (sum(coef_init[ , lambda_idx] != 0) == nCuts_per_group)
    return(names(coef_init[ , lambda_idx][coef_init[ , lambda_idx] != 0]))

  coef_nz <- sort(abs(coef_init[ , lambda_idx][coef_init[ , lambda_idx] != 0]), decreasing = T)
  return(names(coef_nz[c(1 : nCuts_per_group)]))
}



#' Find optimal cut-points where number of cut-points per predictor is fixed
#'
#' @param x input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format (inherit from class "sparseMatrix" as in package Matrix). Requirement: nvars >1; in other words, x should have 2 or more columns.
#' @param y response variable. A matrix with two columns for a Cox family model, and a vector otherwise
#' @param max_nCuts integer, maximum number of cut-points allowed for each predictor
#' @param method specifies biniLasso vs. miniLasso, or "both".
#' @param family any glmnet family option is available.
#' @param lasso_rule Lasso rule for finding optimal lambda; either "1se" (default) or "min".
#' @param lasso_nfolds The number of lambda values - default is 10.
#' @param penalty.factor Separate penalty factors can be applied to each coefficient. This is a number that multiplies lasso lambda to allow differential shrinkage. Can be 0 for some variables, which implies no shrinkage, and that variable is always included in the model. Default is 1 for all variables . Also, any penalty.factor that is set to inf is converted to an exclude, and then internally reset to 1. Note: the penalty factors are internally rescaled to sum to nvars, and the lasso lambda sequence will reflect this change.
#' @param lambda_grid a vector of sequence of lambda values. Can be left NULL.
#' @param cols a vector of numeric column names (characters).
#' @param x_cuts a list of cut-points corresponding to each entry of cols and in the same order.
#'
#' @returns a data.frame in which, each row corresponds to a method, and a column-list containing optimila cut-points per predictor
#' @export
#'
#' @examples
#'
#'
#' simData <- data.frame(t = rexp(1000),
#'                       event = rbinom(1000, 1, 0.6),
#'                       x1 = rnorm(1000),
#'                       x2 = rnorm(1000),
#'                       x3 = rnorm(1000))
#'
#' simData_converted <-
#'   cumBinarizer(data = simData,
#'                cols = c("x1", "x2", "x3"),
#'                method = "quantile",
#'                n_bins = 30)
#'
#' simData_sim_optCuts <-
#'  opt_cuts_finder(x = simData_converted$x,
#'                  y = survival::Surv(simData$t, simData$event),
#'                  method = "biniLasso",
#'                  family = "cox",
#'                  lasso_rule = "1se",
#'                  lasso_nfolds = 10,
#'                  penalty.factor = NULL,
#'                  cols = c("x1", "x2", "x3"),
#'                  x_cuts = simData_converted$x_cuts)
#'
#' opt_fixed_nCuts(data = simData_converted,
#'                 optCuts = .,
#'                 y = Surv(simData_converted$tte, simData_converted$vital_status),
#'                 family = "cox",
#'                 max_nCuts = 2,
#'                 method = "biniLasso",
#'                 cols = c("x1", "x2", "x3"),
#'                 x_cuts = simData_converted$x_cuts)
#'
#'
#'
opt_fixed_nCuts <-
  function (x = gbm_converted_obj$x,
            y = survival::Surv(gbm_data_fnl$tte,
                               gbm_data_fnl$vital_status),
            max_nCuts = 2,
            method = "biniLasso", family = "cox", lasso_rule = "1se",
            lasso_nfolds = 10, penalty.factor = NULL,
            lambda_grid = NULL,
            cols, x_cuts)
  {
    if (method == "both")
      method <- c("biniLasso", "miniLasso")
    if (is.null(penalty.factor)) penalFac_null_flag <- TRUE
    else penalFac_null_flag <- FALSE

    if ("biniLasso" %in% method) {
      cols_selected <- character()
      for (g in unique(cols)) {
        X_mat_tmp <- x[ , grepl(g, colnames(x))]
        cols_selected <- c(cols_selected,
                           cuts_finder(x = X_mat_tmp,
                                       y = y,
                                       nCuts_per_group = max_nCuts,
                                       lambda_grid = lambda_grid,
                                       miniLasso = FALSE))
      }
      x_selected <- x[ , cols_selected]
      if (penalFac_null_flag) penalty.factor <- rep(1, ncol(x_selected))

      bini_cv <- glmnet::cv.glmnet(x = x_selected, y = y, family = family,
                                   nfolds = lasso_nfolds, penalty.factor = penalty.factor)
      bini_fit <- glmnet::glmnet(x = x_selected, y = y, family = family,
                                 lambda = ifelse(lasso_rule == "min", bini_cv$lambda.min,
                                                 bini_cv$lambda.1se),
                                 penalty.factor = penalty.factor)
      x_cuts_bini_opt <- cuts_extractor(glm_fit = bini_fit,
                                        cols = cols, x_cuts = x_cuts)
      x_cuts_bini_opt <- x_cuts_bini_opt[[1]]
      names(x_cuts_bini_opt) <- cols
      x_cuts_bini_opt <- dplyr::tibble(method = "biniLasso",
                                       opt_cuts = list(x_cuts_bini_opt))
    }
    if ("miniLasso" %in% method) {
      cols_selected <- character()
      for (g in unique(cols)) {
        X_mat_tmp <- x[ , grepl(g, colnames(x))]
        cols_selected <- c(cols_selected,
                           cuts_finder(x = X_mat_tmp,
                                       y = y,
                                       nCuts_per_group = max_nCuts,
                                       lambda_grid = lambda_grid,
                                       miniLasso = TRUE))
      }
      x_selected <- x[ , cols_selected]
      if (penalFac_null_flag) penalty.factor <- rep(1, ncol(x_selected))

      ubini_cv <- uniLasso::cv.uniLasso(x = x_selected, y = y, family = family,
                                        nfolds = lasso_nfolds, penalty.factor = penalty.factor)
      ubini_fit <- uniLasso::uniLasso(x = x_selected, y = y, family = family,
                                      penalty.factor = penalty.factor)
      x_cuts_ubini_opt <- cuts_extractor(glm_fit = ubini_fit,
                                         cols = cols, x_cuts = x_cuts,
                                         lambda_opt = ifelse(lasso_rule == "min",
                                                             ubini_cv$lambda.min,
                                                             ubini_cv$lambda.1se))
      x_cuts_ubini_opt <- x_cuts_ubini_opt[[1]]
      names(x_cuts_ubini_opt) <- cols
      x_cuts_ubini_opt <- dplyr::tibble(method = "miniLasso",
                                        opt_cuts = list(x_cuts_ubini_opt))
    }
    if (all(method == "biniLasso"))
      return(x_cuts_bini_opt)
    if (all(method == "miniLasso"))
      return(x_cuts_ubini_opt)
    if (all(method == c("biniLasso", "miniLasso")))
      return(x_cuts_bini_opt %>% dplyr::bind_rows(x_cuts_ubini_opt))
    return(x_cuts_opt)
  }


