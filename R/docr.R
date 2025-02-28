#' CoxPH Modeling
#'
#' @description
#' This function performs CoxPH modeling using the survival::coxph package and
#' return a tidy summary of coefficients
#'
#' @param data_use dataframe
#'
#' @param form right-hand side of the formula with covariates for CoxPH
#'
#' @param risk_var column name of risk data (i.e. lifespan), as a string
#'
#' @param event_var column name of risk data censorship column, as a string.
#' Must be boolean or logical
#'
#' @param value_var column names of metabolite values, as a string. Defaults to
#' \code{norm_abundance}
#'
#' @param n_min minimum number of unique samples required to run coxph.
#' Set to low number to avoid model skipping on low n datasets
#'
#' @importFrom magrittr %>%
#'
#' @returns dataframe of tidy coefficients from coxph modeling
#'
#' @export
coxph_model <- function(data_use,
                        form,
                        risk_var,
                        event_var,
                        value_var = "norm_abundance",
                        n_min = 100) {
  if (!risk_var %in% colnames(data_use) || !is.numeric(data_use[[risk_var]])) {
    stop("risk_var must be a numeric variable column in data_use")
  }
  if (!event_var %in% colnames(data_use)) {
    stop("event_var must be a boolean or binary column in data_use")
  }

  dfx <- data_use %>%
    dplyr::mutate(surv_days.surv = survival::Surv(
      !!rlang::sym(risk_var),
      !!rlang::sym(event_var)
    ))

  if (value_var %in% colnames(dfx)) {
    dfx <- dfx %>%
      tidyr::drop_na(!!rlang::sym(value_var))
  }
  if (nrow(dfx) == 0 || nrow(dfx) < n_min) {
    return(data.frame())
  }
  if (!grepl("~", form)) {
    formula_use <- paste0("surv_days.surv ~ ", form)
  } else {
    formula_use <- paste0("surv_days.surv ", form)
  }

  warns <- NULL
  errs <- NULL
  tryCatch(
    {
      cox_fit <- do.call(
        what = survival::coxph,
        args = list(
          formula = stats::formula(
            x = formula_use,
            env = globalenv()
          ),
          data = dfx
        )
      )

      return_df <- summary(cox_fit)$coefficients %>%
        as.data.frame() %>%
        dplyr::bind_cols(exp(confint(cox_fit))) %>%
        tibble::rownames_to_column("model_term") %>%
        dplyr::mutate(model = formula_use)

      return(return_df)
    },
    error = function(e) {
      errs <<- c(errs, list(e))
    },
    warning = function(w) {
      warns <<- c(warns, list(w))
    }
  )
}


#' Test and Train Data Parsing
#'
#' @description
#' Given a dependent variable, development IDs, and validation IDs,
#' this function returns a list with 4 dataframes corresponding to development
#' and validation dataframes for x and y variables. This function utilizes
#' \code{glmnet} and is designed to work with
#' \code{metstats::metstats_elastic_net} using the development dataframes in
#' the return list
#'
#' @param df a dataframe containing columns described below
#' @param dependent_var column name of regression outcome varibale,
#' as a string
#' @param sample_id_var column name of observation identifiers, as a string,
#' on which \code{dev_ids} and \code{val_ids} will be parsed
#' @param dev_ids character vector corresponding to values in
#' \code{sample_id_var} on which to filter observations for model development
#' @param val_ids character vector corresponding to values in
#' \code{sample_id_var} on which to filter observations for model validation
#' @param key_var column name of identifier variables to be included on the
#' right hand side of the model equation, as a string; defaults to "feature_id"
#' @param value_var column name of values associated with \code{value_var} to
#' be included on the right hand side of the equation, as a string; defaults to
#' "norm_abundance"
#' @param key_var_filter character string on which variables associated with
#' \code{key_var} can be filtered out with grep; i.e. all unknown compounds
#' will be removed when default \code{NULL} is changed to "unk"
#' @param makeX_na_impute logical to indicate if the \code{glmnet::makeX}
#' function should handle missing values (\code{metstats::metstats_elastic_net}
#' will not work with missing values)
#'
#' @importFrom magrittr %>%
#'
#' @returns a list containing dataframes x.dev, y.dev, x.val, y.val for direct
#' input into \code{metstats::metstats_elastic_net}
#'
#' @export
metstats_elastic_data_parse <- function(df,
                                        dependent_var,
                                        sample_id_var,
                                        dev_ids,
                                        val_ids,
                                        key_var = "feature_id",
                                        value_var = "norm_abundance",
                                        key_var_filter = NULL,
                                        makeX_na_impute = TRUE) {
  checkmate::assertString(dependent_var)
  checkmate::assertString(sample_id_var)
  checkmate::assertCharacter(dev_ids)
  checkmate::assertCharacter(val_ids)
  checkmate::assertString(key_var)
  checkmate::assertString(value_var)
  checkmate::assertLogical(makeX_na_impute)

  x_temp <- df %>%
    dplyr::select(
      rlang::sym(value_var), rlang::sym(key_var),
      rlang::sym(dependent_var), rlang::sym(sample_id_var)
    )

  if (!is.null(key_var_filter)) {
    checkmate::assertCharacter(key_var_filter)
    x_temp <- x_temp %>%
      dplyr::filter(!grepl(key_var_filter, !!rlang::sym(key_var)))
  }

  x.dev <- x_temp %>%
    dplyr::filter(!!rlang::sym(sample_id_var) %in% dev_ids) %>%
    tidyr::spread(!!rlang::sym(key_var), !!rlang::sym(value_var)) %>%
    dplyr::select(-rlang::sym(dependent_var)) %>%
    tibble::column_to_rownames(sample_id_var) %>%
    glmnet::makeX(na.impute = makeX_na_impute)
  x.dev <- x.dev[, !apply(x.dev, 2, function(x) all(is.na(x)))]

  y.dev <- x_temp %>%
    dplyr::filter(!!rlang::sym(sample_id_var) %in% dev_ids) %>%
    tidyr::spread(!!rlang::sym(key_var), !!rlang::sym(value_var)) %>%
    dplyr::select(rlang::sym(dependent_var), rlang::sym(sample_id_var)) %>%
    tibble::column_to_rownames(sample_id_var) %>%
    as.matrix()

  x.val <- x_temp %>%
    dplyr::filter(!!rlang::sym(sample_id_var) %in% val_ids) %>%
    tidyr::spread(!!rlang::sym(key_var), !!rlang::sym(value_var)) %>%
    dplyr::select(-rlang::sym(dependent_var)) %>%
    tibble::column_to_rownames(sample_id_var) %>%
    glmnet::makeX(na.impute = makeX_na_impute)
  x.val <- x.val[, !apply(x.val, 2, function(x) all(is.na(x)))]

  y.val <- x_temp %>%
    dplyr::filter(!!rlang::sym(sample_id_var) %in% val_ids) %>%
    tidyr::spread(!!rlang::sym(key_var), !!rlang::sym(value_var)) %>%
    dplyr::select(rlang::sym(dependent_var), rlang::sym(sample_id_var)) %>%
    tibble::column_to_rownames(sample_id_var) %>%
    as.matrix()

  return(list(x.dev = x.dev, y.dev = y.dev, x.val = x.val, y.val = y.val))
}


#' Elastic Net Regularization Regression
#'
#' @description
#' This function trains an elastic net model for multiple alpha levels.
#' There are options to not penalize specific variables (as grepped column
#' names), change the number of iterations (default `its = 100`), and change
#' the percent of development dataset that is used for training versus
#' testing for alpha value (default `train_percent = 0.75`). This function
#' runs \code{glmnet::cv.glmnet}.
#'
#' @param x a dataframe with observations on rows and columns as model
#' ovariates
#' @param y a double or matrix of outcome values with observations
#' corresponding to those in \code{x}
#' @param vars_not_penalized character vector of column names in \code{x} that
#' will recieve no penalty in \code{glmnet::cv.glmnet}
#' @param its number of iterations to run the model; defaults to 100
#' @param train_percent percent of dataset provided on which to determine the
#' appropriate alpha based on lowest MSE; defaults to 0.75
#' @param alpha_levels number of alpha levels to test which will be divided
#' between 0 and 1 proportionally; defaults to (defaults to 5 which corresponds
#' to alpha = 0, 0.25, 0.5, 0.75, and 1)
#' @param key_var column name of identifier variables, as a string;
#' defaults to "feature_id". Used for naming the column in the return dataframe
#' and does not impact modeling
#' @param makeX_na_impute logical to indicate if the \code{glmnet::makeX}
#' function should handle categorical values (convert to binaries) and impute
#' means on missing values
#'
#' @importFrom magrittr %>%
#'
#' @returns a list containing 2 dataframes: \code{Selection} which contains all
#' iteration summary data, and \code{Models} which returns all stored model fits
#'
#' @export
metstats_elastic_net <- function(x,
                                 y,
                                 vars_not_penalized = NULL,
                                 its = 100,
                                 train_percent = 0.75,
                                 alpha_levels = 5,
                                 key_var = "feature_id",
                                 makeX_na_impute = FALSE) {
  checkmate::assertNumber(its)
  checkmate::assertNumber(train_percent)
  checkmate::assertNumber(alpha_levels)
  checkmate::assertLogical(makeX_na_impute)

  selection_df <- data.frame()
  model_list <- list()

  for (j in 1:its) {
    ## Split data into training and testing data sets
    train_rows <- sample(1:nrow(x), train_percent * nrow(x))
    x.train <- x[train_rows, ] %>%
      data.frame(check.names = FALSE)
    x.test <- x[-train_rows, ] %>%
      data.frame(check.names = FALSE)

    ## makeX_na_impute will change categorical variables into binaries,
    ## and will impute missing values
    if (makeX_na_impute) {
      x.train <- glmnet::makeX(x.train, na.impute = makeX_na_impute)
      x.test <- glmnet::makeX(x.test, na.impute = makeX_na_impute)
    } else {
      x.train <- as.matrix(x.train)
      x.test <- as.matrix(x.test)
    }

    y.train <- y[train_rows, ]
    y.test <- y[-train_rows, ]

    ## Manually set penalty factors if provided column name variables indicated
    p.fac <- rep(1, ncol(x.train))
    if (!is.null(vars_not_penalized)) {
      p.fac[which(grepl(vars_not_penalized, colnames(x.train)))] <- 0
      if (j == 1) {
        print(paste0("Not penalized: ", grep(vars_not_penalized, colnames(x.train), value = T)))
      }
    }

    ## Test for various alpha levels (default of 5 yield alphas of 0, 0.25, 0.5, 0.75, 1)
    list.of.fits <- list()
    for (i in 0:(alpha_levels - 1)) {
      fit.name <- paste0("alpha", i / (alpha_levels - 1))

      ## Now fit a model (i.e. optimize lambda) and store it in a list that
      ## uses the variable name we just created as the reference.
      list.of.fits[[fit.name]] <-
        glmnet::cv.glmnet(x.train,
          y.train,
          type.measure = "mse",
          alpha = i / (alpha_levels - 1),
          penalty.factor = p.fac,
          family = "gaussian"
        )
    }

    ## Find which alpha does best
    results <- data.frame()
    for (i in 0:(alpha_levels - 1)) {
      fit.name <- paste0("alpha", i / (alpha_levels - 1))

      ## Use each model to predict 'y' given the test dataset
      predicted <- stats::predict(list.of.fits[[fit.name]],
        s = list.of.fits[[fit.name]]$lambda.1se,
        newx = x.test
      )

      ## Calculate the mean squared error to determine which alpha
      ## performed best
      mse <- mean((y.test - predicted)^2)
      r.squared <- round(summary(lm(predicted ~ y.test))$r.squared * 100, 2)

      ## Store the results
      temp <- data.frame(
        alpha = i / (alpha_levels - 1),
        mse = mse,
        r.squared = r.squared,
        fit.name = fit.name
      )
      results <- dplyr::bind_rows(results, temp)
    }

    # Identify alpha with lowest MSE
    best_alpha_temp <- results %>%
      dplyr::filter(mse == min(mse))

    # If there is a minimum, save that
    if (nrow(best_alpha_temp) == 1) {
      best_alpha <- best_alpha_temp %>%
        dplyr::pull(fit.name)

      # If multiple alphas perform the same, find the mean alpha, and select the
      # closest to the minimum alpha mean (i.e. if 0, 0.25, and 0.)
    } else {
      best_alpha <- best_alpha_temp %>%
        dplyr::mutate(
          alpha_mean = mean(alpha),
          alpha_diff = alpha_mean - alpha
        ) %>%
        dplyr::filter(abs(alpha_diff) == min(abs(alpha_diff))) %>%
        # If models from 2 alphas are equidistant from the minimum alpha mean,
        # randomly select one of the models
        dplyr::slice_sample(n = 1) %>%
        dplyr::pull(fit.name)
    }

    ## View and predictor variables
    coefs <- predict(list.of.fits[[best_alpha]], type = "coef")[1:nrow(predict(list.of.fits[[best_alpha]], type = "coef")), ] %>%
      as.data.frame() %>%
      dplyr::rename(coef = ".") %>%
      dplyr::filter(coef != 0) %>%
      tibble::rownames_to_column(key_var) %>%
      dplyr::mutate(!!rlang::sym(key_var) := gsub("`", "", !!rlang::sym(key_var)),
        alpha = best_alpha,
        iteration = j,
        mse = best_alpha_temp$mse[1],
        r.squared = best_alpha_temp$r.squared[1]
      )

    selection_df <- dplyr::bind_rows(selection_df, coefs)
    model_list[[j]] <- list.of.fits[[best_alpha]]
  }

  return(list(
    Selection = selection_df,
    Models = model_list
  ))
}


#' Linear mixed modeling
#'
#' @description
#' This function trains a linear mixed model through \code{lmer} as specified
#'
#' @param data_use a dataframe with \code{value_var} numeric values
#' corresponding to the feature identified by \code{feature_id_test},
#' containing columns associated with covariates in \code{forms_to_test}
#' @param feature_id_test name of feature being tested. This is not a string
#' for filtering, but is used for print messages during model running. Helpful
#' when model fitting features in parallel
#' @param forms_to_test formulates to test, as character vector. Only include
#' covariates to the right of the tilde
#' @param value_var column names of dependent values, as a string. Defaults to
#' \code{norm_abundance} and is pasted to the left of the tilde in the formula
#' @param n_samples number of unique samples in experiment. Defaults to
#' \code{NULL} to ignore parameter.
#' @param n_sample_threshold percent of \code{n_samples} below which model
#' fitting will skip. Defaults to 0.1 for 10%.
#' @param model_type type of model to run, one of \code{"mixed"} for
#' \code{lmer()} fit, or \code{"simple"} for \code{lm()} fit
#'
#' @importFrom magrittr %>%
#'
#' @returns a dataframe containing the feature_id, model_term, model, coef,
#' se_coef, p_value_coef, t_value, and AIC of the model fit
#'
#' @export
lmer_multi_formula <- function(data_use,
                               feature_id_test,
                               forms_to_test,
                               value_var = "norm_abundance",
                               n_samples = NULL,
                               n_sample_threshold = 0.1,
                               model_type = "mixed") {
  checkmate::assertNumeric(data_use[[value_var]])

  # Initialized empty dataframe
  da <- data.frame()

  # Drop NA values and assign `norm_abundance` to the value variable
  dfx <- data_use %>%
    tidyr::drop_na(!!rlang::sym(value_var)) %>%
    dplyr::mutate(norm_abundance = scale(!!rlang::sym(value_var),
      scale = FALSE, center = TRUE
    ))

  # If there is no data or if data is only present in < n_sample_threshold
  # percent of n_samples, skip modeling
  if (nrow(dfx) == 0) {
    print(paste0(
      "For feature_id ", feature_id_test,
      ", skipping modeling; no data points"
    ))
    return(da)
  }
  if (!is.null(n_samples)) {
    if (nrow(dfx) < n_samples * n_sample_threshold) {
      print(paste0(
        "For feature_id ", feature_id_test,
        paste0(
          ", skipping modeling; data present in < ",
          as.character(n_sample_threshold * 100), " % of samples"
        )
      ))
      return(da)
    }
  }

  # Run analysis through all formulas
  for (f in forms_to_test) {
    # Properly format equation
    if (!grepl("~", f)) {
      f_use <- paste0("norm_abundance ~ ", f)
    } else {
      f_use <- paste0("norm_abundance ", f)
    }

    if (model_type == "mixed") {
      # If only one data point per mouse id, change random intercept to
      # generation_wave
      # Will skip if these columns are not in dataframe (for DO-CR)
      if (all(c("mouse_id", "generation_wave") %in% colnames(dfx))) {
        if (length(unique(dfx$mouse_id)) == nrow(dfx)) {
          f_use <- gsub("mouse_id", "generation_wave", f_use)
          print(paste0(
            "For feature_id ", feature_id_test, ", formula ",
            f_use, ", random intercept changed to generation_wave"
          ))
        }
      }

      # Run linear mixed model on data subset, for each formula
      model1 <- suppressMessages(lmerTest::lmer(formula = f_use, data = dfx))

      # Extract model summary data, including coefficient
      # P-values derived from the lmerTest library
      dfc <- as.data.frame(coef(summary(model1)))
      colnames(dfc) <- c("coef", "se_coef", "df", "t_value", "p_value_coef")
    } else if (model_type == "simple") {
      # Run simple linear model on data subset, for each formula
      model1 <- stats::lm(formula = f_use, data = dfx)

      # Extract model summary data, including coefficient
      dfc <- as.data.frame(coef(summary(model1)))
      colnames(dfc) <- c("coef", "se_coef", "t_value", "p_value_coef")
    } else {
      print("'model_type' must be one of 'mixed' or 'simple'; skipping modeling")
      return(da)
    }

    # Add other information that we want in DA dataframe
    dfc <- dfc %>%
      tibble::rownames_to_column(var = "model_term") %>%
      dplyr::mutate(
        model_term = factor(model_term),
        feature_id = .env$feature_id_test,
        model = .env$f_use,
        AIC = AIC(.env$model1)
      ) %>%
      dplyr::select(
        feature_id, model_term, model,
        coef, se_coef, p_value_coef, t_value, AIC
      )

    # Bind with results from other models
    da <- da %>%
      dplyr::bind_rows(dfc)
  }

  return(da)
}


#' Find outliers from numerical data columns
#'
#' @description
#' This function determines outliers and returns the sample names of outliers
#' based on values in \code{outlier_col_names} that are \code{outlier_sd}
#' (default = 3) standard deviations above or below the mean for all columns
#' provided. The default columns are "PC1", "PC2", and "PC3" from principle
#' component analysis.
#'
#' @param df a dataframe containing, at minimum, columns corresponding to
#' \code{outlier_col_names} and \code{sample_name_col}
#' @param outlier_col_names column name(s) as string or character vector on
#' which to calculate mean and standard deviation
#' @param outlier_sd standard deviation above or below which a sample will
#' be flagged as an outlier
#' @param sample_name_col column name as string of outlier sample names to return
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns a character vector of \code{sample_name_col} values that were
#' identified as outliers
#'
#' @export
find_outliers <- function(df,
                          outlier_col_names = c("PC1", "PC2", "PC3"),
                          outlier_sd = 3,
                          sample_name_col = "sample_name") {
  # Calculate mean and standard deviation on numerical columns
  for (col in outlier_col_names) {
    outlier_col_mean <- mean(df[[col]], na.rm = TRUE)
    outlier_col_sd <- stats::sd(df[[col]], na.rm = TRUE)
    outlier_col_lower <- outlier_col_mean - (outlier_sd * outlier_col_sd)
    outlier_col_upper <- outlier_col_mean + (outlier_sd * outlier_col_sd)

    df <- df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(!!rlang::sym(paste0(col, "_outlier")) :=
        ifelse(!!rlang::sym(col) < outlier_col_lower ||
          !!rlang::sym(col) > outlier_col_upper,
        TRUE,
        FALSE
        )) %>%
      dplyr::ungroup()
  }

  outlier_cols <- grep("_outlier", names(df), value = TRUE)
  outlier_names <- df[rowSums(df[outlier_cols] == TRUE) == length(outlier_cols), sample_name_col] %>%
    tidyr::drop_na() %>%
    dplyr::pull()

  return(outlier_names)
}


#' ANOVA modeling for variance partitioning
#'
#' @description
#' This function runs \code{anova()} on one or more formulas to determine the
#' percent variance explained by each variable.
#'
#' @param data_use a dataframe with \code{value_var} numeric values
#' corresponding to the feature identified by \code{feature_id_test},
#' containing columns associated with covariates in \code{forms_to_test}
#' @param feature_id_test name of feature being tested. This is not a string
#' for filtering, but is used for print messages during model running.
#' Helpful when model fitting features in parallel
#' @param forms_to_test formulates to test, as character vector. Only
#' include covariates to the right of the tilde
#' @param value_var column names of dependent values, as a string. Defaults to
#' \code{norm_abundance} and is pasted to the left of the tilde in the formula
#' @param n_min minimum number of unique samples required to run ANOVA
#' @param model_type type of model to run, one of \code{"mixed"} for
#' \code{lmer()} fit, or \code{"simple"} for \code{lm()} fit
#' @param per_exp_var variable on which to calculate the percent of variance
#' explained, one of \code{"Sum Sq"} for net sum of squares, or
#' \code{"Mean Sq"} for mean sum of squares for each variable (sum of squares
#' divided by the degrees of freedom)
#'
#' @importFrom magrittr %>%
#'
#' @returns a dataframe containing the model terms and percent variance
#' explained from \code{anova()}
#'
#' @export
anova_multi_formula <- function(data_use,
                                feature_id_test,
                                forms_to_test,
                                value_var = "norm_abundance",
                                n_min = 100,
                                model_type = "simple",
                                per_exp_var = "Sum Sq") {
  checkmate::assertNumeric(data_use[[value_var]])

  # Initialized empty dataframe
  da <- data.frame()

  # Drop NA values and assign `norm_abundance` to the value variable
  dfx <- data_use %>%
    tidyr::drop_na(!!rlang::sym(value_var)) %>%
    dplyr::mutate(norm_abundance = scale(!!rlang::sym(value_var), scale = F, center = T))

  # If there is no data or if data is only present in < n_min samples, skip modeling
  if (nrow(dfx) == 0 || nrow(dfx) < n_min) {
    print(paste0("For feature_id ", feature_id_test, ", skipping modeling; not enough data points"))
    return(da)
  }

  # DOCR specific check
  if ("diet" %in% colnames(dfx)) {
    if (length(unique(dfx$diet)) < 5) {
      print(paste0("For feature_id ", feature_id_test, ", skipping modeling; < 5 diets represented"))
      return(da)
    }
  }

  # Run Variance analysis
  for (f in forms_to_test) {
    # Properly format equation
    if (!grepl("~", f)) {
      f_use <- paste0("norm_abundance ~ ", f)
    } else {
      f_use <- paste0("norm_abundance ", f)
    }

    # Run linear model on data subset, and run anova to determine percent
    # of variance explained. The variance can be decomposed on the linear mixed
    # effect model as well, but the residual variance cannot be extracted
    if (model_type == "simple") {
      anova_model <- stats::anova(stats::lm(formula = f_use, data = dfx))
    } else if (model_type == "mixed") {
      anova_model <- stats::anova(lme4::lmer(formula = f_use, data = dfx))
    } else {
      stop("model_type must be one of `simple` for lm() linear model, or `mixed` for lmer() mixed effects model")
    }

    ss <- anova_model[[per_exp_var]]
    percent_explained <- anova_model %>%
      dplyr::bind_cols(PerExp = ss / sum(ss) * 100) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("model_term") %>%
      dplyr::mutate(
        model = .env$f,
        feature_id = .env$feature_id_test
      )

    # Bind with results from other models
    da <- da %>%
      dplyr::bind_rows(percent_explained)
  }

  return(da)
}
