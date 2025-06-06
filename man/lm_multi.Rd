% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regression.R
\name{lm_multi}
\alias{lm_multi}
\title{Linear regression for multiple conditions}
\usage{
lm_multi(feature_id, ref_conds, metadata, df, quant_var = "log2_abundance")
}
\arguments{
\item{feature_id}{groupId(s) from mzrolldb}

\item{ref_conds}{RefCondition(s) in metadata}

\item{metadata}{metadata of experiment}

\item{df}{dataframe in long format}

\item{quant_var}{string of quant variable to be used for lm}
}
\value{
one-row data from linear regression
}
\description{
Linear regression for multiple conditions
}
