
#' fits the models to training data and generates the predictions.
#'
#' @param allmods 
#' @param allLagsT 
#' @param allLagsO 
#' @param results.path 
#'
#' @return
#' @export
#'
#' @examples
models_lags = function(allmods, allLagsT, allLagsO, results.path){
  allfits <- map(allmods, ~gam(as.formula(.x), data=allLagsT, family=nb()))
  # only predict from AIC best model
  AICfits <- map_dbl(allfits, MuMIn::AICc)
  best <- which.min(AICfits)
  form <- Reduce(paste, deparse(allfits[[best]]$formula[3]))
  
 if (grepl("year", form) == TRUE )
{ preds <- predict_wYr(allfits[[best]], allLagsT, allLagsO)
  } else if (grepl("year", form) == FALSE ) {
    preds <- predict_noYr(allfits[[best]], allLagsT, allLagsO)
  }

  return(list(predictions= preds,
              other = list(fittedModels=allfits,
                           AICfits = AICfits,
                           best = best)))
}



