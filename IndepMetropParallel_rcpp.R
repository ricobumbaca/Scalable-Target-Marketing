IndepMetropParallel = function (Data, Draws, Mcmc, Options)
{
# Rico Bumbaca 10/10/2019
  R = dim(Draws)[1]
  B = Mcmc$B2
  keep = Mcmc$keep2
  
  betadrawList = IndepMetropParallel_rcpp_loop(Data, Draws, Mcmc, Options)
  
  if (is.null(betadrawList$betadraw))
  {
    betadraw = NULL
  } else {
    betadraw = betadrawList$betadraw
  }

  if (is.null(betadraw) || dim(betadraw)[3] == 1) betadraw = NULL
  if(Options$DiscardBetas == TRUE)
  {
    betadraw = NULL
  }
  
  return(list(betadraw=betadraw))
}
