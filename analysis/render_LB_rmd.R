## Script to render all .Rmds from LB's DAR analysis
## Using this to render selected .Rmds because analysis_final.rmd errors out
##  on my end for some reason

rmds <- c(
  "analysis/dar-analysis_AB.Rmd",
  "analysis/dar-analysis_AC.Rmd",
  "analysis/dar-analysis_BC.Rmd",
  "analysis/enrichment_AB.Rmd",
  "analysis/enrichment_AC.Rmd",
  "analysis/enrichment_BC.Rmd"
)

workflowr::wflow_build(rmds)
