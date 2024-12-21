mean_ci_boot <- function(boot, lci=0.025, uci=0.975) {
  n<-length(boot[[1]])
  out<-data.frame()
  for (i in 1:n) {
    out<-rbind(out,
               c(mean(unlist(lapply(boot, function(x) x[i]))),
                 quantile(unlist(lapply(boot, function(x) x[i])), probs=c(lci, uci), na.rm=T)))
  }
  names(out)<-c("Mean", "LCI", "UCI")
  rownames(out)<-names(boot[[1]])
  return(out)
}