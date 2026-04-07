library(tidyverse)
#set local wd

localDir = "HMSC_MER"
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir,"models")
resultDir = file.path(localDir,"results")


library(Hmsc)

# model parameters
thin <- 50
samples <- 3000
nChains <- 2
filename=file.path(modelDir, paste0("MER_model_thin_",as.character(thin),"_samples_",as.character(samples),"_chains_",as.character(nChains),'.Rdata'))

coda.object <- convertToCodaObject(models[[1]])


# view trace for Gamma parameter
viewTraceGamma <- function(coda.object){
  out.list <- lapply(coda.object$Gamma, as.data.frame)
  out.list <- lapply(out.list, function(x){
    x$sample <- 1:nrow(x)
    return(x)})
  out.df <- bind_rows(out.list, .id='chain') %>% 
    pivot_longer(cols=ends_with(']'), names_to='parameter', values_to='value') %>% 
    mutate(parameter = gsub('\\[', ', ', parameter), 
           parameter = gsub('\\]$', '', parameter)) %>% 
    separate(parameter, c('estimate', 'parameter', 'species'), sep=', ') %>% 
    select(-species)
  
  ggplot(out.df, aes(x=sample, y=value, colour=chain)) + 
    geom_line() + 
    geom_hline(aes(yintercept=0)) + 
    facet_wrap(vars(parameter), scales='free_y')
}

# view trace for Beta parameter
viewTraceBeta <- function(coda.object, file=NULL, height=10, width=10){
  out.list <- lapply(coda.object$Beta, as.data.frame)
  out.list <- lapply(out.list, function(x){
    x$sample <- 1:nrow(x)
    return(x)})
  out.df <- bind_rows(out.list, .id='chain') %>% 
    pivot_longer(cols=ends_with(']'), names_to='parameter', values_to='value') %>% 
    mutate(parameter = gsub('\\[', ', ', parameter), 
           parameter = gsub('\\]$', '', parameter)) %>% 
    separate(parameter, c('estimate', 'parameter', 'species'), sep=', ')
  
  plot_list <- lapply(split(out.df, out.df$species), function(x){
    sp <- x$species[1]
    ggplot(x, aes(x = sample, y = value, colour = chain)) + 
      geom_line() + 
      geom_hline(aes(yintercept = 0)) + 
      labs(title = sp) + 
      facet_wrap(vars(parameter), scales = 'free_y')
  })
  # View plots
  for (p in plot_list) print(p)
  
  invisible(plot_list)  # return silently if you want to store or use it later
}

# summarise trace for Beta parameter
summaryTraceBeta <- function(coda.object){
  out.list <- lapply(coda.object$Beta, as.data.frame)
  out.list <- lapply(out.list, function(x){
    x$sample <- 1:nrow(x)
    return(x)})
  out.df <- bind_rows(out.list, .id='chain') %>% 
    pivot_longer(cols=ends_with(']'), names_to='parameter', values_to='value') %>% 
    mutate(parameter = gsub('\\[', ', ', parameter), 
           parameter = gsub('\\]$', '', parameter)) %>% 
    separate(parameter, c('estimate', 'parameter', 'species'), sep=', ')
}

# view trace for Sigma parameter
viewTraceSigma <- function(coda.object){
  out.list <- lapply(coda.object$Sigma, as.data.frame)
  out.list <- lapply(out.list, function(x){
    x$sample <- 1:nrow(x)
    return(x)})
  out.df <- bind_rows(out.list, .id='chain') %>% 
    pivot_longer(cols=ends_with(']'), names_to='parameter', values_to='value') %>% 
    mutate(parameter = gsub('\\[', ', ', parameter), 
           parameter = gsub('\\]$', '', parameter)) %>% 
    separate(parameter, c('estimate', 'species'), sep=', ')
  
  ggplot(out.df, aes(x=sample, y=log(value), colour=chain)) + 
    geom_line() + 
    facet_wrap(vars(species), scales='free_y')
}



viewTraceGamma(coda.object)
viewTraceBeta(coda.object)
summaryTraceBeta(coda.object)
#look at each OTU
viewTraceSigma(coda.object)
