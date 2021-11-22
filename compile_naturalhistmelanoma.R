
compile_naturalhist_weibull <- function(prop_adv, mortrates, shapes, subgroup_probs) {
  
  stage_probs <- c(Early=1-prop_adv, Advanced=prop_adv)
  df <- lapply(names(mortrates),
               function(x) {
                 return(data.frame(prop=subgroup_probs*stage_probs[x],
                                   stage=x,
                                   subgroup=names(subgroup_probs),
                                   mortrate=mortrates[x],
                                   shape=shapes[x],
                                   stringsAsFactors=FALSE))
               })
  
  df <- plyr::ldply(df)
  if (round(sum(df$prop),1)!=1) stop('Check that subgroup_probs sum to 1')
  class(df) <- append(class(df), 'naturalhist')
  return(df)
}


create_stageshift_map_melanoma <- function(x) {

  if (!'naturalhist'%in%class(x)) stop('x must have class naturalhist')
  create_stageshift_map <- function(x) {

    if (!'naturalhist'%in%class(x)) stop('x must have class naturalhist')
    subgroups <- as.character(unique(x$subgroup))
    stage_pairs <- sapply(subgroups, function(group, df) which(x$subgroup==group), x)
    rownames(stage_pairs) <- x$stage[stage_pairs[,1]]

    return(stage_pairs)
  }
  stage_pairs <- sapply(subgroups, function(group, df) which(x$subgroup==group), x)
  rownames(stage_pairs) <- x$stage[stage_pairs[,1]]

  return(stage_pairs)
}

