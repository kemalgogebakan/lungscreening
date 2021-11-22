compile_naturalhist_lung <- function(prop_adv, mortrates, subgroup_names, sub_probs_early, sub_probs_adv) {
  stage_probs1<-c(Early=1-prop_adv, Advanced=prop_adv)
  stage_probs2<-c(1-prop_adv, prop_adv)
  propss<-c(1:4)
  propss[1:2]<-sub_probs_early*stage_probs2[1]
  propss[3:4]<-sub_probs_adv*stage_probs2[2]
  df <- lapply(names(stage_probs1),
               function(x) {
                 return(data.frame(prop=propss,
                                   stage=x,
                                   subgroup=names(subgroup_names),
                                   mortrate=mortrates,
                                   stringsAsFactors=FALSE))
               })

 df <- plyr::ldply(df)
 df <- df %>% slice(-c(3:6))
  if (round(sum(df$prop),1)!=1) stop('Check that subgroup_probs sum to 1')
  class(df) <- append(class(df), 'naturalhist')
  return(df)
}