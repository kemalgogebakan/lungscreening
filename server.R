

library(shinydashboard)
library(reshape)  
library(ggplot2)
library(knitr)
library(plyr)
load("ex1.rda")
load("incidence_mel.rda")
load("incidence_risk.rda")
load("othercausemort.rda")
load("agestructure.rda")
library(bcimodel)
#bcimodel is the library for Ruth's model.
source("compile_naturalhistmelanoma.R",local = TRUE)
source("simpolicies_melanoma.R", local = TRUE)
source("simpolicies_melanoma_risk.R", local = TRUE)
source("databases_melanoma.R", local = TRUE)
source("initialize_melanoma.R", local=TRUE)
source("initialize_melanoma_risk.R", local=TRUE)
source("outcomes_melanoma.R", local=TRUE)
source("cantrance_melanoma.R", local=TRUE)
source("earlydetect_melanoma.R", local = TRUE)
source("general_melanoma.R", local = TRUE)
source("initialize_melanoma.R", local = TRUE)
source("outcomes_melanoma.R", local=TRUE)
source("parallelized_melanoma.R", local=TRUE)
source("screentreatlibrary_melanoma.R", local=TRUE)
source("systime_melanoma.R", local=TRUE)
source("treatment_melanoma.R", local=TRUE)



shinyServer <- function(input, output, session){
  
  observe({gen<-input$gender
               if(gen=="male")
                  {updateSliderInput(session, "advst", value=12.8)}
                 else if(gen=="female")
                  {updateSliderInput(session, "advst", value=9.2)}
               })
  
  observe({gen<-input$gender
  if(gen=="male")
  {updateSliderInput(session, "advin", value=6.4)}
  else if(gen=="female")
  {updateSliderInput(session, "advin", value=4.6)}
  })
  
  observe({gen<-input$gender
  if(gen=="male")
  {updateSliderInput(session, "stage1", value=73.1)}
  else if(gen=="female")
  {updateSliderInput(session, "stage1", value=79.3)}
  })
  
  observe({gen<-input$gender
  surv<-input$survival
  if(gen=="male")
  {
    if(surv==5)
    {updateSliderInput(session, "advsur", value=44.8)
      updateSliderInput(session, "earlysur", value=95.6)}
    else if(surv==10)
    {updateSliderInput(session, "advsur", value=40.3)
      updateSliderInput(session, "earlysur", value=92.7)}
  }
    
  else if(gen=="female")
    if(surv==5)
    {updateSliderInput(session, "advsur", value=53.4)
      updateSliderInput(session, "earlysur", value=96.8)}
  else if(surv==10)
  {updateSliderInput(session, "advsur", value=53.7)
    updateSliderInput(session, "earlysur", value=96.3)}
  }
  
  )
  

  
  
  sliderValues <- reactive({
    
    data.frame(
      Parameter = c("Percent advanced stage, before screening",
               "Percent advanced stage, after screening",
               "Percent reduction in advanced stage due to screening campaign or trial"),
      Value = as.character(c(input$advst,
                             input$advin,
                             round(((input$advst-input$advin)/input$advst)*100, digits = 0))),
      stringsAsFactors = FALSE)
    
  })
  
  
  
  sliderFeatures <- reactive({
    
    data.frame(
      Parameter = c("Gender", "Minimum age", "Maximum age","Year of survival statistic (k)", 
                    "Percent surviving k years at advanced stage",
                    "Percent surviving k years at early stage",
                    "Percent presenting in advanced stage, before screening",
                    "Percent presenting in advanced stage due to stage shift of screening","Risk group",
                    "Population size", "Number of simulations"),
      Value = as.character(c(input$gender,input$tfd[1],input$tfd[2],input$survival, input$advsur,
                             input$earlysur,
                             input$advst, input$advin, input$riskstrat, input$sizepop, input$nsim)),
      stringsAsFactors = FALSE)
    
  })
  
  sliderTreatmentsadvbrafp <- reactive({
    
    data.frame(
      "Treatment" = c("Chemotherapy-Immunotherapy",
                      "BRAF Inhibitor",
                      "BRAF-MEK inhibitor Combination",
                      "Anti-PD1",
                      "Anti-CTLA-4",
                      "Anti-PD1-CTLA-4 Combination",
                      "Other(Temozolomide, Corboplatin)"),
     "Historical Treatments" = as.character(c(input$advchemostbp,0,0,0,0,0,0)),
     "Novel Treatments" = as.character(c(0,input$advbrafintbp,
                                      input$advbrafmekintbp,
                                      input$advpd1intbp, input$advctla4intbp,input$advpd1ctla4intbp, input$advotherbp)),
      stringsAsFactors = FALSE)
    
  })
  
  sliderTreatmentsadvbrafn <- reactive({
    
    data.frame(
      "Treatment" = c("Chemotherapy-Immunotherapy",
                      "BRAF Inhibitor",
                      "BRAF-MEK inhibitor Combination",
                      "Anti-PD1",
                      "Anti-CTLA-4",
                      "Anti-PD1-CTLA-4 Combination",
                      "Other(Temozolomide, Corboplatin):"),
      "Historical Treatments" = as.character(c(input$advchemostbn,0,0,0,0,0,0)),
      "Novel Treatments" = as.character(c(0,input$advbrafintbn,
                                      input$advbrafmekintbn,
                                      input$advpd1intbn, input$advctla4intbn,input$advpd1ctla4intbn, input$advotherbn)),
      stringsAsFactors = FALSE)
    
  })
  
  sliderEfficacy <- reactive({
  
  data.frame(
    Stage = c("Advanced",
              "",
              "",
              "",
              "",
              "",
              ""),
    Treatment=c("Chemotherapy-Immunotherapy (BRAF +/-)",
                "BRAF Inhibitor (BRAF +)",
                "BRAF-MEK inhibitor Combination(BRAF +)",
                "Anti-PD1 (BRAF +/-)",
                "Anti-CTLA-4 (BRAF +/-)",
                "Anti-PD1-CTLA-4 Combination (BRAF +/-)",
                "Other(Temozolomide, Corboplatin (BRAF +/-))"),
    "Hazard_Ratio" = as.character(c(1,
                                   0.81,
                                    0.55,
                                    0.50,
                                    0.69,
                                    0.33,
                                    1)),
    "Survival_Improvement(%)"=as.character(c(0,19,45, 50, 31, 67,0)),
    stringsAsFactors = FALSE
  )
  })
  
  sliderPapers <- reactive({
    
    data.frame(
      Treatment = c("Chemotherapy-Immunotherapy",
                    "BRAF Inhibitor",
                    "BRAF-MEK inhibitor Combination",
                    "Anti-PD1",
                    "Anti-CTLA-4",
                    "Anti-PD1-CTLA-4 Combination",
                    "Other(Temozolomide, Corboplatin)"),
      "Reference Study"=c(  "Garbe, C., et al., Systematic review of medical treatment in melanoma: current status and future prospects. Oncologist, 2011. 16(1): p. 5-24.",
                  "Chapman, P.B., et al., Vemurafenib in patients with BRAFV600 mutation-positive metastatic melanoma: final overall survival results of the randomized BRIM-3 study. Ann Oncol, 2017. 28(10): p. 2581-2587.",
                  "Robert, C., et al., Five-Year Outcomes with Dabrafenib plus Trametinib in Metastatic Melanoma. N Engl J Med, 2019. 381(7): p. 626-636b",
                  "Pembrolizumab versus ipilimumab in advanced melanoma (KEYNOTE-006): post-hoc 5-year results from an open-label, multicentre, randomised, controlled, phase 3 study. Lancet Oncol, 2019. 20(9): p. 1239-1251",
                  "Pembrolizumab versus ipilimumab in advanced melanoma (KEYNOTE-006): post-hoc 5-year results from an open-label, multicentre, randomised, controlled, phase 3 study. Lancet Oncol, 2019. 20(9): p. 1239-1251",
                  "Larkin, J., et al., Five-Year Survival with Combined Nivolumab and Ipilimumab in Advanced Melanoma. N Engl J Med, 2019. 381(16): p. 1535-1546.",
                  "Bhatia, S., Tykodi, S. S., & Thompson, J. A. (2009). Treatment of metastatic melanoma: an overview. Oncology (Williston Park, NY), 23(6), 488.")
      ,
                  stringsAsFactors = FALSE
    )
  })
  
  slidersimulation <- reactive({
    ppv<-6.5
    braf<-input$braf1/100
    overd<-input$overd
    advanced1<-input$advst/100
    advanced2<-input$advin/100
    Stage1<-input$stage1/100
    surr<-input$survival
    if(surr==5){survv<-5} else if(surr==10){survv<-10}
    #The input data for the example are pre-loaded into the package in an object called ex1, a list with four elements:  pol, nh, map, and tx.
    melanomUS <- vector('list', length=length(ex1))
    names(melanomUS) <- names(ex1)
    propERpos<-1
    pop_chars <- list(male = data.frame(male = c(0), prop = c(1)))
    #pol:Policy
    melanomUS$pol  <- data.frame(num=c(1:4),
                                    id=c('adv.base', 
                                         'advbase.treatment','adv.shift', 'adv.shift.treatment'),
                                    name=c('M0: Historical Treatments Era without Screening',
                                                 'M1: Novel Treatments Era with No Screening',
                                                 'M2: Historical Treatments Era with Screening',
                                                 'M3: Novel Treatments Era with Screening'),
                                           pairnum=c(NA, NA, c(1,2)),
                                           earlydetHR=c(rep(1, 2), rep(input$advin/input$advst, 2)),
                                           stringsAsFactors=FALSE)
    
    #Calculating baseline survival for early stages
  
    ff = function(x, p=propERpos, surg=100/100, hsurg=1, So=input$earlysur/100) { 
      p*(1-surg)*x + p*surg*x^hsurg- So}
    
    Sb.early <- uniroot(ff, lower=0.5, upper=input$earlysur/100, tol = 0.0001)
    early.mrate.co = cumsurv_to_exprate(Sb.early$root, year=survv)
    round(100*Sb.early$root)
    
    #Calculating baseline survival for advanced stages
    f = function(x, p=propERpos, chm=input$advchemostbp/100,
    hchm=1, So=input$advsur/100) {p*(1-chm)*x + p*chm*x^hchm-So}
    
    Sb.late <- uniroot(f, lower=0.2, upper=input$advsur/100, tol = 0.0001)
    late.mrate.co = cumsurv_to_exprate(Sb.late$root, year=survv)
    round(100*Sb.late$root)
    
    #Natural history (nh) parameters
    melanomUS$nh<-compile_naturalhist(prop_adv = input$advst/100, mortrates=c(Early=early.mrate.co, Advanced=late.mrate.co),subgroup_probs=c(`BRAF+`=braf, `BRAF-`=1-braf))
    #Stage shift map
    melanomUS$map <- create_stageshift_map(melanomUS$nh)
    #Treatment
    melanomUS$tx <- data.frame(expand.grid(txSSid=c("Surgical Excision",
                                                    "Chemotherapy-Immunotherapy",
                                                    "BRAF Inhibitor",
                                                    "BRAF-MEK inhibitor Combination",
                                                    "Anti-PD1",
                                                    "Anti-CTLA-4",
                                                    "Anti-PD1-CTLA-4 Combination",
                                                    "Other(Temozolomide, Corboplatin)"),
                                              SSno=1:nrow(melanomUS$nh)),
                                  stringsAsFactors=FALSE)
    ntreat <- nrow(melanomUS$tx)
    ntx <- length(unique(melanomUS$tx$txSSid))
    # Proportions sum to 1 within stages
    melanomUS$tx <- transform(melanomUS$tx, 
                                 SSid=c(rep('Early.BRAF+', ntx), rep('Early.BRAF-', ntx), 
                                        rep('Advanced.BRAF+', ntx), rep('Advanced.BRAF-', ntx)), 
                                 txSSno=1:ntreat)
    melanomUS$tx <- transform(melanomUS$tx,
                              txHR=c(1, 1, 1, 1, 1, 1, 1,1,
                                     1, 1, 1, 1, 1, 1, 1,1,
                                     1, 1, 0.81, 0.55, 0.50, 0.69, 0.33,1,
                                     1, 1, 1, 1, 0.50, 0.69, 0.33,1))
    advA <- vector('list', length=length(unique(melanomUS$tx$SSid)))
    # Order of SSid: Early Advanced
    # Order of treatments: None, Interferon, Chemo, Checkpoint, Targeted
    advA[[1]] <- c(100, 0, 0, 0, 0, 0, 0, 0)
    advA[[2]] <- c(100, 0, 0, 0, 0, 0, 0, 0)
    advA[[3]] <- c(0, input$advchemostbp, 0,0,0,0,0,0)
    advA[[4]] <- c(0, input$advchemostbn, 0,0,0,0,0,0)
    # After treatment advances:
    advB<-advA
    advA[[3]] <- c(0,0,input$advbrafintbp, input$advbrafmekintbp,
                   input$advpd1intbp, input$advctla4intbp,input$advpd1ctla4intbp,input$advotherbp)
    advA[[4]] <- c(0,0,input$advbrafintbn, input$advbrafmekintbn,
                   input$advpd1intbn, input$advctla4intbn,input$advpd1ctla4intbn,input$advotherbn)
    # Put together complete vectors
    advA.v <- do.call('c', advA)
    advB.v <- do.call('c', advB)
    props <- data.frame(advB.v, advA.v, advB.v, advA.v)
    
    # colSums(props)  
    colnames(props) <- melanomUS$pol$id
    props2 <- props
    colnames(props2) <- melanomUS$pol$name
    toprint <- data.frame(melanomUS$tx, props2, stringsAsFactors=FALSE, check.names=FALSE)
    melanomUS$tx <- data.frame(melanomUS$tx, props, stringsAsFactors=FALSE)
    melanomUS$tx <- data.frame(melanomUS$tx, stringsAsFactors=FALSE)
    startclock <- proc.time()
    
   if (input$gender=="male")
   {
      withProgress(message = 'Simulation in Progress', {
        
          
          observeEvent(input$point, {
            removeUI(
              selector = "div:has(> #txt)"
            )
          })
        
        manuscript_melanoma_male <- simpolicies_mel(melanomUS$pol, melanomUS$nh, melanomUS$tx,overd,ppv,advanced1,advanced2,Stage1,
                                                    incsource='Male',
                                                    mortsource='Male', 
                                                    returnstats=c('mean','lower', 
                                                                  'upper'), 
                                                    popsize =input$sizepop,sims=input$nsim,
                                                    futimes=c(10,25),
                                                    minage=input$tfd[1], maxage=input$tfd[2],
                                                    denom=input$sizepop)
        # Runtime
        
        # Format and save results
        finaltab <- format_bounds_list(manuscript_melanoma_male, 
                                       paren=TRUE, includemean=TRUE, 
                                       digits=c(0,0,0,0,2,0,0,2),
                                       compileall=TRUE)
      })
   
   }
    
  else if (input$gender=="female")
  {
    withProgress(message = 'Simulation in Progress', {
      
        
        observeEvent(input$point, {
          removeUI(
            selector = "div:has(> #txt)"
          )
        })
      manuscript_melanoma_male <- simpolicies_mel(melanomUS$pol, melanomUS$nh, melanomUS$tx,overd,ppv,advanced1,advanced2,Stage1,
                                                  incsource='Female',
                                                  mortsource='Female', 
                                                  returnstats=c('mean','lower', 
                                                                'upper'), 
                                                  popsize =input$sizepop,sims=input$nsim,
                                                  futimes=c(10,25),
                                                  minage=input$tfd[1], maxage=input$tfd[2],
                                                  denom=input$sizepop)
      # Runtime
      
      # Format and save results
      finaltab <- format_bounds_list(manuscript_melanoma_male, 
                                     paren=TRUE, includemean=TRUE, 
                                     digits=c(0,0,0,0,2,0,0,2),
                                     compileall=TRUE)
    })
  } 
  })
  
  
  slidersimulation1 <- reactive({
    ppv<-5
    braf<-input$braf1/100
    overd<-input$overd
    risk<-input$riskstrat
    if (risk==5)
    {
      riskk<-5
      RR<-5.57
      }
   
    else if (risk==10)
    {
      riskk<-10
      RR<-4.08
     }
    
    else if (risk==15) 
    {
      riskk<-15
      RR<-3.26
      }
    
    else if (risk==20)
    {
      riskk<-20
      RR<-2.87
    }
    
   
   Stage1<-input$stage1/100
     advanced1<-input$advst/100
     advanced2<-input$advin/100
     
    surr<-input$survival
    if(surr==5){survv<-5} else if(surr==10){survv<-10}
    
    melanomUS <- vector('list', length=length(ex1))
    names(melanomUS) <- names(ex1)
    propERpos<-1
    pop_chars <- list(male = data.frame(male = c(0), prop = c(1)))

    melanomUS$pol  <- data.frame(num=c(1:4),
                                 id=c('adv.base', 
                                      'advbase.treatment','adv.shift', 'adv.shift.treatment'),
                                 name=c('M0: Historical Treatments Era without Screening',
                                        'M1: Novel Treatments Era with without Screening',
                                        'M2: Historical Treatments Era with Screening',
                                        'M3: Novel Treatments Era with Screening'),
                                 pairnum=c(NA, NA, c(1,2)),
                                 earlydetHR=c(rep(1, 2), rep(input$advin/input$advst, 2)),
                                 stringsAsFactors=FALSE)
    
    #Calculating baseline survival for early stages
    
    ff = function(x, p=propERpos, surg=100/100, hsurg=1, So=input$earlysur/100) { 
      p*(1-surg)*x + p*surg*x^hsurg- So}
    
    Sb.early <- uniroot(ff, lower=0.5, upper=input$earlysur/100, tol = 0.0001)
    early.mrate.co = cumsurv_to_exprate(Sb.early$root, year=survv)
    round(100*Sb.early$root)
    
    #Calculating baseline survival for advanced stages
    f = function(x, p=propERpos, chm=input$advchemostbp/100,
                 hchm=1, So=input$advsur/100) {p*(1-chm)*x + p*chm*x^hchm-So}
    
    Sb.late <- uniroot(f, lower=0.2, upper=input$advsur/100, tol = 0.0001)
    late.mrate.co = cumsurv_to_exprate(Sb.late$root, year=survv)
    round(100*Sb.late$root)
    
    #Natural history (nh) parameters
    melanomUS$nh<-compile_naturalhist(prop_adv = input$advst/100, mortrates=c(Early=early.mrate.co, Advanced=late.mrate.co),subgroup_probs=c(`BRAF+`=braf, `BRAF-`=1-braf))
    #Stage shift map
    melanomUS$map <- create_stageshift_map(melanomUS$nh)
    #Treatment
    melanomUS$tx <- data.frame(expand.grid(txSSid=c("Surgical Excision",
                                                    "Chemotherapy-Immunotherapy",
                                                    "BRAF Inhibitor",
                                                    "BRAF-MEK inhibitor Combination",
                                                    "Anti-PD1",
                                                    "Anti-CTLA-4",
                                                    "Anti-PD1-CTLA-4 Combination",
                                                    "Other(Temozolomide, Corboplatin)"),
                                           SSno=1:nrow(melanomUS$nh)),
                               stringsAsFactors=FALSE)
    ntreat <- nrow(melanomUS$tx)
    ntx <- length(unique(melanomUS$tx$txSSid))
    # Proportions sum to 1 within stages
    melanomUS$tx <- transform(melanomUS$tx, 
                              SSid=c(rep('Early.BRAF+', ntx), rep('Early.BRAF-', ntx), 
                                     rep('Advanced.BRAF+', ntx), rep('Advanced.BRAF-', ntx)), 
                              txSSno=1:ntreat)
    melanomUS$tx <- transform(melanomUS$tx,
                              txHR=c(1, 1, 1, 1, 1, 1, 1,1,
                                     1, 1, 1, 1, 1, 1, 1,1,
                                     1, 1, 0.81, 0.55, 0.50, 0.69, 0.33,1,
                                     1, 1, 1, 1, 0.50, 0.69, 0.33,1))
    advA <- vector('list', length=length(unique(melanomUS$tx$SSid)))
    # Order of SSid: Early Advanced
    # Order of treatments: None, Interferon, Chemo, Checkpoint, Targeted
    advA[[1]] <- c(100, 0, 0, 0, 0, 0, 0, 0)
    advA[[2]] <- c(100, 0, 0, 0, 0, 0, 0, 0)
    advA[[3]] <- c(0, input$advchemostbp, 0,0,0,0,0,0)
    advA[[4]] <- c(0, input$advchemostbn, 0,0,0,0,0,0)
    # After treatment advances:
    advB<-advA
    advA[[3]] <- c(0,0,input$advbrafintbp, input$advbrafmekintbp,
                   input$advpd1intbp, input$advctla4intbp,input$advpd1ctla4intbp,input$advotherbp)
    advA[[4]] <- c(0,0,input$advbrafintbn, input$advbrafmekintbn,
                   input$advpd1intbn, input$advctla4intbn,input$advpd1ctla4intbn,input$advotherbn)
    # Put together complete vectors
    advA.v <- do.call('c', advA)
    advB.v <- do.call('c', advB)
    props <- data.frame(advB.v, advA.v, advB.v, advA.v)
    
    # colSums(props)  
    colnames(props) <- melanomUS$pol$id
    props2 <- props
    colnames(props2) <- melanomUS$pol$name
    toprint <- data.frame(melanomUS$tx, props2, stringsAsFactors=FALSE, check.names=FALSE)
    melanomUS$tx <- data.frame(melanomUS$tx, props, stringsAsFactors=FALSE)
    melanomUS$tx <- data.frame(melanomUS$tx, stringsAsFactors=FALSE)
    startclock <- proc.time()
    spop<-input$sizepop
    spoprisk<-spop*riskk/100
   
    if (input$gender=="male")
    {
      
      withProgress(message = 'Simulation in Progress', {
        
        observeEvent(input$point_risk, {
          removeUI(
            selector = "div:has(> #txt)"
          )
        })
         manuscript_melanoma_male_risk <- simpolicies_mel_risk(melanomUS$pol, melanomUS$nh, melanomUS$tx,overd,ppv,advanced1,advanced2,Stage1,RR,
                                                         incsource='Male',
                                                         mortsource='Male', 
                                                         returnstats=c('mean','lower', 
                                                                       'upper'), 
                                                         popsize=spoprisk, sims=input$nsim,
                                                         futimes=c(10,25),
                                                         minage=input$tfd[1], maxage=input$tfd[2],
                                                         denom=spoprisk)
      })
      
      # Format and save results
      finaltab2 <- format_bounds_list(manuscript_melanoma_male_risk, 
                                      paren=TRUE, includemean=TRUE, 
                                      digits=c(0,0,0,1,0),
                                      compileall=TRUE)
      
    }
    
    else if (input$gender=="female")
    {
      withProgress(message = 'Simulation in Progress', {
        observeEvent(input$point_risk, {
          removeUI(
            selector = "div:has(> #txt)"
          )
        })
       
      
        manuscript_melanoma_female_risk <- simpolicies_mel_risk(melanomUS$pol, melanomUS$nh, melanomUS$tx,overd,ppv,advanced1,advanced2,Stage1,RR,
                                                              incsource='Female',
                                                              mortsource='Female', 
                                                              returnstats=c('mean','lower', 
                                                                            'upper'), 
                                                              popsize=spoprisk, sims=input$nsim,
                                                              futimes=c(10,25),
                                                              minage=input$tfd[1], maxage=input$tfd[2],
                                                              denom=spoprisk)
        
      })
      # Format and save results
      finaltab2 <- format_bounds_list(manuscript_melanoma_female_risk, 
                                      paren=TRUE, includemean=TRUE, 
                                      digits=c(0,0,0,1,0),
                                      compileall=TRUE)
      
    }
    
  })
  
  
  
    

 
    
  
  
  
  # Show the values in an HTML table ----
  output$values <- renderTable({
    sliderValues()
  })
  
  output$features<-renderTable({
    sliderFeatures()
  })
  
  output$treatmentadv<-renderTable({
    sliderTreatmentsAdv()
  })
  
  output$treatmentear<-renderTable({
    sliderTreatmentsEar()
  })
  
  output$efficacy<- renderTable({
    sliderEfficacy()
  })
  
  output$papers<- renderTable({
    sliderPapers()
  })
  
  output$policy<-renderTable({
    slidersimulation()
  })
  
  output$policy_risk<-renderTable({
    slidersimulation1()
  })
}
