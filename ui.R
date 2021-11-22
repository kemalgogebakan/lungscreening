#-------------------------------------------------------------------------
#  This application is governed by the CeCILL-B license. 
#  You can  use, modify and/ or redistribute this code under the terms
#  of the CeCILL license:  http://www.cecill.info/index.en.html
#
#  Marc Lavielle, Inria Saclay
#  May 11th, 2015
#-------------------------------------------------------------------------

library(shinydashboard)
library(bcimodel)
load("ex1.rda")
load("incidence_mel.rda")
load("othercausemort.rda")
source("compile_naturalhistmelanoma.R")
source("simpolicies_melanoma.R")
source("databases_melanoma.R")
source("initialize_melanoma.R")
source("outcomes_melanoma.R")
source("cantrance_melanoma.R")
source("earlydetect_melanoma.R")
source("general_melanoma.R")
source("initialize_melanoma.R")
source("outcomes_melanoma.R")
source("parallelized_melanoma.R")
source("screentreatlibrary_melanoma.R")
source("systime_melanoma.R")
source("treatment_melanoma.R")




sidebar <- dashboardSidebar(
  hr(),
  width=450,
  sidebarMenu(id="tabs",
              menuItem("Introduction", tabName = "introduction"),
              
              menuItem("Inputs to the Microsimulation Model", tabName ="inputs_model"),
              
              menuItem("Population Characteristics", tabName = "disease", 
                       menuSubItem("Age and Gender", tabName = "age", icon = icon("angle-right"))),
              
              menuItem("Disease Characteristics", tabName = "disease2",
                       menuSubItem("BRAF Status", tabName = "braf", icon = icon("angle-right"))),
              
              menuItem("Stage Distribution Before and After Screening", tabName = "stage"),
              
              menuItem("Treatment", tabName = "treatment",
                       menuSubItem("Distribution in Historical Treatment Period", tabName = "stdcare", icon = icon("angle-right")),
                       menuSubItem("Distribution in Novel Treatment", tabName = "intervene", icon = icon("angle-right")),
                       menuSubItem("Efficacy of Each Treatment", tabName = "efficacy", icon = icon("angle-right"))
              ),
              
              menuItem("Mortality", tabName = "mortality",
                       menuSubItem("Melanoma Cancer Survival", tabName = "survival", icon = icon("angle-right")),
                       menuSubItem("Other Cause Mortality", tabName = "other", icon = icon("angle-right"))
              ),
              
              menuItem("Overdiagnosis Rate", tabName = "overdiag_ppv"),
              
              menuItem("Cohort Size and Number of Simulations", tabName = "size"),
  
              menuItem("Incidence and Mortality Outcomes", tabName = "result",
                       menuSubItem("Review inputs", tabName = "review", icon = icon("angle-right")),
                       menuSubItem("Simulation Results", tabName = "point", icon = icon("angle-right"))
                       ),
  
              menuItem("Risk Stratified Screening", tabName = "riskstrat",
                       menuSubItem("Risk Group", tabName = "review_risk", icon = icon("angle-right")),
                       menuSubItem("Simulation Results of Risk Stratified Screening", tabName = "point_risk", icon = icon("angle-right"))
              )
              
          )
)

body <- dashboardBody(
    tabItems(
      tabItem(tabName = "introduction",
              box(width = NULL, status = "primary", solidHeader = TRUE, title="Overview:",
              h5("This interface allows you to model the survival benefit of a risk stratified melanoma 
              screening and/or treatment intervention in a virtual population of your choosing.
              The model posits a simple stage-shift mechanism of screening benefit. 
              For more information about stage-shift models, visit:"),
              tags$a(href="https://cancerpolicy.shinyapps.io/breastcancer/", "https://cancerpolicy.shinyapps.io/breastcancer/"))),
      
      tabItem(tabName = "inputs_model",
              box(width = NULL, status = "primary", solidHeader = TRUE, title="1.Incidence:",
              h5(strong("Clinical Melanoma Incidence by Age, for Year Of Diagnosis 2011-2015")),
              h5(strong("Source:"),"Surveillance, Epidemiology, and End Results (SEER) Program (www.seer.cancer.gov) 
                 SEER*Stat Database: Incidence - SEER 18 Regs Research Data, Nov 2018 Sub (1975-2016), 
                 Linked To County Attributes - Total U.S., 1969-2017 Counties, National Cancer Institute, 
                 DCCPS, Surveillance Research Program, released April 2019, based on the November 2018 submission.")),
              
              box(width = NULL, status = "primary", solidHeader = TRUE, title="2.Default Advanced Stage and Stage I Fraction of Population:",
                h5(strong("Proportion of cases that are advanced stage (stage III and IV) at clinical diagnosis in 2011-2015 by gender:")),
                h5(strong("Source:"),"Surveillance, Epidemiology, and End Results (SEER) Program (www.seer.cancer.gov) 
                 SEER*Stat Database: Incidence - SEER 18 Regs Research Data, 
                 Nov 2018 Sub (1975-2016), 
                 Linked To County Attributes - Total U.S., 1969-2017 Counties, National Cancer Institute, 
                 DCCPS, Surveillance Research Program, released April 2019, based on the November 2018 submission.")),
              
              box(width = NULL, status = "primary", solidHeader = TRUE, title="3.Default and Improved Treatment Distribution of Population:",
                  h5("Model assumes that all patients receive surgical therapy."),
                  h5("For other treatments and proportions see Treatment Section.")),
              
              box(width = NULL, status = "primary", solidHeader = TRUE, title="4.Treatment Efficacies:",
                  h5("See Treatment Section.")),
              
              box(width = NULL, status = "primary", solidHeader = TRUE, title="5.Default Disease Specific Survival Rates by Stage and Gender at diagnosis between 2006-2010:",
                  h5(strong("5 year survival rates:")),
                  h5("For males: Early:95.6%, Advanced:44.8%"),
                  h5("For females: Early:96.8%, Advanced:53.4%"),
                  h5(strong("10 year survival rates:")),
                  h5("For males: Early:92.4%, Advanced:41.8%"),
                  h5("For females: Early:96.1%, Advanced:53.5%"),
                  h5(strong("Source:"),"Surveillance, Epidemiology, and End Results (SEER) Program (www.seer.cancer.gov) 
                 SEER*Stat Database: Incidence - SEER 9 Regs Research Data, Nov 2018 Sub (1975-2016), 
                 Linked To County Attributes - Total U.S., 1969-2017 Counties, National Cancer Institute, 
                 DCCPS, Surveillance Research Program, released April 2019, based on the November 2018 submission.")),
                          
              box(width = NULL, status = "primary", solidHeader = TRUE, title="6.Other cause mortality rates by age:",
                  h5(strong("Source:"),"Human Mortality Database, for Year 2015")),
              
              box(width = NULL, status = "primary", solidHeader = TRUE, title="7.Default Stage Shift From Advanced to Early Stage:",
                  h5("In the baseline analysis, screening reduces the incidence of presentation of advanced stage disease by 50%")),
              
              box(width = NULL, status = "primary", solidHeader = TRUE, title="8.Fraction of Advanced Melanoma Patients with a BRAF V600 Mutaion:",
                  h5("This fraction is 45%, which is assumed to be same for male and female patients"),
                  h5("Reference: Whitman, E. D., Liu, F. X., Cao, X., Diede, S. J., Haiderali, A., & Abernethy, A. P. (2018). Treatment patterns and outcomes for patients with advanced melanoma in US oncology clinical practices. Future Oncology, 15(5), 459-471.")),
              
              box(width = NULL, status = "primary", solidHeader = TRUE, title="6.Overdiagnosis rate:",
                  h5("Enter fraction of Stage 1 patients and overdiagnosis rate which is the rate of increase in Stage 1 patients after screening")),
              
              box(width = NULL, status = "primary", solidHeader = TRUE, title="7.Risk level to be screened:",
                  h5("Select risk group to be screened"),
                  h5("Relative risks of high risk patients are taken from from Williams and colleagues’ 
                     2011 Washington State melanoma case-control study where self-assessed risk factors were collected"),
                  h5("Reference: Williams, L.H., et al., Identifying Persons at Highest Risk of Melanoma Using Self-Assessed Risk Factors. 
                     J Clin Exp Dermatol Res, 2011. 2(6).")),
             
              box(width = NULL, status = "primary", solidHeader = TRUE, title="For Model Customization:",
              h5("Specify inputs using the left panel to navigate."))
              ),
    
             
      
      tabItem(tabName = "age", box(width = NULL, status = "primary", solidHeader = TRUE, title="Define the age range and gender of interest",
                              h5("The model will track outcomes for a cohort of of these ages and gender")),
                        box(width = NULL, status = "primary", solidHeader = TRUE, title="Select Gender:",          
                            radioButtons(inputId="gender", label = "Gender",
                                         c("Male" = "male","Female" = "female"), selected="male")
                        ),
                        box(width = NULL, status = "primary", solidHeader = TRUE, title="Define The Age Range:",         
                            sliderInput("tfd","", value=c(35,75), min=0, max = 100, step=1)
                            )
                        ),
      
      tabItem(tabName="stage",
              box(width = NULL, status = "primary", solidHeader = TRUE, title="Enter the percentage of advanced stage cases under two scenarios:",
                      h5("For each of the standard-of-care and intervention scenarios, this is the percent of cases who are advanced-stage at the time of clinical diagnosis. If the intervention is expected to detect cases early, the percent of cases diagnosed in advanced stage should be lower.")
                  ),
              box(width = NULL, status = "primary", solidHeader = TRUE, title="Percentage of advanced stage melanoma, before screening:",         
                  sliderInput("advst","", value=12.8, min=0, max = 100, step=0.1)
              ),
              box(width = NULL, status = "primary", solidHeader = TRUE, title="Percentage of advanced stage melanoma, after screening:",         
                  sliderInput("advin","", value=6.4, min=0, max = 100, step=0.1)),
                  
              
              box(width = NULL, status = "primary", solidHeader = TRUE, title="Summary of specified stage distributions:", 
            
                tableOutput("values")),
                
              
              ),
      
      
      tabItem(tabName="braf",
              box(width = NULL, status = "primary", solidHeader = TRUE, title="Select the  percentage of advanced stage patients who are BRAF positive"),
              box(width = NULL, status = "primary", solidHeader = TRUE, title="BRAF positive percentage",         
                  sliderInput("braf1","BRAF Positive", value=45, min=0, max = 100, step=1))
            
      ),
    
    tabItem(tabName="stdcare",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Select the treatment distributions under standard of care for early and advanced stages of melanoma",
                h5("The proportion for each melanoma treatment should be entered as summing up to 100.")),
            
            
            box(width = 4, status = "primary", solidHeader = TRUE, title="BRAF-mutant Type Advanced Melanoma",         
                sliderInput("advchemostbp", "Chemotherapy-Immunotherapy:", value=100, min=0, max = 100, step=1)
            ),
            
            box(width = 4, status = "primary",solidHeader = TRUE, title="BRAF-wild Type Advanced Melanoma",        
                sliderInput("advchemostbn", "Chemotherapy-Immunotherapy:", value=100, min=0, max = 100, step=1)
            ),
            
            box(width = 4, status = "primary",solidHeader = TRUE, title="Early Melanoma",        
                sliderInput("earlyst", "Surgical Excision:", value=100, min=0, max = 100, step=1))
            
    ),
    
    tabItem(tabName="intervene",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Select the treatment distributions under advances in care for early and advanced stages of melanoma",
                h5("The proportion for each melanoma treatment should be entered as summing up to 100.")),
            
            box(width = 4, status = "primary", solidHeader = TRUE, title="BRAF-mutant Type Advanced Melanoma",         
                sliderInput("advbrafintbp", "BRAF Inhibitor Therapy:", value=5, min=0, max = 100, step=1),
                sliderInput("advbrafmekintbp", "BRAF-MEK Inhibitor Therapy:", value=31, min=0, max = 100, step=1),
                sliderInput("advpd1intbp", "Anti-PD1 Therapy:", value=30, min=0, max = 100, step=1),
                sliderInput("advctla4intbp", "Anti-CTLA4 Therapy:", value=6, min=0, max = 100, step=1),
                sliderInput("advpd1ctla4intbp", "Anti PD1-CTLA-4 Therapy:", value=25, min=0, max = 100, step=1),
                sliderInput("advotherbp", "Other(Temozolomide, Corboplatin):", value=3, min=0, max = 100, step=1)
            ),
            
            box(width = 4, status = "primary", solidHeader = TRUE, title="BRAF-wild Type Advanced Melanoma",         
                sliderInput("advbrafintbn", "BRAF Inhibitor Therapy:", value=0, min=0, max = 100, step=1),
                sliderInput("advbrafmekintbn", "BRAF-MEK Inhibitor Therapy:", value=0, min=0, max = 100, step=1),
                sliderInput("advpd1intbn", "Anti-PD1 Therapy:", value=58, min=0, max = 100, step=1),
                sliderInput("advctla4intbn", "Anti-CTLA4 Therapy:", value=11, min=0, max = 100, step=1),
                sliderInput("advpd1ctla4intbn", "Anti PD1-CTLA-4 Therapy:", value=27, min=0, max = 100, step=1),
                sliderInput("advotherbn", "Other(Temozolomide, Corboplatin):", value=4, min=0, max = 100, step=1)
            ),
            
            box(width = 4, status = "primary",solidHeader = TRUE, title="Early Melanoma",        
                sliderInput("earlyint", "Surgical Excision:", value=100, min=0, max = 100, step=1))
            
    ),
    
    tabItem(tabName = "efficacy",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Efficacy of treatment is obtained by literature review",
            h5("Hazard ratios of treatments of melanoma are listed according to the results of randomized controlled trials.")
            ),
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Hazard Ratios", 
                tableOutput("efficacy")
            ),
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Research articles presenting efficacy values of treatments:",
                tableOutput("papers")
            )
             
    ),
            
      
    
    tabItem(tabName="survival",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Specify the percent of cases surviving at k years after diagnosis (i.e. baseline survival)",
            selectInput("survival", label = h3("Year of survival statistic"), 
                        choices = c(5, 10), 
                        selected = 5)),
            
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Advanced cases:baseline survival at k years",         
                sliderInput("advsur","", value=44.8, min=0, max = 100, step=0.1)),
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Early cases:baseline survival at k years",         
                sliderInput("earlysur","", value=95.6, min=0, max = 100, step=0.1))
    ),
    
    tabItem(tabName="other",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Other cause mortality",
                h5(strong("Source:"),"Human Mortality Database, for Year 2015"),
                tags$a(href="https://www.mortality.org/", "https://www.mortality.org/"))
    ),
    
   
    
    tabItem(tabName="overdiag_ppv",
            
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Select the percentage of Stage 1 patients before screening and 
                expected  rate of increase in Stage 1 patients (i.e. Overdiagnosis Rate) after the Screening Program"),
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Fraction of Stage 1 Patients",         
                sliderInput("stage1", "Stage1 Fraction", value=73.1, min=0, max = 100, step=1)),
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Overdiagnosis Rate",         
                sliderInput("overd", "Overdiagnosis Rate", value=10, min=0, max = 100, step=1))
            ),
   

  
    tabItem(tabName="size",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Select the size of population and number of simulations"),
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Size of population",         
                sliderInput("sizepop","", value=100000, min=50000, max = 1000000, step=5000)),
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Number of simulations",         
                sliderInput("nsim", "", value=50, min=0, max = 300, step=5))
            ),
    
    tabItem(tabName="review",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Confirmation page of selected parameters for simulation",
            h5("If you want to make changes, please revisit the previous pages for model reparametrization")),
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Population Features:", 
                tableOutput("features")),
            
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Treatment distributions, BRAF-mutant Advanced stage:", 
                tableOutput("treatmentsadvbrafp")),
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Treatment distributions, BRAF-wild Advanced stage:", 
                tableOutput("treatmentsadvbrafn"))
            
    ),
    
    tabItem(tabName="point",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Simulation Results",

            
            tableOutput("policy"))
    ),
    
    tabItem(tabName="review_risk",
            
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Select top x% risk group of the population to be screened",
            selectInput("riskstrat", label = h3("Risk level (%)"), 
                        choices = c(5, 10, 15, 20), 
                        selected = 20))
            ),
    
    tabItem(tabName="point_risk",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Simulation Results of Risk Stratified Screening",
                
                
                tableOutput("policy_risk")))

  )
    
)




dashboardPage(
  dashboardHeader(title = "Melanoma Microsimulation Model of Early Detection and Treatment", titleWidth = 750),
  sidebar,
  body
)