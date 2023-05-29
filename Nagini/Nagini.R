#####################################################
#Title: app per analisi differenziale di reti       #
#Authors: Francesco Chiodo, Pietro Hiram Guzzi      #
#Version: 1.0                                       #
#Last mod: 04/02/2023 10:00                         #
#####################################################

#librerie

library(shiny)
library(shinythemes)
library(xml2)
library(XML)
library(modeldata)
library(DataExplorer)
library(plotly)
library(tidyverse)
library(igraph)
library(glmnet)
library(bslib)
library(readxl)

#####################################
#Funzioni fornite                   #
#####################################

#odata=read_xlsx('/Users/pietrohiramguzzi/Downloads/Data.xlsx',col_types="numeric")

#dataf=data.frame(data)

#Label=read_xlsx('/Users/pietrohiramguzzi/Downloads/Label.xlsx')

#labelframe=data.frame(Label)




chNetModified <- function(X, group, subsampling, R, lambar, parallel = FALSE, nCpus = 4){
  
  
  p = dim(X)[2]
  
  #ind1 = which(group == unique(group)[1])
  #len1=length()
  #ind2 = which(group == unique(group)[2])
  #len2=length(ind2)
  #finqui
  #####subsampling experiment
  res.weight = matrix(1,p,p)
  print('start')
  
  
  
  cc = compute.statis(X, group, parallel, nCpus)
  print('end statis')
  lams.est <- get.knots(cc$tmu, cc$trho)
  lams.edge.temp = lams.est$int.knots
  
  diag(lams.edge.temp)=0
  
  lams.edge  <- (lams.edge.temp >= lambar)*1
  diff.edge.weight = lams.edge * res.weight
  
  lams.edge = (diff.edge.weight >0)*1
  
  lams.diffgene  <- (lams.est$main.knots >= lambar)*1
  
  
  
  Diff.edge.temp = upper.tri(lams.edge)*lams.edge
  Diff.edge=Diff.edge.temp+t(Diff.edge.temp)
  row.names(Diff.edge) = colnames(X)
  colnames(Diff.edge) = colnames(X)
  Degree.Diff = apply(Diff.edge, 1, sum)
  Diff.graph.connected = Diff.edge[Degree.Diff>0, Degree.Diff>0]
  Diff.net.connected  = graph_from_adjacency_matrix(Diff.graph.connected, mode = "undirected", weighted = TRUE, diag = FALSE)
  Diff.net = graph_from_adjacency_matrix(Diff.edge, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  return(list( diff.edge.weight = diff.edge.weight,diff.edge=lams.edge, diff.gene=lams.diffgene,
               Diff.net=Diff.net,Diff.net.connected=Diff.net.connected))
}



lasso.estimator <- function(X, lam1,  parallel = FALSE, nCpus = 4){
  
  print('start lasso')
  n = dim(X)[1]
  p = dim(X)[2]
  
  centered.X = scale(X, scale = FALSE)
  Sigma = cov(centered.X)
  DSigma = diag(Sigma)
  
  lam2 = sqrt(DSigma*log(p)/n)
  beta = array(0, dim=c(p,p,length(lam1)))
  res = array(0,dim=c(n,p,length(lam1)))
  
  wrapper = function(i){
    print('wrapper')
    #fit=glmnet(centered.X[,-i], centered.X[,i], lambda= lam1*lam2[i])
    fit=glmnet(centered.X[,-i], centered.X[,i], lambda= 0.0001)
    print('end glmnet')
    fit$beta=as.matrix(fit$beta)
    if(ncol(fit$beta)<length(lam1)){
      tmp = matrix(0,nrow = nrow(fit$beta),ncol = length(lam1))
      tmp[,1:ncol(fit$beta)]=fit$beta
      tmp[,ncol(fit$beta):length(lam1)] = fit$beta[,ncol(fit$beta)]
      fit$beta = tmp
    }
    
    if(i==1){
      beta[2:p,i,]=fit$beta
    }else if(i==p){
      beta[1:(p-1),i,]=fit$beta
    }else{
      beta[1:(i-1),i,]=fit$beta[1:(i-1),]
      beta[(i+1):p,i,]=fit$beta[i:nrow(fit$beta),]
    }
    res[,i,] = matrix(rep(centered.X[,i],length(lam1)),ncol = length(lam1)) - centered.X[,-i]%*%fit$beta
    out = list(beta = beta[,i,], res = res[,i,])
    return(out)
  }
  
  
  if(parallel){
    fit =mclapply(1:p, wrapper, mc.cores=nCpus)
  }else{
    fit =lapply(1:p, wrapper)
  }
  
  
  for(i in 1:p){
    beta[,i,]=fit[[i]]$beta
    res[,i,]=fit[[i]]$res
  }
  
  
  r.tilde = array(0, dim=c(p,p,length(lam1)))
  r.hat = array(0, dim=c(p,p,length(lam1)))
  
  for(k in seq(length(lam1))){
    r.tilde[,,k] = cov(res[,,k])*(n-1)/(n)
    r.hat[,,k] =r.tilde[,,k] + diag(diag(r.tilde[,,k]))%*%beta[,,k] + t(beta[,,k])%*%diag(diag(r.tilde[,,k]))
  }
  
  out = list(beta = beta, res = res, r.tilde = r.tilde, r.hat = r.hat)
  return(out)
}

## calculate test statistics
compute.statis<- function(X, group, parallel = FALSE, nCpus = 4 ) {
  print('Start Computestatis')
  n <- nrow(X)
  p <- ncol(X)
  
  lam1 = seq(40,4)/20
  nlam1=length(lam1)
  
  stopifnot(length(group) == n)
  
  stopifnot(table(group) > 3) # need at least 2 observations in each class
  
  X1 = X[group==unique(group)[1],]
  X2 = X[group==unique(group)[2],]
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  
  #step 1  calculate the t-statistical of mu
  print('step 1 compute')
  xbar1 <- colMeans( X1)
  xbar2 <- colMeans( X2)
  s1 <- apply( X1, 2, sd)
  s2 <- apply( X2, 2, sd)
  
  tmu <- (xbar1-xbar2) / sqrt(s1^2/n1+ s2^2/n2)
  
  
  #step 2  calculate the t-statistical of rho
  print('step 2 cmpute')
  fit1 = lasso.estimator(X1, lam1, parallel, nCpus)
  fit2 = lasso.estimator(X2, lam1, parallel, nCpus)
  
  WW = array(0, dim=c(p,p,nlam1))
  score = rep(0,nlam1)
  for(k in seq(nlam1)){
    # print(k)
    Dr1 = diag((diag(fit1$r.tilde[,,k]))^(-0.5))
    Dr2 = diag((diag(fit2$r.tilde[,,k]))^(-0.5))
    T1 = Dr1%*%fit1$r.hat[,,k]%*%Dr1
    T2 = Dr2%*%fit2$r.hat[,,k]%*%Dr2
    rho1 = T1*((abs(T1)>=2*sqrt(log(p)/n1))*1)
    rho2 = T2*((abs(T2)>=2*sqrt(log(p)/n2))*1)
    WW[,,k] = (T1-T2)/sqrt((((1-rho1^2)^2)/n1)+(((1-rho2^2)^2)/n2))
    diag(WW[,,k] ) = 0
    WW[,,k] = WW[,,k]*upper.tri(WW[,,k])
    for (l in seq(1,10)){
      score[k] = score[k] + (sum(abs(WW[,,k])>qnorm(1-l*(1-pnorm(sqrt(log(p))))/10))/(l*(p^2-p)*(1-pnorm(sqrt(log(p))))/10)-1)^2
    }
  }
  
  kmin = which(score==min(score))
  
  if(length(kmin)>1){
    k.min=kmin[1]
  }else{k.min=kmin}
  
  W.opt = WW[,,k.min]
  
  W.opt =as.matrix(  W.opt )
  diag(W.opt) = 0
  
  cc= list(tmu=tmu, trho=W.opt)
  
}

# A closed-formed solution is derived to solve the optimization model.
get.knots <- function(w, z) {
  
  p <- length(w)
  #cat("p=",p)
  #cat("nrow(z) =",nrow(z) )
  if (nrow(z) != p) browser()
  stopifnot(nrow(z) == p, ncol(z) == p)
  aw <- abs(w)
  diag(z) <- rep(0, p)
  az <- abs(z)
  
  main.knots <- pmax(aw, (aw + apply(az, 1, max)) / 2)
  int.knots <- matrix(NA, p, p)
  for (j in seq(p)) {
    for (k in seq(p)) {
      if (j == k) next
      int.knots[j, k] <- max(aw[j] - sum(pmax(az[j, -j] - az[j, k], 0)), 0) / 2
    }
  }
  int.knots <- pmin(az, az / 2 + int.knots)
  int.knots <- pmax(int.knots, t(int.knots)) # use max(lam_jk, lam_kj)
  
  return(list(main.knots=main.knots, int.knots=int.knots))
  
}

########################
#Definizione Client    #
########################


# user interface
ui <- fluidPage(theme = shinytheme("cyborg"),
  
  navbarPage("CHnet_myAPP", id = "inTabset",
             
             #################################
             #Prima pagina: seleziona azione #           
             #################################
             
             tabPanel(title = "Selection what to do ", value = "panel1",
                      
                      #####################################
                      #Prima pagina: selezione funzione   #
                      #####################################
                      
                      # Panel per titolo
                      titlePanel("Select a function: "),
                      
                      #Pannello principale
                      mainPanel(
                        
                        #Bottone per richiedere il servizio
                        actionButton('jumpToP2', "CHnet", class = "btn-primary", width = 500),
                        
                        # Divisore
                        tags$hr(),
                        
                        #Bottone per richiedere il servizio
                        actionButton("act1", "servizio2", width = 500),
                        
                        # Divisore
                        tags$hr(),
                        
                        #Bottone per richiedere il servizio
                        actionButton("act2", "servizio3",  width = 500)
                        
                        
                      )
                      
                      ),
             
             #####################################
             #Seconda pagina: carica files e var #
             #####################################
             
             
             tabPanel(title = "Upload files ",  value = "panel2",
                      # Panel per titolo
                      titlePanel("Uploading Files"),
                      
                      # Sidebar senza definizione in e out
                      sidebarLayout(
                        
                        # Sidebar panel 
                        sidebarPanel(
                          
                          # Selezione file
                          fileInput("file1", "Choose first XLS file"),
                          
                          # Divisore
                          tags$hr(),
                          
                          
                          # Selezione file
                          fileInput("file2", "Choose second XLS file"),
                          
                          # Divisore
                          tags$hr(),
                          
                          # Divisore
                          tags$hr(),
                          
                          # Numero di righe da mostrare
                          radioButtons("varBool", "Subsampling",
                                       choices = c("True" = 1,
                                                   "False" = 0),
                                       selected = "head"),
                                       fluidRow(column(3, verbatimTextOutput("varBool"))),
                          
                          #Input per valore numerico R (1- +inf.)
                          numericInput(inputId = "obs",
                                       label = "Number of observations to view (R param.):",
                                       value = 10),
                          
                          #Input per valore numerico lambda (2- +inf.)
                          numericInput(inputId = "obs2",
                                       label = "Lambda value:",
                                       value = 0.1),
                          
                          
                          #Bottone per richiedere il servizio
                          actionButton("jumpToP3", "Procedi col caricamento", class = "btn-primary")
                          
                          
                        ),
                        
                        # Pannello dove viene mostrato l'output
                        mainPanel(
                          
                          # Output
                          tableOutput("contents"),
                          
                          tableOutput("Data1"),
                          
                          tableOutput("Data2"),
                          
                          tableOutput("result"),
                          
                        )
                        
                      )
                      
                      
                      
                      
                      
                      
                      ),
             
             
             tabPanel(title = "Plot results ", value = "panel3"),
             
             mainPanel(
               
               tableOutput("plotRes")
             )
             
  ),
  
)
  
  
  



########################
#Definizione Server    #
########################




# Definizione server
server <- function(input, output, session) {
  
  #Observer per il tasto CHnet
  observeEvent(input$jumpToP2, {
        updateTabsetPanel(session, "inTabset",
            selected = "panel2")
  })
  
  #Observer per il tasto submitData
  observeEvent(input$jumpToP3, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel3")
  })
  
  output$contents <- renderTable({
    
    req(input$file1)
    
    # Valuta errori sulle virgole dei CSV
    tryCatch(
      {
        
        #Primo file
        df1 <- eventReactive(input$file1, {
          dfdata = data.frame(read_excel(input$file1$datapath))
        })
        output$Data1 <- renderTable( head( df1()))
        
        #Secondo file
        df2 <- eventReactive(input$file2, {
          dfdata = data.frame(read_excel(input$file2$datapath))
        })
        output$Data2 <- renderTable( head( df2()))
        
        #Valore di subsampling
        output$boolChose <- renderPrint({ input$varBool })
        
        #Valore di R
        output$valueR <- renderText({ input$obs })
        
        #Valore di lambda
        output$valueLambda <- renderText({ input$obs2 })
        
        #Runno CHNetModified
        output$result <- chNetModified(df1, df2, boolChose, valueR, valueLambda)
        
        output$plotRes <- plot(result$Diff.net)
        
      },
      error = function(e) {
        # Stop in caso di errore
        stop(safeError(e))
      }
    )
    # TODO if size(df1) != size(df2) => break;
    
  })
  
}

#Creazione app runnable
shinyApp(ui = ui, server = server)