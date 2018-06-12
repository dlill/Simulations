library(shiny)
library(miniUI)
library(ggplot2)

ropt_gadget <- function(which_pars_opt) {
  
  ui <- miniPage(
    gadgetTitleBar("Explore the algorithm's behaviour depending on complexes and alphas"),
    miniContentPanel(
      # The brush="brush" argument means we can listen for
      # brush events on the plot using input$brush.
      tableOutput("r0"),
      tableOutput("rkept"),
      tableOutput("ropt"),
      selectInput("whichParOpt", "which pars to compare", which_pars_opt, selected = which_pars_opt, multiple = T),
      sliderInput("value", "value", -4,4, -10, step = 0.001),
      actionButton("save", "Save")
    )
  )
  
  server <- function(input, output, session) {
    r0 <- R_fun(pars_opt = pars_opt,
          perturbation_prediction = perturbation_prediction,
          obs_fun = g,
          p_fun = (p*p_pert),
          pars = pars) %>% local_response_matrix_eq10()
    output$r0 <- renderTable(r0)
    
    rkept <- reactive(r_kept_fun(pars_opt = pars_opt[input$whichParOpt],
                                        perturbation_prediction = perturbation_prediction, 
                                        obs_fun = g,
                                        p_fun = (p * p_pert), 
                                        pars = pars,
                                        alpha = 
                                   -log(.Machine$double.eps) + input$value))
    output$rkept <- renderTable(rkept() %>% {mymat <- .;mymat[!mymat] <- mymat[!mymat] %>% str_replace("FALSE", "0"); mymat})
    
    ropt <- reactive({
      myfits <- mstrust(obj_alpha, 
                        center =  structure(rep(0, length(input$whichParOpt)), names = input$whichParOpt), 
                        studyname = "Fits", 
                        cores = 3, 
                        fits = 3, 
                        sd = 1, 
                        mypars = pars, 
                        perturbation_prediction = perturbation_prediction, 
                        r_kept = rkept(), 
                        p_fun = (p * p_pert),
                        obs_fun = g)
      
      try(myfits %>% as.parframe()) %>% 
      {e <- .; if (inherits(e, "try-error")) print(myfits)}
      
      
      out <- r_alpha_fun(pars_opt = myfits %>% as.parframe() %>% as.parvec() %>% unclass() , 
                  pars = pars, 
                  perturbation_prediction = perturbation_prediction,
                  p_fun = (p*p_pert))
      attr(out, "pars_opt") <- myfits %>% as.parframe() %>% as.parvec() %>% unclass()
    
      out
    })
    
    output$ropt <- renderTable(ropt())
    
    observeEvent(input$save, {
      tpaste0 <- function(...) {
        paste0(format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
               , "_" , ...)
      }
      filename <- paste0("~/Promotion/Projects/MRA/Simulations/Models/CascadeComprehensiveFeedbacks/Exploration/", tpaste0("result"))
      out <- rbind(r0 %>% round(2), rkept() %>% round(2),ropt() %>% round(2)) %>% as.data.frame %>% 
        rbind(0) %>% rbind(0) %>% rbind(0)
      pars_opt <- ropt() %>% attr("pars_opt") %>% {bla <- .; attr(bla, "alpha") <- input$value; bla} %>% {bla <- .; c(bla, rep(0, nrow(out)-length(bla)))}
      write.csv(cbind(out, pars_opt, names(pars_opt)), file = paste0(filename, ".csv"))
      saveRDS(ropt() %>% attr("pars_opt") %>% {bla <- .; attr(bla, "alpha") <- input$value; bla}, paste0(filename, "_pars.rds"))
    })
  
    observeEvent(input$done, stopApp(NULL))  
  }
  
  runGadget(ui, server)
}
