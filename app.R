library(shiny)

#### Front-end 

ui <- fluidPage(
    #theme = "bootswatch-cerulean.css",
    titlePanel("Modelo de Solow"),
    h2("Emmanuel Maruri, Diego López."),
    h3(" Colegio de México."),
    h4("Simulación de las trayectorias de las variables endógenas y simulación de shock exógeno."),
    h4("Puedes ver el código utilizado al final de la página."),
tabsetPanel(
    tabPanel("Trayectorias",
             ### Sidebar layout 1
             sidebarLayout(position = "right",
                           sidebarPanel(
                               tags$img(height=50,width=50,src="colmex.jpg"),
                               tags$hr(),
                               tags$p("Visita nuestros sitios"),
                               tags$hr(),
                               tags$a(href="https://e-maruri.github.io/index.html", "Emmanuel Maruri  "),
                               tags$br(),
                               tags$a(href="https://diego-eco.github.io/", "  Diego López"),
                               h2("Parámetros del modelo"),
                               sliderInput(inputId = "tec", 
                                           label = "Tecnología (A)", 
                                           value = 1, min = 0.5, max = 10),
                               sliderInput(inputId = "k", 
                                           label = "Capital Inicial (k_0)", 
                                           value = 4, min = 0.9, max = 30),
                               sliderInput(inputId = "alpha", 
                                           label = "alpha", 
                                           value = 0.5, min = 0.1, max = 0.9),
                               sliderInput(inputId = "s", 
                                           label = "Ahorro (s)", 
                                           value = 0.3, min = 0.1, max = 0.9),
                               sliderInput(inputId = "delta", 
                                           label = "Tasa de descuento (d)", 
                                           value = 0.1, min = 0.1, max = 0.9)
                           ),
                           mainPanel(
                               h2("Gráficos interactivos"),
                               h4("Función de producción Y=Ak^(alpha)"),
                               splitLayout(plotOutput("primera"),
                                           plotOutput("segunda")),
                               splitLayout(plotOutput("tercera"),
                                           plotOutput("cuarta"))
                               
                           )
             ) ### End Sidebar layout 1
             ), ### End tabPanel Trayectorias
    tabPanel("Shocks",
             ### Sidebar layout 2
             sidebarLayout(position = "right",
                           sidebarPanel(
                               h4("Parámetros de Shock Exógeno"),
                               sliderInput(inputId = "tec_1", 
                                           label = "Tecnología (A)", 
                                           value = 1, min = 0.5, max = 10),
                               sliderInput(inputId = "k_1", 
                                           label = "Capital Inicial (k_0)", 
                                           value = 9, min = 0.9, max = 30),
                               sliderInput(inputId = "alpha_1", 
                                           label = "alpha", 
                                           value = 0.5, min = 0.1, max = 0.9),
                               sliderInput(inputId = "s_1", 
                                           label = "Ahorro (s)", 
                                           value = 0.3, min = 0.1, max = 0.9),
                               sliderInput(inputId = "delta_1", 
                                           label = "Tasa de descuento (d)", 
                                           value = 0.1, min = 0.1, max = 0.9)
                           ),
                           mainPanel(
                               h2("Shock simulado en el periodo 10"),
                               p("Pruebe partir del nivel de estado estacionario, k*=9, y*=3, 
                        con una tasa de ahorro del 3%, supongamos que esta sube exógenamente a 4%."),
                               plotOutput("quinta")
                           )
             ) #End sidebarLayout 2
             ) # End tabPanel Shocks
        ) # End TabsetPanel    
        ) # End fluidPage

#### End of Front-end

# Definimos una función que crea una matriz con la simulación de 100 periodos
# Notar que está fuera del Back-end (server) porque sólo la necesitamos definir en
# el ambiente global del servidor de R.

solow_matrix <- function(A=1,k_0=4, s=0.3, alpha=0.5, delta=0.1, n=100) {
    df <- data.frame(matrix(0, nrow = n, ncol = 7))
    colnames(df) <- c("Periodo","k", "y", "c", "inv", "dk", "Dk")
    df$Periodo <- 1:n
    
    # Parámetros 
    y = (A*k_0)^alpha
    c = y - s*y 
    inv = s*y
    depre = delta*k_0
    cambiok = s*y - delta*k_0
    
    # Los ponemos como 1a observación (periodo)
    df$k[1] = k_0
    df$y[1] = y
    df$c[1] = c
    df$inv[1] = inv 
    df$dk[1] = depre
    df$Dk[1] = cambiok
    
    # Iteramos
    for (j in 2:n) {
        df$k[j] <- df$k[j-1] + df$Dk[j-1]
        df$y[j] <- df$k[j]^(0.5)
        df$c[j] <- df$y[j]*(1-s)
        df$inv[j] <- s*df$y[j] 
        df$dk[j] <- delta*df$k[j]
        df$Dk[j] <- df$inv[j] - delta*df$k[j]
    }
    return(df)
}

#### Back-end 

server <- function(input, output) {
    library(tidyverse)
    ### Gráfica Consumo (K)
    output$primera <- renderPlot({
        datos <- solow_matrix(A=input$tec,k_0=input$k, s=input$s, alpha=input$alpha, delta=input$delta)
        plot(datos$Periodo,datos$k,type = "l",main = "Trayectoria del Capital (K)",
             xlab="Periodo", ylab="Capital (k)",lwd=3,col=2)
        #abline(h=9, col="blue")
    })
    ### Gráfica Consumo (C)
    output$segunda <- renderPlot({
        datos2 <- solow_matrix(A=input$tec,k_0=input$k, s=input$s, alpha=input$alpha, delta=input$delta)
        plot(datos2$Periodo,datos2$c,type = "l",main = "Trayectoria del Consumo (c)",
             xlab="Periodo", ylab="Consumo (c)",lwd=3,col=3)
    })
    ### Gráfica Inversión (I)
    output$tercera <- renderPlot({
        datos3 <- solow_matrix(A=input$tec,k_0=input$k, s=input$s, alpha=input$alpha, delta=input$delta)
    plot(datos3$Periodo,datos3$inv,type = "l",main = "Trayectoria de la Inversión (I)",
         xlab="Periodo", ylab="Inversión (I)",lwd=3,col=4)
    })
    ### Gráfica Producto (Y)
    output$cuarta <- renderPlot({
        datos4 <- solow_matrix(A=input$tec,k_0=input$k, s=input$s, alpha=input$alpha, delta=input$delta)
        plot(datos4$Periodo,datos4$y,type = "l",main = "Trayectoria del Producto (Y)",
             xlab="Periodo", ylab="Producto (Y)",lwd=3,col=6)
    })
    ### Gráfica Shocks
    output$quinta <- renderPlot({
        ee <- solow_matrix(k_0 = 9, s = .3) #Estado Estacionario
        yy <- solow_matrix(A=input$tec_1,k_0=input$k_1, s=input$s_1, alpha=input$alpha_1, delta=input$delta_1)
        zz <- rbind.data.frame(ee[1:10,], yy)
        
        zz$Periodo[1:110] <- 1:110
        
        zz %>% ggplot(aes(Periodo)) + 
            geom_line(aes(y=c, col="Consumo (C)")) + 
            geom_line(aes(y=inv, col="Inversión (I)")) + 
            geom_line(aes(y=y, col="Producto (Y)")) + 
            geom_hline(yintercept = 2.1, linetype = 2) +
            geom_vline(xintercept = 10, linetype = 2) + 
            labs(x = "Periodo", y = "Variables", 
                 title = "", 
                 subtitle = (""), 
                 caption = "") +
            theme_classic() +
            theme(legend.background = element_rect(fill = "transparent") ) + 
            scale_color_brewer(name= NULL,
                               palette = "Dark2")
    })
}

shinyApp(ui = ui, server = server)