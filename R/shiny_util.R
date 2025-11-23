
plot_ui <- function(id, width=800, height=600, margin_controls=TRUE, ...) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::div(
            style="display: grid; grid-template-columns: auto auto auto auto 1fr; gap: 1em;",
            shiny::numericInput(ns("width"), "Plot width", width, min=100, max=10000, step=50, width="10em"),
            shiny::numericInput(ns("height"), "Plot height", height, min=100, max=10000, step=50, width="10em"),
            if (margin_controls) shiny::numericInput(ns("margin_left"), "Left margin (0-1)", NA, min=0,max=1,step=0.1, width="10em") else shiny::div(),
            if (margin_controls) shiny::numericInput(ns("margin_right"), "Right margin (0-1)", NA, min=0,max=1,step=0.1, width="10em") else shiny::div(),
            shiny::div(shiny::tags$label("Download"), shiny::tags$br(),
                shiny::downloadButton(ns("pdf"), "PDF"),
                shiny::downloadButton(ns("eps"), "EPS"),
                shiny::downloadButton(ns("svg"), "SVG"),
                shiny::downloadButton(ns("png"), "PNG")
            )
        ),
        shiny::plotOutput(ns("plot"), width="auto", height="auto", ...)
    )   
}

plot_server <- function(id, callback, dlname="plot", dpi=96) { moduleServer(id, function(input, output, session) {
    do_plot <- function() {
        p <- callback()
        if (!is.null(input$margin_left) || !is.null(input$margin_right)) {
            p <- plot_ensure_margin(p, left=input$margin_left, right=input$margin_right)
            plot(p)
        } else {
            print(p)
        }
    }
    
    output$plot <- shiny::renderPlot(
        withProgress({ do_plot() }, message="Plotting", value=NA),
        width=\() input$width,
        height=\() input$height,
        res=dpi
    )
    
    output$pdf <- shiny::downloadHandler(
        paste0(dlname,".pdf"),
        function(filename) {
            pdf(filename, width=input$width/dpi, height=input$height/dpi)
            do_plot()
            dev.off()
        }
    )
    
    output$eps <- shiny::downloadHandler(
        paste0(dlname,".eps"),
        function(filename) {
            postscript(filename, width=input$width/dpi, height=input$height/dpi,
                       paper="special", onefile=FALSE, horizontal=FALSE)
            do_plot()
            dev.off()
        }
    )
    
    output$svg <- shiny::downloadHandler(
        paste0(dlname,".svg"),
        function(filename) {
            svg(filename, width=input$width/dpi, height=input$height/dpi)
            do_plot()
            dev.off()
        }
    )
    
    output$png <- shiny::downloadHandler(
        paste0(dlname,".png"),
        function(filename) {
            png(filename, width=input$width*2, height=input$height*2, res=dpi*2)
            do_plot()
            dev.off()
        }
    )
})}
