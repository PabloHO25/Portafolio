#install.packages(c("shiny", "DT", "foreign", "ggplot2", "cluster"))
install.packages("shinylive")


library(shiny)
library(DT)
library(foreign)
library(ggplot2)
options(shiny.maxRequestSize = 100 * 1024^2)

# Las métricas no alteran el objetivo del recocido: sigue siendo la pureza.
calcular_metricas <- function(X, modelo, clases, semilla) {
  grupos <- modelo$cluster; n <- nrow(X); k <- nrow(modelo$centers)
  m <- table(grupos, clases); a <- rowSums(m); b <- colSums(m)
  c2 <- function(z) z * (z - 1) / 2
  pares <- c2(n); acuerdos <- sum(c2(m)); sa <- sum(c2(a)); sb <- sum(c2(b))
  rand <- (pares - sa - sb + 2 * acuerdos) / pares
  dispersion <- vapply(seq_len(k), function(j) {
    delta <- sweep(X[grupos == j, , drop = FALSE], 2, modelo$centers[j, ], "-")
    mean(sqrt(rowSums(delta^2)))
  }, numeric(1))
  dist_centros <- as.matrix(dist(modelo$centers))
  razones <- outer(dispersion, dispersion, "+") / dist_centros
  diag(razones) <- -Inf
  db <- if (any(dist_centros[row(dist_centros) != col(dist_centros)] == 0)) NA_real_ else mean(apply(razones, 1, max))
  # Silueta exacta con todos los registros y distancia euclídea.
  gs <- as.integer(factor(grupos))
  sil <- if (k >= 2 && k < n) mean(cluster::silhouette(gs, dist(X))[, "sil_width"]) else NA_real_
  data.frame(Tipo = c("Externa", "Interna", "Externa", "Interna"),
             Medida = c("Pureza (%)", "Silueta media", "Índice Rand", "Davies–Bouldin"),
             Valor = c(100 * sum(apply(m, 1, max)) / n, sil, rand, db),
             Preferencia = c(rep("Mayor es mejor", 3), "Menor es mejor"), Registros = rep(n, 4),
             Nota = c("Usa las clases de Result", "Todos los registros; distancia euclídea",
                      "Acuerdo entre pares; sin ajuste por azar", "NA si hay centroides coincidentes"), check.names = FALSE)
  
}

ui <- fluidPage(
  titlePanel("Selección de variables mediante recocido simulado"),
  sidebarLayout(
    sidebarPanel(
      fileInput("archivo", "Dataset ARFF", accept = ".arff"),
      numericInput("temperatura", "Temperatura inicial (> 1)", 10000, min = 1.01),
      numericInput("alpha", "Factor de enfriamiento (0 < α < 1)", 0.99, min = 0.001, max = 0.9999, step = 0.001),
      helpText("Ejemplo: α = 0.99 reduce la temperatura un 1 % por iteración. La temperatura final es 1."),
      numericInput("clusters", "Número de clusters", 2, min = 2, step = 1),
      numericInput("semilla", "Semilla", 150, min = 0, max = 1000000, step = 1),
      actionButton("ejecutar", "Ejecutar recocido", class = "btn-primary"), br(), br(),
      downloadButton("descargar", "Descargar historial CSV"),
      helpText("Objetivo: maximizar la pureza respecto a Result. Puede tardar varios minutos. Las etiquetas se usan para seleccionar variables; esta pureza no es una evaluación independiente de generalización.")
    ),
    mainPanel(
      verbatimTextOutput("resumen"),
      tabsetPanel(
        tabPanel("Convergencia", plotOutput("convergencia", height = "430px")),
        tabPanel("Gráfico de clusters", plotOutput("clusters_plot", height = "480px")),
        tabPanel("Métricas", br(),
                 helpText("Medidas de la mejor solución. Las externas comparan clusters con Result; las internas usan las variables seleccionadas y estandarizadas. NA indica una medida no definida. Todas las medidas usan todos los registros. La silueta calcula las distancias entre todos los pares."),
                 DTOutput("metricas"), br(), downloadButton("descargar_metricas", "Descargar métricas CSV")),
        tabPanel("Historial", DTOutput("historial")),
        tabPanel("Variables seleccionadas", DTOutput("variables")),
        tabPanel("Clusters y clases", DTOutput("contingencia"))
      )
    )
  )
)

server <- function(input, output, session) {
  resultado <- reactiveVal(NULL)
  observeEvent(input$archivo, { resultado(NULL) })
  observeEvent(input$ejecutar, {
    req(input$archivo)
    tryCatch({
      t0 <- input$temperatura; alpha <- input$alpha; k <- input$clusters; semilla <- input$semilla
      if (!is.finite(t0) || t0 <= 1) stop("La temperatura inicial debe ser mayor que 1.")
      if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("El factor de enfriamiento debe estar entre 0 y 1.")
      if (!is.finite(k) || k < 2 || k != floor(k)) stop("El número de clusters debe ser un entero mayor o igual a 2.")
      if (!is.finite(semilla) || semilla < 0 || semilla > 1000000 || semilla != floor(semilla)) stop("La semilla debe ser un entero entre 0 y 1000000.")
      n_iter <- ceiling(log(1 / t0) / log(alpha))
      if (!is.finite(n_iter) || n_iter > 20000) stop("Esta configuración supera 20000 iteraciones; disminuye α o la temperatura inicial.")
      df <- as.data.frame(read.arff(input$archivo$datapath))
      if (!"Result" %in% names(df)) stop("El archivo debe contener la columna Result.")
      y <- df$Result; X <- df[, names(df) != "Result", drop = FALSE]
      if (!ncol(X)) stop("No hay variables predictoras.")
      X[] <- lapply(X, function(x) suppressWarnings(as.numeric(as.character(x))))
      if (anyNA(y) || anyNA(X) || any(!is.finite(as.matrix(X)))) stop("Hay valores faltantes, infinitos o predictores no numéricos. Limpia el dataset antes de cargarlo.")
      variables_validas <- vapply(X, function(x) length(unique(x)) > 1, logical(1))
      omitidas <- names(X)[!variables_validas]; X <- X[, variables_validas, drop = FALSE]
      if (!ncol(X)) stop("Todas las variables son constantes.")
      if (nrow(X) <= k) stop("El número de filas debe ser mayor que el número de clusters.")
      X <- scale(X); p <- ncol(X); cache <- new.env(parent = emptyenv())
      evaluar <- function(mask) {
        clave <- paste(mask, collapse = "")
        if (exists(clave, envir = cache, inherits = FALSE)) return(get(clave, envir = cache))
        xs <- X[, mask == 1, drop = FALSE]
        if (nrow(unique(xs)) < k) return(NULL)
        # Restaurar el estado aleatorio para no alterar la generación de vecinos.
        rng <- .Random.seed
        on.exit(assign(".Random.seed", rng, envir = .GlobalEnv))
        set.seed(semilla)
        modelo <- kmeans(xs, centers = k, nstart = 25, iter.max = 100)
        matriz <- table(Cluster = modelo$cluster, Clase_real = y)
        res <- list(score = 100 * sum(apply(matriz, 1, max)) / sum(matriz), matriz = matriz)
        assign(clave, res, envir = cache); res
      }
      set.seed(semilla)
      mask <- sample(0:1, p, replace = TRUE)
      if (!any(mask)) mask[1] <- 1L
      actual <- evaluar(mask)
      if (is.null(actual)) { mask <- rep(1L, p); actual <- evaluar(mask) }
      if (is.null(actual)) stop("No hay suficientes observaciones distintas para formar los clusters.")
      score <- best <- actual$score; best_mask <- mask; temp <- t0
      fila <- function(i, t, aceptada) data.frame(Iteracion = i, Temperatura = t, Pureza_actual = score,
                                                  Mejor_pureza = best, Variables = sum(mask), Aceptada = aceptada, Mask = paste(mask, collapse = ","))
      historial <- vector("list", n_iter + 1L); historial[[1]] <- fila(0, temp, NA)
      inicio <- proc.time()[["elapsed"]]
      withProgress(message = "Ejecutando recocido simulado", value = 0, {
        for (i in seq_len(n_iter)) {
          candidata <- mask
          posiciones <- sample.int(p, sample.int(min(3L, p), 1))
          candidata[posiciones] <- 1L - candidata[posiciones]
          aceptada <- FALSE
          if (any(candidata)) {
            evaluacion <- evaluar(candidata)
            if (!is.null(evaluacion)) {
              diferencia <- evaluacion$score - score
              aceptada <- diferencia >= 0 || runif(1) < exp(diferencia / temp)
              if (aceptada) { mask <- candidata; score <- evaluacion$score }
              if (score > best) { best <- score; best_mask <- mask }
            }
          }
          # Temperatura registrada: la utilizada para decidir la aceptación.
          historial[[i + 1L]] <- fila(i, temp, aceptada)
          temp <- temp * alpha
          if (i %% 10L == 0L || i == n_iter) setProgress(value = i / n_iter, detail = paste(i, "de", n_iter))
        }
      })
      xs_mejor <- X[, best_mask == 1, drop = FALSE]
      set.seed(semilla)
      modelo_mejor <- kmeans(xs_mejor, centers = k, nstart = 25, iter.max = 100)
      if (ncol(xs_mejor) >= 2) {
        pca <- prcomp(xs_mejor, center = FALSE, scale. = FALSE)
        coordenadas <- pca$x[, 1:2, drop = FALSE]
        varianza <- 100 * pca$sdev^2 / sum(pca$sdev^2)
        ejes <- sprintf("CP%d (%.1f %%)", 1:2, varianza[1:2])
      } else {
        coordenadas <- cbind(xs_mejor[, 1], rep(0, nrow(xs_mejor)))
        ejes <- c(paste0(colnames(xs_mejor)[1], " (estandarizada)"), "")
      }
      grafico_clusters <- data.frame(X = coordenadas[, 1], Y = coordenadas[, 2], Cluster = factor(modelo_mejor$cluster))
      resultado(list(metricas = calcular_metricas(xs_mejor, modelo_mejor, y, semilla), grafico_clusters = grafico_clusters, ejes = ejes, historial = do.call(rbind, historial), best = best, mask = best_mask,
                     nombres = colnames(X), omitidas = omitidas, matriz = evaluar(best_mask)$matriz,
                     n = nrow(X), iteraciones = n_iter, segundos = proc.time()[["elapsed"]] - inicio,
                     t0 = t0, alpha = alpha, final = temp, k = k, semilla = semilla))
      showNotification("Ejecución terminada.", type = "message")
    }, error = function(e) showNotification(conditionMessage(e), type = "error", duration = NULL))
  })
  output$resumen <- renderText({
    r <- resultado()
    if (is.null(r)) return("Carga tu archivo ARFF y pulsa Ejecutar recocido.")
    paste0("Mejor pureza: ", round(r$best, 4), " % | Variables: ", sum(r$mask), "/", length(r$mask),
           "\nFilas: ", r$n, " | Iteraciones: ", r$iteraciones, " | Tiempo: ", round(r$segundos, 1), " s",
           "\nTemperatura inicial: ", r$t0, " | α: ", r$alpha, " | Temperatura final: ", signif(r$final, 5),
           "\nClusters: ", r$k, " | Semilla: ", r$semilla,
           if (length(r$omitidas)) paste0("\nVariables constantes excluidas: ", paste(r$omitidas, collapse = ", ")) else "")
  })
  output$convergencia <- renderPlot({
    r <- resultado(); req(r); h <- r$historial
    ggplot(h, aes(x = Iteracion)) +
      geom_line(aes(y = Pureza_actual, color = "Solución actual"), linewidth = 0.45, alpha = 0.65) +
      geom_line(aes(y = Mejor_pureza, color = "Mejor encontrada"), linewidth = 1) +
      scale_color_manual(values = c("Solución actual" = "#8B98A5", "Mejor encontrada" = "#087E8B")) +
      labs(title = "Convergencia del recocido simulado", x = "Iteración", y = "Pureza (%)", color = NULL) +
      theme_minimal(base_size = 13) + theme(legend.position = "bottom")
  })
  output$clusters_plot <- renderPlot({
    r <- resultado(); req(r)
    ggplot(r$grafico_clusters, aes(X, Y, color = Cluster)) +
      geom_point(size = 1.5, alpha = 0.6) +
      labs(title = "Clusters de la mejor solución", x = r$ejes[1], y = r$ejes[2],
           subtitle = if (sum(r$mask) >= 2) "Proyección PCA de las variables seleccionadas" else "Solución con una variable seleccionada") +
      theme_minimal(base_size = 13) + theme(legend.position = "bottom")
  })
  output$metricas <- renderDT({
    r <- resultado(); req(r)
    formatRound(datatable(r$metricas, rownames = FALSE,
                          options = list(dom = "t", scrollX = TRUE, pageLength = 10)), "Valor", digits = 4)
  })
  output$descargar_metricas <- downloadHandler(
    filename = function() "metricas_recocido_simulado.csv",
    content = function(file) { req(resultado()); write.csv(resultado()$metricas, file, row.names = FALSE) }
  )
  output$historial <- renderDT({
    req(resultado())
    formatRound(datatable(resultado()$historial, rownames = FALSE, filter = "top",
                          options = list(pageLength = 10, scrollX = TRUE, order = list(list(0, "asc")))),
                c("Temperatura", "Pureza_actual", "Mejor_pureza"), digits = 4)
  })
  output$variables <- renderDT({
    r <- resultado(); req(r)
    datatable(data.frame(Variable = r$nombres, Seleccionada = ifelse(r$mask == 1, "Sí", "No")),
              rownames = FALSE, filter = "top", options = list(pageLength = 15))
  })
  output$contingencia <- renderDT({
    r <- resultado(); req(r); m <- r$matriz
    d <- data.frame(Cluster = rownames(m), as.data.frame.matrix(m), check.names = FALSE)
    d$Pureza_cluster <- 100 * apply(m, 1, max) / rowSums(m)
    formatRound(datatable(d, rownames = FALSE, options = list(dom = "t", scrollX = TRUE)), "Pureza_cluster", 4)
  })
  output$descargar <- downloadHandler(
    filename = function() "historial_recocido_simulado.csv",
    content = function(file) { req(resultado()); write.csv(resultado()$historial, file, row.names = FALSE) }
  )
}
shinyApp(ui, server)
