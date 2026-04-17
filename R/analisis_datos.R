#' Asistente Interactivo de Análisis de Datos
#'
#'
#' Esta es la función principal y punto de entrada del paquete. Actúa como un consultor
#' automatizado que guía al usuario desde la carga de datos hasta la generación del código final,
#' seleccionando la metodología adecuada según la naturaleza de la variable objetivo.
#'
#' @section ¿Cómo funciona este Asistente paso a paso?:
#' \enumerate{
#'   \item **Escaneo y Validación:**
#'         Analiza la base de datos para asegurar que la variable objetivo existe.
#'         Además, realiza un test interno de autocorrelación para verificar si tiene sentido aplicar modelos de Series de Tiempo (evitando análisis espurios).
#'   \item **Menú Interactivo:**
#'         Despliega opciones en la consola basándose en los datos detectados, preguntando al usuario qué enfoque desea (Regresión Clásica, Modelos Generalizados o Series de Tiempo).
#'   \item **Generación de Estrategias (El "Motor"):**
#'         \itemize{
#'            \item \emph{Para Regresión y GLM:} Crea dos versiones del código.
#'            \enumerate{
#'               \item **Modelo Full:** Incluye todas las variables disponibles.
#'               \item **Modelo Optimizado:** Ejecuta internamente un algoritmo "Stepwise" basado en AIC para limpiar el modelo y quedarse solo con los predictores relevantes.
#'            }
#'            \item \emph{Para Series de Tiempo:} Ejecuta una "Búsqueda en Rejilla" (Grid Search) probando combinaciones de parámetros (p, d, q). Identifica automáticamente si la serie necesita diferenciación (ARIMA), si es puramente autoregresiva (AR), de media móvil (MA) o mixta (ARMA), y genera el código para la función específica.
#'         }
#' }
#'
#' @section Guía de la Gráfica Exploratoria Generada:
#' Dependiendo de la opción elegida, la función genera una visualización preliminar diferente para entender los datos antes de modelar:
#' \describe{
#'   \item{\strong{Escenario 1: Regresión Lineal (Histograma + Gauss)}}{
#'     Muestra un histograma azul de tus datos superpuesto con una curva roja (Campana de Gauss teórica).
#'     \emph{Objetivo:} Evaluar visualmente si la variable respuesta se parece a una distribución Normal (supuesto clave de la regresión lineal).
#'   }
#'   \item{\strong{Escenario 2: GLM (Barras o Histograma)}}{
#'     Detecta el tipo de variable:
#'     \itemize{
#'       \item *Si es Binaria (0/1):* Gráfico de barras naranja. Ayuda a ver si las clases están desbalanceadas (ej. muchos "0" y pocos "1").
#'       \item *Si es Conteo:* Histograma verde. Muestra la frecuencia de eventos (típico de Poisson).
#'     }
#'   }
#'   \item{\strong{Escenario 3: Serie de Tiempo (Línea de Tendencia)}}{
#'     Grafica la evolución de la variable a lo largo del índice temporal.
#'     \emph{Objetivo:} Identificar tendencias (subidas/bajadas a largo plazo) o ciclos repetitivos antes de ajustar el modelo.
#'   }
#' }
#'
#' @param data Data frame o base de datos que contiene las variables.
#' @param variable_objetivo Cadena de texto con el nombre exacto de la columna que deseas analizar (ej. "Ventas").
#'
#' @return No devuelve un objeto al entorno (return invisible). Su función es:
#' \enumerate{
#'   \item Imprimir en consola el **código R listo para copiar y pegar**.
#'   \item Generar el gráfico exploratorio en el panel de gráficos.
#' }
#'
#' @examples
#' # Supongamos que tienes una base llamada 'rendimiento_maiz_ton_ha'
#' # maiz<- analisis_datos(datos_rendimiento_maiz, "rendimiento_maiz_ton_ha")
#' # 1. Se abrirá el menú.
#' # 2. Eliges opción 1 (Regresión).
#' # 3. Copias el código que aparece en consola para ajustar tu modelo.
#'
#' @importFrom stats lm glm step cor shapiro.test arima gaussian binomial poisson acf
#' @importFrom utils menu
#' @import ggplot2
#' @export
analisis_datos <- function(data, variable_objetivo) {

  # --- Verificación de Dependencias ---
  # Se asegura de que ggplot2 esté instalado para poder generar las gráficas.
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Instala 'ggplot2'.")

  # --- 1. Validaciones e Información de la Base ---

  # Verificamos si la columna objetivo realmente existe en el data frame proporcionado.
  if (!(variable_objetivo %in% names(data))) stop("❌ Error: La variable objetivo no existe en la base de datos.")

  # Extracción de la variable respuesta (y) y lista de variables predictoras (vars_x).
  y <- data[[variable_objetivo]]
  vars_x <- setdiff(names(data), variable_objetivo)
  n <- nrow(data)

  # Obtenemos el nombre del objeto 'data' tal como lo escribió el usuario (ej. "mi_base")
  # para usarlo textualmente en la generación automática del código sugerido.
  nombre_data <- deparse(substitute(data))
  # --- 2A. Detección Automática y Recomendación ---

  # --- CORRECCIÓN: Definimos p_externas aquí ---
  p_externas <- length(vars_x)
  # ---------------------------------------------

  # 1. Detección de Series de Tiempo
  col_fecha <- names(data)[sapply(data, inherits, what = c("Date", "POSIXct", "POSIXt"))]
  tiene_fecha <- length(col_fecha) == 1
  hay_autocorr <- FALSE

  if (tiene_fecha) {
    # Ordenar por fecha y probar autocorrelación
    data <- data[order(data[[col_fecha]]), ]
    y_ord <- data[[variable_objetivo]]
    acf_obj <- try(stats::acf(y_ord, plot = FALSE), silent = TRUE)
    if (!inherits(acf_obj, "try-error")) {
      limite_acf <- 2 / sqrt(length(y_ord))
      hay_autocorr <- any(abs(acf_obj$acf[-1]) > limite_acf)
    }
  }
  es_ts <- tiene_fecha && hay_autocorr

  # --- 2. Detección Automática y Menú Interactivo (CORREGIDO) ---

  # --- 2A. Definición de Variables y Detección Inteligente de Fechas ---
  p_externas <- length(vars_x) # Definimos esto primero para evitar errores
  col_fecha <- NULL

  # 1. Buscar si ya existe columna fecha nativa
  cols_nativas <- names(data)[sapply(data, inherits, what = c("Date", "POSIXct", "POSIXt"))]

  if (length(cols_nativas) > 0) {
    col_fecha <- cols_nativas[1]
  } else {
    # 2. Intentar convertir columnas de texto a fecha
    formatos_prueba <- c("%d/%m/%Y", "%Y-%m-%d", "%m/%d/%Y", "%d-%m-%Y", "%Y/%m/%d")
    cols_texto <- names(data)[sapply(data, function(x) is.character(x) || is.factor(x))]

    for (col in cols_texto) {
      vals <- as.character(data[[col]])
      exito <- FALSE
      for (fmt in formatos_prueba) {
        fechas_test <- try(as.Date(vals, format = fmt), silent = TRUE)
        # Criterio: Si no da error y tiene pocos valores nulos (<20%), es fecha válida
        if (!inherits(fechas_test, "try-error") && sum(is.na(fechas_test)) < (0.2 * length(vals))) {
          data[[col]] <- fechas_test # ¡Auto-conversión!
          col_fecha <- col
          exito <- TRUE
          cat(sprintf("ℹ️  Auto-corrección: Columna '%s' detectada como FECHA (formato '%s').\n", col, fmt))
          break
        }
      }
      if (exito) break
    }
  }

  tiene_fecha <- !is.null(col_fecha)
  hay_autocorr <- FALSE

  # 3. Verificar Autocorrelación (solo si hay fecha)
  if (tiene_fecha) {
    data <- data[order(data[[col_fecha]]), ] # Ordenar cronológicamente
    y_ord <- data[[variable_objetivo]]
    acf_obj <- try(stats::acf(y_ord, plot = FALSE), silent = TRUE)
    if (!inherits(acf_obj, "try-error")) {
      limite_acf <- 2 / sqrt(length(y_ord))
      hay_autocorr <- any(abs(acf_obj$acf[-1]) > limite_acf)
    }
  }
  es_ts <- tiene_fecha && hay_autocorr

  # --- 2B. Lógica de Recomendación ---
  y_temp <- data[[variable_objetivo]]

  es_binaria <- all(y_temp %in% c(0, 1), na.rm=TRUE) || (is.factor(y_temp) && length(levels(y_temp))==2)
  es_conteo <- is.numeric(y_temp) && all(y_temp == floor(y_temp), na.rm=TRUE) && !es_binaria && all(y_temp >= 0, na.rm=TRUE)

  opcion_sugerida <- "1"
  razon_sugerencia <- "Variable continua (Regresión Lineal)."

  if (es_ts) {
    opcion_sugerida <- "3"
    razon_sugerencia <- sprintf("Serie de Tiempo (Fecha detectada: '%s' + Autocorrelación).", col_fecha)
  } else if (es_binaria) {
    opcion_sugerida <- "2"
    razon_sugerencia <- "Variable Binaria (Logística)."
  } else if (es_conteo) {
    opcion_sugerida <- "2"
    razon_sugerencia <- "Variable de Conteo (Poisson)."
  }

  # --- 2C. Menú Visual ---
  cat("\nAsistente de Análisis de Datos\n")
  cat("===========================================\n")
  cat(sprintf("   Base de datos: %s\n", nombre_data))
  cat(sprintf("   Objetivo: '%s'\n", variable_objetivo))
  cat(sprintf("   Observaciones: %d | Predictoras: %d\n\n", n, p_externas))

  cat("RECOMENDACIÓN AUTOMÁTICA:\n")
  cat(sprintf("   -> Sugerencia: Opción %s\n", opcion_sugerida))
  cat(sprintf("   -> Motivo: %s\n\n", razon_sugerencia))

  cat("Opciones Disponibles:\n")

  # Marcas visuales para guiar al usuario
  mark1 <- if(opcion_sugerida == "1") " ★ (Recomendado)" else ""
  mark2 <- if(opcion_sugerida == "2") " ★ (Recomendado)" else ""
  mark3 <- if(opcion_sugerida == "3") " ★ (Recomendado)" else ""

  cat(sprintf("   1. Regresión Lineal%s\n", mark1))
  cat(sprintf("   2. Modelos GLM (Logística/Poisson)%s\n", mark2))
  cat(sprintf("   3. Series de Tiempo%s\n", mark3))

  # AQUÍ ESTABA EL ERROR: Ahora sí definimos 'seleccion'
  seleccion <- readline(prompt = "Selecciona una opción (1-3) [Enter para auto]: ")

  if (seleccion == "") {
    seleccion <- opcion_sugerida
    cat(sprintf("-> Selección automática aplicada: %s\n", seleccion))
  }

  # --- 3. Procesamiento según Selección ---
  # Inicialización de variables para almacenar el código generado y explicaciones.
  codigo_full <- ""
  codigo_opt <- ""
  explicacion <- ""
  p_graph <- NULL

  # ==============================================================================
  # OPCIÓN 1: REGRESIÓN LINEAL
  # Estrategia: Ajustar mínimos cuadrados ordinarios (RLS/RLM).
  # ==============================================================================
  if (seleccion == "1") {
    cat("\n Analizando para Regresión Lineal...\n")

    # A) Código Modelo Completo
    # Genera el string del código para un modelo que incluye TODAS las variables (formula = y ~ .).
    codigo_full <- paste0(
      sprintf("modelo_rlm_full <- rlm(%s ~ ., data = %s)\n", variable_objetivo, nombre_data),
      "modelo_rlm_full\n",
      "plot(modelo_rlm_full)"
    )

    # B) Código Modelo Optimizado (Stepwise AIC)
    # Ejecuta internamente un 'stepwise selection' para encontrar qué variables
    # realmente aportan valor al modelo y descartar las redundantes.
    formula_opt <- "."
    vars_seleccionadas <- "Todas"

    if (length(vars_x) > 0) {
      # Ajuste temporal para evaluación (no se muestra al usuario).
      mod_temp <- stats::lm(as.formula(paste(variable_objetivo, "~ .")), data = data)
      # Selección de variables basada en AIC (trace=0 silencia la salida técnica).
      mod_step <- stats::step(mod_temp, trace = 0)
      vars_opt <- names(coef(mod_step))[-1] # Extraemos nombres de coeficientes (excluyendo intercepto).

      # Si se logró reducir variables, actualizamos la fórmula optimizada.
      if (length(vars_opt) > 0 && length(vars_opt) < length(vars_x)) {
        formula_opt <- paste(vars_opt, collapse = " + ")
        vars_seleccionadas <- paste(vars_opt, collapse = ", ")
      }
    }

    # Genera el string del código para el modelo optimizado.
    codigo_opt <- paste0(
      sprintf("modelo_rlm_opt <- rlm(%s ~ %s, data = %s)\n", variable_objetivo, formula_opt, nombre_data),
      "modelo_rlm_opt\n",
      "plot(modelo_rlm_opt)"
    )

    # C) Gráfica Exploratoria con NOTAS
    # Crea un histograma de la variable objetivo y superpone una curva normal teórica
    # para enseñar al usuario sobre el supuesto de normalidad visualmente.
    media_y <- mean(y, na.rm=TRUE); sd_y <- sd(y, na.rm=TRUE)
    p_graph <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[variable_objetivo]])) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30, fill = "#3498DB", color = "white", alpha = 0.7) +
      ggplot2::stat_function(fun = dnorm, args = list(mean = media_y, sd = sd_y), color = "red", size = 1.2) +
      ggplot2::labs(
        title = paste("Distribución de", variable_objetivo),
        subtitle = "Histograma (Datos) vs Campana de Gauss (Teórico)",
        x = variable_objetivo, y = "Densidad",
        caption = "NOTA:\n• La línea ROJA representa la normalidad perfecta.\n• Si las barras azules se alejan mucho de la forma de campana roja,\n  los supuestos de normalidad podrían no cumplirse (considera transformar variables)."
      ) + ggplot2::theme_minimal()
  }

  # ==============================================================================
  # OPCIÓN 2: GLM (Logística / Poisson)
  # Estrategia: Detectar automáticamente si es binaria (0/1) o conteo (enteros positivos).
  # ==============================================================================
  else if (seleccion == "2") {
    cat("\n Analizando para GLM...\n")

    # Detección heurística del tipo de variable respuesta
    es_binaria <- all(y %in% c(0, 1), na.rm=TRUE) || (is.factor(y) && length(levels(y))==2)
    es_conteo <- is.numeric(y) && all(y == floor(y)) && !es_binaria

    # Configuración de parámetros según el tipo detectado
    if (es_binaria) {
      fun_glm <- "logistico"; nombre_mod <- "modelo_log"; familia <- "binomial"
      cat("   -> Detectada Variable Binaria. Sugerencia: Regresión Logística.\n")
      nota_grafica <- "NOTA:\n• Revisa el balance entre las barras (0 y 1).\n• Si una barra es diminuta comparada con la otra (Desbalance),\n  el modelo tendrá dificultades para predecir la clase minoritaria."
    } else if (es_conteo) {
      fun_glm <- "reg_poisson"; nombre_mod <- "modelo_pois"; familia <- "poisson"
      cat("   -> Detectada Variable de Conteo. Sugerencia: Regresión Poisson.\n")
      nota_grafica <- "NOTA:\n• En Poisson, es común ver barras altas a la izquierda (valores bajos)\n  que descienden hacia la derecha. Esto es normal y esperado."
    } else {
      # Fallback a gaussiana si no es clara la detección
      fun_glm <- "glm_gen"; nombre_mod <- "modelo_glm"; familia <- "gaussian"
      cat("⚠️ Tipo ambiguo. Usando configuración genérica.\n")
      nota_grafica <- "NOTA: Distribución de la variable respuesta."
    }

    # A) Código Modelo Completo
    codigo_full <- paste0(
      sprintf("%s_full <- %s(%s ~ ., data = %s)\n", nombre_mod, fun_glm, variable_objetivo, nombre_data),
      sprintf("%s_full\n", nombre_mod),
      sprintf("plot(%s_full)", nombre_mod)
    )

    # B) Código Modelo Optimizado
    # Intenta realizar una selección de variables stepwise para GLM.
    formula_opt <- "."
    if (length(vars_x) > 0) {
      try({
        mod_temp <- stats::glm(as.formula(paste(variable_objetivo, "~.")), data = data, family = familia)
        mod_step <- stats::step(mod_temp, trace = 0)
        vars_opt <- names(coef(mod_step))[-1]
        if (length(vars_opt) > 0 && length(vars_opt) < length(vars_x)) {
          formula_opt <- paste(vars_opt, collapse = " + ")
        }
      }, silent = TRUE)
    }

    codigo_opt <- paste0(
      sprintf("%s_opt <- %s(%s ~ %s, data = %s)\n", nombre_mod, fun_glm, variable_objetivo, formula_opt, nombre_data),
      sprintf("%s_opt\n", nombre_mod),
      sprintf("plot(%s_opt)", nombre_mod)
    )

    # C) Gráfica
    # Gráfico de barras para binaria o histograma para conteo/continuo
    p_graph <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[variable_objetivo]])) +
      (if(es_binaria) ggplot2::geom_bar(fill = "orange", alpha=0.8) else ggplot2::geom_histogram(fill="#2ECC71", bins=20, color="white")) +
      ggplot2::labs(
        title = "Distribución de la Variable Respuesta",
        caption = nota_grafica
      ) + ggplot2::theme_minimal()
  }

  # ==============================================================================
  # OPCIÓN 3: SERIES DE TIEMPO (Sin Alertas)
  # Estrategia: Grid Search (fuerza bruta controlada) para hallar (p, d, q) óptimos por AIC.
  # ==============================================================================
  else if (seleccion == "3") {
    cat("\n Análisis de Serie de Tiempo (AR / MA / ARMA / ARIMA)...\n")

    if (p_externas > 0) {
      cat("   Los modelos AR/MA/ARMA/ARIMA son univariados.\n")
    }


    mejor_aic <- Inf
    mejor_orden <- c(0, 0, 0) # Orden inicial (p, d, q)

    # Grid Search con silenciador de warnings (suppressWarnings)
    # Se iteran combinaciones de p (AR), d (Diferencias) y q (MA).
    for (d_test in 0:1) {
      for (p_test in 0:3) {
        for (q_test in 0:3) {
          if (p_test==0 && q_test==0 && d_test==0) next # Saltar modelo nulo
          try({
            # Ajuste de prueba. suppressWarnings() es clave aquí para evitar que el usuario
            # se asuste con alertas de "posible no convergencia" durante la búsqueda automática.
            mod <- suppressWarnings(stats::arima(y, order = c(p_test, d_test, q_test), method="CSS-ML"))

            # Guardamos el modelo si tiene el menor AIC encontrado hasta ahora
            if (!is.null(mod$aic) && mod$aic < mejor_aic) {
              mejor_aic <- mod$aic
              mejor_orden <- c(p_test, d_test, q_test)
            }
          }, silent = TRUE)
        }
      }
    }

    p_opt <- mejor_orden[1]; d_opt <- mejor_orden[2]; q_opt <- mejor_orden[3]

    # Selección de Función de la paquetería según el orden óptimo encontrado
    if (d_opt > 0) {
      # Si d > 0, requiere diferenciación -> ARIMA
      funcion_chapi <- "modelo_arima"; nombre_modelo <- "modelo_ARIMAchapi"
      args <- sprintf("p = %d, d = %d, q = %d", p_opt, d_opt, q_opt)
    } else {
      # Si d = 0, es estacionaria (teóricamente). Verificamos componentes AR y MA.
      if (p_opt > 0 && q_opt == 0) {
        funcion_chapi <- "modelo_ar"; nombre_modelo <- "modelo_AR"; args <- sprintf("p = %d", p_opt)
      } else if (p_opt == 0 && q_opt > 0) {
        funcion_chapi <- "modelo_ma"; nombre_modelo <- "modelo_MA"; args <- sprintf("q = %d, include_mean = TRUE", q_opt)
      } else {
        funcion_chapi <- "modelo_arma"; nombre_modelo <- "modelo_ARMA"; args <- sprintf("p = %d, q = %d, include_mean = TRUE", p_opt, q_opt)
      }
    }

    explicacion <- sprintf("AIC óptimo encontrado con orden (%d, %d, %d).", p_opt, d_opt, q_opt)

    # Código Único: En series de tiempo univariadas no seleccionamos variables,
    # sino que seleccionamos el orden del modelo, por eso solo hay un bloque de código.
    codigo_salida <- paste0(
      sprintf("ts_data <- %s$%s\n", nombre_data, variable_objetivo),
      sprintf("%s <- %s(ts_data, %s)\n", nombre_modelo, funcion_chapi, args),
      sprintf("%s\n", nombre_modelo),
      sprintf("plot(%s)", nombre_modelo)
    )

    # Gráfica de Serie de Tiempo con Notas
    p_graph <- ggplot2::ggplot(data, ggplot2::aes(x = 1:n, y = .data[[variable_objetivo]])) +
      ggplot2::geom_line(color = "#8E44AD") +
      ggplot2::labs(
        title = paste("Evolución Temporal:", variable_objetivo),
        x = "Tiempo / Índice", y = "Valor",
        caption = "NOTA:\n• Observa la gráfica: ¿Los datos suben o bajan con el tiempo (Tendencia)?\n• ¿Se repiten patrones regularmente (Ciclos)?\n• El modelo sugerido abajo ha sido seleccionado para manejar estas características."
      ) + ggplot2::theme_minimal()

    codigo_full <- codigo_salida # Para unificar la salida final
    codigo_opt <- "" # No aplica segunda opción en TS univariado
  } else {
    stop("❌ Opción inválida.")
  }

  # --- 4. Salida Final ---
  # Presentación de resultados al usuario
  cat("\n✅ ANÁLISIS COMPLETADO.\n")
  if (explicacion != "") cat(paste("ℹ️  NOTA:", explicacion, "\n"))
  cat("--------------------------------------------------\n")

  # Imprimir Opción 1 (Código sugerido)
  if (seleccion == "3") {
    cat(" Código Sugerido (Mejor Ajuste AIC):\n")
    cat("--------------------------------------------------\n")
    cat(codigo_full, "\n")
    cat("--------------------------------------------------\n")
  } else {
    cat(" OPCIÓN A: Modelo con TODAS las variables\n")
    cat("--------------------------------------------------\n")
    cat(codigo_full, "\n")
    cat("--------------------------------------------------\n\n")

    cat(" OPCIÓN B: Modelo OPTIMIZADO (Variables Seleccionadas)\n")
    cat("   (Elimina variables que aportan ruido o son redundantes)\n")
    cat("--------------------------------------------------\n")
    cat(codigo_opt, "\n")
    cat("--------------------------------------------------\n")
  }

  # Mostrar gráfico generado
  print(p_graph)
}


