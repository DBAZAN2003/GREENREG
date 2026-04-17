#' Modelo de Media Móvil (MA)
#'
#' Esta función ajusta un modelo de series de tiempo donde el valor actual ($Y_t$)
#' se explica mediante una combinación lineal de los errores (shocks) pasados
#' ($\epsilon_{t-1}, \epsilon_{t-2}, ...$). Es ideal para modelar fenómenos con
#' efectos de memoria corta o impactos transitorios.
#'
#' @details
#' El flujo interno de la función se divide en cinco etapas fundamentales:
#' \enumerate{
#'   \item **Estimación (MLE):** Utiliza Máxima Verosimilitud para estimar los
#'         parámetros $\theta$ (Theta). A diferencia del AR, el MA requiere
#'         algoritmos iterativos ya que los errores no son directamente observables.
#'   \item **Análisis de Invertibilidad:** Calcula las raíces inversas del polinomio
#'         característico. Un modelo MA debe ser "Invertible" para que el error
#'         pueda expresarse como una función convergente de los datos observados.
#'   \item **Identificación Estructural:** Evalúa la ACF (Autocorrelación). En un
#'         proceso MA(q), la ACF debe "cortarse" abruptamente después del lag $q$.
#'   \item **Auditoría de Residuos:** Verifica que el error remanente sea Ruido Blanco
#'         mediante la prueba de Ljung-Box y la distribución Gaussiana (Shapiro-Wilk).
#'   \item **Traducción Didáctica:** Genera un reporte que explica la persistencia
#'         de los choques y valida si el modelo es único y matemáticamente sólido.
#' }
#'
#' @section Guía de las 9 Gráficas de Diagnóstico (`plot`):
#' Al ejecutar \code{plot(modelo)}, se despliega una auditoría de "shocks" temporales:
#' \describe{
#'   \item{\strong{1. Inspección Visual}}{Muestra la serie bruta con su media histórica. El modelo MA supone que la serie siempre es atraída hacia este nivel central.}
#'   \item{\strong{2. Ajuste Temporal}}{Superpone el modelo (azul) sobre la realidad. Permite ver cómo el modelo "suaviza" los picos de la serie.}
#'   \item{\strong{3. Precisión (Real vs Ajustado)}}{Gráfica de dispersión para evaluar la correlación. Los puntos deben seguir la diagonal roja.}
#'   \item{\strong{4. Círculo Unitario}}{Prueba de Invertibilidad. Los rombos azules deben estar dentro del círculo para asegurar que el modelo sea único.}
#'   \item{\strong{5. Residuos en el Tiempo}}{Evaluación de aleatoriedad. No debe haber patrones; los errores deben ser impredecibles.}
#'   \item{\strong{6. ACF de Residuos}}{Prueba crítica de Ruido Blanco. Si hay barras largas (números rojos), el orden \code{q} es insuficiente.}
#'   \item{\strong{7. PACF de Residuos}}{Complemento de la ACF para verificar que no existan dependencias directas remanentes.}
#'   \item{\strong{8. Normalidad (Q-Q Plot)}}{Valida que los choques sigan una campana de Gauss, asegurando la validez de los intervalos de confianza.}
#'   \item{\strong{9. Estacionariedad Móvil}}{Confirma que la media y varianza no se desplazan. Vital para asegurar que el modelo MA sea aplicable.}
#' }
#'
#' @param x Vector numérico o serie de tiempo (\code{ts}). No debe contener valores \code{NA}.
#' @param q Orden del modelo (número de choques pasados a considerar). Por defecto es 1.
#' @param include_mean Lógico. ¿Incluir el promedio constante de la serie? Por defecto \code{TRUE}.
#'
#' @return Un objeto de clase \code{"ma_greenreg"} con coeficientes $\theta$,
#'         pruebas de invertibilidad, diagnósticos de residuos y dataframes para
#'         visualización con \code{ggplot2}.
#' @examples
#' data("datos_anomalia_temperatura")
#' modelo_MA <- modelo_ma(datos_anomalia_temperatura$anomalia_c, q = 1, include_mean = TRUE)
#' modelo_MA
#' plot(modelo_MA)
#'
#' @importFrom stats arima coef residuals fitted Box.test pnorm qnorm shapiro.test acf pacf logLik filter sd
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
#' @export
modelo_ma <- function(x, q = 1, include_mean = TRUE) {

  # --- 1. Validación y Ajuste ---
  if (!is.numeric(x)) stop("El argumento 'x' debe ser numérico.")

  # Ajuste usando ARIMA(0,0,q) que equivale a MA(q)
  modelo_base <- stats::arima(x, order = c(0, 0, q), include.mean = include_mean)

  # --- 2. Coeficientes y Estadísticas ---
  coefs <- modelo_base$coef
  se <- sqrt(diag(modelo_base$var.coef))
  z_val <- coefs / se
  p_val <- 2 * (1 - stats::pnorm(abs(z_val)))

  tabla_coefs <- cbind(Estimate = coefs, `Std. Error` = se, `z value` = z_val, `Pr(>|z|)` = p_val)

  # --- 3. Análisis de Invertibilidad (Raíces Unitarias) ---
  # La invertibilidad en MA es equivalente a la estabilidad en AR
  # Extraemos solo los coeficientes MA (ignoramos intercepto)
  ma_coefs <- coefs[grepl("^ma", names(coefs))]

  if (length(ma_coefs) > 0) {
    # Polinomio característico: 1 + theta1*z + theta2*z^2 ...
    poly_roots <- polyroot(c(1, ma_coefs))
    inverse_roots <- 1 / poly_roots
    modulos <- Mod(inverse_roots)
    # Regla: Invertible si TODOS los módulos son < 1
    es_invertible <- all(modulos < 1)
  } else {
    inverse_roots <- NULL
    es_invertible <- TRUE
  }

  # --- 4. Diagnósticos de Residuos ---
  residuos <- stats::residuals(modelo_base)

  # Ljung-Box (Ruido Blanco)
  lb_test <- stats::Box.test(residuos, type = "Ljung-Box", lag = q + 5)

  # Shapiro-Wilk (Normalidad)
  if (length(residuos) >= 3 && length(residuos) <= 5000) {
    shapiro_test <- stats::shapiro.test(residuos)
  } else {
    shapiro_test <- list(p.value = NA)
  }

  # --- 5. Datos para Gráficas ---
  ajustados <- x - residuos
  sigma <- sqrt(modelo_base$sigma2)
  # Intervalos de confianza 95%
  lower_ci <- ajustados - 1.96 * sigma
  upper_ci <- ajustados + 1.96 * sigma

  df_plot <- data.frame(
    Tiempo = 1:length(x),
    Observado = as.numeric(x),
    Ajustado = as.numeric(ajustados),
    Lower = as.numeric(lower_ci),
    Upper = as.numeric(upper_ci),
    Residuos = as.numeric(residuos)
  )

  # Datos Círculo Unitario
  df_raices <- NULL
  if (!is.null(inverse_roots)) {
    df_raices <- data.frame(Real = Re(inverse_roots), Imaginario = Im(inverse_roots))
  }

  # --- 6. Objeto Final ---
  resultado <- list(
    x_original = x,
    modelo = modelo_base,
    orden_q = q,
    coeficientes_tabla = tabla_coefs,
    invertibilidad = list(es_invertible = es_invertible, raices = inverse_roots),
    diagnosticos = list(lb = lb_test, shapiro = shapiro_test, aic = modelo_base$aic),
    data_plot = df_plot,
    data_raices = df_raices
  )

  class(resultado) <- "ma_greenreg"
  return(resultado)
}

#' Impresión Académica para Modelos de Media Móvil (MA)
#' @export
print.ma_greenreg <- function(x, ...) {
  cat("\n==========================================================\n")
  cat("      MODELO DE MEDIA MÓVIL MA(", x$orden_q, ") \n")
  cat("==========================================================\n\n")

  # --- 1. MODELO Y ECUACIÓN ---
  cat("--- 1. MODELO Y ECUACIÓN ---\n")
  theta <- round(x$modelo$coef[grepl("^ma", names(x$modelo$coef))], 3)
  mu <- round(x$modelo$coef["intercept"], 3)

  # Construcción visual de la ecuación de choques
  eq <- paste0("Y(t) = ", mu, " + ε(t)")
  if(length(theta) > 0) {
    for(i in 1:length(theta)) {
      eq <- paste0(eq, " + (", theta[i], " * ε(t-", i, "))")
    }
  }

  cat("Ecuación de Medias Móviles Estimada:\n", eq, "\n")
  cat("AIC (Criterio de Akaike):", round(x$diagnosticos$aic, 2), "\n")
  cat(" NOTA: El modelo MA(", x$orden_q, ") describe cómo los choques (ε) del pasado afectan al presente.\n\n")

  # --- 2. COEFICIENTES ESTIMADOS ---
  cat("--- 2. COEFICIENTES ESTIMADOS (Parámetros Theta) ---\n")
  stats::printCoefmat(x$coeficientes_tabla, digits = 4, signif.stars = TRUE)
  cat("\n GUÍA: Los coeficientes Theta miden la persistencia de los choques externos.\n\n")

  # --- 3. IDENTIFICACIÓN DEL MODELO (ORDEN q) ---
  cat("--- 3. IDENTIFICACIÓN DEL MODELO (Orden q) ---\n")
  cat("• Orden seleccionado: q =", x$orden_q, "\n")
  cat("• Criterio ACF: La autocorrelación debe cortarse abruptamente después del lag", x$orden_q, ".\n")
  cat("• Criterio PACF: Los coeficientes parciales deben decaer exponencialmente hacia cero.\n")
  cat(" NOTA: El modelo MA es el 'espejo' del AR. Aquí la ACF es la que manda para elegir q.\n\n")

  # --- 4. VERIFICACIÓN DE SUPUESTOS ---
  cat("--- 4. VERIFICACIÓN DE SUPUESTOS ---\n")

  # A) Estacionariedad e Invertibilidad
  cat("A) Estacionariedad e Invertibilidad:\n")
  cat("   • Estacionariedad: Un modelo MA puro es SIEMPRE estacionario por definición. ✅\n")
  if (!is.null(x$invertibilidad$raices)) {
    modulos <- Mod(x$invertibilidad$raices)
    cat("   • Módulos de raíces inversas (Invertibilidad):", paste(round(modulos, 3), collapse = ", "), "\n")
    if(x$invertibilidad$es_invertible) {
      cat("     -> ✅ OK: El modelo es Invertible (Raíces < 1). Los errores son estimables.\n")
    } else {
      cat("     -> ⚠️ ALERTA: Modelo NO Invertible. Los choques pasados tienen un peso infinito.\n")
    }
  }

  # B) Errores ~ Ruido Blanco
  cat("\nB) Diagnóstico de Errores (ε ~ RB):\n")
  # Media Cero
  mean_res <- mean(x$data_plot$Residuos, na.rm = TRUE)
  cat(sprintf("   • Media de los residuos: %.4f (Esperado: 0)\n", mean_res))

  # No Autocorrelación (Ljung-Box)
  lb_p <- x$diagnosticos$lb$p.value
  cat(sprintf("   • No Autocorrelación (Ljung-Box p-valor): %.4f\n", lb_p))
  if (lb_p > 0.05) {
    cat("     -> ✅ OK: Los residuos son ruido blanco. Toda la señal fue capturada.\n")
  } else {
    cat("     -> ⚠️ ALERTA: Queda estructura en los residuos. Considera aumentar el orden q.\n")
  }

  # Varianza y Normalidad
  sh_p <- x$diagnosticos$shapiro$p.value
  cat(sprintf("   • Normalidad (Shapiro p-valor): %s\n",
              if(is.na(sh_p)) "N/D" else round(sh_p, 4)))
  cat("     ->", if(!is.na(sh_p) && sh_p > 0.05) "✅ OK: Errores normales." else "⚠️ Errores no normales.\n")

  # --- 5. RESUMEN DEL MODELO ---
  cat("\n--- 5. RESUMEN DEL MODELO ---\n")
  if(x$invertibilidad$es_invertible && lb_p > 0.05) {
    cat("✅ RESULTADO FINAL: El modelo MA(", x$orden_q, ") es VÁLIDO y ÚNICO.\n", sep="")
  } else {
    cat("⚠️ RECOMENDACIÓN: Revisa el orden q o la invertibilidad para asegurar un ajuste correcto.\n")
  }

  cat("\n Usa plot(modelo) para visualizar la serie y la ACF de los residuos.\n")
  cat("==========================================================\n")
}



#' Gráficas Didácticas para Modelo MA
#' @export
plot.ma_greenreg <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Instala ggplot2.")

  df <- x$data_plot

  # TEMA ESTÁNDAR
  mi_tema <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, color = "#2C3E50"),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#7F8C8D"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 10, face = "italic", color = "#555555", margin = ggplot2::margin(t = 10))
    )

  cat("Generando gráficas de diagnóstico MA... (Presiona [Enter] para avanzar)\n")

  # --- GRÁFICA 1: SERIE ORIGINAL (VALORES VS TIEMPO) ---

  # Calculamos la media para la línea de referencia
  media_serie <- mean(df$Observado, na.rm = TRUE)

  g1 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Observado)) +
    # Línea de la serie original
    ggplot2::geom_line(color = "#34495E", size = 0.7) +
    # Puntos sutiles
    ggplot2::geom_point(color = "#34495E", alpha = 0.3, size = 1) +
    # Línea de la media histórica (Referencia fundamental en MA)
    ggplot2::geom_hline(yintercept = media_serie, linetype = "dashed",
                        color = "#2980B9", size = 0.8) +
    ggplot2::labs(
      title = "1. Comportamiento de la Serie Original",
      subtitle = "Visualización de la serie temporal bruta",
      x = "Tiempo (Periodos)",
      y = "Valor Observado",
      caption = paste0("INTERPRETACIÓN:\n",
                       "🔵 LÍNEA AZUL: Es el valor promedio histórico.\n",
                       "✅ MEDIA MÓVIL: El modelo MA supone que la serie siempre regresa a esta media.\n",
                       " NOTA: Si observas una tendencia marcada (subida o bajada constante),\n",
                       "   un modelo MA puro podría no ser suficiente sin diferenciar primero.")
    ) +
    mi_tema

  print(g1)
  readline(prompt = "Gráfica 1: Inspección de la serie pura > ")

  # --- GRÁFICA 2: AJUSTE DEL MODELO (Real vs Ajustado) ---

  g2 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo)) +
    # Intervalo de confianza (Sombra azul)
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper),
                         fill = "#2980B9", alpha = 0.15) +
    # Serie Observada (Realidad)
    ggplot2::geom_line(ggplot2::aes(y = Observado, color = "Observado"), size = 0.7) +
    # Serie Ajustada (Modelo MA)
    ggplot2::geom_line(ggplot2::aes(y = Ajustado, color = "Modelo MA"), size = 0.8) +
    # Configuración de colores
    ggplot2::scale_color_manual(values = c("Observado" = "#34495E", "Modelo MA" = "#2980B9")) +
    ggplot2::labs(
      title = paste0("2. Ajuste del Modelo MA(", x$orden_q, ")"),
      subtitle = "Serie Observada vs. Valores Ajustados (Fitted Values)",
      x = "Tiempo (Periodos)",
      y = "Valor",
      color = "Leyenda:",
      caption = paste0("INTERPRETACIÓN:\n",
                       "✅ BIEN: La línea azul debe capturar la tendencia central de la oscura.\n",
                       "🔵 SOMBRA AZUL: Margen de error esperado (IC 95%).\n",
                       " NOTA: El modelo MA suele ser más 'suave' que la serie original,\n",
                       "   ya que filtra el ruido aleatorio para quedarse con el promedio móvil.")
    ) +
    mi_tema

  print(g2)
  readline(prompt = "Gráfica 2: Comparación del modelo contra los datos reales > ")

  # --- GRÁFICA 3: PRECISIÓN DEL MODELO (REAL VS AJUSTADO) ---

  # Calculamos la correlación para el reporte visual
  cor_ma <- stats::cor(df$Observado, df$Ajustado, use = "complete.obs")

  g3 <- ggplot2::ggplot(df, ggplot2::aes(x = Ajustado, y = Observado)) +
    # Puntos de dispersión (Color azul MA)
    ggplot2::geom_point(color = "#2980B9", alpha = 0.5, size = 2) +
    # Línea de Identidad (El ideal de perfección)
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", size = 1) +
    # Línea de tendencia real con fórmula explícita para evitar avisos
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#2C3E50", fill = "#BDC3C7", alpha = 0.2) +
    ggplot2::labs(
      title = "3. Precisión del Modelo (Predicho vs Real)",
      subtitle = paste0("Correlación de Pearson: ", round(cor_ma, 3)),
      x = "Valor Ajustado (Modelo MA)",
      y = "Valor Observado (Realidad)",
      caption = paste0("INTERPRETACIÓN:\n",
                       "🔴 LÍNEA ROJA: Diagonal de 45° (Donde el modelo es perfecto).\n",
                       "✅ BIEN: Los puntos deben agruparse cerca de la línea roja.\n",
                       "⚠️ NOTA: Si la nube de puntos es muy ancha, el modelo MA tiene poca precisión.\n",
                       " CLAVE: Un valor de correlación cercano a 1 indica un ajuste excelente.")
    ) +
    mi_tema

  print(g3)
  readline(prompt = "Gráfica 3: Análisis de correlación y precisión. > ")

  # --- GRÁFICA 4: INVERTIBILIDAD (CÍRCULO UNITARIO) ---

  if (!is.null(x$data_raices)) {
    # 1. Crear el círculo de referencia
    theta <- seq(0, 2*pi, length.out = 200)
    df_circle <- data.frame(x = cos(theta), y = sin(theta))

    # 2. Construcción de la Gráfica con Estilo MA (Azul)
    g4 <- ggplot2::ggplot() +
      # Círculo unitario punteado
      ggplot2::geom_path(data = df_circle, ggplot2::aes(x, y),
                         color = "#7F8C8D", linetype = "dashed", size = 0.8) +
      # Ejes cartesianos
      ggplot2::geom_vline(xintercept = 0, color = "#BDC3C7") +
      ggplot2::geom_hline(yintercept = 0, color = "#BDC3C7") +
      # Área de Invertibilidad (Verde tenue)
      ggplot2::geom_polygon(data = df_circle, ggplot2::aes(x, y),
                            fill = "#27AE60", alpha = 0.05) +
      # Raíces inversas (Puntos azules para MA)
      ggplot2::geom_point(data = x$data_raices,
                          ggplot2::aes(x = Real, y = Imaginario),
                          color = "#2980B9", size = 4, shape = 18) +
      # Ajuste de proporciones
      ggplot2::coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
      ggplot2::labs(
        title = "4. Invertibilidad (Círculo Unitario)",
        subtitle = paste0("Raíces Inversas del Polinomio MA(", x$orden_q, ")"),
        x = "Parte Real",
        y = "Parte Imaginaria",
        caption = paste0("INTERPRETACIÓN:\n",
                         "✅ INVERTIBLE: Todas las raíces (rombos azules) están DENTRO del círculo.\n",
                         "⚠️ ALERTA: Si una raíz toca o sale del círculo, el modelo no es único.\n",
                         " NOTA: La invertibilidad asegura que los errores del modelo sean estimables\n",
                         "   a partir de los datos observados de forma convergente.")
      ) +
      mi_tema

    print(g4)

  } else {
    cat("ℹ️ El modelo MA(0) no tiene raíces; es ruido blanco y siempre es invertible.\n")
  }
  readline(prompt = "Gráfica 4: Análisis de raíces del proceso de media móvil > ")

  # --- GRÁFICA 5: RESIDUOS EN EL TIEMPO (ε vs t) ---

  # Calculamos la desviación estándar para las bandas de control
  sd_res_ma <- sd(df$Residuos, na.rm = TRUE)

  g5 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Residuos)) +
    # Banda de confianza (±2 desviaciones estándar) en color Morado/Azul
    ggplot2::geom_ribbon(ggplot2::aes(ymin = -2*sd_res_ma, ymax = 2*sd_res_ma),
                         fill = "#8E44AD", alpha = 0.1) +
    # Línea de los residuos
    ggplot2::geom_line(color = "#2C3E50", size = 0.5) +
    # Puntos para identificar errores individuales
    ggplot2::geom_point(color = "#8E44AD", alpha = 0.4, size = 1) +
    # Línea central (Media Cero)
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
    # Umbrales críticos
    ggplot2::geom_hline(yintercept = c(-2*sd_res_ma, 2*sd_res_ma),
                        linetype = "dotted", color = "#E74C3C", alpha = 0.6) +
    ggplot2::labs(
      title = "5. Residuos en el Tiempo",
      subtitle = "Análisis de Aleatoriedad (ε ~ Ruido Blanco)",
      x = "Tiempo (Periodos)",
      y = "Error (Residuo)",
      caption = paste0("INTERPRETACIÓN:\n",
                       "✅ BIEN: Los puntos deben oscilar al azar alrededor de la línea negra.\n",
                       "⚠️ BANDAS ROJAS: El 95% de los errores deben estar dentro de estos límites.\n",
                       "❌ MAL: Si ves ciclos o tendencias, el orden 'q' es insuficiente.\n",
                       " CLAVE: Si la variabilidad aumenta con el tiempo, el modelo es inestable.")
    ) +
    mi_tema

  print(g5)
  readline(prompt = "Gráfica 5: Verificación de la independencia temporal de los errores. > ")

  # --- GRÁFICA 6: AUTOCORRELACIÓN DE RESIDUOS (ACF) ---

  # 1. Cálculo de la ACF de los residuos (máximo 20 lags)
  acf_res_ma <- stats::acf(df$Residuos, plot = FALSE, lag.max = 20)

  # Creamos el dataframe omitiendo el Lag 0 (que siempre es 1)
  df_acf_ma <- data.frame(
    Lag = as.numeric(acf_res_ma$lag[-1]),
    ACF = as.numeric(acf_res_ma$acf[-1])
  )

  # Límite de significancia estadística (95%)
  n_obs_ma <- nrow(df)
  limite_ma <- 1.96 / sqrt(n_obs_ma)

  # Lógica de etiquetas rojas para alertar sobre fallos en lags específicos
  df_acf_ma$Etiqueta <- ifelse(abs(df_acf_ma$ACF) > limite_ma, as.character(df_acf_ma$Lag), "")
  df_acf_ma$vjust_pos <- ifelse(df_acf_ma$ACF > 0, -0.5, 1.5)

  # 2. Construcción de la Gráfica
  g6 <- ggplot2::ggplot(df_acf_ma, ggplot2::aes(x = Lag, y = ACF)) +
    # Barras de correlación (Color oscuro para contraste)
    ggplot2::geom_col(fill = "#34495E", width = 0.4) +
    # Líneas de significancia (Azules para MA)
    ggplot2::geom_hline(yintercept = c(limite_ma, -limite_ma),
                        linetype = "dashed", color = "#2980B9", size = 0.8) +
    # Línea base en cero
    ggplot2::geom_hline(yintercept = 0, color = "black") +
    # Etiquetas de error (Números Rojos)
    ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                       color = "#E74C3C", fontface = "bold", size = 3.5) +
    # Escala del eje X
    ggplot2::scale_x_continuous(breaks = 1:20) +
    ggplot2::labs(
      title = "6. Autocorrelación de Residuos (ACF)",
      subtitle = "¿Quedó memoria de choques sin capturar?",
      x = "Lag (Retraso temporal)",
      y = "Coeficiente ACF",
      caption = paste0("INTERPRETACIÓN:\n",
                       "✅ BIEN: Si todas las barras están dentro de las líneas azules.\n",
                       "⚠️ NÚMEROS ROJOS: Indican que el error en ese lag NO es aleatorio.\n",
                       "❌ MAL: Si el Lag 1 es rojo, tu modelo MA(", x$orden_q, ") es insuficiente.\n",
                       " CLAVE: Los residuos deben ser 'Ruido Blanco' (puro azar).")
    ) +
    mi_tema

  print(g6)
  readline(prompt = "Gráfica 6: Diagnóstico de independencia de los errores. > ")

  # --- GRÁFICA 7: AUTOCORRELACIÓN PARCIAL (PACF) DE RESIDUOS ---

  # 1. Cálculo de la PACF de los residuos
  pacf_res_ma <- stats::pacf(df$Residuos, plot = FALSE, lag.max = 20)

  df_pacf_ma <- data.frame(
    Lag = as.numeric(pacf_res_ma$lag),
    PACF = as.numeric(pacf_res_ma$acf)
  )

  # Límite de significancia (Mismo que en ACF)
  n_obs_ma <- nrow(df)
  limite_ma <- 1.96 / sqrt(n_obs_ma)

  # Lógica de etiquetas rojas para fallas
  df_pacf_ma$Etiqueta <- ifelse(abs(df_pacf_ma$PACF) > limite_ma, as.character(df_pacf_ma$Lag), "")
  df_pacf_ma$vjust_pos <- ifelse(df_pacf_ma$PACF > 0, -0.5, 1.5)

  # 2. Construcción de la Gráfica con Estilo MA (Azul/Naranja)
  g7 <- ggplot2::ggplot(df_pacf_ma, ggplot2::aes(x = Lag, y = PACF)) +
    # Barras de la PACF (Color Azul MA)
    ggplot2::geom_col(fill = "#2980B9", width = 0.4) +
    # Líneas de significancia (Cambiamos a Naranja para diferenciar de ACF)
    ggplot2::geom_hline(yintercept = c(limite_ma, -limite_ma),
                        linetype = "dashed", color = "#E67E22", size = 0.8) +
    # Línea base
    ggplot2::geom_hline(yintercept = 0, color = "black") +
    # Etiquetas de error (Chivatos en Rojo)
    ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                       color = "#E74C3C", fontface = "bold", size = 3.5) +
    # Escala del eje X
    ggplot2::scale_x_continuous(breaks = 1:20) +
    ggplot2::labs(
      title = "7. Autocorrelación Parcial (PACF)",
      subtitle = "Identificación de dependencias directas en residuos",
      x = "Lag (Retraso temporal)",
      y = "Coeficiente PACF",
      caption = paste0("INTERPRETACIÓN:\n",
                       "✅ BIEN: Ninguna barra sobresale de las líneas naranjas.\n",
                       "⚠️ NÚMEROS ROJOS: Indican que ese retraso específico aún tiene información.\n",
                       " NOTA: En un modelo MA bien ajustado, tanto ACF como PACF de residuos \n",
                       "   deberían verse como Ruido Blanco (barras cortas).")
    ) +
    mi_tema

  print(g7)
  readline(prompt = "Gráfica 7: Verificación de la limpieza estructural del error > ")

  # --- GRÁFICA 8: NORMALIDAD DE RESIDUOS (Q-Q PLOT) ---

  g8 <- ggplot2::ggplot(df, ggplot2::aes(sample = Residuos)) +
    # Puntos de los cuantiles observados (Color Naranja para MA)
    ggplot2::stat_qq(color = "#E67E22", alpha = 0.6, size = 2) +
    # Línea de referencia teórica (Normalidad perfecta)
    ggplot2::stat_qq_line(color = "#2C3E50", linetype = "dashed", size = 1) +
    ggplot2::labs(
      title = "8. Normalidad de los Errores (Q-Q Plot)",
      subtitle = "Comparación de Cuantiles de Residuos vs. Normal Teórica",
      x = "Cuantiles Teóricos (Distribución Normal)",
      y = "Cuantiles Observados (Residuos MA)",
      caption = paste0("INTERPRETACIÓN:\n",
                       "✅ BIEN: Los puntos deben seguir la línea diagonal punteada.\n",
                       "⚠️ EXTREMOS: Si los puntos se alejan en las puntas, hay valores atípicos.\n",
                       " NOTA: La normalidad valida que el error típico (Sigma) sea constante.\n",
                       " CLAVE: Si los puntos forman una curva, la serie tiene asimetría.")
    ) +
    mi_tema

  print(g8)
  readline(prompt = "Gráfica 8: Evaluación de la distribución gaussiana del error > ")

  # --- GRÁFICA 9: DIAGNÓSTICO DE ESTACIONARIEDAD (MÓVIL) ---

  # 1. Definimos la ventana (10% de los datos o mínimo 5)
  ventana_ma <- max(5, round(nrow(df) * 0.1))

  df_est_ma <- df

  # 2. CÁLCULOS (Es vital forzar as.numeric para ggplot2)
  df_est_ma$Media_Movil <- as.numeric(stats::filter(df$Observado, rep(1/ventana_ma, ventana_ma), sides = 2))

  df_est_ma$SD_Movil <- sapply(1:nrow(df), function(i) {
    idx <- max(1, i-ventana_ma):min(nrow(df), i+ventana_ma)
    sd(df$Observado[idx], na.rm = TRUE)
  })

  # 3. FILTRADO (Eliminamos los NAs antes de entrar a ggplot)
  df_final_ma <- df_est_ma[!is.na(df_est_ma$Media_Movil), ]

  # 4. CONSTRUCCIÓN DE LA GRÁFICA
  g9 <- ggplot2::ggplot(df_final_ma, ggplot2::aes(x = Tiempo)) +
    # Serie original completa de fondo (usando df original para no perder datos visuales)
    ggplot2::geom_line(data = df, ggplot2::aes(y = Observado), color = "gray80", size = 0.5) +

    # Franja de Varianza Móvil (Azul MA)
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Media_Movil - SD_Movil,
                                      ymax = Media_Movil + SD_Movil),
                         fill = "#2980B9", alpha = 0.15) +

    # Línea de la Media Móvil (Nivel local)
    ggplot2::geom_line(ggplot2::aes(y = Media_Movil), color = "#2980B9", size = 1.2) +

    # Línea de la Media Global (Referencia fija en Rojo)
    ggplot2::geom_hline(yintercept = mean(df$Observado, na.rm = TRUE),
                        linetype = "dashed", color = "#E74C3C", size = 0.8) +

    ggplot2::labs(
      title = "9. Diagnóstico de Estacionariedad",
      subtitle = paste0("Análisis de Media y Varianza Móvil (Ventana k = ", ventana_ma, ")"),
      x = "Tiempo (Periodos)",
      y = "Valor de la Serie",
      caption = paste0("INTERPRETACIÓN:\n",
                       "🔵 LÍNEA AZUL: Media local móvil. Debe ser horizontal y estable.\n",
                       "🔴 LÍNEA ROJA: Media global (El objetivo del modelo MA).\n",
                       "✅ BIEN: Si la franja azul mantiene su ancho y nivel constante.\n",
                       " NOTA: Se han filtrado los extremos para evitar avisos de valores nulos.")
    ) +
    mi_tema

  print(g9)
  cat("\n Gráfica 9: Verificación final de estabilidad estadística.\n")
}


