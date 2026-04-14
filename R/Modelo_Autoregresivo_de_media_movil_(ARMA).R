#' Modelo Autorregresivo de Media Móvil (ARMA)
#'
#' Esta función ajusta un modelo híbrido ARMA(p, q) que combina la dependencia
#' de los valores pasados (AutoRegresivo) con la influencia de los errores
#' pasados (Media Móvil). Es la herramienta estándar para modelar series de
#' tiempo estacionarias con estructuras de dependencia complejas.
#'
#' @details
#' El flujo interno de la función realiza una auditoría estructural en cinco pasos:
#' \enumerate{
#'   \item **Estimación Conjunta:** Utiliza Máxima Verosimilitud (MLE) para encontrar
#'         simultáneamente los coeficientes \eqn{\phi} (AR) y \eqn{\theta} (MA).
#'   \item **Estabilidad (Parte AR):** Verifica que el componente autorregresivo no sea
#'         explosivo. Las raíces inversas del polinomio AR deben estar dentro del círculo.
#'   \item **Invertibilidad (Parte MA):** Verifica que el componente de media móvil sea
#'         unívoco. Las raíces inversas del polinomio MA deben estar dentro del círculo.
#'   \item **Diagnóstico de Residuos:** Evalúa si el error resultante es Ruido Blanco
#'         mediante el test de Ljung-Box, asegurando que no quede señal sin capturar.
#'   \item **Criterio de Información:** Proporciona el AIC para facilitar la comparación
#'         entre diferentes combinaciones de (p, q), buscando el equilibrio entre
#'         complejidad y precisión.
#' }
#'
#' @section Guía de las 9 Gráficas de Diagnóstico (`plot`):
#' Al ejecutar \code{plot(modelo)}, se despliega una secuencia de validación mixta:
#' \describe{
#'   \item{\strong{1. Inspección de Serie}}{Muestra la serie bruta con su media. Fundamental para confirmar visualmente la estacionariedad.}
#'   \item{\strong{2. Ajuste Mixto}}{Superpone el modelo (verde) sobre la realidad. Permite evaluar cómo interactúan la inercia y los shocks.}
#'   \item{\strong{3. Residuos en el Tiempo}}{Evaluación de la aleatoriedad del error. Los puntos deben formar una nube sin patrones.}
#'   \item{\strong{4. Círculo Unitario Dual}}{Gráfica crítica exclusiva de ARMA. Muestra simultáneamente la Estabilidad (Puntos Rojos) e Invertibilidad (Rombos Azules).}
#'   \item{\strong{5. ACF de Residuos}}{Busca autocorrelaciones remanentes. Las barras deben ser mínimas para confirmar que el modelo es óptimo.}
#'   \item{\strong{6. PACF de Residuos}}{Ayuda a identificar si el orden autorregresivo (p) es suficiente.}
#'   \item{\strong{7. Normalidad (Q-Q Plot)}}{Valida que los errores sigan una campana de Gauss, requisito para la validez de los intervalos de confianza.}
#'   \item{\strong{8. Precisión (Real vs Ajustado)}}{Muestra la correlación directa. Un ARMA exitoso suele tener una alineación estrecha sobre la diagonal.}
#'   \item{\strong{9. Estacionariedad Móvil}}{Prueba final de estabilidad. Una franja verde horizontal confirma que el modelo no se desplaza en el tiempo.}
#' }
#'
#' @param x Vector numérico o serie de tiempo (\code{ts}).
#' @param p Orden Autorregresivo (pasados de la variable).
#' @param q Orden de Media Móvil (pasados del error).
#' @param include_mean Lógico. ¿Incluir constante/intercepto? Por defecto \code{TRUE}.
#'
#' @return Un objeto de clase \code{"arma_greenreg"} que contiene los coeficientes
#'         mixtos, el estatus de estabilidad/invertibilidad, diagnósticos de residuos
#'         y coordenadas para el círculo unitario doble.
#'
#' @examples
#' # data("datos_anomalia_temperatura")
#' # ts_data <- datos_anomalia_temperatura$anomalia_c
#' # modelo <- modelo_arma(ts_data, p = 1, q = 1)
#' # modelo
#' # plot(modelo)
#'
#' @importFrom stats arima coef residuals fitted Box.test pnorm qnorm shapiro.test acf pacf filter sd
#' @import ggplot2
#' @export
modelo_arma <- function(x, p = 1, q = 1, include_mean = TRUE) {

  # --- 1. Validación y Ajuste ---
  if (!is.numeric(x)) stop("El argumento 'x' debe ser numérico.")

  # Ajuste ARIMA(p, 0, q)
  modelo_base <- stats::arima(x, order = c(p, 0, q), include.mean = include_mean)

  # --- 2. Coeficientes y Estadísticas ---
  coefs <- modelo_base$coef
  se <- sqrt(diag(modelo_base$var.coef))
  z_val <- coefs / se
  p_val <- 2 * (1 - stats::pnorm(abs(z_val)))

  tabla_coefs <- cbind(Estimate = coefs, `Std. Error` = se, `z value` = z_val, `Pr(>|z|)` = p_val)

  # --- 3. Análisis de Raíces (Doble Verificación) ---

  # A) Estabilidad (Parte AR) - Raíces en ROJO
  ar_coefs <- coefs[grepl("^ar", names(coefs))]
  raices_ar_inv <- NULL
  es_estable <- TRUE

  if (length(ar_coefs) > 0) {
    poly_roots_ar <- polyroot(c(1, -ar_coefs))
    raices_ar_inv <- 1 / poly_roots_ar
    es_estable <- all(Mod(raices_ar_inv) < 1)
  }

  # B) Invertibilidad (Parte MA) - Raíces en AZUL
  ma_coefs <- coefs[grepl("^ma", names(coefs))]
  raices_ma_inv <- NULL
  es_invertible <- TRUE

  if (length(ma_coefs) > 0) {
    poly_roots_ma <- polyroot(c(1, ma_coefs))
    raices_ma_inv <- 1 / poly_roots_ma
    es_invertible <- all(Mod(raices_ma_inv) < 1)
  }

  # --- 4. Diagnósticos de Residuos ---
  residuos <- stats::residuals(modelo_base)

  # Ljung-Box (Ruido Blanco)
  lb_test <- stats::Box.test(residuos, type = "Ljung-Box", lag = p + q + 5)

  # Shapiro-Wilk (Normalidad)
  if (length(residuos) >= 3 && length(residuos) <= 5000) {
    shapiro_test <- stats::shapiro.test(residuos)
  } else {
    shapiro_test <- list(p.value = NA)
  }

  # --- 5. Datos para Gráficas ---
  ajustados <- x - residuos
  sigma <- sqrt(modelo_base$sigma2)
  # Intervalos 95%
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

  # Datos combinados para círculo unitario
  df_raices <- data.frame(Real = numeric(), Imaginario = numeric(), Tipo = character())
  if (!is.null(raices_ar_inv)) {
    df_raices <- rbind(df_raices, data.frame(Real = Re(raices_ar_inv), Imaginario = Im(raices_ar_inv), Tipo = "AR (Estabilidad)"))
  }
  if (!is.null(raices_ma_inv)) {
    df_raices <- rbind(df_raices, data.frame(Real = Re(raices_ma_inv), Imaginario = Im(raices_ma_inv), Tipo = "MA (Invertibilidad)"))
  }

  # --- 6. Objeto Final ---
  resultado <- list(
    modelo = modelo_base,
    ordenes = c(p = p, q = q),
    coeficientes_tabla = tabla_coefs,
    estatus = list(estable = es_estable, invertible = es_invertible),
    diagnosticos = list(lb = lb_test, shapiro = shapiro_test, aic = modelo_base$aic),
    data_plot = df_plot,
    data_raices = df_raices
  )

  class(resultado) <- "arma_greenreg"
  return(resultado)
}

#' Impresión para Modelo ARMA
#' @export
print.arma_greenreg <- function(x, ...) {
  cat("\n==========================================================\n")
  cat(" Modelo Mixto ARMA(", x$ordenes['p'], ",", x$ordenes['q'], ")  \n")
  cat("==========================================================\n\n")

  cat("--- 1. COEFICIENTES ESTIMADOS ---\n")
  stats::printCoefmat(x$coeficientes_tabla, digits = 4, signif.stars = TRUE)
  cat("\nAIC:", round(x$diagnosticos$aic, 2), "(Menor es mejor)\n\n")

  cat("--- 2. VERIFICACIÓN ESTRUCTURAL (Raíces Unitarias) ---\n")

  # Check AR
  if (x$ordenes['p'] > 0) {
    if (x$estatus$estable) cat("✅ Parte AR: ESTABLE (Raíces dentro del círculo).\n")
    else cat("❌ Parte AR: INESTABLE (Explosivo).\n")
  } else { cat("ℹ️ Parte AR: No aplica (p=0).\n") }

  # Check MA
  if (x$ordenes['q'] > 0) {
    if (x$estatus$invertible) cat("✅ Parte MA: INVERTIBLE (Modelo único).\n")
    else cat("❌ Parte MA: NO INVERTIBLE (Problemas de identificación).\n")
  } else { cat("ℹ️ Parte MA: No aplica (q=0).\n") }
  cat("\n")

  cat("--- 3. DIAGNÓSTICO DE RESIDUOS ---\n")
  lb_p <- x$diagnosticos$lb$p.value
  cat(sprintf("• Independencia (Ljung-Box): p-valor = %.4f\n", lb_p))

  if (lb_p > 0.05) {
    cat("✅ Ruido Blanco: El modelo capturó la dinámica temporal.\n")
  } else {
    cat("⚠️ Autocorrelación: El modelo ARMA actual es insuficiente.\n")
  }

  cat("\n Graficas de modelo series de tiempo modelo (ARMA).\n")
  cat("==========================================================\n")
}

#' Gráficas para Modelo ARMA
#' @export
plot.arma_greenreg <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Instala ggplot2.")

  df <- x$data_plot
  p  <- x$ordenes['p'] # Extraemos p del objeto x
  q  <- x$ordenes['q'] # Extraemos q del objeto x

  # TEMA ESTÁNDAR
  mi_tema <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, color = "#2C3E50"),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#7F8C8D"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 10, face = "italic", color = "#555555", margin = ggplot2::margin(t = 10))
    )

  cat("Generando gráficas ARMA... (Presiona [Enter] para avanzar)\n")

  # --- GRÁFICA 1: SERIE ORIGINAL (VALORES VS TIEMPO) ---

  media_global <- mean(df$Observado, na.rm = TRUE)

  g1 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Observado)) +
    # Línea de la serie
    ggplot2::geom_line(color = "#34495E", size = 0.7) +
    # Puntos de datos
    ggplot2::geom_point(color = "#34495E", alpha = 0.3, size = 1) +
    # Línea de la media (Referencia de Estacionariedad)
    ggplot2::geom_hline(yintercept = media_global, linetype = "dashed",
                        color = "#16A085", size = 0.8) +
    ggplot2::labs(
      title = "1. Comportamiento de la Serie Original",
      subtitle = "Visualización de la serie temporal antes del ajuste mixto",
      x = "Tiempo (Periodos)",
      y = "Valor Observado",
      caption = paste0("INTERPRETACIÓN:\n",
                       "🟢 LÍNEA VERDE: Promedio histórico de la serie.\n",
                       "✅ OBJETIVO: El modelo ARMA asume que la serie fluctúa alrededor\n",
                       "   de este nivel constante (Estacionariedad en media).\n",
                       " NOTA: Si ves una 'escalera' o tendencia, el modelo ARMA(p,q)\n",
                       "   podría requerir una diferencia previa (Modelo ARIMA).")
    ) +
    mi_tema

  print(g1)
  readline(prompt = "Gráfica 1: Inspección de los datos brutos) > ")


# --- GRÁFICA 2: AJUSTE DEL MODELO ARMA (Real vs Ajustado) ---

g2 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo)) +
  # Sombra de Confianza (Intervalo del 95%)
  ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper),
                       fill = "#16A085", alpha = 0.15) +
  # Serie Observada (Datos reales)
  ggplot2::geom_line(ggplot2::aes(y = Observado, color = "Observado"), size = 0.7) +
  # Serie Ajustada (Predicción del modelo ARMA)
  ggplot2::geom_line(ggplot2::aes(y = Ajustado, color = "Modelo ARMA"), size = 0.8) +
  # Escala de colores personalizada
  ggplot2::scale_color_manual(values = c("Observado" = "#34495E", "Modelo ARMA" = "#16A085")) +
  ggplot2::labs(
    title = paste0("2. Ajuste del Modelo ARMA(", p, ",", q, ")"),
    subtitle = "Comparación de valores reales vs. estimación mixta",
    x = "Tiempo (Periodos)",
    y = "Valor",
    color = "Referencia:",
    caption = paste0("INTERPRETACIÓN:\n",
                     "✅ BIEN: La línea verde debe 'abrazar' los movimientos de la oscura.\n",
                     "🟢 SOMBRA: Zona de confianza. El 95% de los datos deben caer aquí.\n",
                     " NOTA: El ARMA es más potente que el AR o MA por separado,\n",
                     "   logrando un equilibrio entre tendencia y choques aleatorios.")
  ) +
  mi_tema

print(g2)
readline(prompt = "Gráfica 2: Evaluación visual del ajuste mixto > ")

# --- GRÁFICA 3: RESIDUOS (ERRORES) VS TIEMPO ---

# Calculamos la desviación estándar de los residuos para las bandas de control
sd_res_arma <- sd(df$Residuos, na.rm = TRUE)

g3 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Residuos)) +
  # Banda de confianza (±2 SD) - Representa el 95% de probabilidad
  ggplot2::geom_ribbon(ggplot2::aes(ymin = -2*sd_res_arma, ymax = 2*sd_res_arma),
                       fill = "#16A085", alpha = 0.1) +
  # Línea de los residuos en el tiempo
  ggplot2::geom_line(color = "#2C3E50", size = 0.5) +
  # Puntos para identificar shocks específicos
  ggplot2::geom_point(color = "#16A085", alpha = 0.4, size = 1) +
  # Línea central en Cero (Media esperada)
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
  # Líneas críticas de control
  ggplot2::geom_hline(yintercept = c(-2*sd_res_arma, 2*sd_res_arma),
                      linetype = "dotted", color = "#E74C3C", alpha = 0.7) +
  ggplot2::labs(
    title = "3. Análisis de Residuos (ε vs t)",
    subtitle = "Verificación de Aleatoriedad y Media Cero",
    x = "Tiempo (Periodos)",
    y = "Residuo (Error)",
    caption = paste0("INTERPRETACIÓN:\n",
                     "✅ BIEN: Los errores deben ser una 'nube' aleatoria centrada en cero.\n",
                     "⚠️ ROJO: Puntos fuera de estas líneas son shocks no explicados (Outliers).\n",
                     "❌ MAL: Si ves ciclos, 'olas' o tendencias, el modelo no es suficiente.\n",
                     " NOTA: El ARMA busca que el error de hoy no dependa del error de ayer.")
  ) +
  mi_tema

print(g3)
readline(prompt = "Gráfica 3: Inspección de la calidad del error residual > ")

# --- GRÁFICA 4: ESTRUCTURA (CÍRCULO UNITARIO AR Y MA) ---

if (nrow(x$data_raices) > 0) {
  # 1. Crear el círculo de referencia
  theta_circ <- seq(0, 2*pi, length.out = 200)
  df_circle_arma <- data.frame(x = cos(theta_circ), y = sin(theta_circ))

  # 2. Construcción de la Gráfica Dual
  g4 <- ggplot2::ggplot() +
    # Círculo unitario de referencia
    ggplot2::geom_path(data = df_circle_arma, ggplot2::aes(x, y),
                       color = "#7F8C8D", linetype = "dashed", size = 0.8) +
    # Área de validez (Sombreado muy tenue)
    ggplot2::geom_polygon(data = df_circle_arma, ggplot2::aes(x, y),
                          fill = "#BDC3C7", alpha = 0.1) +
    # Ejes cartesianos
    ggplot2::geom_vline(xintercept = 0, color = "gray80") +
    ggplot2::geom_hline(yintercept = 0, color = "gray80") +
    # Graficar Raíces AR (Rojas) y MA (Azules)
    ggplot2::geom_point(data = x$data_raices,
                        ggplot2::aes(x = Real, y = Imaginario, color = Tipo, shape = Tipo),
                        size = 4, stroke = 1.2) +
    # Definición de colores y formas para distinguir AR de MA
    ggplot2::scale_color_manual(values = c("AR (Estabilidad)" = "#E74C3C",
                                           "MA (Invertibilidad)" = "#2980B9")) +
    ggplot2::scale_shape_manual(values = c("AR (Estabilidad)" = 16,
                                           "MA (Invertibilidad)" = 18)) +
    # Ajuste de proporciones
    ggplot2::coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
    ggplot2::labs(
      title = "4. Estructura de Raíces (Círculo Unitario)",
      subtitle = paste0("Diagnóstico Dual: ARMA(", p, ",", q, ")"),
      x = "Parte Real",
      y = "Parte Imaginaria",
      color = "Componente:", shape = "Componente:",
      caption = paste0("INTERPRETACIÓN:\n",
                       "🔴 PUNTOS ROJOS (AR): Si están dentro, el modelo es ESTABLE.\n",
                       "🔵 ROMBOS AZULES (MA): Si están dentro, el modelo es INVERTIBLE.\n",
                       "⚠️ CRÍTICO: Si algún punto toca o sale del círculo gris,\n",
                       "   los resultados estadísticos NO son confiables.")
    ) +
    mi_tema + ggplot2::theme(legend.position = "right")

  print(g4)

} else {
  cat("ℹ️ El modelo ARMA(0,0) no tiene raíces; es ruido blanco.\n")
}
readline(prompt = "Gráfica 4: Verificación de estabilidad e invertibilidad mixta. > ")

# --- GRÁFICA 5: AUTOCORRELACIÓN DE RESIDUOS (ACF) ---

# 1. Cálculo de la ACF de los residuos (excluyendo Lag 0)
acf_res_arma <- stats::acf(df$Residuos, plot = FALSE, lag.max = 20)
df_acf_arma <- data.frame(
  Lag = as.numeric(acf_res_arma$lag[-1]),
  ACF = as.numeric(acf_res_arma$acf[-1])
)

# Límite de significancia estadística (95%)
limite_arma <- 1.96 / sqrt(nrow(df))

# Lógica de detección: Marcar en rojo los fallos
df_acf_arma$Etiqueta <- ifelse(abs(df_acf_arma$ACF) > limite_arma, as.character(df_acf_arma$Lag), "")
df_acf_arma$vjust_pos <- ifelse(df_acf_arma$ACF > 0, -0.5, 1.5)

# 2. Construcción de la Gráfica
g5 <- ggplot2::ggplot(df_acf_arma, ggplot2::aes(x = Lag, y = ACF)) +
  # Barras de correlación
  ggplot2::geom_col(fill = "#34495E", width = 0.4) +
  # Líneas de significancia (Color Verde Esmeralda para ARMA)
  ggplot2::geom_hline(yintercept = c(limite_arma, -limite_arma),
                      linetype = "dashed", color = "#16A085", size = 0.8) +
  # Línea base
  ggplot2::geom_hline(yintercept = 0, color = "black") +
  # Números rojos de alerta
  ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                     color = "#E74C3C", fontface = "bold", size = 3.5) +
  ggplot2::scale_x_continuous(breaks = 1:20) +
  ggplot2::labs(
    title = "5. Autocorrelación de Residuos (ACF)",
    subtitle = "¿Quedó estructura temporal en el error?",
    x = "Lag (Retraso)",
    y = "Coeficiente ACF",
    caption = paste0("INTERPRETACIÓN:\n",
                     "✅ BIEN: Si todas las barras están 'atrapadas' en las líneas verdes.\n",
                     "⚠️ NÚMEROS ROJOS: Indican qué lags NO son ruido blanco.\n",
                     " NOTA: El ARMA busca 'limpiar' la ACF. Si el Lag 1 o 2 es rojo,\n",
                     "   el modelo mixto aún no es óptimo.")
  ) +
  mi_tema

print(g5)
readline(prompt = "Gráfica 5: Diagnóstico de independencia de los errores mixtos > ")

# --- GRÁFICA 6: AUTOCORRELACIÓN PARCIAL (PACF) DE RESIDUOS ---

# 1. Cálculo de la PACF de los residuos (Máximo 20 lags)
pacf_res_arma <- stats::pacf(df$Residuos, plot = FALSE, lag.max = 20)

df_pacf_arma <- data.frame(
  Lag = as.numeric(pacf_res_arma$lag),
  PACF = as.numeric(pacf_res_arma$acf)
)

# Límite de significancia (Mismo que en ACF)
limite_arma <- 1.96 / sqrt(nrow(df))

# Lógica de detección: Etiquetas rojas para los fallos directos
df_pacf_arma$Etiqueta <- ifelse(abs(df_pacf_arma$PACF) > limite_arma, as.character(df_pacf_arma$Lag), "")
df_pacf_arma$vjust_pos <- ifelse(df_pacf_arma$PACF > 0, -0.5, 1.5)

# 2. Construcción de la Gráfica con Estilo ARMA (Verde/Naranja)
g6 <- ggplot2::ggplot(df_pacf_arma, ggplot2::aes(x = Lag, y = PACF)) +
  # Barras de la PACF (Color Verde Esmeralda)
  ggplot2::geom_col(fill = "#16A085", width = 0.4) +
  # Líneas de significancia (Naranja para contraste preventivo)
  ggplot2::geom_hline(yintercept = c(limite_arma, -limite_arma),
                      linetype = "dashed", color = "#E67E22", size = 0.8) +
  # Línea base en cero
  ggplot2::geom_hline(yintercept = 0, color = "black") +
  # Etiquetas de error (Chivatos en Rojo)
  ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                     color = "#E74C3C", fontface = "bold", size = 3.5) +
  # Escala del eje X
  ggplot2::scale_x_continuous(breaks = 1:20) +
  ggplot2::labs(
    title = "6. Autocorrelación Parcial (PACF)",
    subtitle = "Identificación de dependencias directas en errores mixtos",
    x = "Lag (Retraso)",
    y = "Coeficiente PACF",
    caption = paste0("INTERPRETACIÓN:\n",
                     "✅ BIEN: Ninguna barra sobresale de las líneas naranjas.\n",
                     "⚠️ NÚMEROS ROJOS: Indican que ese lag específico aún tiene 'señal'.\n",
                     " NOTA: Si el Lag 1 es significativo, el componente AR(", p, ") \n",
                     "   no está capturando bien la inercia inmediata.")
  ) +
  mi_tema

 print(g6)
 readline(prompt = "Gráfica 6: Auditoría final de la estructura del error > ")


 # --- GRÁFICA 7: NORMALIDAD DE RESIDUOS (Q-Q PLOT) ---

 g7 <- ggplot2::ggplot(df, ggplot2::aes(sample = Residuos)) +
   # Puntos de los cuantiles observados (Color Verde ARMA)
   ggplot2::stat_qq(color = "#16A085", alpha = 0.6, size = 2) +
   # Línea de referencia teórica
   ggplot2::stat_qq_line(color = "#2C3E50", linetype = "dashed", size = 1) +
   ggplot2::labs(
     title = "7. Normalidad de los Errores (Q-Q Plot)",
     subtitle = "Comparación de Cuantiles de Residuos Mixtos vs. Normal Teórica",
     x = "Cuantiles Teóricos (Normal)",
     y = "Cuantiles Observados (Residuos ARMA)",
     caption = paste0("INTERPRETACIÓN:\n",
                      "✅ BIEN: Los puntos deben seguir la línea diagonal punteada.\n",
                      "⚠️ COLAS: Si los extremos se curvan, hay eventos extremos (outliers).\n",
                      " NOTA: La normalidad asegura que el error típico sea consistente.\n",
                      " CLAVE: Si los puntos se alejan mucho, los p-valores del modelo son dudosos.")
   ) +
   mi_tema

 print(g7)
 readline(prompt = "Gráfica 7: Verificación de la distribución gaussiana del error > ")

 # --- GRÁFICA 8: PRECISIÓN DEL MODELO (REAL VS AJUSTADO) ---

 # Calculamos la correlación de Pearson para el reporte
 cor_arma <- stats::cor(df$Observado, df$Ajustado, use = "complete.obs")

 g8 <- ggplot2::ggplot(df, ggplot2::aes(x = Ajustado, y = Observado)) +
   # Puntos de dispersión (Color Verde ARMA)
   ggplot2::geom_point(color = "#16A085", alpha = 0.5, size = 2) +
   # Línea de Identidad (Puntos donde Real = Modelo)
   ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", size = 1) +
   # Línea de tendencia con fórmula explícita
   ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#2C3E50", fill = "#BDC3C7", alpha = 0.2) +
   ggplot2::labs(
     title = "8. Precisión del Modelo (Predicho vs Real)",
     subtitle = paste0("Correlación de Pearson: ", round(cor_arma, 3)),
     x = "Valor Ajustado (Modelo ARMA)",
     y = "Valor Observado (Realidad)",
     caption = paste0("INTERPRETACIÓN:\n",
                      "🔴 LÍNEA ROJA: El ideal donde el modelo predice exactamente la realidad.\n",
                      "✅ BIEN: Los puntos verdes deben seguir la trayectoria de la línea roja.\n",
                      "⚠️ DISPERSIÓN: Si los puntos están muy lejos de la línea, el modelo es pobre.\n",
                      " CLAVE: Un ARMA exitoso debe tener una correlación superior a 0.7.")
   ) +
   mi_tema

 print(g8)
 readline(prompt = "Gráfica 8: Análisis de correlación y capacidad predictiva. > ")


 # --- GRÁFICA 9: DIAGNÓSTICO DE ESTACIONARIEDAD (MÓVIL) ---

 # 1. Definir ventana (10% de los datos)
 ventana_arma <- max(5, round(nrow(df) * 0.1))

 df_est_arma <- df

 # 2. CÁLCULOS (Primero calculamos para poder filtrar después)
 # Media móvil simple
 df_est_arma$Media_Movil <- as.numeric(stats::filter(df$Observado, rep(1/ventana_arma, ventana_arma), sides = 2))

 # Desviación estándar móvil
 df_est_arma$SD_Movil <- sapply(1:nrow(df), function(i) {
   idx <- max(1, i-ventana_arma):min(nrow(df), i+ventana_arma)
   sd(df$Observado[idx], na.rm = TRUE)
 })

 # 3. FILTRO (Ahora sí, eliminamos los NAs generados por el filtro en los extremos)
 df_final_grafica <- df_est_arma[!is.na(df_est_arma$Media_Movil), ]

 # 4. CONSTRUCCIÓN DE LA GRÁFICA (Usando el DF filtrado)
 g9 <- ggplot2::ggplot(df_final_grafica, ggplot2::aes(x = Tiempo)) +
   # Serie original (Fondo gris)
   ggplot2::geom_line(data = df, ggplot2::aes(y = Observado), color = "gray85", size = 0.5) +
   # Franja de Varianza Móvil (Verde ARMA)
   ggplot2::geom_ribbon(ggplot2::aes(ymin = Media_Movil - SD_Movil,
                                     ymax = Media_Movil + SD_Movil),
                        fill = "#16A085", alpha = 0.15) +
   # Línea de la Media Móvil
   ggplot2::geom_line(ggplot2::aes(y = Media_Movil), color = "#16A085", size = 1.2) +
   # Línea de la Media Global (Referencia fija)
   ggplot2::geom_hline(yintercept = mean(df$Observado, na.rm = TRUE),
                       linetype = "dashed", color = "#E74C3C", size = 0.8) +
   ggplot2::labs(
     title = "9. Diagnóstico de Estacionariedad",
     subtitle = paste0("Análisis de Media y Varianza Móvil (Ventana k = ", ventana_arma, ")"),
     x = "Tiempo (Periodos)",
     y = "Valor de la Serie",
     caption = paste0("INTERPRETACIÓN:\n",
                      "🟢 LÍNEA VERDE: Media local. Debe ser horizontal y estable.\n",
                      "🔴 LÍNEA ROJA: Promedio global. La base del modelo.\n",
                      "✅ BIEN: Si la franja verde mantiene un ancho constante.\n",
                      " NOTA: Se han omitido los extremos para evitar errores de cálculo móvil.")
   ) +
   mi_tema

 print(g9)
 cat("Recuerda que un buen modelo ARMA es el equilibrio entre parsimonia y precisión.\n")
  cat("✅ Diagnóstico ARMA completado.\n")

}
