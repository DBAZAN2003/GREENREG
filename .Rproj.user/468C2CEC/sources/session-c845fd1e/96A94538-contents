#' Regresión Lineal Múltiple (RLM)
#'
#' Esta función ajusta un modelo de regresión lineal con múltiples variables predictoras
#' ($y = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + ...$) mediante Mínimos Cuadrados Ordinarios (OLS).
#' Incluye un robusto motor de diagnóstico de Multicolinealidad (VIF) y pruebas de supuestos
#' con recomendaciones de acción para el investigador.
#'
#' @details
#' El flujo interno de la función se divide en seis etapas:
#' \enumerate{
#'   \item **Validación de Datos:** Asegura que la fórmula sea válida y que el dataframe no contenga errores estructurales.
#'   \item **Ajuste (Engine):** Utiliza \code{stats::lm}. En RLM, el foco principal es el **R² Ajustado**,
#'         que penaliza la inclusión de variables innecesarias que no aportan información real.
#'   \item **Detección de Multicolinealidad (VIF):** Calcula manualmente el Factor de Inflación de
#'         la Varianza ($VIF = 1 / (1 - R^2_j)$) para cada predictor. Esto identifica si las variables
#'         X están correlacionadas entre sí, lo cual infla los errores estándar y vuelve inestables los coeficientes.
#'   \item **Auditoría de Supuestos:** Ejecuta de forma nativa pruebas de Normalidad (Shapiro),
#'         Homocedasticidad (Breusch-Pagan), Independencia (Durbin-Watson) y Linealidad (Rainbow Test).
#'   \item **Traducción Automática:** Evalúa los estadísticos y genera "Notas" (✅, ⚠️, ❌)
#'         y sugerencias de mejora (ej. transformar variables o eliminar predictores redundantes).
#'   \item **Preparación de Diagnóstico:** Construye un dataframe optimizado para \code{ggplot2} que
#'         incluye residuos estandarizados, distancias de Cook y apalancamiento (leverage).
#' }
#'
#' @section Guía de las 10 Gráficas de Diagnóstico (`plot`):
#' Al ejecutar \code{plot(modelo)}, se despliega una secuencia educativa de 10 visualizaciones:
#' \describe{
#'   \item{\strong{1. Distribución de Variables (Y, Xi)}}{Muestra histogramas y densidades de todas las variables. Permite detectar sesgos o datos atípicos antes del ajuste.}
#'   \item{\strong{2. Relación de Variables (Y vs Xi)}}{Scatterplots individuales de la respuesta contra cada predictor con una curva LOESS para verificar la linealidad visualmente.}
#'   \item{\strong{3. Ajuste (Observados vs Ajustados)}}{Muestra la precisión global. Los puntos deben seguir la diagonal de 45°.}
#'   \item{\strong{4. Q-Q Plot (Normalidad)}}{Valida si los residuos siguen una distribución normal. Crucial para la validez de los P-valores.}
#'   \item{\strong{5. Homocedasticidad}}{Residuos vs. Predichos. Busca una "nube de estrellas" sin forma de embudo.}
#'   \item{\strong{6. Independencia}}{Residuos vs. Orden de recolección. Detecta patrones temporales o espaciales en los errores.}
#'   \item{\strong{7. Predicción vs Realidad (Identidad)}}{Versión simplificada del ajuste para evaluar la correlación de Pearson.}
#'   \item{\strong{8. Multicolinealidad (VIF)}}{Gráfica de barras exclusiva de RLM. Las barras rojas (VIF > 10) indican variables que deben ser eliminadas por redundancia.}
#'   \item{\strong{9. Distancia de Cook}}{Identifica observaciones influyentes. Barras que superan el umbral $4/n$ indican datos que deforman los resultados.}
#'   \item{\strong{10. Leverage vs Residuos}}{Cruza la "rareza" de X (Leverage) con el error en Y (Residuo). Detecta puntos que actúan como palancas sobre la recta.}
#' }
#'
#' @param formula Fórmula modelo.
#' @param data Dataframe con los datos.
#' @param nivel_confianza Nivel de confianza para intervalos (default 0.95).
#'
#' @return Un objeto de clase \code{"rlm_greenreg"} que contiene:
#' \itemize{
#'   \item \code{modelo}: Objeto \code{lm} original.
#'   \item \code{variables}: Lista con los nombres de Y y los predictores X.
#'   \item \code{coeficientes}: Matriz con estimaciones, p-valores e intervalos.
#'   \item \code{vif}: Vector con los valores de inflación de varianza por variable.
#'   \item \code{supuestos}: Resultados de los tests de Shapiro, BP, DW y Rainbow.
#'   \item \code{notas}: Reporte interpretado con semáforos de advertencia.
#'   \item \code{data_plot}: Dataframe listo para visualización avanzada.
#' }
#' @examples
#' # data(datos_rendimiento_maiz)
#' # modelo_rlm <- rlm(rendimiento_maiz_ton_ha ~., data = datos_rendimiento_maiz)
#' # modelo_rlm       # Ver reporte
#' # plot(modelo_rlm) # Ver gráficas incluyendo VIF
#'
#' @importFrom stats lm anova shapiro.test residuals fitted coef confint sd median cor pf model.matrix cooks.distance hatvalues rstandard
#' @importFrom lmtest bptest dwtest raintest
#' @import ggplot2
#' @export
rlm <- function(formula, data, nivel_confianza = 0.95) {

  # --- 1. Validación y Ajuste ---
  if (!inherits(formula, "formula")) stop("El argumento 'formula' debe ser una fórmula (ej. y ~ x1 + x2).")
  if (!is.data.frame(data)) stop("El argumento 'data' debe ser un data.frame.")

  modelo_base <- stats::lm(formula, data = data)
  resumen <- summary(modelo_base)

  # Identificar variables
  var_y <- names(modelo_base$model)[1]
  vars_x <- names(modelo_base$model)[-1] # Todas las predictoras

  # --- 2. Estadísticos Clave ---
  # En RLM, el R2 Ajustado es el rey
  r2_adj <- resumen$adj.r.squared
  sigma <- resumen$sigma
  f_statistic <- resumen$fstatistic

  # P-valor global
  if (!is.null(f_statistic)) {
    p_valor_modelo <- pf(f_statistic[1], f_statistic[2], f_statistic[3], lower.tail = FALSE)
  } else {
    p_valor_modelo <- NA
  }

  coefs <- resumen$coefficients
  intervalos <- stats::confint(modelo_base, level = nivel_confianza)
  tabla_coefs <- cbind(coefs, intervalos)
  tabla_anova <- stats::anova(modelo_base)

  # --- 3. Diagnósticos de Supuestos ---
  residuos <- stats::residuals(modelo_base)

  # A) Normalidad (Shapiro)
  if (length(residuos) >= 3 && length(residuos) <= 5000) {
    p_norm <- stats::shapiro.test(residuos)$p.value
  } else {
    p_norm <- NA
  }

  # B) Homocedasticidad (Breusch-Pagan)
  if (requireNamespace("lmtest", quietly = TRUE)) {
    p_homo <- lmtest::bptest(modelo_base)$p.value
  } else {
    p_homo <- NA
  }

  # C) Independencia (Durbin-Watson)
  if (requireNamespace("lmtest", quietly = TRUE)) {
    dw_res <- lmtest::dwtest(modelo_base)
    stat_dw <- dw_res$statistic
  } else {
    stat_dw <- NA
  }

  # D) MULTICOLINEALIDAD (VIF) - Manual robusto para no depender de 'car'
  # VIF_j = 1 / (1 - R2_j)
  # Calculamos el VIF solo si hay más de 1 predictor
  vif_vals <- numeric()
  if (length(vars_x) > 1) {
    tryCatch({
      # Matriz de diseño sin intercepto
      X_mat <- stats::model.matrix(modelo_base)[, -1, drop = FALSE]
      if(ncol(X_mat) > 1) {
        for(i in 1:ncol(X_mat)) {
          # Regresión de X_i contra los otros X
          y_aux <- X_mat[, i]
          x_aux <- X_mat[, -i, drop = FALSE]
          r2_aux <- summary(lm(y_aux ~ x_aux))$r.squared
          vif_vals[colnames(X_mat)[i]] <- 1 / (1 - r2_aux)
        }
      }
    }, error = function(e) {
      warning("No se pudo calcular el VIF (posiblemente variables categóricas complejas).")
    })
  }
  # E) LINEALIDAD (Rainbow Test)
  if (requireNamespace("lmtest", quietly = TRUE)) {
    test_lineal <- lmtest::raintest(modelo_base)
    p_lineal <- test_lineal$p.value
  } else {
    p_lineal <- NA
  }
  # --- 4. Notas (Interpretación) ---

  # R2 Ajustado
  nota_r2 <- paste0("El modelo explica el ", round(r2_adj * 100, 2),
                    "% de la variabilidad (R² Ajustado). ",
                    ifelse(r2_adj > 0.7, "Ajuste BUENO.",
                           ifelse(r2_adj > 0.4, "Ajuste MODERADO.", "Ajuste BAJO.")))

  # Nota Global
  nota_global <- if (!is.na(p_valor_modelo) && p_valor_modelo < 0.05) {
    "✅ El modelo en conjunto es SIGNIFICATIVO (al menos una variable ayuda a predecir)."
  } else {
    "⚠️ El modelo en conjunto NO es significativo (ninguna variable parece servir)."
  }

  # Semáforos
  nota_norm <- if (!is.na(p_norm) && p_norm > 0.05) "✅ OK: Residuos normales." else "⚠️ ALERTA: Residuos NO normales."
  nota_homo <- if (!is.na(p_homo) && p_homo > 0.05) "✅ OK: Varianza constante." else "⚠️ ALERTA: Heterocedasticidad detectada."
  nota_indep <- if (!is.na(stat_dw) && stat_dw > 1.5 && stat_dw < 2.5) "✅ OK: Independencia." else "⚠️ ALERTA: Posible autocorrelación."
  nota_lineal <- if (!is.na(p_lineal) && p_lineal > 0.05) {
    "✅ OK: Relación lineal adecuada."
  } else {
    "⚠️ ALERTA: Posible No-Linealidad. 💡 SUGERENCIA: Prueba transformar variables o agregar términos cuadráticos."
  }
  # Semáforo VIF
  msg_vif <- "No aplica (1 variable)."
  if (length(vif_vals) > 0) {
    max_vif <- max(vif_vals)
    if (max_vif < 5) {
      msg_vif <- "✅ OK: No hay problemas graves de Multicolinealidad (VIF < 5)."
    } else if (max_vif < 10) {
      msg_vif <- "⚠️ PRECAUCIÓN: Hay cierta correlación entre variables (5 < VIF < 10)."
    } else {
      msg_vif <- "❌ PELIGRO: Multicolinealidad severa (VIF > 10). Tus variables se estorban entre sí."
    }
  }

  # --- 5. Objeto Final ---
  resultado <- list(
    modelo = modelo_base,
    variables = list(y = var_y, x = vars_x),
    coeficientes = tabla_coefs,
    anova = tabla_anova,
    vif = vif_vals,
    resumen_estadistico = list(R2_adj = r2_adj, P_valor = p_valor_modelo),
    supuestos = list(norm = p_norm, homo = p_homo, dw = stat_dw, lin = p_lineal),
    notas = list(
      r2 = nota_r2,
      global = nota_global,
      supuestos = c(nota_norm, nota_homo, nota_indep, msg_vif)
    ),
    data_plot = data.frame(
      obs_y = data[[var_y]],
      pred = fitted(modelo_base),
      resid = residuals(modelo_base),
      index = 1:length(residuos)
    )
  )

  class(resultado) <- "rlm_greenreg"
  return(resultado)
}


#' Reporte de Regresión Lineal Múltiple
#' @description Imprime un resumen detallado del modelo RLM.
#' @param x Objeto de clase rlm_greenreg.
#' @param ... Argumentos adicionales.
#' @export
print.rlm_greenreg <- function(x, ...) {

  # 0. Preparación de la Ecuación Dinámica
  coefs_vals <- x$coeficientes[, 1]
  nombres_x <- x$variables$x

  # Construcción de la cadena: Y = b0 + b1*X1 + b2*X2...
  ecuacion <- paste(round(coefs_vals[1], 4)) # Intercepto
  for(i in 1:length(nombres_x)) {
    val <- round(coefs_vals[i+1], 4)
    signo <- if(val >= 0) " + " else " - "
    ecuacion <- paste0(ecuacion, signo, abs(val), " * ", nombres_x[i])
  }

  cat("\n==================================================================\n")
  cat("       REGRESIÓN LINEAL MÚLTIPLE (RLM)        \n")
  cat("==================================================================\n\n")

  # --- 1. MODELO MATEMÁTICO ESTIMADO ---
  cat("--- 1. MODELO MATEMÁTICO ESTIMADO ---\n")
  cat("• Variable Dependiente (Y): ", x$variables$y, "\n")
  cat("• Predictores (X):           ", paste(x$variables$x, collapse = ", "), "\n")
  cat("• Ecuación Estimada:         ", x$variables$y, " =", ecuacion, "\n")
  cat("• Ajuste (R² Ajustado):      ", round(x$resumen_estadistico$R2_adj * 100, 2), "%\n\n")

  # --- 2. COEFICIENTES ESTIMADOS ---
  cat("--- 2. COEFICIENTES ESTIMADOS ---\n")
  # Imprimimos las primeras 4 columnas (Est, Err, t, p)
  printCoefmat(x$coeficientes[, 1:4], digits = 4, signif.stars = TRUE, has.Pvalue = TRUE)

  # Si existe VIF, lo mostramos aquí mismo para contextualizar los coeficientes
  if (length(x$vif) > 0) {
    cat("\nMulticolinealidad (VIF):\n")
    print(round(x$vif, 2))
  }
  cat("\n")

  # --- 3. INTERVALOS DE CONFIANZA ---
  cat("--- 3. INTERVALOS DE CONFIANZA ---\n")
  # Columnas 5 y 6 del objeto coeficientes
  print(round(x$coeficientes[, 5:6], 4))
  cat("\n")

  # --- 4. ANÁLISIS DE VARIANZA (ANOVA) ---
  cat("--- 4. ANÁLISIS DE VARIANZA (ANOVA) ---\n")
  print(x$anova)
  cat("\n")

  # --- 5. VERIFICACIÓN DE SUPUESTOS ---
  cat("--- 5. VERIFICACIÓN DE SUPUESTOS ---\n")

  limpiar_txt <- function(txt) {
    if(is.null(txt) || is.na(txt)) return("No disponible")
    res <- strsplit(txt, ": ")[[1]]
    if(length(res) > 1) return(res[2]) else return(res[1])
  }

  # 1. Linealidad (Nuevo)
  cat(sprintf("• Linealidad (Rainbow Test):   [p-valor: %-7.4f] -> %s\n",
              x$supuestos$lin, limpiar_txt(x$notas$supuestos[which(grepl("lineal", x$notas$supuestos, ignore.case=T))])))

  # 2. Normalidad
  cat(sprintf("• Normalidad (Shapiro-Wilk):   [p-valor: %-7.4f] -> %s\n",
              x$supuestos$norm, limpiar_txt(x$notas$supuestos[1])))

  # 3. Homocedasticidad
  cat(sprintf("• Homocedasticidad (B-Pagan):  [p-valor: %-7.4f] -> %s\n",
              x$supuestos$homo, limpiar_txt(x$notas$supuestos[2])))

  # 4. Independencia
  cat(sprintf("• Independencia (D-Watson):    [Estad.: %-8.2f] -> %s\n",
              x$supuestos$dw, limpiar_txt(x$notas$supuestos[3])))

  # 5. VIF (Solo si hay más de 1 X)
  if (length(x$variables$x) > 1) {
    cat(sprintf("• Multicolinealidad (VIF):     %s\n",
                x$notas$supuestos[length(x$notas$supuestos)]))
  }

  # --- 6. RESUMEN DEL MODELO ---
  cat("\n--- 6. RESUMEN DEL MODELO ---\n")
  cat("GLOBAL:    ", x$notas$global, "\n")
  cat("AJUSTE:    ", x$notas$r2, "\n")

  # Verificación rápida de fallos
  hay_error <- any(grepl("❌|⚠️", x$notas$supuestos))
  cat("VALIDEZ:   ", if(hay_error) "❗ PRECAUCIÓN: Se detectaron fallas en los supuestos fundamentales."
      else "El modelo cumple con los supuestos básicos.", "\n")

  cat("\nUse 'plot(modelo)' para ver diagnósticos visuales.\n")
  cat("==================================================================\n\n")
}



#' Gráficas para RLM
#' @export
plot.rlm_greenreg <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Instala ggplot2.")

  df <- x$data_plot
  var_y <- x$variables$y

  # Tema base
  # Temas comunes para consistencia visual
  mi_tema <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, color = "#2C3E50"),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#34495E"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 10, face = "italic", color = "#555555", margin = ggplot2::margin(t = 10)),
      axis.title = ggplot2::element_text(face = "bold")
    )


  cat("Generando gráficas... (Presiona [Enter] para avanzar)\n")

  # --- GRÁFICA 1: HISTOGRAMAS (Y, Xi) (DISTRIBUCIÓN) ---
  # Paso 1: Preparar datos en formato largo para ggplot2
  # Extraemos Y y todas las X del modelo original
  datos_raw <- x$modelo$model
  df_dist <- data.frame(
    Valor = as.numeric(unlist(datos_raw)),
    Variable = rep(names(datos_raw), each = nrow(datos_raw))
  )

  g1 <- ggplot2::ggplot(df_dist, ggplot2::aes(x = Valor, fill = Variable)) +
    # Histograma con densidad superpuesta para cada variable
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = 15, color = "white", alpha = 0.7) +
    ggplot2::geom_density(color = "#E74C3C", size = 1) +
    # Dividir en páneles (facetas)
    ggplot2::facet_wrap(~Variable, scales = "free") +
    # Estética y Colores
    ggplot2::scale_fill_viridis_d(option = "mako", begin = 0.3, end = 0.8) +
    ggplot2::labs(
      title = "1. Distribución de Variables (Y, Xi)",
      subtitle = "¿Cómo se distribuyen los datos originales antes del ajuste?",
      x = "Valor observado",
      y = "Densidad / Frecuencia",
      caption = paste0("INTERPRETACIÓN:\n",
                       "✅ BIEN: Las variables muestran distribuciones equilibradas (campanas).\n",
                       "⚠️ ALERTA: Si una variable está muy sesgada (cola larga), podría afectar los coeficientes.\n",
                       " NOTA: La línea roja ayuda a identificar la forma de la distribución.")
    ) +
    mi_tema +
    ggplot2::theme(legend.position = "none")
  print(g1)
  readline(prompt = "Gráfica 1/10 (Distribución Y, Xi) > ")


  # --- GRÁFICA 2: RELACIÓN DE VARIABLES (Y vs Xi) ---
  datos_raw <- x$modelo$model
  var_y_name <- names(datos_raw)[1]

  # Creamos un set de datos donde Y se repite para cada X
  df_relacion <- data.frame()
  for(col in names(datos_raw)[-1]) {
    temp <- data.frame(
      Y = datos_raw[[1]],
      X_val = datos_raw[[col]],
      X_name = col
    )
    df_relacion <- rbind(df_relacion, temp)
  }
  g2 <- ggplot2::ggplot(df_relacion, ggplot2::aes(x = X_val, y = Y)) +
    # Puntos con el estilo de tu tema original
    ggplot2::geom_point(color = "#2C3E50", size = 2, alpha = 0.5) +
    # Añadimos una línea de tendencia suave (LOESS) para ver la forma de la relación
    ggplot2::geom_smooth(method = "loess", color = "#3498DB", fill = "#3498DB", alpha = 0.2, size = 1) +
    # Dividir por cada variable X
    ggplot2::facet_wrap(~X_name, scales = "free_x") +
    ggplot2::labs(
      title = "2. Relación de Variables (Y frente a Predictoras)",
      subtitle = paste("Visualizando la tendencia individual de cada X respecto a", var_y_name),
      x = "Valor de la Variable Predictora (Xi)",
      y = paste("Variable Respuesta (", var_y_name, ")"),
      caption = paste0("INTERPRETACIÓN:\n",
                       "✅ BIEN: Se observa una tendencia clara (lineal o curva suave) en los puntos.\n",
                       "❌ MAL: Nubes de puntos totalmente dispersas sugieren que esa variable no aporta al modelo.\n",
                       " NOTA: La línea azul ayuda a identificar si la relación es realmente lineal.")
    ) +
    mi_tema
  print(g2)
  readline(prompt = "Gráfica 2/10 (Relación Y vs Xi) > ")


  # --- GRÁFICA 3: AJUSTE (OBSERVADOS VS AJUSTADOS) ---
  # Comparamos qué tanto se aleja la realidad (Y) de la predicción del modelo (Ŷ)

  # Calculamos la correlación para mostrarla en el subtítulo
  cor_ajuste <- stats::cor(df$obs_y, df$pred)

  g3 <- ggplot2::ggplot(df, ggplot2::aes(x = obs_y, y = pred)) +
    # Línea de identidad (diagonal de 45°) - El ideal
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", size = 1) +
    # Puntos observados vs predichos
    ggplot2::geom_point(color = "#2C3E50", size = 3, alpha = 0.6) +
    # Añadimos una tendencia suave para ver desviaciones sistemáticas
    ggplot2::geom_smooth(method = "lm", color = "#3498DB", fill = "#3498DB", alpha = 0.1, linetype = "dotted") +
    ggplot2::labs(
      title = "3. Ajuste: Valores Observados vs. Ajustados",
      subtitle = paste("Correlación de ajuste (Pearson):", round(cor_ajuste, 4)),
      x = paste("Valores Reales (", x$variables$y, ")"),
      y = "Valores Predichos por el Modelo (Ŷ)",
      caption = paste0("INTERPRETACIÓN:\n",
                       "✅ BIEN: Los puntos se agrupan estrechamente sobre la línea roja diagonal.\n",
                       "❌ MAL: Nube muy dispersa o con forma curvada. El modelo no está capturando la realidad.\n",
                       " CLAVE: La línea roja representa la perfección (Predicho = Real).")
    ) +
    mi_tema

  print(g3)
  readline(prompt = "Gráfica 3/10 (Observados vs Ajustados) > ")


  # 4. Normalidad QQ
  g4 <- ggplot2::ggplot(df, ggplot2::aes(sample = resid)) +
    ggplot2::geom_qq(color = "#3498DB") +
    ggplot2::geom_qq_line(color = "#E74C3C", linetype = "dashed") +
    ggplot2::labs(
      title = "4. Gráfico Q-Q (Validación de Normalidad)",
      subtitle = "Los puntos deben seguir la línea diagonal sólida.",
      x = "Cuantiles Teóricos",
      y = "Cuantiles de los Residuos",
      caption = "INTERPRETACIÓN:\n✅ BIEN: Los puntos están pegados a la línea (como hormigas en fila).\n⚠️ CUIDADO: Si los puntos de los extremos se curvan o alejan mucho,\ntus pruebas de significancia (P-Values) podrían ser poco fiables."
    ) + mi_tema
  print(g4); readline(prompt = "Grafica 4/10 (NORMALIDAD) > ")

  # 5. Homocedasticidad
  g5 <- ggplot2::ggplot(df, ggplot2::aes(x = pred, y = resid)) +
    ggplot2::geom_point(color = "#8E44AD", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = "5. Homocedasticidad (Varianza Constante)",
      subtitle = "Busca una 'nube aleatoria' (cielo estrellado).",
      x = "Valores Predichos",
      y = "Errores (Residuos)",
      caption = "INTERPRETACIÓN:\n✅ BIEN: Una nube rectangular sin forma definida.\n❌ MAL (Embudo): Si los puntos se abren como un abanico <, la varianza cambia (error grave).\n❌ MAL (Curva): Si ves una forma de 'U', te falta elevar una variable al cuadrado."
    ) + mi_tema
  print(g5); readline(prompt = "Grafica 5/10 (HOMOCEDASTICIDAD) > ")

  # 6. Independencia
  g6 <- ggplot2::ggplot(df, ggplot2::aes(x = index, y = resid)) +
    ggplot2::geom_line(color = "#95A5A6") +
    ggplot2::geom_point(color = "#D35400") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = "6. Independencia de los Errores",
      subtitle = "¿Ves algún patrón repetitivo u olas en los datos?",
      caption = "INTERPRETACIÓN:\n✅ BIEN: Ruido puro, picos arriba y abajo sin orden aparente.\n❌ MAL: Si ves olas suaves, zig-zags constantes o grupos, \nsignifica que el pasado afecta al futuro (falta una variable temporal)."
    ) + mi_tema
  print(g6); readline(prompt = "Grafica 6/10 (INDEPENDENCIA) > ")

  # 7. Observados vs Predichos
  correlacion <- cor(df$obs_y, df$pred)
  g7 <- ggplot2::ggplot(df, ggplot2::aes(x = obs_y, y = pred)) +
    ggplot2::geom_point(color = "#2980B9", alpha = 0.6, size = 2.5) +
    # Agregamos la línea diagonal perfecta (Ideal)
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
    ggplot2::labs(
      title = "7. Predicción vs Realidad (Precisión)",
      subtitle = "Lo ideal es que los puntos caigan sobre la línea roja discontinua.",
      x = "Valor Real (Observado)",
      y = "Valor Predicho por el Modelo",
      caption = "INTERPRETACIÓN:\n✅ BIEN: Los puntos forman una línea diagonal estrecha (45 grados).\n❌ MAL: Si ves una nube redonda dispersa, el modelo no está prediciendo nada útil."
    ) + mi_tema
  print(g7); readline(prompt = "Grafica 7/10(OBSERVADOS VS PREDICHOS)  > ")


  # 8. VIF (Multicolinealidad)
  if (length(x$vif) > 0) {
    df_vif <- data.frame(Variable = names(x$vif), VIF = x$vif)
    g8 <- ggplot2::ggplot(df_vif, ggplot2::aes(x = Variable, y = VIF, fill = VIF > 10)) +
      ggplot2::geom_col() +
      ggplot2::geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
      ggplot2::scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#2ECC71"), guide = "none") +
      ggplot2::labs(
        title = "8. Multicolinealidad (Inflación de Varianza)",
        subtitle = "Variables que dicen lo mismo (Redundancia).",
        y = "Valor VIF",
        x = "Variables Predictoras",
        caption = "INTERPRETACIÓN:\n✅ BIEN (Verde): La variable aporta información única.\n❌ MAL (Rojo): VIF > 10. Esta variable es redundante (es combinación de otras).\nACCIÓN: Elimina las variables con barra roja para estabilizar el modelo."
      ) + mi_tema
    print(g8); readline(prompt = "Grafica 8/10 (MULTICOLINEALIDAD) > ")
  } else {
    cat("ℹ️ Gráfica de VIF no generada (solo hay 1 variable predictora).\n")
  }
  # --- 9. Distancia de Cook (Valores Influyentes) ---
  cooksd <- stats::cooks.distance(x$modelo)
  n_datos <- length(cooksd)
  umbral_cook <- 4 / n_datos  # Criterio estándar: 4/n

  df_cook <- data.frame(Index = 1:n_datos, Cook = cooksd)

  # Identificamos los puntos peligrosos para ponerles etiqueta
  df_cook$Label <- ifelse(df_cook$Cook > umbral_cook, as.character(df_cook$Index), "")

  g9 <- ggplot2::ggplot(df_cook, ggplot2::aes(x = Index, y = Cook)) +
    ggplot2::geom_bar(stat = "identity", fill = "#8E44AD", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = umbral_cook, color = "red", linetype = "dashed") +
    ggplot2::geom_text(ggplot2::aes(label = Label), vjust = -0.5, size = 3) +
    ggplot2::labs(
      title = "9. Valores Influyentes (Distancia de Cook)",
      subtitle = "¿Algún dato está forzando al modelo a cambiar?",
      x = "Número de Observación (Fila)",
      y = "Influencia (Cook)",
      caption = "INTERPRETACIÓN:\n⚠️ ALERTA: Las barras que cruzan la línea roja son datos sospechosos.\nACCIÓN: Revisa esas filas en tu Excel. ¿Son errores de dedo o casos especiales?\nSi son errores, elimínalos y corre el modelo de nuevo."
    ) + mi_tema

  print(g9)
  readline(prompt = "Grafica 9/10 (INFLUENCIA) > ")

  # --- 10. Residuos vs Apalancamiento (Leverage) ---
  lev <- stats::hatvalues(x$modelo)
  res_std <- stats::rstandard(x$modelo)
  n_datos <- length(lev)
  p_params <- length(stats::coef(x$modelo)) # Número de coeficientes (incluye intercepto)

  # Umbral de leverage para RLM: 2 * (p / n)
  umbral_lev <- 2 * (p_params / n_datos)

  df_diag <- data.frame(
    Index = 1:n_datos,
    Leverage = lev,
    ResidStd = res_std
  )

  # Etiquetar si el leverage es alto (X extremo) O si el residuo es grande (Outlier en Y)
  df_diag$Label <- ifelse(df_diag$Leverage > umbral_lev | abs(df_diag$ResidStd) > 2,
                          as.character(df_diag$Index), "")

  g10 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = Leverage, y = ResidStd)) +
    # El tamaño del punto ayuda a ver el peso del leverage
    ggplot2::geom_point(ggplot2::aes(size = Leverage), color = "#C0392B", alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    # Líneas de advertencia para Outliers (Residuos > 2 o < -2)
    ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "red", alpha = 0.5) +
    # Línea de advertencia para Leverage alto
    ggplot2::geom_vline(xintercept = umbral_lev, linetype = "dashed", color = "blue", alpha = 0.5) +
    # Etiquetas de los números de fila (usando vjust para que no tapen el punto)
    ggplot2::geom_text(ggplot2::aes(label = Label), vjust = -1, size = 3, fontface = "bold") +
    ggplot2::labs(
      title = "10. Diagnóstico de Influencia (Leverage vs Residuos)",
      subtitle = "Identificación de puntos con combinaciones de variables extremas",
      x = "Apalancamiento (Leverage)",
      y = "Residuos Estandarizados",
      caption = paste0("INTERPRETACIÓN:\n",
                       "⚠️ EJE X: Puntos a la derecha de la línea azul tienen valores de X inusuales.\n",
                       "⚠️ EJE Y: Puntos fuera de las líneas rojas (2, -2) son Outliers.\n",
                       " PELIGRO: Si un punto está muy a la derecha Y muy lejos del 0 en Y, es una 'palanca' que deforma el modelo.")
    ) +
    mi_tema +
    ggplot2::theme(legend.position = "none")

  print(g10)
  readline(prompt = "Grafica 10/10 (LEVERAGE) > ")

  cat("✅ Gráficas completadas.\n")
}
