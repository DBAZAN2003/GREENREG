#' Regresión Lineal Simple (RLS)
#'
#' Esta función ajusta un modelo de regresión lineal ($y = \beta_0 + \beta_1 x$) mediante
#' Mínimos Cuadrados Ordinarios (OLS). A diferencia de `lm()` estándar, esta función
#' devuelve un objeto enriquecido con interpretaciones en lenguaje natural, verificaciones
#' automáticas de supuestos y sugerencias para el investigador.
#'
#' @details
#' El flujo interno de la función se divide en cinco etapas críticas:
#' \enumerate{
#'   \item **Preparación:** Valida que los datos sean coherentes y extrae las variables de la fórmula.
#'   \item **Estimación:** Calcula los coeficientes, el error estándar residual ($\sigma$) y la significancia global (F-test).
#'   \item **Inferencia:** Genera intervalos de confianza para los parámetros ($\beta$) al nivel especificado.
#'   \item **Auditoría de Supuestos:** Ejecuta pruebas de Shapiro-Wilk (Normalidad), Breusch-Pagan (Homocedasticidad)
#'         y Durbin-Watson (Independencia).
#'   \item **Traducción:** Convierte los estadísticos abstractos en "Notas" de texto que recomiendan acciones
#'         específicas (ej. aplicar logaritmos o revisar datos atípicos).
#' }
#'
#' @section Guía de las 9 Gráficas de Diagnóstico (`plot`):
#' Al ejecutar \code{plot(modelo)}, se despliega una secuencia educativa:
#' \describe{
#'   \item{\strong{1. Histogramas Exploratorios}}{Compara la distribución de X e Y. Útil para detectar sesgos iniciales en las variables originales.}
#'   \item{\strong{2. Relación Original (Scatterplot)}}{Muestra la nube de puntos pura. Ayuda a confirmar visualmente si la relación parece lineal o curva antes de ver la recta.}
#'   \item{\strong{3. Recta de Regresión}}{Visualiza el ajuste final con bandas de confianza. Permite ver qué tan bien la recta representa la tendencia de los datos.}
#'   \item{\strong{4. Q-Q Plot (Normalidad)}}{Grafica los cuantiles de los residuos contra una normal teórica. Si los puntos se alejan de la diagonal, los p-valores del modelo podrían no ser válidos.}
#'   \item{\strong{5. Homocedasticidad}}{Residuos vs. Predichos. Si la dispersión aumenta (forma de abanico), el modelo pierde precisión en valores altos.}
#'   \item{\strong{6. Independencia}}{Residuos vs. Orden de datos. Vital para detectar si el error de hoy depende del error de ayer (Autocorrelación).}
#'   \item{\strong{7. Precisión (Predicho vs Real)}}{Muestra la correlación directa. En un modelo perfecto, todos los puntos caerían sobre la línea de identidad (45°).}
#'   \item{\strong{8. Distancia de Cook}}{Identifica observaciones "influyentes". Una barra alta indica un dato que, si se elimina, cambiaría drásticamente la pendiente del modelo.}
#'   \item{\strong{9. Leverage vs Residuos}}{Diferencia entre 'Outliers' (error en Y) y 'Apalancamiento' (valor extremo en X). Identifica puntos que "jalonean" la recta injustificadamente.}
#' }
#'
#'
#' @param formula Una fórmula de R estándar (ej. `rendimiento ~ fertilizante`).
#' @param data El data.frame que contiene las variables.
#' @param nivel_confianza Nivel de confianza para los intervalos (por defecto 0.95).
#'
#' @return Un objeto de clase \code{"rls_greenreg"} que es una lista con:
#' \itemize{
#'   \item \code{modelo}: El objeto \code{lm} original.
#'   \item \code{variables}: Nombres de X e Y.
#'   \item \code{coeficientes}: Matriz con Estimaciones, Errores y P-valores.
#'   \item \code{supuestos}: Resultados de los tests de Shapiro, BP y DW.
#'   \item \code{notas}: Diccionario de interpretaciones y sugerencias de mejora.
#'   \item \code{data_plot}: Dataframe listo para ggplot2 con residuos estandarizados e índices de influencia.
#' }
#'
#'
#' @examples
#' # data(maiz)
#' # modelo <- rls(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = datos_rendimiento_maiz)
#' # modelo       # Imprime el reporte
#' # plot(modelo) # Genera las 9 gráficas
#'
#' @importFrom stats lm anova shapiro.test residuals fitted coef confint sd median cor pf cooks.distance hatvalues rstandard
#' @importFrom lmtest bptest dwtest
#' @import ggplot2
#' @export
rls <- function(formula, data, nivel_confianza = 0.95) {

  # --- 1. Validación y Preparación ---
  if (!inherits(formula, "formula")) stop("El argumento 'formula' debe ser una fórmula (ej. y ~ x).")
  if (!is.data.frame(data)) stop("El argumento 'data' debe ser un data.frame.")

  # Ajuste del modelo base (Usamos el motor robusto de R)
  modelo_base <- stats::lm(formula, data = data)
  resumen <- summary(modelo_base)

  # Nombres de variables
  var_y <- names(modelo_base$model)[1]
  var_x <- names(modelo_base$model)[2]

  # --- 2. Extracción de Estadísticos Clave ---
  r2 <- resumen$r.squared
  r2_adj <- resumen$adj.r.squared
  sigma <- resumen$sigma # Error estándar residual
  f_statistic <- resumen$fstatistic

  # Calculo manual del p-value del modelo global (ANOVA)
  if (!is.null(f_statistic)) {
    p_valor_modelo <- pf(f_statistic[1], f_statistic[2], f_statistic[3], lower.tail = FALSE)
  } else {
    p_valor_modelo <- NA
  }

  # Coeficientes e Intervalos
  coefs <- resumen$coefficients
  intervalos <- stats::confint(modelo_base, level = nivel_confianza)
  # Combinar coeficientes con sus intervalos
  tabla_coefs <- cbind(coefs, intervalos)

  # Tabla ANOVA
  tabla_anova <- stats::anova(modelo_base)

  # --- 3. Pruebas de Diagnóstico (Supuestos) ---
  residuos <- stats::residuals(modelo_base)

  # A) Normalidad (Shapiro-Wilk)
  # Nota: Shapiro solo funciona bien entre 3 y 5000 datos.
  if (length(residuos) >= 3 && length(residuos) <= 5000) {
    test_shapiro <- stats::shapiro.test(residuos)
    p_norm <- test_shapiro$p.value
  } else {
    p_norm <- NA
  }

  # B) Homocedasticidad (Breusch-Pagan)
  if (requireNamespace("lmtest", quietly = TRUE)) {
    test_bp <- lmtest::bptest(modelo_base)
    p_homo <- test_bp$p.value
  } else {
    p_homo <- NA
    warning("Instale el paquete 'lmtest' para la prueba de Homocedasticidad.")
  }

  # C) Independencia (Durbin-Watson)
  if (requireNamespace("lmtest", quietly = TRUE)) {
    test_dw <- lmtest::dwtest(modelo_base)
    stat_dw <- test_dw$statistic
    p_indep <- test_dw$p.value
  } else {
    stat_dw <- NA; p_indep <- NA
  }

  # --- 4. Generación de Notas (CON RECOMENDACIONES) ---

  # Interpretación de R2
  nota_r2 <- paste0("El modelo explica el ", round(r2 * 100, 2),
                    "% de la variabilidad de '", var_y, "'. ",
                    ifelse(r2 > 0.7, "Es un ajuste ALTO (Bueno).",
                           ifelse(r2 > 0.4, "Es un ajuste MODERADO.", "Es un ajuste BAJO.")))

  # Interpretación de la Pendiente (Beta 1) + Recomendación
  beta1 <- coefs[2, 1]
  pval_beta1 <- coefs[2, 4]
  nota_beta1 <- if (pval_beta1 < 0.05) {
    paste0("La relación es SIGNIFICATIVA. Por cada aumento de 1 unidad en '", var_x,
           "', la variable '", var_y, "' cambia en promedio ", round(beta1, 4), " unidades.")
  } else {
    paste0("La variable '", var_x, "' NO tiene un efecto estadísticamente significativo. ",
           "💡 SUGERENCIA: Intenta recolectar más datos, revisa si la relación es curva (no lineal) o consulta a un experto para agregar más variables explicativas.")
  }

  # Interpretación de Supuestos con Recomendaciones de Acción
  nota_norm <- if (!is.na(p_norm) && p_norm > 0.05) {
    "✅ OK: Los residuos parecen normales (p > 0.05)."
  } else {
    "⚠️ ALERTA: Los residuos NO son normales. 💡 ACCIÓN: Los p-values pueden ser engañosos. Prueba transformar Y (ej. log(y)) o usa modelos no paramétricos."
  }

  nota_homo <- if (!is.na(p_homo) && p_homo > 0.05) {
    "✅ OK: Varianza constante (Homocedasticidad)."
  } else {
    "⚠️ ALERTA: Heterocedasticidad detectada. 💡 ACCIÓN: La precisión es inestable. Se recomienda usar errores estándar robustos o un modelo GLM."
  }

  nota_indep <- if (!is.na(stat_dw) && stat_dw > 1.5 && stat_dw < 2.5) {
    "✅ OK: Residuos independientes (DW cercano a 2)."
  } else {
    "⚠️ ALERTA: Posible autocorrelación. 💡 ACCIÓN: Si tus datos son temporales, este modelo es inválido. Consulta a un experto en Series de Tiempo."
  }

  # Nota de aviso general si algo falla
  fallo_critico <- (!is.na(p_norm) && p_norm <= 0.05) || (!is.na(p_homo) && p_homo <= 0.05) || (!is.na(stat_dw) && (stat_dw <= 1.5 | stat_dw >= 2.5))
  nota_experto <- if(fallo_critico) "❗ RECOMENDACIÓN FINAL: Se detectaron fallas en los supuestos. Los resultados podrían no ser confiables; consulta a un analista de datos." else "✨ El modelo cumple con los supuestos básicos."

  # --- 5. Construcción del Objeto Final ---
  resultado <- list(
    modelo = modelo_base,
    variables = list(y = var_y, x = var_x),
    info_datos = list(n = length(residuos), media_y = mean(data[[var_y]], na.rm=T), sd_y = sd(data[[var_y]], na.rm=T)),
    anova = tabla_anova,
    coeficientes = tabla_coefs,
    resumen_estadistico = list(R2 = r2, R2_adj = r2_adj, Sigma = sigma, F_stat = f_statistic),
    supuestos = list(
      shapiro_p = p_norm,
      bp_p = p_homo,
      dw_stat = stat_dw
    ),
    notas = list(
      r2 = nota_r2,
      impacto = nota_beta1,
      supuestos = c(nota_norm, nota_homo, nota_indep),
      alerta_experto = nota_experto
    ),
    data_plot = data.frame(
      obs_y = data[[var_y]],
      obs_x = data[[var_x]],
      pred = fitted(modelo_base),
      resid = residuals(modelo_base),
      index = 1:length(residuos)
    )
  )

  class(resultado) <- "rls_greenreg"
  return(resultado)
}

#' Impresión de resultados rls_greenreg
#' @description Método optimizado para imprimir el reporte de la regresión.
#' @param x Objeto de clase rls_greenreg.
#' @param ... Argumentos adicionales.
#' @export
print.rls_greenreg <- function(x, ...) {

  # 0. Preparación de la Ecuación
  b0 <- round(x$coeficientes[1, 1], 4)
  b1 <- round(x$coeficientes[2, 1], 4)
  # Creamos el signo dinámico con espacios para legibilidad
  signo_txt <- if(b1 >= 0) " + " else " - "

  cat("\n==================================================================\n")
  cat("       REGRESIÓN LINEAL SIMPLE        \n")
  cat("==================================================================\n\n")

  # --- 1. MODELO Y ECUACIÓN ---
  cat("\n--- 1. MODELO MATEMÁTICO ESTIMADO ---\n")
  cat("• Variable Dependiente (Y): ", x$variables$y, "\n")
  cat("• Variable Independiente (X):", x$variables$x, "\n")
  cat("• Ecuación de la Recta:  ", x$variables$y, " =", b0, signo_txt, abs(b1), "*", x$variables$x, "\n")
  cat("• Bondad de Ajuste (R²):  ", round(x$resumen_estadistico$R2 * 100, 2), "%\n\n")

  # --- 2. COEFICIENTES ESTIMADOS ---
  cat("\n--- 2. COEFICIENTES ESTIMADOS ---\n")
  printCoefmat(x$coeficientes[, 1:4],
               digits = 4,
               signif.stars = TRUE)


  # --- 3. INTERVALOS DE CONFIANZA ---
  cat("\n--- 3. INTERVALOS DE CONFIANZA ---\n")
  print(x$coeficientes[, 5:6])


  # --- 4. ANÁLISIS DE VARIANZA (ANOVA) ---
  cat("\n--- 4. ANÁLISIS DE VARIANZA (ANOVA) ---\n")
  # Imprimimos la tabla ANOVA estándar de R
  print(x$anova)
  cat("\n")


  # --- 5. VERIFICACIÓN DE SUPUESTOS ---
  cat("\n--- 5. VERIFICACIÓN DE SUPUESTOS ---\n")

  # Formateo con alineación fija para que parezca una tabla profesional
  cat(sprintf("• Normalidad (Shapiro-Wilk):   [p-valor: %-7.4f] -> %s\n",
              x$supuestos$shapiro_p, x$notas$supuestos[1]))

  cat(sprintf("• Homocedasticidad (B-Pagan):  [p-valor: %-7.4f] -> %s\n",
              x$supuestos$bp_p, x$notas$supuestos[2]))

  cat(sprintf("• Independencia (D-Watson):    [Estad.: %-8.2f] -> %s\n",
              x$supuestos$dw_stat, x$notas$supuestos[3]))


  # --- 6. RESUMEN DEL MODELO (Interpretación Final) ---
  cat("\n--- 6. RESUMEN DEL MODELO ---\n")
  cat("IMPACTO:  ", x$notas$impacto, "\n")
  cat("AJUSTE:   ", x$notas$r2, "\n")
  cat("VALIDEZ:  ", x$notas$alerta_experto, "\n")

  cat("\n Use 'plot(modelo)' para ver diagnósticos visuales.\n")
  cat("==================================================================\n\n")
}
  #' Gráficas de diagnóstico para rls_greenreg
  #' @description Genera una secuencia de 8 gráficas para validar el modelo.
  #' @param x Objeto de clase rls_greenreg.
  #' @param ... Argumentos adicionales.
  #' @export
plot.rls_greenreg <- function(x, ...) {
  # Asegurar ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("El paquete ggplot2 es necesario. Instálalo con install.packages('ggplot2')")
  }

  df <- x$data_plot
  var_y <- x$variables$y
  var_x <- x$variables$x

  # Tema visual consistente
  mi_tema <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, color = "#2C3E50"),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#34495E"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 10, face = "italic", color = "#555555", margin = ggplot2::margin(t = 10)),
      axis.title = ggplot2::element_text(face = "bold")
    )

  cat("Generando gráficas... (Presione [Enter] en la consola para avanzar)\n")

  # ==============================================================================
  # --- GRÁFICA 1: HISTOGRAMAS EXPLORATORIOS (X, Y) ---
  # ==============================================================================

  # Paso 1: Preparación de datos (Formato Largo)
  # Creamos un dataframe temporal donde las variables X e Y están en una sola columna
  # y otra columna indica a qué variable pertenece el dato.
  df_long <- data.frame(
    Valor = c(df$obs_x, df$obs_y),
    Variable = c(rep(paste("Variable X:", var_x), nrow(df)),
                 rep(paste("Variable Y:", var_y), nrow(df)))
  )

  # Paso 2: Creación de la gráfica combinada con ggplot2
  g1 <- ggplot2::ggplot(df_long, ggplot2::aes(x = Valor, fill = Variable)) +
    # Histograma base
    ggplot2::geom_histogram(bins = 15, color = "white", alpha = 0.8, show.legend = FALSE) +
    # Separar en dos páneles (uno para X, otro para Y)
    ggplot2::facet_wrap(~Variable, scales = "free") + # scales="free" ajusta el eje X a cada variable
    # Paleta de colores manual (Azul para X, Verde para Y)
    ggplot2::scale_fill_manual(values = c("#3498DB", "#16A085")) +
    # Etiquetas y Notas
    ggplot2::labs(
      title = "1. Exploración Inicial de Datos (Histogramas)",
      subtitle = "¿Cómo se distribuyen las variables originales antes del modelo?",
      x = "Valor de la Variable",
      y = "Frecuencia Absoluta",
      caption = paste0("INTERPRETACIÓN:\n",
                       "✅ BIEN: Las barras forman una figura aproximadamente simétrica (campana).\n",
                       "⚠️ ALERTA: Si un histograma tiene una 'cola' muy larga hacia la derecha o izquierda,\n",
                       "   la variable está sesgada. Esto puede afectar la precisión de la recta de regresión.")
    ) +
    # Aplicación estricta de 'mi_tema'
    mi_tema
  print(g1)
  readline(prompt = "Gráfica 1/9 (Histogramas X,Y) > Presione [Enter] para continuar...")

  # --- GRÁFICA 2: RELACIÓN ORIGINAL (SCATTERPLOT) ---
  # El objetivo es ver la nube de puntos pura para detectar la forma de la relación.
  g2 <- ggplot2::ggplot(df, ggplot2::aes(x = obs_x, y = obs_y)) +
    # Puntos con el color original de tu tema
    ggplot2::geom_point(color = "#2C3E50", size = 3, alpha = 0.6) +
    ggplot2::labs(
      title = "2. Relación Original (Scatterplot)",
      subtitle = paste("Explorando la correlación visual entre", var_x, "y", var_y),
      x = var_x,
      y = var_y,
      caption = paste0("INTERPRETACIÓN:\n",
                       "✅ BIEN: Los puntos forman una hilera o elipse clara (ascendente o descendente).\n",
                       "❌ MAL: Los puntos parecen una nube circular o ruidosa. Indica que no hay relación lineal.")
    ) +
    mi_tema

  print(g2)
  readline(prompt = "Gráfica 2/9 (Relación Original) > ")

  # --- GRÁFICA 3: RECTA DE REGRESIÓN (AJUSTE) ---
  # Aquí combinamos los puntos reales con la línea estimada por el modelo.
  g3 <- ggplot2::ggplot(df, ggplot2::aes(x = obs_x, y = obs_y)) +
    # Puntos originales (con transparencia para ver solapamiento)
    ggplot2::geom_point(color = "#2C3E50", size = 3, alpha = 0.4) +
    # La recta de regresión con su intervalo de confianza (sombra)
    ggplot2::geom_smooth(method = "lm", color = "#E74C3C", fill = "#E74C3C", alpha = 0.2, size = 1.2) +
    ggplot2::labs(
      title = "3. Ajuste del Modelo (Recta de Regresión)",
      subtitle = paste("Modelo Estimado:", var_y, "~", var_x),
      x = var_x,
      y = var_y,
      caption = paste0("INTERPRETACIÓN:\n",
                       "✅ BIEN: La línea roja cruza por el centro de la nube de puntos.\n",
                       " NOTA: La banda sombreada es el error de estimación. Si es muy ancha,\n",
                       "   el modelo tiene mucha incertidumbre en esa zona.")
    ) +
    mi_tema

  print(g3)
  readline(prompt = "Gráfica 3/9 (Ajuste del Modelo) > ")

  # --- GRÁFICA 4: NORMALIDAD (QQ-PLOT) ---
  g4 <- ggplot2::ggplot(df, ggplot2::aes(sample = resid)) +
    ggplot2::geom_qq(color = "#3498DB", size = 2) +
    ggplot2::geom_qq_line(color = "#E74C3C", size = 1.2, linetype = "dashed") +
    ggplot2::labs(
      title = "4. Validación de Normalidad (Q-Q Plot)",
      subtitle = "Comparación contra una distribución normal ideal",
      x = "Cuantiles Teóricos", y = "Cuantiles Observados",
      caption = "INTERPRETACIÓN:\n✅ BIEN: Los puntos azules están pegados a la línea roja discontinua.\n❌ MAL: Si los puntos en los extremos se curvan o alejan mucho de la línea."
    ) + mi_tema
  print(g4)
  readline(prompt = "Gráfica 4/9 (Normalidad) > ")

  # --- GRÁFICA 5: HOMOCEDASTICIDAD ---
  g5 <- ggplot2::ggplot(df, ggplot2::aes(x = pred, y = resid)) +
    ggplot2::geom_point(color = "#8E44AD", size = 2.5, alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::labs(
      title = "5. Homocedasticidad (Varianza Constante)",
      subtitle = "¿El error es igual para valores chicos y grandes?",
      x = "Valores Predichos", y = "Residuos (Errores)",
      caption = "INTERPRETACIÓN:\n✅ BIEN: Nube dispersa sin forma (como cielo estrellado).\n❌ MAL: Forma de 'Embudo' (<) o abanico. Significa que el error crece con el valor."
    ) + mi_tema
  print(g5)
  readline(prompt = "Gráfica 5/9 Homocedasticidad > ")

  # --- GRÁFICA 6: INDEPENDENCIA ---
  g6 <- ggplot2::ggplot(df, ggplot2::aes(x = index, y = resid)) +
    ggplot2::geom_line(color = "#95A5A6", alpha = 0.5) +
    ggplot2::geom_point(color = "#D35400", size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = "6. Independencia (Residuos vs Orden)",
      subtitle = "Detectando patrones en el tiempo de recolección",
      x = "Orden de Recolección", y = "Residuos",
      caption = "INTERPRETACIÓN:\n✅ BIEN: Ruido aleatorio (picos arriba y abajo sin orden).\n❌ MAL: Si ves ondas suaves o tendencias, los datos no son independientes."
    ) + mi_tema
  print(g6)
  readline(prompt = "Gráfica 6/9 Independencia > ")

  # --- GRÁFICA 7: PREDICCIÓN VS REALIDAD ---
  correlacion <- cor(df$obs_y, df$pred)
  g7 <- ggplot2::ggplot(df, ggplot2::aes(x = obs_y, y = pred)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", size=1) +
    ggplot2::geom_point(color = "#27AE60", size = 3, alpha = 0.7) +
    ggplot2::labs(
      title = "7. Precisión: Predicción vs Realidad",
      subtitle = paste("Correlación:", round(correlacion, 3)),
      x = "Valor Real", y = "Valor Predicho",
      caption = "INTERPRETACIÓN:\n✅ BIEN: Los puntos caen sobre la línea gris punteada (Identidad).\n❌ MAL: Nube dispersa y redonda (el modelo no predice bien)."
    ) + mi_tema
  print(g7)
  readline(prompt = "Gráfica 7/9 > ")

  # --- GRÁFICA 8: DISTANCIA DE COOK (NUEVA) ---
  # Cálculo de influencia
  cooksd <- stats::cooks.distance(x$modelo)
  n_datos <- length(cooksd)
  umbral_cook <- 4 / n_datos
  df_cook <- data.frame(Index = 1:n_datos, Cook = cooksd)
  # Etiquetar solo los peligrosos
  df_cook$Label <- ifelse(df_cook$Cook > umbral_cook, as.character(df_cook$Index), "")

  g8 <- ggplot2::ggplot(df_cook, ggplot2::aes(x = Index, y = Cook)) +
    ggplot2::geom_bar(stat = "identity", fill = "#E67E22", alpha = 0.8) +
    ggplot2::geom_hline(yintercept = umbral_cook, color = "red", linetype = "dashed") +
    ggplot2::geom_text(ggplot2::aes(label = Label), vjust = -0.5, size = 3) +
    ggplot2::labs(
      title = "8. Valores Influyentes (Distancia de Cook)",
      subtitle = "Detección de datos atípicos peligrosos",
      x = "Número de Observación", y = "Influencia",
      caption = "INTERPRETACIÓN:\n⚠️ ALERTA: Las barras que cruzan la línea roja son datos sospechosos.\nACCIÓN: Revisa esas filas en tu Excel, podrían ser errores de captura."
    ) + mi_tema
  print(g8)

  readline(prompt = "Gráfica 8/9 Influyentes > ")

  # --- GRÁFICA 9: RESIDUOS VS APALANCAMIENTO ---
  # Cálculo de leverage y residuos estandarizados
  lev <- stats::hatvalues(x$modelo)
  res_std <- stats::rstandard(x$modelo)
  n_datos <- length(lev)

  # Umbral típico de leverage: 2 * (p/n). Para RLS p = 2 (interceptor y pendiente)
  umbral_lev <- 4 / n_datos

  df_diag <- data.frame(
    index = 1:n_datos,
    leverage = lev,
    resid_std = res_std
  )

  # Etiquetar solo si el leverage es alto O el residuo es muy grande (Outlier)
  df_diag$Label <- ifelse(df_diag$leverage > umbral_lev | abs(df_diag$resid_std) > 2,
                          as.character(df_diag$index), "")

  g9 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = leverage, y = resid_std)) +
    ggplot2::geom_point(ggplot2::aes(size = leverage), color = "#C0392B", alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "red", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = umbral_lev, linetype = "dashed", color = "blue", alpha = 0.5) +
    # Etiquetas de los puntos
    ggrepel::geom_text_repel(ggplot2::aes(label = Label), size = 3.5, fontface = "bold") +
    ggplot2::labs(
      title = "9. Diagnóstico de Influencia (Leverage vs Residuos)",
      subtitle = "Identificando puntos que 'jalonean' la recta",
      x = "Apalancamiento (Leverage)",
      y = "Residuos Estandarizados",
      caption = paste0("INTERPRETACIÓN:\n",
                       "⚠️ EJE X: Puntos a la derecha de la línea azul tienen valores de X muy raros.\n",
                       "⚠️ EJE Y: Puntos fuera de las líneas rojas punteadas son Outliers (errores en Y).\n",
                       " CLAVE: Si un punto está muy a la derecha Y muy arriba/abajo, está deformando tu modelo.")
    ) +
    mi_tema +
    ggplot2::theme(legend.position = "none")

  print(g9)

  cat("✅ Secuencia de gráficas completada.\n")
 }

