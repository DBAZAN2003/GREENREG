#' Regresión Logística Binaria (MLG Logit)
#'
#' Esta función ajusta un Modelo Lineal Generalizado (GLM) de la familia binomial para
#' predecir la probabilidad de un evento dicotómico (Sí/No, 0/1). A diferencia de la
#' regresión lineal, el motor de cálculo estima los coeficientes mediante Máxima
#' Verosimilitud (MLE) y ofrece un diagnóstico especializado en clasificación.
#'
#' @details
#' El flujo interno de la función se divide en cinco etapas fundamentales:
#' \enumerate{
#'   \item **Ajuste Logit:** Transforma la variable respuesta mediante la función de enlace
#'         \eqn{logit(p) = ln(p / (1-p))}, permitiendo que la probabilidad se mantenga entre 0 y 1.
#'   \item **Conversión de Impacto:** Los coeficientes estimados se transforman automáticamente
#'         a **Odds Ratios (OR)** mediante \eqn{e^\beta}. Esto permite concluir cuánto
#'         multiplica cada variable la posibilidad de que ocurra el evento.
#'   \item **Evaluación de Bondad (McFadden):** Calcula el Pseudo-\eqn{R^2} comparando la
#'         verosimilitud del modelo ajustado contra un modelo nulo (sin predictores).
#'   \item **Métricas de Clasificación:** Genera la Matriz de Confusión basada en el \code{umbral}
#'         especificado y calcula la Precisión Global (Accuracy).
#'   \item **Diagnóstico de Influencia:** Evalúa si existen observaciones que deforman
#'         las probabilidades mediante el cálculo de la Distancia de Cook y el Apalancamiento (Leverage).
#' }
#'
#' @section Guía de las 9 Gráficas de Diagnóstico (`plot`):
#' Al ejecutar \code{plot(modelo)}, se despliega una secuencia de validación predictiva:
#' \describe{
#'   \item{\strong{1. Capacidad de Separación}}{Muestra si el modelo realmente distingue entre grupos. Buscamos "montañas" de densidad que no se traslapen.}
#'   \item{\strong{2. Matriz de Confusión Visual}}{Mapa de calor que resume aciertos (diagonal) y errores (falsos positivos/negativos).}
#'   \item{\strong{3. Odds Ratios (Forest Plot)}}{Visualiza la importancia de las variables. Los puntos a la derecha del 1 son factores de riesgo; a la izquierda, protectores.}
#'   \item{\strong{4. Curva de Calibración}}{Compara la probabilidad que "dice" el modelo contra la frecuencia real observada. Debe seguir la diagonal.}
#'   \item{\strong{5. Curva ROC (AUC)}}{Mide la potencia global de diagnóstico. Un AUC de 1.0 es una clasificación perfecta; 0.5 es puro azar.}
#'   \item{\strong{6. Distancia de Cook}}{Detecta casos atípicos en la clasificación que tienen un peso excesivo en los coeficientes.}
#'   \item{\strong{7. Trade-off de Umbrales}}{Permite ver cómo cambia la Sensibilidad y Especificidad al mover el punto de corte.}
#'   \item{\strong{8. Residuos de Pearson}}{Busca patrones sistemáticos en el error. Idealmente, la línea roja debe mantenerse cerca del cero.}
#'   \item{\strong{9. Leverage vs Residuos}}{Identifica observaciones con combinaciones de variables X muy raras que "jalonean" las probabilidades.}
#' }
#'
#' @param formula Una fórmula de R (ej. \code{presencia_plaga ~ humedad + altitud}).
#' @param data Un \code{data.frame} con las observaciones.
#' @param umbral Punto de corte para clasificar como "Evento (1)" (por defecto 0.5).
#'
#' @return Un objeto de clase \code{"logistico_greenreg"} que contiene:
#' \itemize{
#'   \item \code{modelo}: Objeto \code{glm} original.
#'   \item \code{odds_ratios}: Tabla con OR, intervalos de confianza y p-valores.
#'   \item \code{metricas}: Lista con Pseudo-R2, AIC y Accuracy.
#'   \item \code{confusion}: Tabla cruzada de Real vs. Predicho.
#'   \item \code{predicciones}: Dataframe con probabilidades y métricas de diagnóstico (leverage/cook).
#'   \item \code{notas}: Reporte interpretado con semáforos de validación técnica.
#' }
#'
#' @examples
#' # data("datos_roya_cafe")
#' # modelo_log <- logistico(presencia_roya ~ humedad_relativa_pct + altitud_msnm + manejo_sombra, data = datos_roya_cafe)
#' # modelo_log
#' # plot(modelo_log)
#'
#' @importFrom stats glm binomial fitted coef confint formula na.omit predict model.frame logLik aggregate pchisq cooks.distance hatvalues rstandard
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
#' @export
# Asegúrate de usar esta versión actualizada de tu función principal
logistico <- function(formula, data, umbral = 0.5) {

  # --- 1. Preparación y Ajuste ---
  # Dejamos que glm maneje el punto (.) y los NAs automáticamente
  modelo_base <- stats::glm(formula, data = data, family = stats::binomial(link = "logit"), na.action = stats::na.omit)
  resumen <- summary(modelo_base)

  # Extraemos los datos EXACTOS que usó el modelo (sin NAs)
  data_model <- modelo_base$model

  # Identificar la variable respuesta (siempre es la 1ra columna del model frame)
  var_y_name <- names(data_model)[1]
  y_real <- modelo_base$y # vector 0/1 real usado en el ajuste

  # --- 2. Interpretación de Coeficientes (Odds Ratios) ---
  coefs_log <- coef(modelo_base)
  coefs_or <- exp(coefs_log)
  cis_or <- exp(stats::confint(modelo_base)) # Intervalos de confianza para OR

  tabla_or <- cbind(OR = coefs_or, `2.5 %` = cis_or[,1], `97.5 %` = cis_or[,2], `P-valor` = resumen$coefficients[,4])

  # --- 3. Métricas de Bondad de Ajuste ---
  # Pseudo R2 de McFadden: 1 - (LogLik_Modelo / LogLik_Nulo)
  ll_modelo <- as.numeric(logLik(modelo_base))

  # Ajustamos modelo nulo sobre LOS MISMOS DATOS exactos (data_model)
  modelo_nulo <- stats::glm(as.formula(paste(var_y_name, "~ 1")), data = data_model, family = binomial)
  ll_nulo <- as.numeric(logLik(modelo_nulo))

  pseudo_r2 <- 1 - (ll_modelo / ll_nulo)

  # Matriz de Confusión y Accuracy
  probs_pred <- fitted(modelo_base)
  clase_pred <- ifelse(probs_pred > umbral, 1, 0)

  tabla_confusion <- table(Real = y_real, Predicho = clase_pred)
  accuracy <- sum(diag(tabla_confusion)) / sum(tabla_confusion)

  # --- 4. Notas Didácticas ---
  nota_r2 <- paste0("Pseudo R² (McFadden) = ", round(pseudo_r2, 3),
                    ". (Valores entre 0.2 y 0.4 se consideran excelentes en logística).")

  nota_acc <- paste0("Precisión Global (Accuracy): ", round(accuracy * 100, 1),
                     "%. El modelo clasificó correctamente este porcentaje de casos.")

  # Mensaje sobre Odds Ratios significativos
  significativos <- rownames(tabla_or)[tabla_or[, "P-valor"] < 0.05 & rownames(tabla_or) != "(Intercept)"]
  if (length(significativos) > 0) {
    ejemplo_var <- significativos[1]
    or_val <- tabla_or[ejemplo_var, "OR"]
    direccion <- ifelse(or_val > 1, "AUMENTA", "DISMINUYE")
    nota_interp <- paste0("Ejemplo de interpretación: La variable '", ejemplo_var,
                          "' es significativa. Su OR es ", round(or_val, 2),
                          ", lo que significa que ", direccion, " la probabilidad del evento.")
  } else {
    nota_interp <- "Ninguna variable parece tener un efecto estadísticamente significativo."
  }
  # --- 4.5 Semáforo de Validación (Diagnóstico Técnico) ---
  # A) Significancia del Modelo (Likelihood Ratio Test)
  p_modelo <- 1 - pchisq(modelo_base$null.deviance - modelo_base$deviance,
                         df = modelo_base$df.null - modelo_base$df.residual)

  nota_signif <- if(p_modelo < 0.05) "✅ OK: El modelo es significativamente mejor que el nulo." else "⚠️ ALERTA: El modelo no aporta información relevante (p > 0.05)."

  # B) Multicolinealidad (VIF)
  vif_vals <- if(requireNamespace("car", quietly = TRUE)) car::vif(modelo_base) else NA
  nota_vif <- if(any(vif_vals > 10, na.rm = TRUE)) "⚠️ ALERTA: Multicolinealidad alta (VIF > 10). Revisa las variables redundantes." else "✅ OK: No hay multicolinealidad severa."

  # C) Influencia (Puntos Cook > Umbral)
  cooksd <- stats::cooks.distance(modelo_base)
  umbral_cook <- 4 / length(cooksd)
  hay_influencia <- any(cooksd > umbral_cook)
  nota_influencia <- if(hay_influencia) "⚠️ ALERTA: Se detectaron casos influyentes (Distancia de Cook). Revisa los casos atípicos." else "✅ OK: No hay puntos influyentes críticos."

  # --- NUEVO: D) Leverage (Apalancamiento) ---
  leverage <- stats::hatvalues(modelo_base) # Calculamos h
  num_vars <- length(coefs_log) # número de parámetros (incluyendo intercepto)
  n_obs <- length(y_real)
  umbral_leverage <- 2 * num_vars / n_obs # Umbral común: 2*k/n
  hay_leverage_alto <- any(leverage > umbral_leverage)
  nota_leverage <- if(hay_leverage_alto) "⚠️ ALERTA: Hay observaciones con Leverage alto (puntos atípicos en X). Usa plot() para ver el gráfico." else "✅ OK: No hay puntos con apalancamiento extremo."

  # --- 5. Construcción del Objeto Final ---
  # Consolidamos TODO en una sola lista
  resultado <- list(
    formula = formula, # <--- Guardamos la fórmula (necesario para el print sección 1)
    modelo = modelo_base,
    variables = list(y = var_y_name),
    coeficientes = resumen$coefficients,
    odds_ratios = tabla_or,
    metricas = list(pseudo_r2 = pseudo_r2, aic = modelo_base$aic, accuracy = accuracy),
    confusion = tabla_confusion,
    predicciones = data.frame(
      obs = y_real,
      prob = probs_pred,
      clasif = clase_pred,
      leverage = leverage # <--- Guardamos leverage para el plot
    ),
    diag_leverage = list(umbral = umbral_leverage, hay_alto = hay_leverage_alto), # Info diagnóstica
    notas = list(
      r2 = nota_r2,
      acc = nota_acc,
      interp = nota_interp
    ),
    # Aquí unificamos el semáforo que creaste arriba, incluyendo el nuevo de leverage
    supuestos = c(nota_signif, nota_vif, nota_influencia, nota_leverage)
  )

  class(resultado) <- "logistico_greenreg"
  return(resultado)
}



#' Impresión para Regresión Logística
#' @export
print.logistico_greenreg <- function(x, ...) {
  cat("\n==========================================================\n")
  cat("               REGRESIÓN LOGÍSTICA \n")
  cat("==========================================================\n\n")

  # --- 1. MODELO Y ECUACIÓN ---
  cat("--- 1. MODELO Y ECUACIÓN ---\n")
  cat("Fórmula aplicada:", deparse(x$formula), "\n")
  # Extracción de coeficientes para mostrar la forma lineal del logit
  b <- round(coef(x$modelo), 3)
  vars <- names(b)
  ecuacion <- paste0("logit(p) = ", b[1])
  if(length(b) > 1) {
    for(i in 2:length(b)) {
      ecuacion <- paste0(ecuacion, " + (", b[i], " * ", vars[i], ")")
    }
  }
  cat("Ecuación estimada:\n", ecuacion, "\n\n")

  # --- 2. COEFICIENTES ESTIMADOS ---
  cat("--- 2. COEFICIENTES ESTIMADOS ---\n")
  cat("Estimaciones en escala Logit (log-odds):\n")
  stats::printCoefmat(x$coeficientes, digits = 4, signif.stars = TRUE)
  cat("\nInterpretación en Odds Ratios (OR):\n")
  print(round(x$odds_ratios, 4))
  cat(" GUÍA: OR > 1 Aumenta prob. | OR < 1 Disminuye prob.\n\n")

  # --- 3. EVALUACIÓN DEL MODELO ---
  cat("--- 3. EVALUACIÓN DEL MODELO ---\n")
  # Prueba de Razón de Verosimilitud (LRT)
  chi_cuad <- x$modelo$null.deviance - x$modelo$deviance
  df_chi <- x$modelo$df.null - x$modelo$df.residual
  p_val_lrt <- 1 - pchisq(chi_cuad, df_chi)

  cat("Prueba de Razón de Verosimilitud (LRT):\n")
  cat(" • Chi-cuadrado:", round(chi_cuad, 3), "| GL:", df_chi, "| p-valor:", format.pval(p_val_lrt), "\n")
  cat(" • Resultado:", if(p_val_lrt < 0.05) "Modelo Globalmente Significativo ✅" else "Modelo No Significativo ⚠️", "\n\n")

  # --- 4. BONDAD DE AJUSTE ---
  cat("--- 4. BONDAD DE AJUSTE ---\n")
  cat("• Log-Verosimilitud (LogLik):", round(logLik(x$modelo), 3), "\n")
  cat("• AIC (Criterio de Akaike): ", round(x$metricas$aic, 3), "\n")
  cat("•", x$notas$r2, "\n")
  cat("• Prueba de Wald: Revisar asteriscos en la sección 2 (Significancia individual).\n\n")

  # --- 5. MATRIZ DE CONFUSIÓN ---
  cat("--- 5. MATRIZ DE CONFUSIÓN (Real vs Predicho) ---\n")
  print(x$confusion)
  cat("\n INTERPRETACIÓN:\n")
  cat("Buscamos una diagonal principal (Real y Predicho coinciden) con valores altos.\n")
  cat("Los valores fuera de la diagonal representan Falsos Positivos y Falsos Negativos.\n")
  cat("•", x$notas$acc, "\n\n")

  # --- 6. VERIFICACIÓN DE SUPUESTOS ---
  cat("--- 6. VERIFICACIÓN DE SUPUESTOS ---\n")
  cat("• Linealidad en el Logit: Supuesta (Usar plot(x) para verificar gráficos de residuos).\n")
  cat("• Independencia: Observaciones independientes asumidas por diseño.\n")
  # Mostrar las alertas que calculamos en el semáforo (VIF, Influencia, y el NUEVO LEVERAGE)
  for (nota in x$supuestos) {
    cat(" ", nota, "\n")
  }
  cat("\n")

  # --- 7. RESUMEN DEL MODELO ---
  cat("--- 7. RESUMEN DEL MODELO ---\n")
  cat("Conclusión rápida:", x$notas$interp, "\n")

  if(any(grepl("⚠️", x$supuestos))) {
    cat("\n❗ ALERTA: El modelo presenta advertencias técnicas (Multicolinealidad, Puntos influyentes o Leverage).\n")
    cat(" Se recomienda usar plot(nombre_del_objeto) para ver los gráficos de diagnóstico.\n")
  }

  cat("\n==========================================================\n")
}

#' Generación de Gráficas para Regresión Logística
#' @export
plot.logistico_greenreg <- function(x, ...) {
  # Verificar dependencias básicas
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Necesitas instalar 'ggplot2'.")

  # --- CORRECCIÓN 1: Usar el nombre correcto de la lista ---
  df_preds <- x$predicciones  # Antes decía x$preds (que no existe)
  # ---------------------------------------------------------

  modelo   <- x$modelo

  # Usar el umbral si existe, sino 0.5
  umbral   <- if(!is.null(x$umbral)) x$umbral else 0.5

  # Tema visual
  mi_tema <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, color = "#2C3E50"),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#7F8C8D"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 10, face = "italic", color = "#555555", margin = ggplot2::margin(t = 10))
    )

  cat("Generando gráficas de diagnóstico logístico... (Presiona [Enter] para avanzar)\n")

  # --- GRÁFICA 1: SEPARACIÓN DE CLASES (DENSIDAD) ---
  g1 <- ggplot2::ggplot(df_preds, ggplot2::aes(x = prob, fill = factor(obs))) +
    ggplot2::geom_density(alpha = 0.6, color = "white") +
    ggplot2::scale_fill_manual(values = c("0" = "#3498DB", "1" = "#E74C3C"), name = "Clase Real") +
    ggplot2::labs(
      title = "1. Capacidad de Separación (Discriminación)",
      subtitle = "Distribución de las probabilidades predichas para cada grupo",
      x = "Probabilidad Predicha por el Modelo (0 a 1)",
      y = "Densidad",
      caption = "INTERPRETACIÓN:\n✅ BIEN: Dos montañas separadas. La azul a la izquierda (cerca de 0) y la roja a la derecha (cerca de 1).\n❌ MAL: Si las montañas están una encima de la otra (solapadas), el modelo está confundido."
    ) + mi_tema
  print(g1)
  readline(prompt = "Gráfica 1/9 > ")

  # --- GRÁFICA 2: MATRIZ DE CONFUSIÓN VISUAL ---
  pred_clase <- ifelse(df_preds$prob > umbral, 1, 0)
  tabla <- table(Real = df_preds$obs, Predicho = pred_clase)
  df_matriz <- as.data.frame(tabla)

  g2 <- ggplot2::ggplot(df_matriz, ggplot2::aes(x = Predicho, y = Real, fill = Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Freq), size = 8, color = "white", fontface = "bold") +
    ggplot2::scale_fill_gradient(low = "#95A5A6", high = "#2C3E50") +
    ggplot2::labs(
      title = paste("2. Matriz de Confusión (Umbral:", umbral, ")"),
      subtitle = "Conteo de Aciertos y Errores",
      caption = "INTERPRETACIÓN:\n✅ DIAGONAL (0-0 y 1-1): Son los Aciertos correctos.\n❌ CRUZADA (0-1 y 1-0): Son los Errores (Falsos Positivos y Falsos Negativos).\nObjetivo: Que los cuadros oscuros estén en la diagonal."
    ) + mi_tema
  print(g2)
  readline(prompt = "Gráfica 2/9 > ")

  # --- GRÁFICA 3: ODDS RATIOS (FOREST PLOT MEJORADO) ---
  # --- CORRECCIÓN 2: Verificar que 'broom' esté instalado ---
  if (requireNamespace("broom", quietly = TRUE)) {
    coefs <- broom::tidy(modelo, conf.int = TRUE, exponentiate = TRUE)
    coefs <- coefs[coefs$term != "(Intercept)", ]

    coefs$Tipo <- "No Significativo"
    coefs$Tipo[coefs$p.value < 0.05 & coefs$estimate > 1] <- "Riesgo (Aumenta Prob)"
    coefs$Tipo[coefs$p.value < 0.05 & coefs$estimate < 1] <- "Protector (Baja Prob)"

    g3 <- ggplot2::ggplot(coefs, ggplot2::aes(x = estimate, y = reorder(term, estimate), color = Tipo)) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      ggplot2::geom_pointrange(ggplot2::aes(xmin = conf.low, xmax = conf.high), size = 0.8) +
      ggplot2::scale_color_manual(values = c("Riesgo (Aumenta Prob)" = "#E74C3C",
                                             "Protector (Baja Prob)" = "#27AE60",
                                             "No Significativo" = "#95A5A6")) +
      ggplot2::labs(
        title = "3. Importancia de Variables (Odds Ratios)",
        subtitle = "Impacto de cada variable en la probabilidad del evento",
        x = "Odds Ratio (Escala Logarítmica)", y = "Variables",
        caption = "INTERPRETACIÓN:\n➡️ Derecha del 1: AUMENTA la probabilidad del evento (Riesgo).\n⬅️ Izquierda del 1: DISMINUYE la probabilidad (Protección).\n⚠️ GRIS (Cruza el 1): La variable NO sirve, estadísticamente no afecta."
      ) + mi_tema
    print(g3)
  } else {
    cat("⚠️ Instala el paquete 'broom' para ver la gráfica de Odds Ratios: install.packages('broom')\n")
  }
  readline(prompt = "Gráfica 3/9 > ")

  # --- GRÁFICA 4: CURVA DE CALIBRACIÓN ---
  df_preds$bin <- cut(df_preds$prob, breaks = seq(0, 1, 0.1), include.lowest = TRUE)

  # Agregamos na.rm=TRUE para evitar errores con valores faltantes
  calibracion <- aggregate(obs ~ bin, data = df_preds, mean, na.rm=TRUE)
  prob_medias <- aggregate(prob ~ bin, data = df_preds, mean, na.rm=TRUE)
  calibracion$prob_media <- prob_medias$prob

  g4 <- ggplot2::ggplot(calibracion, ggplot2::aes(x = prob_media, y = obs)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::geom_line(color = "#8E44AD", size = 1) +
    ggplot2::geom_point(color = "#8E44AD", size = 3) +
    ggplot2::xlim(0,1) + ggplot2::ylim(0,1) +
    ggplot2::labs(
      title = "4. Calibración del Modelo",
      subtitle = "¿Dice la verdad el modelo sobre sus probabilidades?",
      x = "Probabilidad que dijo el modelo",
      y = "Porcentaje Real de Eventos",
      caption = "INTERPRETACIÓN:\n✅ BIEN: Los puntos morados siguen la línea gris punteada.\n❌ MAL: Si la curva está muy por debajo de la línea, el modelo es 'exagerado'."
    ) + mi_tema
  print(g4)
  readline(prompt = "Gráfica 4/9 > ")

  # --- GRÁFICA 5: CURVA ROC ---
  if (requireNamespace("pROC", quietly = TRUE)) {
    roc_obj <- pROC::roc(df_preds$obs, df_preds$prob, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    roc_df <- data.frame(esp = roc_obj$specificities, sens = roc_obj$sensitivities)

    g5 <- ggplot2::ggplot(roc_df, ggplot2::aes(x = 1 - esp, y = sens)) +
      ggplot2::geom_line(color = "#D35400", size = 1.2) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      ggplot2::labs(
        title = paste0("5. Curva ROC (AUC = ", round(auc_val, 2), ")"),
        subtitle = "Capacidad global de diagnóstico",
        x = "Tasa de Falsos Positivos", y = "Sensibilidad (Tasa de Aciertos)",
        caption = "INTERPRETACIÓN:\n✅ BIEN: La curva naranja se acerca a la esquina superior izquierda.\n❌ MAL: Si la curva está cerca de la diagonal gris, el modelo es puro azar."
      ) + mi_tema
    print(g5)
  } else {
    cat("⚠️ Instala el paquete 'pROC' para ver la Curva ROC: install.packages('pROC')\n")
  }
  readline(prompt = "Gráfica 5/9 > ")

  # --- GRÁFICA 6: VALORES INFLUYENTES (COOK) ---
  cooksd <- stats::cooks.distance(modelo)
  n_datos <- length(cooksd)
  umbral_cook <- 4 / n_datos
  df_cook <- data.frame(Index = 1:n_datos, Cook = cooksd)
  df_cook$Label <- ifelse(df_cook$Cook > umbral_cook, as.character(df_cook$Index), "")

  g6 <- ggplot2::ggplot(df_cook, ggplot2::aes(x = Index, y = Cook)) +
    ggplot2::geom_bar(stat = "identity", fill = "#34495E", alpha = 0.8) +
    ggplot2::geom_hline(yintercept = umbral_cook, color = "red", linetype = "dashed") +
    ggplot2::geom_text(ggplot2::aes(label = Label), vjust = -0.5, size = 3) +
    ggplot2::labs(
      title = "6. Datos Influyentes (Distancia de Cook)",
      subtitle = "Detección de casos que distorsionan el modelo",
      x = "Número de Caso", y = "Influencia",
      caption = "INTERPRETACIÓN:\n⚠️ ALERTA: Las barras altas que cruzan la línea roja son casos atípicos.\nRevísalos en tu base de datos."
    ) + mi_tema
  print(g6)
  readline(prompt = "Gráfica 6/9 > ")

  # --- NUEVA GRÁFICA: BARRIDO DE UMBRALES (Trade-off) ---
  thresholds <- seq(0, 1, length = 100)
  sens <- sapply(thresholds, function(t) mean(df_preds$prob[df_preds$obs == 1] >= t, na.rm=TRUE))
  spec <- sapply(thresholds, function(t) mean(df_preds$prob[df_preds$obs == 0] < t, na.rm=TRUE))
  df_sweep <- data.frame(thresholds, sens, spec)

  g7 <- ggplot2::ggplot(df_sweep) +
    ggplot2::geom_line(ggplot2::aes(x = thresholds, y = sens, color = "Sensibilidad"), size = 1) +
    ggplot2::geom_line(ggplot2::aes(x = thresholds, y = spec, color = "Especificidad"), size = 1) +
    ggplot2::geom_vline(xintercept = umbral, linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("Sensibilidad" = "#E74C3C", "Especificidad" = "#2980B9")) +
    ggplot2::labs(
      title = "7. Trade-off: Umbral de Clasificación",
      subtitle = "Cómo el umbral afecta los errores",
      x = "Umbral", y = "Tasa", color = "Métrica",
      caption = "INTERPRETACIÓN:\nEl punto de cruce es el umbral 'balanceado'.\nSi necesitas detectar todos los casos, mueve el umbral a la izquierda (más sensibilidad)."
    ) + mi_tema
  print(g7)

  readline(prompt = "Gráfica 7/9 > ")

  # --- NUEVA GRÁFICA: RESIDUOS DE PEARSON ---
  # pearson_resid = (obs - prob) / sqrt(prob * (1 - prob))
  df_preds$pearson_res <- (df_preds$obs - df_preds$prob) / sqrt(df_preds$prob * (1 - df_preds$prob))

  g8 <- ggplot2::ggplot(df_preds, ggplot2::aes(x = prob, y = pearson_res)) +
    ggplot2::geom_point(alpha = 0.5, color = "#2C3E50") +
    ggplot2::geom_smooth(method = "loess", color = "#E74C3C", se = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = "8. Análisis de Residuos (Pearson)",
      subtitle = "¿Existen patrones ocultos en el error?",
      x = "Probabilidad Predicha", y = "Residuo de Pearson",
      caption = "INTERPRETACIÓN:\nSi la línea roja se aleja del 0, el modelo tiene problemas estructurales."
    ) + mi_tema
  print(g8)
  readline(prompt = "Gráfica 8/9 > ")
  # --- GRÁFICA 9: RESIDUOS VS APALANCAMIENTO (Estilo RLS) ---

  # 1. Extracción de métricas de diagnóstico del modelo GLM
  lev <- stats::hatvalues(x$modelo)
  res_std <- stats::rstandard(x$modelo) # Residuos de Pearson estandarizados
  n_datos <- length(lev)
  p_vars <- length(coef(x$modelo)) # Número de parámetros k (incluye intercepto)

  # Umbral típico de leverage: 2 * (k/n)
  umbral_lev <- (2 * p_vars) / n_datos

  # 2. Creación del dataframe para ggplot
  df_diag <- data.frame(
    index = 1:n_datos,
    leverage = lev,
    resid_std = res_std
  )

  # 3. Lógica de Etiquetado (Solo puntos críticos)
  # Etiquetar si el leverage supera el umbral O si el residuo es un outlier (>|2|)
  df_diag$Label <- ifelse(df_diag$leverage > umbral_lev | abs(df_diag$resid_std) > 2,
                          as.character(df_diag$index), "")

  # 4. Construcción de la Gráfica con ggplot2
  g9 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = leverage, y = resid_std)) +
    # Puntos con tamaño basado en su leverage
    ggplot2::geom_point(ggplot2::aes(size = leverage), color = "#C0392B", alpha = 0.6) +
    # Líneas de referencia
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "red", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = umbral_lev, linetype = "dashed", color = "blue", alpha = 0.5) +
    # Etiquetas inteligentes para evitar solapamiento
    ggrepel::geom_text_repel(ggplot2::aes(label = Label), size = 3.5, fontface = "bold") +
    # Etiquetas y títulos
    ggplot2::labs(
      title = "9. Diagnóstico de Influencia (Leverage vs Residuos)",
      subtitle = "Identificando puntos que 'jalonean' el ajuste logístico",
      x = "Apalancamiento (Leverage)",
      y = "Residuos Estandarizados",
      caption = paste0("INTERPRETACIÓN:\n",
                       "⚠️ EJE X: Puntos a la derecha de la línea azul tienen valores de X muy atípicos.\n",
                       "⚠️ EJE Y: Puntos fuera de las líneas rojas punteadas son Outliers (errores en clasificación).\n",
                       " CLAVE: Puntos arriba/abajo a la derecha están deformando las probabilidades del modelo.")
    ) +
    mi_tema +
    ggplot2::theme(legend.position = "none")
  print(g9)

  cat("✅ Secuencia de gráficas completada. Se ha generado la gráfica 9 de influencia.\n")
}

