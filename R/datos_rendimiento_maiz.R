#' Rendimiento de maíz bajo distintos niveles de fertilización
#'
#' Conjunto de datos que presenta el rendimiento del maíz en función
#' de diferentes variables climáticas y de manejo del cultivo.
#' Este dataset es útil para ejemplificar modelos de regresión lineal simple y múltiple.
#'
#' @format Un data frame con 250 observaciones y 4 variables:
#' \describe{
#'   \item{rendimiento_maiz_ton_ha}{Rendimiento del cultivo de maíz (ton/ha)}
#'   \item{precipitacion_mm}{Precipitación registrada durante el ciclo agrícola (mm)}
#'   \item{temperatura_c}{Temperatura media (°C)}
#'   \item{fertilizante_kg_ha}{Cantidad de fertilizante aplicado (kg/ha)}
#' }
#'
#' @details Los datos fueron generados de forma controlada, tomando como referencia
#' distribuciones y rangos reales reportados por SIAP e INEGI. Esta base se diseñó
#' específicamente para probar las funciones de regresión del paquete \pkg{CHAPIREG}.
#'
#' @examples
#' # Análisis completo con regresión múltiple
#' modelo_completo <- rlm(rendimiento_maiz_ton_ha ~ precipitacion_mm + temperatura_c + fertilizante_kg_ha,
#'                        data = rendimiento_maiz)
#'
#' # Análisis simple para ver solo el efecto de la lluvia
#' modelo_simple <- rls(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = rendimiento_maiz)
#'
#' # Aplicar Análisis de Componentes Principales (ACP)
#' predictores <- rendimiento_maiz[, c("precipitacion_mm", "temperatura_c", "fertilizante_kg_ha")]
#' modelo_acp(predictores, escala = TRUE)
#'
#' @source Elaboración propia con fines didácticos tomando referencia de datos del SIAP e INEGI.
"datos_rendimiento_maiz"
