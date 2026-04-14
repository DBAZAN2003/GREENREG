#' Presencia de roya del café en diferentes localidades
#'
#' Base de datos simulada que presenta la presencia/ausencia de roya del café
#' junto con variables ambientales y de manejo agronómico relevantes.
#' Es útil para ejercicios de regresión logística.
#'
#' @format Un data frame con 300 observaciones y 4 variables:
#' \describe{
#'   \item{presencia_roya}{Presencia (1) o ausencia (0) de roya del café}
#'   \item{humedad_relativa_pct}{Humedad relativa promedio (%)}
#'   \item{altitud_msnm}{Altitud de la finca (m s.n.m.)}
#'   \item{manejo_sombra}{Porcentaje de cobertura de sombra en la plantación}
#' }
#'
#' @details Los datos fueron generados de forma controlada utilizando parámetros
#' similares a los observados por SIAP y estudios agronómicos del INIFAP.
#' Esta base fue creada para ser utilizada en funciones de regresión logística del paquete \pkg{CHAPIREG}.
#'
#' @examples
#' # Predecir la presencia o ausencia de roya
#' modelo_roya <- logistico(presencia_roya ~ humedad_relativa_pct + altitud_msnm + manejo_sombra,
#'                          data = roya_cafe)
#'
#' @source Elaboración propia con fines didácticos basada en información de SIAP e INIFAP.
"datos_roya_cafe"
