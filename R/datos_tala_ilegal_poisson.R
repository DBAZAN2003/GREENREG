#' Registros de tala ilegal en distintas regiones
#'
#' Conjunto de datos que contiene el número de eventos de tala ilegal observados
#' junto con variables explicativas relacionadas con la densidad poblacional,
#' superficie forestal y esfuerzo de vigilancia. Ideal para ejemplificar la aplicación
#' de modelos de regresión Poisson.
#'
#' @format Un data frame con 200 observaciones y 4 variables:
#' \describe{
#'   \item{numero_talas}{Número de eventos de tala ilegal observados}
#'   \item{densidad_poblacion}{Densidad poblacional de la región (hab/km²)}
#'   \item{superficie_forestal_ha}{Superficie forestal en hectáreas}
#'   \item{dias_vigilancia}{Número de días de vigilancia reportados}
#' }
#'
#' @details Los datos fueron generados de manera simulada con base en valores reportados
#' por CONAFOR y SEMARNAT. Esto permite realizar ejercicios reproducibles y seguros
#' sin exponer datos sensibles o confidenciales.
#'
#' @examples
#' # Modelar el número de talas ilegales
#' modelo_conteo <- poisson(numero_talas ~ densidad_poblacion + superficie_forestal_ha + dias_vigilancia,
#'                          data = datos_tala_ilegal_poisson)
#'
#' @source Elaboración propia con fines didácticos tomando referencia de CONAFOR e INEGI.
"datos_tala_ilegal_poisson"

