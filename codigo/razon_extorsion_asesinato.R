# Cálculos complementarios del artículo con los conteos observados y = (124, 188, 44, 144)

y <- c(124, 188, 44, 144); alpha <- 0.02; z <- qnorm(1 - alpha / 2)

# Extorsión frente a asesinato, theta_4 / theta_1: intervalo score condicional
s <- y[4] + y[1]; p <- y[4] / s
ext <- (p + z^2 / (2 * s) + c(-1, 1) * z * sqrt(p * (1 - p) / s + z^2 / (4 * s^2))) / (1 + z^2 / s)
cat("Razón estimada extorsión / asesinato:", round(y[4] / y[1], 4), "\n")
cat("Intervalo score condicional del 98 %:", round(ext / (1 - ext), 4), "\n")

# Umbrales de gravedad 1 / tau_j con los límites del intervalo score condicional
L <- c(tau_1 = 0.5044, tau_3 = 0.1589, tau_4 = 0.5924); U <- c(0.8626, 0.3447, 0.9903)
print(round(rbind("1 / límite inferior" = 1 / L, "1 / límite superior" = 1 / U), 2))

# Potencia exacta de la UMP no aleatorizada para la extorsión (tau_4 = 0.70, psi_4 = 0.68) según n
potencia <- function(n, tau = 0.70, psi = 0.28 + 0.40) {
  s <- 0:n
  sum(dbinom(s, n, psi) * pbinom(qbinom(alpha, s, 0.5) - 1, s, tau / (1 + tau)))
}
cat("Potencia con n = 500:", round(potencia(500), 4), "; con n = 1000:", round(potencia(1000), 4), "\n")

# Intervalos score expresados frente a los 188 hurtos observados
print(round(rbind("mínimo frente a 188 hurtos" = L * y[2], "máximo frente a 188 hurtos" = U * y[2])))
