# Numerales 8 y 9: pruebas de hipótesis, potencia, errores tipo I y II y prueba UMP (alpha = 0.02)
# H0j: tau_j >= 1  vs  H1j: tau_j < 1, j = 1, 3, 4. Global: el hurto es el delito más frecuente.

if (!file.exists("fase2_delta_intervalos.R"))
  stop("Fije como directorio de trabajo la carpeta codigo (setwd).")
source("fase2_delta_intervalos.R", encoding = "UTF-8")

tau0 <- 1
pi0  <- tau0 / (1 + tau0)
z1   <- qnorm(1 - alpha)
pruebas <- c("Wald_original", "Wald_log", "UMP_aleatorizada", "UMP_no_aleatorizada")

# phi(y, s): probabilidad de rechazar H0 con Y_j = y y S_j = Y_j + Y_2 = s.
# Si Y_2 = 0, tau_hat no existe y no se rechaza.

# Numeral 8. Pruebas de Wald con el error estándar del método delta
phi_wald <- function(y, s, escala = c("original", "log")) {
  escala <- match.arg(escala); m <- s - y
  Z <- if (escala == "original") (y / m - tau0) / sqrt(y * s / m^3)
       else (log(y / m) - log(tau0)) / sqrt(1 / y + 1 / m)
  Z[y == 0 & m > 0] <- -Inf
  as.numeric(m > 0 & Z < -z1)
}

# Numeral 9. UMP condicional: Y_j | S_j = s ~ Bin(s, pi_j); rechaza si y < c(s)
# (la versión aleatorizada rechaza con probabilidad gamma(s) si y = c(s))
corte_ump <- function(s) {
  c_s <- qbinom(alpha, s, pi0)
  gam <- (alpha - pbinom(c_s - 1, s, pi0)) / dbinom(c_s, s, pi0)
  list(c = c_s, gamma = gam)
}
phi_ump <- function(y, s, aleatorizada = TRUE) {
  cg <- corte_ump(s)
  (y < cg$c) + aleatorizada * cg$gamma * (y == cg$c)
}
phi <- function(prueba, y, s) switch(prueba,
  Wald_original       = phi_wald(y, s, "original"),
  Wald_log            = phi_wald(y, s, "log"),
  UMP_aleatorizada    = phi_ump(y, s, TRUE),
  UMP_no_aleatorizada = phi_ump(y, s, FALSE))

# Potencia exacta: S_j ~ Bin(n, psi_j), psi_j = theta_j + theta_2, y Y_j | S_j ~ Bin(S_j, tau/(1 + tau))
rejilla_s <- data.frame(s = rep(0:n, 0:n + 1), y = sequence(0:n + 1) - 1)
PHI <- sapply(pruebas, function(p) phi(p, rejilla_s$y, rejilla_s$s))

potencia_exacta <- function(tau, psi) {
  pr <- dbinom(rejilla_s$s, n, psi) * dbinom(rejilla_s$y, rejilla_s$s, tau / (1 + tau))
  colSums(pr * PHI)
}

# Tamaño condicional de la UMP: alpha (aleatorizada) y a lo sumo alpha (no aleatorizada)
tam_cond <- sapply(1:n, function(s) {
  y <- 0:s; b <- dbinom(y, s, pi0)
  c(sum(b * phi_ump(y, s, TRUE)), sum(b * phi_ump(y, s, FALSE)))
})
stopifnot(max(abs(tam_cond[1, ] - alpha)) < 1e-12, all(tam_cond[2, ] <= alpha + 1e-12))

# Razón de verosimilitudes decreciente en y (Neyman-Pearson y Karlin-Rubin)
for (pi1 in c(0.05, 0.2, 0.35, 0.49)) for (s in c(10, 150, 320))
  stopifnot(all(diff(dbinom(0:s, s, pi1, log = TRUE) - dbinom(0:s, s, pi0, log = TRUE)) < 0))

# Dualidad: la UMP no aleatorizada rechaza si la cota superior exacta de tau_j es menor que 1
ok_m <- rejilla_s$s > rejilla_s$y
cota_sup <- with(rejilla_s[ok_m, ], (y + 1) * qf(1 - alpha, 2 * (y + 1), 2 * (s - y)) / (s - y))
stopifnot(all((cota_sup < tau0) == (PHI[ok_m, "UMP_no_aleatorizada"] == 1)))

# Pruebas con los datos
y_obs <- as.vector(Y[idx_tau]); m_obs <- as.vector(Y[2]); s_obs <- y_obs + m_obs
cg_obs <- corte_ump(s_obs)
prueba_datos <- data.frame(
  parametro = nombres_tau, delito = delitos[idx_tau],
  y_j = as.vector(y_obs), y_2 = as.vector(m_obs), s = as.vector(s_obs),
  tau_hat = as.vector(tau_hat), pi_hat = as.vector(y_obs / s_obs),
  Z_wald_original = as.vector((tau_hat - tau0) / ee_delta),
  Z_wald_log = as.vector(log(tau_hat / tau0) / ee_log),
  p_wald_original = as.vector(pnorm((tau_hat - tau0) / ee_delta)),
  p_wald_log = as.vector(pnorm(log(tau_hat / tau0) / ee_log)),
  c_ump = cg_obs$c, gamma_ump = cg_obs$gamma,
  p_ump_condicional = as.vector(pbinom(y_obs, s_obs, pi0)),
  row.names = NULL)
prueba_datos$decision <- t(sapply(seq_along(nombres_tau), function(l)
  sapply(pruebas, function(p) phi(p, y_obs[l], s_obs[l]))))
stopifnot(all(prueba_datos$decision %in% c(0, 1)))
print(prueba_datos, digits = 4)

# Prueba global (intersección-unión): se rechaza si se rechazan las tres individuales
decision_global <- apply(prueba_datos$decision, 2, min)
print(decision_global)

# Errores tipo I y II, estimados con theta_hat y reales con theta generador
psi_hat <- setNames(as.vector(theta_hat_MV[idx_tau] + theta_hat_MV[2]), nombres_tau)
psi_gen <- setNames(as.vector(theta_verdadero[idx_tau] + theta_verdadero[2]), nombres_tau)

errores <- do.call(rbind, lapply(seq_along(nombres_tau), function(l) {
  a_hat <- potencia_exacta(tau0, psi_hat[l])
  a_gen <- potencia_exacta(tau0, psi_gen[l])
  p_hat <- potencia_exacta(tau_hat[l], psi_hat[l])
  p_gen <- potencia_exacta(tau_verdadero[l], psi_gen[l])
  data.frame(parametro = nombres_tau[l], prueba = pruebas,
             tipo_I_psi_hat = a_hat, tipo_I_psi_gen = a_gen,
             potencia_estimada = p_hat, tipo_II_estimado = 1 - p_hat,
             potencia_real = p_gen, tipo_II_real = 1 - p_gen, row.names = NULL)
}))
print(errores, digits = 4)

# Tamaño: supremo en la frontera tau = 1, recorriendo el parámetro de estorbo psi
psi_rejilla <- seq(0.02, 0.98, by = 0.02)
tam_psi <- t(sapply(psi_rejilla, function(p) potencia_exacta(tau0, p)))
tamano <- data.frame(prueba = pruebas, tamano_sup_psi = apply(tam_psi, 2, max),
                     psi_del_sup = psi_rejilla[apply(tam_psi, 2, which.max)],
                     tamano_min_psi = apply(tam_psi, 2, min), row.names = NULL)
stopifnot(max(abs(tam_psi[, "UMP_aleatorizada"] - alpha)) < 1e-10,
          all(tam_psi[, "UMP_no_aleatorizada"] <= alpha + 1e-12))
print(tamano, digits = 4)

# Funciones de potencia con psi = psi_hat_j
tau_rejilla <- lapply(seq_along(nombres_tau), function(l)
  seq(0.5 * min(tau_verdadero[l], tau_hat[l]), 1.25, length.out = 150))
curvas <- lapply(seq_along(nombres_tau), function(l)
  t(sapply(tau_rejilla[[l]], function(t) potencia_exacta(t, psi_hat[l]))))
monotonia <- sapply(curvas, function(C) apply(C, 2, function(v) all(diff(v) <= 1e-12)))
dimnames(monotonia) <- list(pruebas, nombres_tau)
print(monotonia)
stopifnot(all(monotonia[c("UMP_aleatorizada", "UMP_no_aleatorizada"), ]))
ventaja_sobre_UMP <- sapply(seq_along(nombres_tau), function(l) {
  h1 <- tau_rejilla[[l]] < tau0
  apply(curvas[[l]][h1, , drop = FALSE] - curvas[[l]][h1, "UMP_aleatorizada"], 2, max)
})
colnames(ventaja_sobre_UMP) <- nombres_tau
print(signif(ventaja_sobre_UMP, 3))

# Monte Carlo: error tipo I con el MV restringido a H0, potencia con theta_hat y con theta generador
set.seed(20264)
M_ph <- 20000
theta_tilde <- lapply(idx_tau, function(j) {
  th <- theta_hat_MV; th[c(j, 2)] <- (Y[j] + Y[2]) / (2 * n); setNames(as.vector(th), delitos) })
names(theta_tilde) <- nombres_tau
escenarios <- c(list(theta_hat = theta_hat_MV, theta_generador = theta_verdadero),
                setNames(theta_tilde, paste0("theta_tilde_", nombres_tau)))
mc_pruebas <- do.call(rbind, lapply(names(escenarios), function(e) {
  Ym <- rmultinom(M_ph, n, escenarios[[e]])
  D <- lapply(seq_along(nombres_tau), function(l)
    sapply(pruebas, function(p) phi(p, Ym[idx_tau[l], ], Ym[idx_tau[l], ] + Ym[2, ])))
  rbind(data.frame(escenario = e, parametro = rep(nombres_tau, each = length(pruebas)),
                   prueba = pruebas, tasa_rechazo = unlist(lapply(D, colMeans)), row.names = NULL),
        data.frame(escenario = e, parametro = "global_IUT", prueba = pruebas,
                   tasa_rechazo = colMeans(D[[1]] * D[[2]] * D[[3]]), row.names = NULL))
}))
mc_pruebas$error_mc <- sqrt(mc_pruebas$tasa_rechazo * (1 - mc_pruebas$tasa_rechazo) / M_ph)
print(mc_pruebas, digits = 4)

# Contraste exacto vs Monte Carlo
cruce <- merge(subset(mc_pruebas, parametro != "global_IUT"),
  rbind(transform(errores[, c("parametro", "prueba")], escenario = "theta_hat", exacto = errores$potencia_estimada),
        transform(errores[, c("parametro", "prueba")], escenario = "theta_generador", exacto = errores$potencia_real)))
cruce_tilde <- subset(mc_pruebas, grepl("theta_tilde", escenario) &
                        parametro == sub("theta_tilde_", "", escenario))
cruce_tilde$exacto <- errores$tipo_I_psi_hat[match(paste(cruce_tilde$parametro, cruce_tilde$prueba),
                                                   paste(errores$parametro, errores$prueba))]
cruce <- rbind(cruce, cruce_tilde[, names(cruce)])
cruce$z_dif <- (cruce$tasa_rechazo - cruce$exacto) / pmax(cruce$error_mc, 1 / M_ph)
stopifnot(all(abs(cruce$z_dif) < 4.5))

write.csv(prueba_datos, "../resultados/n8_prueba_con_datos.csv", row.names = FALSE)
write.csv(errores, "../resultados/n8_errores_tipo_I_II.csv", row.names = FALSE)
write.csv(tamano, "../resultados/n8_tamano_sup_psi.csv", row.names = FALSE)
write.csv(mc_pruebas, "../resultados/n8_monte_carlo_potencia.csv", row.names = FALSE)

col_p <- c(Wald_original = "grey45", Wald_log = colores[["Bayes_media"]],
           UMP_aleatorizada = colores[["MV"]], UMP_no_aleatorizada = colores[["Bayes_mediana"]])
lty_p <- c(Wald_original = 3, Wald_log = 4, UMP_aleatorizada = 1, UMP_no_aleatorizada = 2)
etq_p <- c("Wald escala original", "Wald escala log", "UMP aleatorizada", "UMP no aleatorizada")

png("../resultados/n8_potencia.png", width = 1500, height = 540, res = 130)
par(mfrow = c(1, 3), mar = c(4.5, 4.2, 3, 1))
for (l in seq_along(nombres_tau)) {
  matplot(tau_rejilla[[l]], curvas[[l]], type = "l", lty = lty_p, col = col_p, lwd = 2,
          ylim = c(0, 1), xlab = expression(tau), ylab = expression(beta(tau)),
          main = sprintf("%s: potencia, psi = %.3f", nombres_tau[l], psi_hat[l]))
  abline(h = alpha, col = "grey60"); abline(v = tau0, lty = 2)
  abline(v = tau_hat[l], lwd = 1.5); abline(v = tau_verdadero[l], lty = 4, col = "grey40", lwd = 2)
  if (l == 1) legend("topright", bty = "n", cex = 0.78, legend = c(etq_p, "tau_hat", "tau generador"),
                     col = c(col_p, "black", "grey40"), lty = c(lty_p, 1, 4), lwd = 2)
}
invisible(dev.off())

png("../resultados/n8_tamano_vs_psi.png", width = 1100, height = 560, res = 130)
par(mar = c(4.5, 4.5, 3, 1))
matplot(psi_rejilla, tam_psi, type = "l", lty = lty_p, col = col_p, lwd = 2,
        ylim = c(0, max(tam_psi) * 1.05), xlab = expression(psi[j] == theta[j] + theta[2]),
        ylab = "P(rechazar H0 | tau = 1)", main = "Error tipo I en la frontera de H0 según el estorbo psi")
abline(h = alpha, col = "grey60")
abline(v = psi_hat, lty = 3, col = "grey50")
legend("topright", bty = "n", cex = 0.8, legend = c(etq_p, "psi_hat (j = 1, 3, 4)"),
       col = c(col_p, "grey50"), lty = c(lty_p, 3), lwd = c(2, 2, 2, 2, 1))
invisible(dev.off())

resultados_fase3 <- list(tau0 = tau0, alpha = alpha, z1 = z1, prueba_datos = prueba_datos,
                         decision_global = decision_global, errores = errores, tamano = tamano,
                         mc_pruebas = mc_pruebas, cruce_exacto_mc = cruce, M_ph = M_ph,
                         ventaja_sobre_UMP = ventaja_sobre_UMP,
                         psi_hat = psi_hat, psi_gen = psi_gen)
saveRDS(resultados_fase3, "../resultados/resultados_fase3.rds")
