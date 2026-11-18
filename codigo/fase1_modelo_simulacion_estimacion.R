# Numerales 1 a 3: modelo, simulación y estimación de tau_j = theta_j / theta_2

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
semilla <- 2048
set.seed(semilla)
dir.create("../resultados", showWarnings = FALSE)

# Numeral 1. Modelo estadístico
n <- 500
k <- 4
delitos <- c("asesinato", "hurto_mano_armada", "secuestro", "extorsion")
idx_tau <- c(1, 3, 4)
nombres_tau <- c("tau_1", "tau_3", "tau_4")

# Numeral 2. Parámetros generadores y simulación
theta_verdadero <- c(theta_1 = 0.24, theta_2 = 0.40, theta_3 = 0.08, theta_4 = 0.28)

stopifnot(all(theta_verdadero > 0),
          isTRUE(all.equal(sum(theta_verdadero), 1)))
stopifnot(which.max(theta_verdadero) == 2,
          which.min(theta_verdadero) == 3,
          abs(theta_verdadero[1] - theta_verdadero[4]) <= 0.05)

tau_verdadero <- setNames(theta_verdadero[idx_tau] / theta_verdadero[2], nombres_tau)

parametros_simulacion <- list(semilla = semilla, n = n, delitos = delitos,
                              theta = theta_verdadero, tau = tau_verdadero)
print(parametros_simulacion)

write.csv(data.frame(parametro = c(names(theta_verdadero), nombres_tau, "n", "semilla"),
                     delito    = c(delitos, paste(delitos[idx_tau], "/", delitos[2]), "", ""),
                     valor     = c(theta_verdadero, tau_verdadero, n, semilla)),
          "../resultados/parametros_simulacion.csv", row.names = FALSE)

# X_i ~ Categórica(theta) iid, i = 1, ..., n
casos_simulados <- sample(1:k, size = n, replace = TRUE, prob = theta_verdadero)

muestra <- data.frame(
  caso   = 1:n,
  codigo = casos_simulados,
  delito = factor(delitos[casos_simulados], levels = delitos)
)
stopifnot(nrow(muestra) == n)
head(muestra, 10)
write.csv(muestra, "../resultados/muestra_simulada_fase1.csv", row.names = FALSE)

# Estadístico suficiente: Y_j = sum_i I(X_i = j)
Y <- setNames(tabulate(muestra$codigo, nbins = k), delitos)
stopifnot(sum(Y) == n)
print(Y)
proporciones_muestrales <- Y / n

# Numeral 3. Estimadores
a0 <- rep(1, k)   # previa Dirichlet(1, 1, 1, 1)

estimar_tau <- function(Ymat, a0 = rep(1, nrow(as.matrix(Ymat)))) {
  Ymat <- as.matrix(Ymat)
  if (nrow(Ymat) != length(a0)) Ymat <- t(Ymat)
  A  <- Ymat + a0                                     # posterior Dirichlet(a0 + y)
  y2 <- Ymat[2, ]; a2 <- A[2, ]
  num <- t(Ymat[idx_tau, , drop = FALSE])
  aj  <- t(A[idx_tau, , drop = FALSE])
  mv  <- num / y2                                     # MV: Y_j / Y_2
  med <- aj / (a2 - 1); med[a2 <= 1] <- NA            # media posterior
  mna <- (aj / a2) * qf(0.5, 2 * aj, 2 * a2)          # mediana posterior
  mod <- (aj - 1) / (a2 + 1)                          # moda posterior
  out <- list(MV = mv, Bayes_media = med, Bayes_mediana = mna, Bayes_moda = mod)
  lapply(out, function(M) { colnames(M) <- nombres_tau; M })
}

# Máxima verosimilitud y momentos para theta
theta_hat_MV <- setNames(Y / n, paste0("theta_hat_", 1:k))
theta_hat_MM <- setNames(Y / n, paste0("theta_hat_", 1:k))
stopifnot(all.equal(theta_hat_MV, theta_hat_MM))

log_verosimilitud <- function(theta, y) sum(y * log(theta))
logv <- c(nucleo_en_theta_hat       = log_verosimilitud(theta_hat_MV, Y),
          nucleo_en_theta_verdadero = log_verosimilitud(theta_verdadero, Y),
          completa_en_theta_hat       = dmultinom(Y, prob = theta_hat_MV, log = TRUE),
          completa_en_theta_verdadero = dmultinom(Y, prob = theta_verdadero, log = TRUE))
print(logv)
stopifnot(logv[1] >= logv[2])

# Comprobación numérica del máximo, parametrizando con eta = log(tau)
neg_logv_eta <- function(eta) {
  th <- exp(c(eta[1], 0, eta[2:3])); th <- th / sum(th)
  -log_verosimilitud(th, Y)
}
opt <- optim(c(0, 0, 0), neg_logv_eta, method = "BFGS", control = list(reltol = 1e-14))
theta_optim <- exp(c(opt$par[1], 0, opt$par[2:3])); theta_optim <- theta_optim / sum(theta_optim)
stopifnot(max(abs(theta_optim - theta_hat_MV)) < 1e-5)
tau_optim <- setNames(exp(opt$par), nombres_tau)

# Bayes para theta: posterior Dirichlet(a0 + Y), marginales Beta(a_j, sum(a) - a_j)
a_post <- a0 + Y
theta_Bayes <- rbind(
  media   = a_post / sum(a_post),
  mediana = qbeta(0.5, a_post, sum(a_post) - a_post),
  moda    = (a_post - 1) / (sum(a_post) - k)
)
colnames(theta_Bayes) <- paste0("theta_", 1:k)

stopifnot(Y[2] > 0)
est_obs <- estimar_tau(Y, a0)
tau_hat <- est_obs$MV[1, ]
stopifnot(max(abs(tau_hat - tau_optim)) < 1e-4)
tabla_estimadores_tau <- rbind(verdadero = tau_verdadero,
                               t(sapply(est_obs, function(M) M[1, ])))

# Verificación de los estimadores de Bayes simulando la posterior
M_post <- 1e5
G <- sapply(a_post, function(a) rgamma(M_post, shape = a))
tau_post <- G[, idx_tau] / G[, 2]
colnames(tau_post) <- nombres_tau
moda_kernel <- apply(tau_post, 2, function(t) { d <- density(t, n = 2048); d$x[which.max(d$y)] })
verif_Bayes <- rbind(media_formula   = est_obs$Bayes_media[1, ],
                     media_simulada  = colMeans(tau_post),
                     mediana_formula = est_obs$Bayes_mediana[1, ],
                     mediana_simulada = apply(tau_post, 2, median),
                     moda_formula    = est_obs$Bayes_moda[1, ],
                     moda_kernel     = moda_kernel)

resumen_theta <- data.frame(
  delito          = delitos,
  theta_verdadero = as.vector(theta_verdadero),
  Y_j             = as.vector(Y),
  theta_hat_MV_MM = as.vector(theta_hat_MV),
  theta_Bayes_media = as.vector(theta_Bayes["media", ]),
  diferencia_MV   = as.vector(theta_hat_MV - theta_verdadero)
)
resumen_tau <- data.frame(
  parametro     = nombres_tau,
  cociente      = paste(delitos[idx_tau], "/", delitos[2]),
  formula_MV    = c("Y_1/Y_2", "Y_3/Y_2", "Y_4/Y_2"),
  tau_verdadero = as.vector(tau_verdadero),
  tau_hat_MV    = as.vector(tau_hat),
  diferencia    = as.vector(tau_hat - tau_verdadero)
)
print(resumen_theta, digits = 4)
print(resumen_tau, digits = 4)
print(round(theta_Bayes, 4))
print(round(tabla_estimadores_tau, 4))
print(round(verif_Bayes, 4))

resultados_fase1 <- list(parametros = parametros_simulacion, Y = Y,
                         theta_hat_MV = theta_hat_MV, theta_Bayes = theta_Bayes,
                         tau_hat = tau_hat, estimadores_tau = tabla_estimadores_tau,
                         log_verosimilitud = logv, a0 = a0)
saveRDS(resultados_fase1, "../resultados/resultados_fase1.rds")
