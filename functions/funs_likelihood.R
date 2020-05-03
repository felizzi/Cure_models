##############################################################
# Creates likelihod function for estimating cure fractions   #
##############################################################

# Likelihood function for uninformed approach (calc cure fractions) -------
ll.mix <- function(table, parms, time_col, cen_col, rate_col, obj) {
  # Specify containers
  time_v <- as.numeric(table[, time_col])
  cen_v  <- as.numeric(table[, cen_col])
  rate_v <- as.numeric(table[, rate_col])
  pi_    <- parms[length(parms)]

  f_t <- do.call(obj$dfns$d, args = c(
    as.list(parms[1:length(parms) - 1]),
    list(x = time_v, log = FALSE)
  ))
  S_t <- do.call(obj$dfns$p, args = c(
    as.list(parms[1:length(parms) - 1]),
    list(q = time_v, lower.tail = FALSE)
  ))
  
  ## Formulae (see paper)
  h_ <- rate_v + ((1 - pi_) * f_t) / (pi_ + (1 - pi_) * S_t)
  s_ <- pi_ + (1 - pi_) * S_t
  ret_val <- sum(cen_v * log(h_) + log(s_))
  
  ## Return object
  return(-ret_val)
}


# Likelihood function for informed approach (cure fractions as inp --------
ll.mix_f <- function(table, parms, pi_, time_col, cen_col, rate_col, obj) {
  # Specify containers
  time_v <- as.numeric(table[, time_col])
  cen_v  <- as.numeric(table[, cen_col])
  rate_v <- as.numeric(table[, rate_col])
  # pi_  <- parms[length(parms)];

  # shape <- parms[1];  scale <- parms[2]; pi_ <- parms[3];
  f_t <- do.call(obj$dfns$d, args = c(
    as.list(parms[1:length(parms)]),
    list(x = time_v, log = FALSE)
  ))
  S_t <- do.call(obj$dfns$p, args = c(
    as.list(parms[1:length(parms)]),
    list(q = time_v, lower.tail = FALSE)
  ))
  
  ## Formulae (see paper)
  h_ <- rate_v + ((1 - pi_) * f_t) / (pi_ + (1 - pi_) * S_t)
  s_ <- pi_ + (1 - pi_) * S_t
  ret_val <- sum(cen_v * log(h_) + log(s_))
  
  ## Return object
  return(-ret_val)
}
