# FolkDevilR — Simulate the Rise and Fall of Moral Panics
# ---------------------------------------------------------------------------
# This script provides two interchangeable engines:
#   1. A compartmental S‑P‑R model (Susceptible‑Panicked‑Recovered).
#   2. An agent‑based model (ABM) running on either a well‑mixed contact
#      process or a static igraph network.
# Core parameters:
#   beta, gamma, amp, me_intensity / me_prop / me_boost, align, pos_exp,
#   shock_prob, shock_size.
# ---------------------------------------------------------------------------

quiet_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install it first.", pkg), call. = FALSE)
  }
}

# ===========================================================================
# 1. COMPARTMENTAL ENGINE ----------------------------------------------------
# ===========================================================================

simulate_panic <- function(steps = 100, N = 10000, P0 = 10,
                           beta = 0.3, gamma = 0.1, amp = 1,
                           me_intensity = 0, align = 1, pos_exp = 0,
                           shock_prob = 0.01, shock_size = 50, seed = NULL) {
  quiet_require("tibble"); quiet_require("stats")
  if (!is.null(seed)) set.seed(seed)
  
  S <- N - P0; P <- P0; R <- 0
  S_vec <- P_vec <- R_vec <- integer(steps + 1)
  shock_vec <- logical(steps + 1)
  S_vec[1] <- S; P_vec[1] <- P; R_vec[1] <- R
  
  beta_eff <- beta * amp * (1 + me_intensity) * align
  
  for (t in seq_len(steps)) {
    trans_prob  <- 1 - exp(-beta_eff * P / N)
    transmissions <- stats::rbinom(1, S, trans_prob)
    gamma_eff   <- pmin(1, gamma + pos_exp)
    recoveries  <- stats::rbinom(1, P, gamma_eff)
    shock_now   <- stats::runif(1) < shock_prob
    newP        <- if (shock_now) shock_size else 0
    
    S <- S - transmissions
    P <- P + transmissions + newP - recoveries
    R <- R + recoveries
    
    S_vec[t + 1]   <- S
    P_vec[t + 1]   <- P
    R_vec[t + 1]   <- R
    shock_vec[t + 1] <- shock_now
  }
  
  tibble::tibble(t = 0:steps, S = S_vec, P = P_vec, R = R_vec, shock = shock_vec)
}

plot_panic <- function(sim, facet = FALSE) {
  quiet_require("ggplot2"); quiet_require("tidyr")
  sim_long <- tidyr::pivot_longer(sim, c(S, P, R), names_to = "state")
  g <- ggplot2::ggplot(sim_long, ggplot2::aes(t, value, colour = state)) +
    ggplot2::geom_line(size = 1) + ggplot2::theme_minimal() +
    ggplot2::labs(x = "Time", y = "Count", colour = "Compartment",
                  title = "FolkDevilR — S‑P‑R trajectory")
  if (facet) g <- g + ggplot2::facet_wrap(~state, scales = "free_y")
  g
}

run_grid <- function(param_grid, seed = NULL) {
  quiet_require("tibble"); quiet_require("purrr")
  df <- tibble::as_tibble(param_grid) %>%
    dplyr::mutate(.row = dplyr::row_number(),
                  me_intensity = dplyr::coalesce(me_intensity, 0),
                  align = dplyr::coalesce(align, 1),
                  pos_exp = dplyr::coalesce(pos_exp, 0))
  
  purrr::pmap_dfr(df, function(.row, steps, N, P0, beta, gamma, amp,
                               me_intensity, align, pos_exp,
                               shock_prob, shock_size, ...) {
    sim <- simulate_panic(steps, N, P0, beta, gamma, amp,
                          me_intensity, align, pos_exp,
                          shock_prob, shock_size,
                          seed = if (!is.null(seed)) seed + .row - 1 else NULL)
    tibble::tibble(sim = list(sim), .row = .row, steps = steps, N = N, P0 = P0,
                   beta = beta, gamma = gamma, amp = amp, me_intensity = me_intensity,
                   align = align, pos_exp = pos_exp,
                   shock_prob = shock_prob, shock_size = shock_size)
  })
}

# ===========================================================================
# 2. AGENT‑BASED ENGINE ------------------------------------------------------
# ===========================================================================

simulate_panic_abm <- function(steps = 100, N = 1000, P0 = 5,
                               beta = 0.2, gamma = 0.05, amp = 1,
                               me_prop = 0.05, me_boost = 0.5,
                               align = 1, pos_exp = 0,
                               avg_deg = 15, net = NULL, dynamic = TRUE,
                               shock_prob = 0.01, shock_size = 20,
                               seed = NULL) {
  quiet_require("tibble"); if (!dynamic || !is.null(net)) quiet_require("igraph")
  if (!is.null(seed)) set.seed(seed)
  
  if (is.null(net) && !dynamic) {
    p_edge <- avg_deg / (N - 1)
    net <- igraph::sample_gnp(N, p_edge, directed = FALSE)
  }
  if (!is.null(net)) {stopifnot(igraph::vcount(net) == N); dynamic <- FALSE}
  
  state <- factor(rep("S", N), levels = c("S", "P", "R"))
  state[sample.int(N, P0)] <- "P"
  me_flag <- rep(FALSE, N)
  if (me_prop > 0) me_flag[sample.int(N, ceiling(N * me_prop))] <- TRUE
  
  S_vec <- P_vec <- R_vec <- integer(steps + 1)
  shock_vec <- logical(steps + 1)
  S_vec[1] <- N - P0; P_vec[1] <- P0; R_vec[1] <- 0
  
  for (t in seq_len(steps)) {
    contacts <- vector("list", N)
    if (dynamic) {
      for (i in seq_len(N)) {
        k <- stats::rpois(1, avg_deg)
        if (k > 0) contacts[[i]] <- sample.int(N, k, FALSE)
      }
    } else {
      contacts <- igraph::adjacent_vertices(net, seq_len(N))
    }
    
    for (i in seq_len(N)) if (state[i] == "P" && length(contacts[[i]]) > 0) {
      eff_beta <- beta * amp * align * (1 + if (me_flag[i]) me_boost else 0)
      trg <- contacts[[i]]; if (inherits(trg, "igraph.vs")) trg <- as.integer(trg)
      sus <- trg[state[trg] == "S"]
      if (length(sus)) state[sus[stats::runif(length(sus)) < eff_beta]] <- "P"
    }
    
    gamma_eff <- pmin(1, gamma + pos_exp)
    recov <- which(state == "P" & stats::runif(N) < gamma_eff)
    state[recov] <- "R"
    
    shock_now <- stats::runif(1) < shock_prob
    if (shock_now) {
      sus <- which(state == "S"); if (length(sus)) state[sample(sus, min(shock_size, length(sus)))] <- "P"
    }
    
    S_vec[t + 1] <- sum(state == "S")
    P_vec[t + 1] <- sum(state == "P")
    R_vec[t + 1] <- sum(state == "R")
    shock_vec[t + 1] <- shock_now
  }
  
  history <- tibble::tibble(t = 0:steps, S = S_vec, P = P_vec, R = R_vec, shock = shock_vec)
  list(history = history, final_states = state)
}

plot_panic_abm <- function(obj, facet = FALSE) {
  plot_panic(obj$history, facet)
}

# # ===========================================================================
# # 3. USAGE EXAMPLES ----------------------------------------------------------
# # ===========================================================================
# # Remove the leading "#" to run, snippets assume you have already `source()`d this file.
# 
# ## ---- Simple S‑P‑R ---------------------------------------------------------
# sim1 <- simulate_panic(beta = 0.2, gamma = 0.08, steps = 120,
#                        N = 50000, P0 = 10, seed = 1)
# plot_panic(sim1)
# 
# ## ---- Parameter sweep with run_grid ---------------------------------------
# grid <- expand.grid(beta = c(0.1, 0.2, 0.3),
#                     amp  = c(1, 1.5),
#                     align = c(1, 1.3),
#                     gamma = 0.06,
#                     steps = 100,
#                     N     = 30000,
#                     me_intensity = 0,
#                     pos_exp = 0,
#                     P0    = 15,
#                     shock_prob = 0.01,
#                     shock_size = 80)
# sims <- run_grid(grid, seed = 100)
# 
# library(dplyr); library(purrr)
# summary <- sims %>%
#   mutate(peak = map_dbl(sim, ~max(.x$P))) %>%
#   select(beta, amp, align, peak)
# print(summary)
# 
# ## ---- Well‑mixed ABM -------------------------------------------------------
# abm_mix <- simulate_panic_abm(steps = 300, N = 2000, P0 = 10,
#                              beta = 0.18, gamma = 0.04,
#                              amp = 1.3,
#                              me_prop = 0.05, me_boost = 1,
#                              align = 1.3, pos_exp = 0.01,
#                              avg_deg = 25, dynamic = TRUE,
#                              shock_prob = 0.02, shock_size = 40,
#                              seed = 42)
# plot_panic_abm(abm_mix)
# 
# ## ---- ABM on a small‑world network ----------------------------------------
# 
# library(igraph)
# g_sw <- sample_smallworld(dim = 1, size = 3000, nei = 8, p = 0.05)
# abm_sw <- simulate_panic_abm(net = g_sw, N = vcount(g_sw), dynamic = FALSE,
#                              steps = 400, beta = 0.15, gamma = 0.04,
#                              amp = 1.2, me_prop = 0.05, me_boost = 0.8,
#                              align = 1.3, pos_exp = 0.01,
#                              shock_prob = 0.02, shock_size = 40,
#                              seed = 99)
# plot_panic_abm(abm_sw)









