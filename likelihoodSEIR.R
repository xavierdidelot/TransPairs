likelihoodSEIR = function (tij,ti,tj) {
  ti=ti-tij
  tj=tj-tij
  tij=0

  a=365/5 #rate of coalescence, in unit 1/year (1/Ne*g)
  g1=365/4 #inverse of mean latency period, in unit 1/year (1/g1years)
  g2=365/6 #inverse of mean infectious period, in unit 1/year (1/(g2days/365))

  m=min(ti,tj)
  # Following equation is from Maple
  l = g1 * (g2 ^ 2) * a * (0.4e1 * exp((a * m + g1 * tj + g1 * m + 2 * g2 * tij)) *
  (g1 ^ 2) * g2 + exp((a * tij + g2 * m + g2 * tij + g1 * tj + g1 * m)) *
  (a ^ 2) * g2 + exp((a * tij + g2 * m + g1 * tij + g2 * tj + g1 * m)) *
  (g1 ^ 3) - exp((a * tij + g2 * tij + 2 * g1 * m + g2 * tj)) * (a ^ 2) *
  g2 + exp((a * tij + 2 * g2 * m + g1 * tij + g1 * tj)) * g1 * (a ^ 2) -
  exp((a * m + g2 * tj + g1 * m + g2 * tij + g1 * tij)) * (g1 ^ 3) -
  exp((a * tij + 2 * g2 * m + g1 * tij + g1 * tj)) * g1 * a * g2 + a *
  exp((a * m + g1 * tj + g1 * m + 2 * g2 * tij)) * (g1 ^ 2) + 0.2e1 * a *
  exp((a * m + g1 * tj + g1 * m + 2 * g2 * tij)) * (g2 ^ 2) + (a ^ 2) *
  exp((a * m + g1 * tj + g1 * m + 2 * g2 * tij)) * g1 - 0.2e1 *
  exp((a * m + g1 * tj + g1 * m + 2 * g2 * tij)) * g1 * (g2 ^ 2) - 0.2e1 * a *
  exp((a * m + g2 * tj + g1 * m + g2 * tij + g1 * tij)) * (g1 ^ 2) - (a ^ 2) *
  exp((a * m + g1 * tj + g1 * m + 2 * g2 * tij)) * g2 + 0.5e1 * a *
  exp((a * m + g2 * tj + g1 * m + g2 * tij + g1 * tij)) * g1 * g2 - 0.7e1 *
  exp((a * tij + g2 * tij + 2 * g1 * m + g2 * tj)) * a * g1 * g2 - 0.3e1 *
  a * exp((a * m + g1 * tj + g1 * m + 2 * g2 * tij)) * g1 * g2 + 0.4e1 *
  exp((a * tij + g2 * m + g2 * tij + g1 * tj + g1 * m)) * a * g1 * g2 + 0.2e1 *
  exp((a * tij + g2 * m + g1 * tij + g2 * tj + g1 * m)) * g1 * a * g2 - (a ^ 2) *
  exp((a * m + g2 * tj + g1 * m + g2 * tij + g1 * tij)) * g1 - 0.3e1 * a *
  exp((a * m + g2 * tj + g1 * m + g2 * tij + g1 * tij)) * (g2 ^ 2) +
  exp((a * tij + 2 * g2 * m + g1 * tij + g1 * tj)) * (g1 ^ 2) * g2 -
  exp((a * tij + g2 * m + g1 * tij + g2 * tj + g1 * m)) * g1 * (a ^ 2) -
  exp((a * tij + 2 * g2 * m + g1 * tij + g1 * tj)) * (g1 ^ 2) * a + (a ^ 2) *
  exp((a * m + g2 * tj + g1 * m + g2 * tij + g1 * tij)) * g2 - 0.2e1 *
  exp((a * m + g1 * tj + g1 * m + 2 * g2 * tij)) * (g1 ^ 3) + 0.2e1 *
  exp((a * m + g2 * tj + g1 * m + g2 * tij + g1 * tij)) * (g2 ^ 3) - 0.2e1 *
  exp((a * tij + g2 * tij + 2 * g1 * m + g2 * tj)) * (g2 ^ 3) + 0.2e1 *
  exp((a * tij + g2 * m + g2 * tij + g1 * tj + g1 * m)) * (g1 ^ 3) + 0.2e1 *
  exp((a * tij + g2 * m + g2 * tij + g1 * tj + g1 * m)) * (g2 ^ 2) * g1 - 0.5e1 *
  exp((a * tij + g2 * m + g2 * tij + g1 * tj + g1 * m)) * (g1 ^ 2) * g2 - 0.2e1 *
  exp((a * tij + g2 * m + g2 * tij + g1 * tj + g1 * m)) * a * (g2 ^ 2) - 0.2e1 *
  exp((a * tij + g2 * m + g2 * tij + g1 * tj + g1 * m)) * (a ^ 2) * g1 - 0.2e1 *
  exp((a * tij + g2 * m + g1 * tij + g2 * tj + g1 * m)) * (g1 ^ 2) * g2 + 0.3e1 *
  exp((a * tij + g2 * tij + 2 * g1 * m + g2 * tj)) * a * (g2 ^ 2) + 0.2e1 *
  exp((a * tij + g2 * tij + 2 * g1 * m + g2 * tj)) * (a ^ 2) * g1 + 0.5e1 * exp((a *
  tij + g2 * tij + 2 * g1 * m + g2 * tj)) * (g2 ^ 2) * g1 + 0.2e1 * exp((a * tij + g2 *
  tij + 2 * g1 * m + g2 * tj)) * a * (g1 ^ 2) - 0.2e1 * exp((a * tij + g2 * tij + 2 *
  g1 * m + g2 * tj)) * (g1 ^ 2) * g2 - 0.5e1 * exp((a * m + g2 * tj + g1 * m + g2 *
  tij + g1 * tij)) * g1 * (g2 ^ 2) + 0.4e1 * exp((a * m + g2 * tj + g1 * m + g2 * tij +
  g1 * tij)) * (g1 ^ 2) * g2) * exp((-a * m - g2 * ti - g1 * tj - g2 * tj - g1 * m)) /
  (a ^ 3 * g1 ^ 2 - a * g1 ^ 4 - 4 * g1 ^ 3 * g2 ^ 2 + g1 ^ 4 * g2 + 5 * g2 ^ 3 *
  g1 ^ 2 + a ^ 3 * g2 ^ 2 - 3 * a ^ 2 * g2 ^ 3 + 2 * g2 ^ 4 * a - 2 * g2 ^ 4 *
  g1 - 3 * a ^ 2 * g1 ^ 2 * g2 + 4 * a * g1 ^ 3 * g2 - 3 * g1 ^ 2 * g2 ^ 2 * a - 2 *
  a ^ 3 * g1 * g2 + 6 * a ^ 2 * g1 * g2 ^ 2 - 2 * a * g2 ^ 3 * g1)
  if (is.nan(l)) l=0
  return(l)
}

