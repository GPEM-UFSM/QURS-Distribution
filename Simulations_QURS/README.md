# QURS distribution: Monte Carlo simulation scenarios

The simulation evaluates maximum likelihood estimation of the quantile
parameter `mu_tau` in the Quantile Unit Rational Sine (QURS) distribution.
The probability level `tau` is fixed within each configuration.

## Simulation design

| Scenario | True `mu_tau` | Configuration IDs | Configurations | Valid replications |
|:---------|--------------:|:------------------|---------------:|-------------------:|
| C1 | 0.10 | 1–20 | 20 | 100,000 |
| C2 | 0.30 | 21–40 | 20 | 100,000 |
| C3 | 0.50 | 41–60 | 20 | 100,000 |
| C4 | 0.70 | 61–80 | 20 | 100,000 |
| C5 | 0.90 | 81–100 | 20 | 100,000 |
| **Total** | | | **100** | **500,000** |

Within each scenario, all four sample sizes (`n = 50, 100, 200, 400`)
were crossed with all five quantile levels (`tau = 0.10, 0.25, 0.50,
0.75, 0.90`). Thus, each combination of scenario, sample size, and
quantile level defines one configuration. The IDs within each scenario
are ordered by increasing `tau`, with `n` varying fastest.

The article reports scenario C3. The remaining scenarios are reported in
the Supplementary Material.

## Data generation and estimation

- For each replication, an independent sample of size `n` was drawn by
  inverse transform sampling at the specified `mu_tau` and `tau`.
- Maximum likelihood estimation was performed with
  `Rsolnp::gosolnp()` (Rsolnp 2.0.1). The optimizer estimates `mu_tau`
  directly, and the power parameter is obtained from
  `alpha = log((2 / pi) * asin(tau / (2 - tau))) / log(mu_tau)`.
- The optimizer used `n.restarts = 5`, `n.sim = 30`, and a numerical
  boundary margin of `1e-8`. Its control settings were
  `outer.iter = 400`, `inner.iter = 800`, `tol = 1e-8`, and
  `delta = 1e-7`. The base random seed was `20261016`.
- The simulation retained 5,000 valid replications per configuration.
  Acceptance required valid estimation and inference, including finite
  positive standard errors and valid Wald and logit-scale confidence
  intervals. The aggregate run log records 500,000 attempts and 500,000
  accepted replications, with zero discards.

## Reported measures

For each configuration, the study summarizes the mean, bias, empirical
standard deviation, and mean squared error of the estimates of `mu_tau`,
as well as empirical coverage of the 95% confidence intervals. The
simulation summary records coverage for both the Wald interval on the
original parameter scale and the interval constructed on the logit scale.

## Reproduction

The simulation depends on R and the `Rsolnp` package. The archived
study used Rsolnp 2.0.1. Install the package with
`install.packages("Rsolnp")` and use the archived package version for
exact reproduction where possible.

| Script | Purpose |
|:-------|:--------|
| `01_FUNCOES_QURS_GOSOLNP.R` | Density, CDF, quantiles, simulation, maximum likelihood estimation, and confidence intervals |
| `03_SIMULACAO_QURS_GOSOLNP_SELECAO_CENARIOS.R` | Scenario grid, seeded replications, validation, checkpoints, and CSV summaries |
| `05_TABELAS_GRAFICOS_MONTE_CARLO_QURS_GOSOLNP.R` | CSV and LaTeX tables, boxplots, and density plots from completed results |

Run the following commands from this folder. The first runs a small test
with 20 valid replications at C3, `tau = 0.50`, and `n = 50`. The second
generates tables and figures for that test.

```sh
Rscript -e 'options(QURS.mode="piloto", QURS.scenarios="C3", QURS.tau=0.50, QURS.n=50, QURS.parallel=FALSE); source("03_SIMULACAO_QURS_GOSOLNP_SELECAO_CENARIOS.R")'
Rscript 05_TABELAS_GRAFICOS_MONTE_CARLO_QURS_GOSOLNP.R RESULTADOS_MONTE_CARLO/piloto_gosolnp_v1_validas_R00020
```

Run the full 100-configuration study and generate its tables and figures
with the following commands. The full run can take several hours.

```sh
Rscript 03_SIMULACAO_QURS_GOSOLNP_SELECAO_CENARIOS.R
Rscript 05_TABELAS_GRAFICOS_MONTE_CARLO_QURS_GOSOLNP.R
```

Rerun the same simulation command to resume an interrupted run. Use a
fresh output directory when changing settings that affect the checkpoint
signature, including the helper script or the Rsolnp version. Set the
R option `QURS.output_root` to direct a run to another directory. The archived
execution settings show a selection starting at C5. The aggregate results
cover all five scenarios.

For the standardized density plots, each configuration is centered at
its empirical mean and scaled by its empirical standard deviation. The
summary tables and file names retain the archived schema for
compatibility with the study outputs.

The completed run is stored in
`RESULTADOS_MONTE_CARLO/final_gosolnp_v1_validas_R05000/`:

| File | Contents |
|:-----|:---------|
| `Grade_Completa_QURS.csv` | Full grid and configuration IDs |
| `Controle_Global_QURS.csv` | Attempts, accepted replications, discards, and completion status |
| `Resumo_Monte_Carlo_QURS.csv` | Monte Carlo summaries by configuration |
| `Configuracao_Execucao.txt` | Recorded run settings and package versions |

