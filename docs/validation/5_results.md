# Risultati della fase 5 — Port del motore statistico (Caso A)

Data: 2026-11-22. Previsioni e soglie: `docs/validation/5_predictions.md` versione 1, committata
in `c9e7669` **prima** dei run. Dopo i run le soglie non sono state toccate.

**Esito.** Il run completo `c1` passa **tutte** le 71 soglie di §1–§5. Lo scarto più grande vale
1.75σ (in unità della soglia, cioè 3σ = limite): nessun fallimento isolato fra 3σ e 4σ, quindi
non serve il run diagnostico con il seme `20260908`. Il port C++ del Caso A riproduce la
calibrazione del toy numero per numero entro l'errore statistico. Il target di tempo D35 resta
mancato (§6).

Etichette: **(a)** fatto stabilito; **(b)** stima con ipotesi esplicite; **(c)** verificato con
Monte Carlo; **(d)** ipotesi da mettere alla prova.

## Run

Build `release` a `c9e7669` (`git describe` v0.6.0-13-gc9e7669), seme 20260907,
σ_Ep = 0.10 MeV, σ_θ = 0.08 rad, griglie 500 × 500, livelli 0.68 / 0.90 / 0.95. Ogni comando è
passato da `scripts/run_capped.sh -m 5G`, uno alla volta, dopo `scripts/data_root.sh check`.
Gli output stanno sull'SSD, in `stats/<calibration_id>/calibration.json`.

| `calibration_id` | Opzioni | RSS massimo | Tempo |
|---|---|---|---|
| `pilot` | `--pilot` | 16 MB | 4.1 s |
| `c1` | — | 19 MB | 39.6 s |

L'RSS e il tempo trascorso vengono da `/usr/bin/time -v` (a). `run_capped.sh` arriva al
programma con una catena di `exec`, quindi l'RSS misurato è quella di `riptide_stats_calibrate`.
Il pilota rispetta il criterio (< 500 MB, < 30 s) e il run completo è partito. La stima era
~45 s e < 100 MB **(b)**: il tempo torna, la memoria è cinque volte sotto la stima.

Le soglie sono state valutate con uno script fuori dal repo che legge `calibration.json` e
applica le formule di `5_predictions.md` così come sono scritte. Nelle tabelle "scarto" è
3·|valore − centro|/semiampiezza: 3.00 è il bordo della soglia.

## 1. Bias — checklist con la verità dal prior

20 000 eventi, verità ~ U(10, 20), prior piatto su [10, 20] e prior sin 2θ_p.

| Quantità | Valore | Soglia | Scarto | Esito |
|---|---|---|---|---|
| `bias_mean` (MeV) | +0.00006 ± 0.0134 | 0 ± 0.040 | 0.00 | pass **(b)** |
| `pull_mean` | −0.0019 | ±0.023 | 0.25 | pass **(b)** |
| `pull_width` | 1.0011 | [0.983, 1.017] | 0.20 | pass **(b)** |
| `coverage` 0.68 | 0.6798 | 0.68 ± 0.0099 | 0.06 | pass **(b)** |
| `coverage` 0.90 | 0.9002 | 0.90 ± 0.0064 | 0.09 | pass **(b)** |
| `coverage` 0.95 | 0.9526 | 0.95 ± 0.0046 | 1.70 | pass **(b)** |

Toy: −0.000 / 1.001 e 0.680 / 0.898 / 0.950. La calibrazione esatta prevista dalla costruzione
è confermata **(c)**.

**Bias per bin (c).** Il prior attira le stime verso il centro, come nel toy.

| Bin (centro, MeV) | Eventi | Bias (MeV) | Toy | Soglia | Scarto | Esito |
|---|---|---|---|---|---|---|
| 11.0 | 3938 | +1.665 | +1.63 | [+1.497, +1.763] | 0.80 | pass |
| 13.0 | 3974 | +0.856 | +0.83 | [+0.698, +0.962] | 0.60 | pass |
| 15.0 | 4065 | +0.068 | +0.11 | [−0.021, +0.241] | 0.96 | pass |
| 17.0 | 3967 | −0.758 | −0.70 | [−0.832, −0.568] | 1.32 | pass |
| 19.0 | 4056 | −1.782 | −1.83 | [−1.961, −1.699] | 1.09 | pass |

## 2. Risoluzione e pull — checklist con la verità ristretta

20 000 eventi con la verità ~ U(11, 19). Sono un campione nuovo, con una sezione di id propria
(1). `5_predictions.md` dice "gli stessi 20 000": la frase è imprecisa, ma le soglie non
cambiano, perché trattano già il campione come indipendente da quello del toy **(a)**. `rms_error` = 1.685 MeV, `mean_std` =
1.769 MeV (contro 1.891 e 1.728 con la verità dal prior): lontano dai bordi il prior restringe
l'errore più di quanto dica la std del posterior, da cui la sovracopertura.

| Quantità | Valore | Toy | Soglia | Scarto | Esito |
|---|---|---|---|---|---|
| `pull_mean` | +0.0266 | +0.031 | [0.003, 0.059] | 0.48 | pass **(c)** |
| `pull_width` | 0.9259 | 0.929 | [0.909, 0.949] | 0.47 | pass **(c)** |
| `coverage` 0.68 | 0.7545 | 0.755 | [0.742, 0.768] | 0.13 | pass **(c)** |
| `coverage` 0.90 | 0.9443 | 0.943 | [0.936, 0.950] | 0.54 | pass **(c)** |
| `coverage` 0.95 | 0.9775 | 0.975 | [0.970, 0.980] | 1.50 | pass **(c)** |

## 3. Coverage e contrazione con N

Posterior combinati, verità ~ U(11, 19), una per esperimento.

| N / M | `rms_error` (MeV) | `std_sqrt_N` (MeV) | `pull_mean` | `pull_width` |
|---|---|---|---|---|
| 10 / 400 | 0.496 (toy 0.45) | 1.361 ± 0.031 (toy 1.38) | +0.019 (toy +0.15) | 1.040 (toy 0.96) |
| 100 / 200 | 0.102 (toy 0.094) | 0.945 ± 0.015 (toy 0.98) | +0.005 (toy +0.07) | 1.043 (toy 0.94) |
| 1000 / 40 | 0.0274 (toy 0.029) | 0.925 ± 0.017 (toy 0.93) | +0.096 (toy −0.00) | 0.969 (toy 0.97) |

Confronto con il toy **(c)**, scarto in unità della soglia:

| N | `std_sqrt_N` | `rms_error` | `pull_mean` | `pull_width` | Esito |
|---|---|---|---|---|---|
| 10 | 0.41 | 1.74 | 1.75 | 1.48 | pass |
| 100 | 1.51 | 1.07 | 0.62 | 1.36 | pass |
| 1000 | 0.19 | 0.36 | 0.44 | 0.01 | pass |

**Coverage (d).** Ogni livello cade entro 3 errori binomiali del nominale con M esperimenti:

| N / M | 0.68 | 0.90 | 0.95 | Esito |
|---|---|---|---|---|
| 100 / 200 | 0.665 | 0.890 | 0.950 | pass |
| 1000 / 40 | 0.775 | 0.925 | 0.950 | pass |

Con M = 40 l'errore binomiale a 0.68 vale 0.074: la soglia è larga e il test è debole.

**Contrazione 1/√N (d).** `std_sqrt_N` scende da 1.361 (N = 10) a 0.945 (N = 100), poi cambia
del 2.0% fino a 0.925 (N = 1000), sotto il 10% previsto: pass. A N = 10 il prior su [10, 20]
taglia ancora la verosimiglianza e la std è sopra il regime asintotico; da N = 100 si è già nel
regime 1/√N.

## 4. Robustezza al prior su θ_p

20 esperimenti di N = 1000 eventi a E_n = 15 MeV, gli stessi eventi con i due prior.

| Quantità | Valore | Toy | Soglia | Scarto | Esito |
|---|---|---|---|---|---|
| `flat.bias_mean` (MeV) | −0.0973 ± 0.0065 | −0.101 ± 0.009 | [−0.134, −0.068] | 0.34 | pass **(c)** |
| `flat.pull_mean` | −3.75 | −3.9 | [−4.98, −2.82] | 0.41 | pass **(c)** |
| `physical.bias_mean` (MeV) | −0.0040 | −0.006 ± 0.008 | 0 ± 0.019 | 0.62 | pass **(b)** |
| `physical.pull_mean` | −0.134 | — | ±0.64 | 0.63 | pass **(b)** |

Con il prior piatto su θ_p la coverage dei 20 esperimenti è 0.00 / 0.05 / 0.15: il bias di
~0.1 MeV vale ~3.7 std del posterior combinato (0.026 MeV). Il prior fisico sin 2θ_p è quindi
necessario, come nel toy **(c)**.

## 5. Robustezza al prior su E_n

Verità ~ U(10, 20) per esperimento; prior di riferimento piatto su [10, 20], alternativo
piatto su [10, 30] con lo stesso passo (999 punti).

| N / M | `shift_mean` | Toy | Soglia | Scarto | `shift_max_abs` | `std_ratio_mean` |
|---|---|---|---|---|---|---|
| 1 / 400 | +0.930 ± 0.071 | +0.79 | [+0.484, +1.096] | 1.37 | 12.8 | 1.60 |
| 3 / 400 | +0.351 | +0.36 | [+0.111, +0.609] | 0.10 | — | — |
| 10 / 400 | +0.105 | +0.09 | [−0.015, +0.195] | 0.43 | — | — |
| 30 / 200 | +0.030 | +0.05 | [−0.007, +0.107] | 1.06 | — | — |
| 100 / 200 | +0.013 | +0.03 | [−0.010, +0.070] | 1.24 | — | — |

Tutti i confronti con il toy passano **(c)**. `shift_mean` decresce con N in ogni passo
**(d)**: pass. A N = 1 lo spostamento medio vale quasi una std, con punte di 12.8 std negli
eventi con verosimiglianza larga: con un solo evento il prior conta. Da N = 30 lo spostamento è
sotto 0.05 std.

**Coverage e pull.** Riferimento **(b)**, entro 3 errori binomiali del nominale con M
esperimenti, e `pull_mean` entro 0 ± 3·`pull_width`/√M; alternativa a 0.68 **(d)**:

| N / M | Rif. 0.68 | Rif. 0.90 | Rif. 0.95 | Rif. `pull_mean` | Alt. 0.68 | Esito |
|---|---|---|---|---|---|---|
| 1 / 400 | 0.670 | 0.913 | 0.963 | −0.034 | 0.670 | pass |
| 3 / 400 | 0.698 | 0.920 | 0.968 | −0.009 | 0.695 | pass |
| 10 / 400 | 0.720 | 0.915 | 0.965 | +0.055 | 0.710 | pass |
| 30 / 200 | 0.715 | 0.930 | 0.970 | +0.079 | 0.715 | pass |
| 100 / 200 | 0.695 | 0.890 | 0.935 | +0.009 | 0.695 | pass |

La coverage dell'alternativa a 0.68 resta fra 0.67 e 0.72, compatibile con il toy (0.66–0.70):
un prior più largo del vero sposta le stime a N piccolo senza rompere la coverage.

## 6. Tempi

**Benchmark della verosimiglianza (c)** (`riptide_stats_bench`, `release`, thread singolo,
1000 eventi, migliore di 3, griglia 500 × 500, seme 20260907):

| Versione | 1 evento | 1000 eventi | Scarto dalla precedente |
|---|---|---|---|
| v0, semplice | 3.79 ms | 1.183 s | — |
| vA, 1/σ precalcolati, logsumexp senza controlli, buffer riusati | 3.22 ms | 1.140 s | max 1.8·10⁻¹² assoluto, 4.4·10⁻¹⁶ relativo |
| vB, finestra su θ_p con C = 50 (D35) | 0.47 ms | 0.174 s | identica a vA bit per bit |
| vC, salto delle celle sotto max − C | — | 0.176 s | identica a vB, scartata |

**Target D35 mancato (a).** 0.174 s per 1000 eventi contro ≤ 0.10 s, ma 4.1 volte meglio dei
0.71 s del toy. Dopo la finestra il tempo è dominato da `exp`: ~3.7 ns per chiamata scalare di
glibc, sui termini rimasti nella finestra. Scendere ancora richiederebbe un `exp` vettoriale o
approssimato, che cambia i numeri, oppure il parallelismo sugli eventi, che D35 ammette solo con
una decisione separata. Per la fase 6 non serve: il run `c1` completo, 175 600 posterior, dura
40 s.

**Tempi per sezione di `c1`** (stampati dal run, senza soglia):

| Sezione | Tempo |
|---|---|
| checklist `truth_from_prior` | 3.66 s |
| checklist `truth_narrow` | 3.53 s |
| contrazione N = 10 / 100 / 1000 | 0.69 / 3.41 / 6.83 s |
| prior su θ_p | 6.84 s |
| prior su E_n, N = 1 / 3 / 10 / 30 / 100 | 0.19 / 0.56 / 1.85 / 2.77 / 9.20 s |

La checklist ristretta estrae un campione nuovo della stessa dimensione, quindi il suo costo è lo
stesso della prima.

## Riepilogo

| Sezione | Soglie | Passate |
|---|---|---|
| §1 checklist dal prior | 11 | 11 |
| §2 checklist ristretta | 5 | 5 |
| §3 contrazione | 20 | 20 |
| §4 prior su θ_p | 4 | 4 |
| §5 prior su E_n | 31 | 31 |
| **Totale** | **71** | **71** |

Le "cinquanta" soglie stimate in `5_predictions.md` sono in realtà 71, contando una soglia per
livello e per riga. Con 71 prove indipendenti a 3σ la probabilità di almeno un falso fallimento
era ~17% **(b)**, e non si è presentato.
