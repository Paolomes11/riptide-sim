# Previsioni della fase 5 — scritte prima dei run

Versione 1 (2026-11-22). Va committata prima di qualsiasi run di `riptide_stats_calibrate`. I
risultati (`5_results.md`) citeranno l'hash di quel commit. Dopo i run le soglie non si toccano.

La fase 5 non legge dati Geant4. `riptide_stats_calibrate` estrae pseudo-esperienze dal modello
del toy (D35), le passa al motore del Caso A portato in C++ (D32, D33) e scrive
`stats/<calibration_id>/calibration.json`. Le chiavi del file sono descritte in
`docs/architecture/data_layout.md`. I numeri di confronto vengono dalla calibrazione del toy
(`docs/toy_model_findings.md` §4) e dai suoi script, che seguono lo stesso protocollo: livelli,
domini, N, M ed energie sono gli stessi.

Etichette: **(a)** fatto stabilito (sorgenti, dati, costruzione del codice); **(b)** stima con
ipotesi esplicite; **(c)** verificato con Monte Carlo; **(d)** ipotesi da mettere alla prova.

## Run previsti

Si usa la build `release` del commit di questo file. Ogni comando gira sotto
`scripts/run_capped.sh -m 5G`, uno alla volta, dopo `scripts/data_root.sh check`. Il seme è
quello del progetto, `20260907`. Le risoluzioni sono quelle del toy: σ_Ep = 0.10 MeV,
σ_θ = 0.08 rad.

```bash
scripts/run_capped.sh -m 5G -- build/release/analysis/apps/riptide_stats_calibrate \
  --calibration-id <id> [--pilot]
```

| `calibration_id` | Opzioni | Ruolo |
|---|---|---|
| `pilot` | `--pilot` | misura di tempo e RSS prima del run completo |
| `c1` | — | calibrazione completa, quella confrontata con le soglie |

**Costo (b).** In `release` la verosimiglianza di un evento costa ~0.17 ms sulla griglia
500 × 500 e ~0.35 ms su quella larga (999 × 500), secondo il benchmark della finestra D35. Il run
completo valuta in tutto:
- 40 000 eventi nelle due checklist;
- 64 000 eventi nella contrazione;
- 40 000 eventi nel test del prior su θ_p, cioè 20 000 eventi per due modelli;
- 31 600 eventi nel test del prior su E_n, una volta sulla griglia di default e una su quella
  larga.

La stima è ~45 s a thread singolo; il pilota costa circa un decimo. **Memoria (b):** la matrice
più grande ha 1024 × 999 double (8 MB), le tabelle cinematiche dei tre modelli occupano 8 MB e i
campioni meno di 1 MB. L'RSS atteso è < 100 MB.

**Criterio per il run completo.** Il run completo parte se il pilota finisce con RSS < 500 MB e
dura meno di 30 s. Altrimenti si torna al piano.

**Prove multiple.** Le soglie sono una cinquantina, quasi tutte a 3σ. Anche con un port corretto
la probabilità che almeno una fallisca per fluttuazione è ~10%. Una soglia fallita resta fallita
in `5_results.md`. Un fallimento isolato fra 3σ e 4σ è un candidato fluttuazione: come diagnosi,
e non in sostituzione, si può ripetere il run con il seme `20260908` e riportare entrambi gli
esiti. Due fallimenti nella stessa sezione, oppure uno oltre 4σ, indicano un difetto da cercare
nel codice.

**Errori di confronto con il toy.** Quando una quantità si confronta con un numero del toy,
anche quel numero ha il suo errore statistico, con lo stesso N e M ma estrazioni indipendenti.
La tolleranza vale quindi 3·√2·σ, con σ stimato dal nostro run. A questa si aggiunge, dove il
toy riporta due cifre decimali, metà dell'ultima cifra. La precisione float32 del toy
(contro i double del port, D32) e le masse uguali di n e p (scarto relativo ~5·10⁻⁷ sulla
cinematica, D32) restano sotto questi errori **(b)**.

## 1. Checklist con la verità dal prior (`checklist[truth_from_prior]`)

Sono 20 000 eventi, con E_n vera ~ U(10, 20) e un posterior per evento con prior piatto su
[10, 20] e prior sin 2θ_p.

**Calibrazione esatta (b).** Il generatore e il modello coincidono, a parte la discretizzazione
delle griglie e il prior sin 2θ_p, che è la densità di θ_p nel limite non relativistico (D32). Se
la verità viene dal prior, il posterior è quindi calibrato per costruzione: media degli errori
nulla, pull di media 0 e scarto quadratico medio 1, coverage pari al livello nominale. Il toy ha
dato −0.000 / 1.001 e 0.680 / 0.898 / 0.950.

| Quantità | Soglia |
|---|---|
| `pull_mean` | 0 ± 3.3/√20 000 = ±0.023 |
| `pull_width` | 1 ± 3.3/√40 000 = [0.983, 1.017] |
| `coverage` 0.68 | 0.68 ± 0.0099 (3 errori binomiali) |
| `coverage` 0.90 | 0.90 ± 0.0064 |
| `coverage` 0.95 | 0.95 ± 0.0046 |
| `bias_mean` | 0 ± 3·`bias_error` |

**Bias per bin (c).** I bin sono 5, uguali sull'intervallo delle verità, come in numpy. Il toy ha
dato +1.63, +0.83, +0.11, −0.70, −1.83 MeV, perché il prior attira le stime verso il centro. Per
ogni bin la soglia è |bias − toy| ≤ 3·√2·`rms_error`/√`counts` + 0.005 MeV. Lo scarto quadratico
medio di tutto il campione sovrastima l'errore di un bin, quindi la soglia è conservativa.

## 2. Checklist con la verità ristretta (`checklist[truth_narrow]`)

Gli eventi sono gli stessi 20 000, con la verità ~ U(11, 19) e lo stesso posterior. La verità
lontana dai bordi del prior porta a una sovracopertura. Toy: pull +0.031 / 0.929, coverage
0.755 / 0.943 / 0.975 **(c)**. Le soglie sono 3·√2·σ attorno al toy:

| Quantità | σ | Soglia |
|---|---|---|
| `pull_mean` | `pull_width`/√N ≈ 0.0066 | [0.003, 0.059] |
| `pull_width` | `pull_width`/√(2N) ≈ 0.0046 | [0.909, 0.949] |
| `coverage` 0.68 | √(p(1−p)/N) ≈ 0.0030 | [0.742, 0.768] |
| `coverage` 0.90 | ≈ 0.0016 | [0.936, 0.950] |
| `coverage` 0.95 | ≈ 0.0011 | [0.970, 0.980] |

## 3. Contrazione con N (`contraction`)

Posterior combinati con la verità ~ U(11, 19), una per esperimento. Valori del toy **(c)**:

| N / M | `rms_error` (MeV) | `std_sqrt_N` (MeV) | `pull_mean` | `pull_width` |
|---|---|---|---|---|
| 10 / 400 | 0.45 | 1.38 | +0.15 | 0.96 |
| 100 / 200 | 0.094 | 0.98 | +0.07 | 0.94 |
| 1000 / 40 | 0.029 | 0.93 | −0.00 | 0.97 |

Soglie, riga per riga, con w = `pull_width` del nostro run:
- `std_sqrt_N`: |Δ| ≤ 3·√2·`std_sqrt_N_error` + 0.005;
- `rms_error`: |Δ| ≤ 3·`rms_error`/√M + mezza cifra del toy. Questo vale 3·√2·σ con
  σ = rms/√(2M);
- `pull_mean`: |Δ| ≤ 3·√2·w/√M + 0.005;
- `pull_width`: |Δ| ≤ 3·w/√M + 0.005;
- `coverage` a N = 100 e 1000 **(d)**: il toy non ha numeri per questi casi. Ci si aspetta che
  ogni livello cada entro 3 errori binomiali del nominale con M esperimenti, perché a N grande
  il prior su E_n non conta più.

**Contrazione 1/√N (d):** `std_sqrt_N` decresce da N = 10 a N = 100 e cambia meno del 10% da
N = 100 a N = 1000, come nel toy (0.98 → 0.93).

## 4. Prior su θ_p (`theta_prior`)

Sono 20 esperimenti di N = 1000 eventi a E_n = 15 MeV, con gli stessi eventi valutati con il
prior fisico sin 2θ_p e con quello piatto. Il toy ha dato bias −0.101 ± 0.009 MeV e pull −3.9 con
il prior piatto, e bias −0.006 ± 0.008 MeV con quello fisico **(c)**.

| Quantità | Soglia |
|---|---|
| `flat.bias_mean` | \|Δ dal toy\| ≤ 3·√(0.009² + `bias_error`²) |
| `flat.pull_mean` | \|Δ da −3.9\| ≤ 3·√2·`pull_width`/√20 + 0.05 |
| `physical.bias_mean` | 0 ± 3·`bias_error` **(b)**: con il prior corretto, a 15 MeV e N = 1000, il prior su E_n è ininfluente |
| `physical.pull_mean` | 0 ± 3·`pull_width`/√20 **(b)** |

## 5. Prior su E_n (`energy_prior`)

La verità è ~ U(10, 20) per esperimento. Il prior di riferimento è piatto su [10, 20]; quello
alternativo è piatto su [10, 30], con lo stesso passo (`sameStepDomain`). Toy, `shift_mean` alle
righe N/M = 1/400, 3/400, 10/400, 30/200, 100/200: +0.79, +0.36, +0.09, +0.05, +0.03, con
coverage 68% fra 0.66 e 0.70 **(c)**.

| Quantità | Soglia |
|---|---|
| `shift_mean` | \|Δ dal toy\| ≤ 3·√2·`shift_error` + 0.005, per riga |
| `shift_mean` | decrescente con N **(d)** |
| `reference.coverage` (tutti i livelli) | entro 3 errori binomiali del nominale con M esperimenti **(b)**: verità dal prior, calibrazione esatta |
| `reference.pull_mean` | 0 ± 3·`pull_width`/√M **(b)** |
| `alternative.coverage` 0.68 | 0.68 ± 3·√(0.68·0.32/M) **(d)** |

## 6. Tempi

I tempi per sezione stampati dal run entrano in `5_results.md` senza soglia. Il target D35 sulla
verosimiglianza (≤ 0.10 s per 1000 eventi) è già mancato nel benchmark: 0.174 s con la finestra
su θ_p (addendum D35). `5_results.md` riporta quella misura e il motivo dello scarto.
