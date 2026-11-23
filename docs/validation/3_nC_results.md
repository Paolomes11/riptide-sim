# Risultati — pull del n–C non elastico lungo +X (run gemello)

Data: 2026-11-23. Previsioni e regole di lettura: `docs/validation/3_nC_predictions.md`
versione 1, committata in `1c8d070` **prima** del run. Fra le previsioni e il run non è entrato
nessun commit.

Etichette: **(a)** fatto stabilito; **(b)** stima con ipotesi esplicite; **(c)** verificato con
Monte Carlo; **(d)** ipotesi da mettere alla prova.

**Esito:** pull del n–C non elastico contro `jeff33` nel run gemello **+1.79**, sotto la soglia
di 3.0. Per la regola scritta prima del run è una fluttuazione e il punto aperto della fase 3
si chiude, senza D36 **(c)**. Le verifiche `c1` passano tutte le soglie della fase 3.

## Run

| Run | Eventi | Tempo | RSS massimo | Su disco | Byte/evento |
|---|---|---|---|---|---|
| `3_ff_mono_14MeV_x_1e7_s20260908` | 10⁷ | 159.3 s | 178 MB | 1.10 GB | 110.1 |

- Build `release` a `1c8d070` (`git describe` v0.7.0-6-g1c8d070), Geant4 11.4.1, `QGSP_BIC_HP`,
  G4NDL 4.7.1, 8 thread, seme 20260908, `scripts/run_capped.sh -m 5G`, dopo
  `scripts/data_root.sh check`. Comandi come nella §2 delle previsioni.
- Tempo da `run.json`, RSS da `/usr/bin/time -v`, byte dagli `events_t*.root` (c).
- `riptide_sim_check --check-id c1 --bins-u 50 --bins-z 50`: 48 s, RSS 128 MB (c).
- Contro la stima (b): circa 2.8 minuti, meno di 0.2 GB, 50 s, 1.1 GB. Tutto confermato.

## 1. Verifiche `c1` — soglie della fase 3

**Oracoli deterministici (a).** Tutti PASS, con le soglie di `3_results.md` §1:
- `must_be_zero.total` = 0 e `deposits.must_be_zero.total` = 0;
- `deposits.skipped` e `pstar.skipped` sono `null`;
- corda 60 mm entro 10⁻¹⁴, sd 0;
- bilancio dei protoni contenuti fra −5.3·10⁻¹⁵ e 9.3·10⁻⁴ MeV;
- `pstar.samples` + `below_table` = 1 069 645 + 79 = 1 069 724, uguale ai protoni del bilancio;
- `pstar.sha256` è quello atteso;
- `Ep_residuals.max_abs` = 1.8·10⁻¹⁴ MeV;
- `observed_non_hp` = 0;
- offset HP relativo da 5.36581 a 5.48614·10⁻⁴, come nel run originale.

**Test statistici contro `jeff33` (c).** Soglie lungo +X: p > 1.000·10⁻³ e |pull| < 3.291.

| Test | Run originale (seme 20260907) | Gemello (seme 20260908) |
|---|---|---|
| spettro u, χ²/49 (p) | 49.57 (0.450) | 58.78 (0.160) |
| profondità s = z/L_i, χ²/49 (p) | 43.99 (0.676) | 40.42 (0.804) |
| frazioni, χ²/4 (p) | 16.22 (0.00274) | 4.21 (0.379) |
| pull di Legendre a₁…a₆ | 0.96, −1.01, −1.95, −0.36, 0.70, 2.00 | 1.19, −1.79, 1.38, 0.96, 2.08, −2.16 |
| pull dell'accettanza | 1.435 | −1.595 |

**Tutti PASS.** a₁ = −1.413·10⁻² ± 5.7·10⁻⁴, a 25σ come nel run originale.

**PSTAR (c).** Tutti PASS:
- luce, residuo relativo medio −0.0619%, nella banda [−0.070%, −0.055%];
- luce per bin, massimo |run − riferimento| = 0.00010% su 42 bin, sotto 0.005%;
- range, residuo relativo medio −0.1093%, nella banda [−0.125%, −0.035%].

**Previsioni (b), non soglie:** confermate.
- 6 363 772 raggi sul blocco contro circa 6.37·10⁶;
- 2 828 699 prime interazioni n–p, n–C el. e n–C non el., contro circa 2.83·10⁶.
- n–C non elastiche: 664 820 contro 663 409 ± 787. Le previsioni davano 663 889, l'attesa del
  run originale, ma l'attesa scala con i raggi sul blocco, che qui sono 4 608 in meno.

## 2. Test principale — (c)

Pull per categoria contro `jeff33`, nell'ordine [sconosciuta, n–p, n–C el., n–C non el.,
altro]:

| Run | Pull |
|---|---|
| originale | −1.59, −0.19, −0.48, **+3.73**, −1.51 |
| gemello | −1.21, +0.69, −0.36, **+1.79**, −0.50 |

|pull| < 3.0, quindi per la regola delle previsioni è una **fluttuazione**.

Pull combinato richiesto dalla regola: (3.73 + 1.79)/√2 = **3.90**. Va letto con la sua
distorsione. Il run originale è entrato nello studio proprio perché era il pull più grande fra
circa dieci, quindi il suo +3.73 è gonfiato dalla selezione. L'unica prova non distorta è il
gemello: +1.79, cioè p = 0.073 a due code (0.036 a una coda, nel verso dell'eccesso) **(c)**.

## 3. Frazione condizionata e diagnostico — attese, non soglie (c)

Frazione di n–C non elastico fra le prime interazioni n–p, n–C el. e n–C non el. Attesa
`jeff33` 0.234671, uguale per tutti i run a 14 MeV (§1 delle previsioni).

| Campione a 14 MeV (indipendenti) | Interazioni | Frazione | Pull |
|---|---|---|---|
| `3_ff_mono_14MeV_x_1e7` (+X) | 2 831 300 | 0.23552 | +3.37 |
| `3_ff_mono_14MeV_x_1e7_s20260908` (+X, gemello) | 2 828 699 | 0.23503 | +1.41 |
| `3_ff_mono_14MeV_isotropic_1e7` | 608 232 | 0.23411 | −1.03 |
| `2_mono_14MeV_1e7` (fascio a matita, cubo da 100 mm) | 6 241 491 | 0.23433 | −2.03 |
| tutti e quattro insieme | 12 509 722 | 0.23475 | +0.63 |

- Il gemello ha di nuovo un eccesso, ma più piccolo: +0.15% in relativo, contro +0.36%.
- I quattro campioni insieme sono in accordo con `jeff33`.
- Fra loro la dispersione è grande: χ² = 18.1 su 3 gradi di libertà attorno alla media comune. È
  dominata dal run originale, che è il campione selezionato. Senza il run originale, i tre
  campioni restanti danno χ² = 7.2 su 3 gradi di libertà contro `jeff33` (p ≈ 0.07).

**Diagnostico per profondità e per facce sul gemello.** È la stessa macro, fuori dal repo, del
§1 delle previsioni:
- **profondità:** χ² = 2.15 su 5 gradi di libertà, nessun andamento;
- **facce:** χ² = 3.02 su 4 gradi di libertà, pull +0.76 / −0.74 / −0.22 / +0.67 / −1.18 dalla
  faccia verso il centro.

Il pull +2.02 del run originale nel bin più vicino alle facce non si ripete (+0.76). Nessuna
struttura geometrica.

## 4. Conclusione e conseguenze

- **Chiuso (c):** il +3.73σ lungo +X non si riproduce con un seme indipendente, e la regola
  scritta prima del run lo classifica come fluttuazione. Non si apre D36.
- **Resta, sotto la sensibilità di questi run (d):** un eccesso di n–C non elastico lungo +X di
  0.1–0.2% in relativo. Il gemello da solo lo vede a 1.4–1.8σ, e non è escluso. Non tocca la
  fase 6: il taglio `single_np` toglie gli eventi n–C, e il Caso A non usa le frazioni.
- **Da rivedere nella fase 9,** quando gli eventi n–C entreranno nel modello. Se le frazioni
  servono a quel livello di precisione, occorre un run dedicato con più semi e lo stesso
  confronto fra configurazioni.
