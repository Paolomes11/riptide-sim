# Risultati della fase 1b — verifiche fisiche e prestazioni

Data: 2026-11-19. Previsioni e soglie: `docs/validation/1b_predictions.md` e
`docs/validation/predictions/*_1e7.json`, committate in `6de0124` **prima** dei run. Correzione
dopo il primo esito: D17, con le previsioni del rilancio committate in `907f0be` prima del
rilancio.

**Esito: tutte le verifiche passano (`c2`).** Il primo esito (`c1`) aveva due fallimenti nel run
mono 14 MeV, entrambi spiegati da un difetto delle attese e non della simulazione (D17).

## Run

Build `release`, 8 thread, seme 20260907, `scripts/run_capped.sh -m 5G`, un run alla volta.
Simulazione al commit `c8c5e3b`; dati sull'SSD in `sim/<run>` e `validation/<run>/<id>`.

| Run | Config | Eventi | Tempo | RSS massimo | Su disco | Byte/evento |
|---|---|---|---|---|---|---|
| `1b_mono_14MeV_1e7` | `config/mono_14MeV.json` | 10⁷ | 253 s | 143 MiB | 543 MB | 54.3 |
| `1b_flat_10_20MeV_1e7` | `config/flat_10_20MeV.json` | 10⁷ | 266 s | 141 MiB | 617 MB | 61.7 |
| `1b_flat_19p5_20MeV_1e7` | `config/flat_19p5_20MeV.json` | 10⁷ | 265 s | 145 MiB | 580 MB | 58.0 |

Tempo da `run.json` (inizializzazione compresa), RSS da `/usr/bin/time -v`, byte dagli
`events_t*.root` (c). Circa 3.8·10⁴ eventi/s, cioè 4.7·10³ per thread: in linea con la stima
di circa 4 minuti del pilota (§5 delle previsioni). La stima per il controllo dello spazio
libero (`kBytesPerEvent`) passa da 1 KiB provvisorio a 128 B, circa il doppio del massimo
misurato.

`riptide_sim_check` (`--bins-u 50 --bins-z 50`): 29–44 s e RSS al più 164 MiB a run (c).

## Ordine dei fatti

1. `c1`, codice al commit `c8c5e3b` (`checks_version` 4): piatto e bordo PASS; mono FAIL su
   offset relativo massimo e χ² delle frazioni.
2. Cause trovate sui dati esistenti, senza nuovi run (D17):
   - offset: un residuo di arrotondamento di 2·10⁻¹⁵ MeV su un evento con E_p ≈ 9 eV;
   - frazioni: il break-up n–d su H-2 è "other" nella simulazione, ma mancava dalle attese.
3. Correzione delle attese (H-2 da G4NDL, `2486cc8`), oracolo dell'offset al netto
   dell'arrotondamento (`00092fd`), D17 e previsioni del rilancio (`907f0be`). Soglie, α e range
   invariati.
4. `c2`, stessi eventi, codice al commit `907f0be` (`checks_version` 5): tutto PASS, come previsto
   in D17.

## 1. Oracoli deterministici — (a)

| Verifica | Soglia | mono | piatto | bordo |
|---|---|---|---|---|
| `must_be_zero.total` | = 0 | 0 | 0 | 0 |
| `Ep_residuals.max_abs` [MeV] | < 10⁻⁹·E_n,max | 1.95·10⁻¹⁴ | 2.84·10⁻¹⁴ | 3.20·10⁻¹⁴ |
| `model_mix.Ep_residuals.max_abs` [MeV] | < 10⁻¹¹·E_n,max | 0 | 4.36·10⁻¹³ | 4.86·10⁻¹³ |
| eventi non HP nel mono | = 0 | 0 | — | — |
| `np_model_ambiguous` | riportato | 0 | 0 | 0 |

`expected_rare.np_first_without_proton` (n–d con deutone di rinculo): 260 / 242 / 170.

**Offset relativo dalla cinematica fisica** (D14), contro il range previsto con margine 10⁻⁹:

| Run | Range previsto | `c1`: relative_min / max | `c2`: rounding_free min / max | Esito c1 → c2 |
|---|---|---|---|---|
| mono | 5.36581425 – 5.48613868 (·10⁻⁴) | 5.36581436 / **5.48614088** | 5.36581436 / 5.48613626 | FAIL → PASS |
| piatto | 5.33221807 – 5.50321193 | 5.33281530 / 5.50303407 | 5.33281530 / 5.50303403 | PASS → PASS |
| bordo | 5.33221807 – 5.50321193 | 5.33229844 / 5.50313958 | 5.33229844 / 5.50313941 | PASS → PASS |

## 2. Test statistici contro `jeff33` — (c)

Bonferroni con α_family = 0.01 a run: α per test 1.11·10⁻³ (mono, 3 χ² + 6 pull, soglia
|pull| < 3.261) e 9.09·10⁻⁴ (piatto e bordo, 4 χ² + 7 pull, soglia 3.317).

| Verifica | mono `c1` | mono `c2` | piatto `c2` | bordo `c2` |
|---|---|---|---|---|
| spettro u, χ²/ndf (p) | 54.98/49 (0.258) | 54.98/49 (0.258) | 39.60/49 (0.829) | 58.35/49 (0.169) |
| profondità z, χ²/ndf (p) | 59.75/49 (0.140) | 59.83/49 (0.138) | 50.47/49 (0.415) | 64.78/49 (0.065) |
| frazioni, χ²/ndf (p) | **28.34/4 (1.1·10⁻⁵)** | 6.33/4 (0.176) | 1.82/4 (0.769) | 0.19/4 (0.996) |
| pull a₁…a₆ | +0.43 −0.25 −1.83 +1.41 −2.79 +0.06 | uguali | −0.13 +0.13 +0.24 +0.15 +0.08 −0.17 | −1.37 −0.28 −0.56 +1.60 −0.38 +1.14 |
| `model_mix` pull | — | — | +0.71 (40 837 contro 40 753.5) | −0.53 (819 885 contro 820 163.7) |
| `model_mix` χ²/ndf (p) | — | — | 7.41/10 (0.686) | 5.51/10 (0.855) |

In `c1` piatto e bordo passavano già: frazioni 9.95/4 (p 0.041) e 6.67/4 (p 0.154). Spettro u,
Legendre e `model_mix` non dipendono da H-2 e sono identici fra `c1` e `c2`.

**Pull delle frazioni** (unknown, n–p, n–C el., n–C non el., other):

| Run | `c1` | `c2` | Previsto in D17 per `c2` (d) |
|---|---|---|---|
| mono | +0.89 +0.64 +0.19 −2.28 **+4.80** | +1.01 +0.49 +0.22 −2.25 +1.07 | +1.0 +0.5 +0.2 −2.3 +1.1 |
| piatto | −0.12 −1.07 +0.85 +0.29 +2.91 | −0.01 −1.21 +0.89 +0.32 +0.12 | 0.0 −1.2 +0.9 +0.3 +0.1 |
| bordo | −0.04 +0.41 −0.34 −0.01 +2.54 | +0.05 +0.29 −0.30 +0.01 +0.24 | 0.0 +0.3 −0.3 0.0 +0.2 |

L'eccesso di "other" in `c1` aveva lo stesso segno e grandezza in tutti e tre i run (~+64
eventi), come atteso per un difetto sistematico delle attese e non per una fluttuazione.

## 3. Confronto con `endfb80` — attese, non soglie (c)

`c2`: frazioni χ² 305/4, 304/4, 293/4, con n–C non elastico +14.9σ…+16.2σ ed n–C elastico
−10.3σ…−12.0σ, come previsto (differenza fra le valutazioni del carbonio). Spettro u: χ² 64.1,
43.4, 69.5 su 49 gradi; pull di Legendre al più |3.19| (a₃, mono), contro 1.35σ previsti come
differenza media: il resto è la fluttuazione comune con `jeff33` (a₃ = −1.83 contro jeff33).
Si conferma il limite di D13: a 10⁷ eventi la simulazione non distingue le due revisioni di
H-1, e la fase 1b verifica l'implementazione, non i dati n–p.

## Regressione delle verifiche v7 (`c3`, 2026-11-20)

`riptide_sim_check` al commit `09bef22` (`checks_version` 7, build `release`) sugli stessi tre
run, con `--check-id c3 --bins-u 50 --bins-z 50` e `scripts/run_capped.sh -m 5G`. Previsione
scritta prima: stessi numeri di `c2`, `deposits` saltata perché i run hanno schema 1.

**Esito: nessuna regressione.** Per tutti e tre i run (a):
- i sei CSV (`spectrum_u`, `depth_z`, `fractions`, `legendre`, `model_mix`, `residuals_Ep`)
  sono identici byte per byte a quelli di `c2`;
- in `checks.json` gli unici campi diversi sono `check_id`, `checks_version` (5 → 7),
  `versions.git`, `versions.riptide` (0.2.0 → 0.3.0) e `versions.schema` (1 → 2). Inoltre
  compare `results.deposits` = `{"skipped": "schema version 1 has no deposit columns"}`;
- nessun CSV di fase 2 è stato scritto.

`versions.schema` non era nella previsione. È il `kSchemaVersion` compilato nel programma di
verifica, non lo schema del run letto, e vale 2 da D19. Lo schema del run è quello di
`run.json`, e qui è 1, come mostra `deposits.skipped`. La descrizione in `data_layout.md` è
stata corretta.

Tempo 33–47 s, RSS massimo 164 MiB (c), come in `c2`.

## Punti aperti

- **Confronto n–p indipendente** (D13): JEFF-3.3 ed ENDF/B-VIII.0 per H-1 sono la stessa
  analisi; resta da fare con una valutazione o misura indipendente.
- **H-2** è verificato solo per coerenza con G4NDL (D17), non contro una valutazione
  indipendente.
- **C-13** (D15): residuo non modellato di circa 3·10⁻⁴ sul non elastico del carbonio a 14 MeV,
  sotto la sensibilità di questi run (≈ 0.4σ).
- **Bordo a 20 MeV** (D16): sopra 19.5 MeV lo stato finale elastico è in parte CHIPS; per la
  fase 2 il bordo della sorgente va trattato come sistematica o evitato.
- **Massa HP di H-1** (D14): E_p più alto di circa 5.4·10⁻⁴ in relativo; correzione, modello o
  sistematica da decidere prima dell'analisi della fase 2.
