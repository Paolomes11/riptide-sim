# Risultati della fase 2 — deposito di energia e quenching di Birks

Data: 2026-11-21. Previsioni e soglie: `docs/validation/2_predictions.md` versione 2, committata
in `85283ba` **prima** dei run di produzione. La versione 1 (`83b959d`) è stata corretta dopo il
pilota, che non è un run di validazione (D22).

**Esito (`c1`):**
- **tre run da 10⁷:** tutte le soglie passano. Oracoli, statistiche contro `jeff33`, bande e
  bin di PSTAR.
- **run di sensibilità a kB:** passano tutte le soglie tranne una, «identici», per due campi in
  virgola mobile: `contained_balance_MeV.mean` e `pstar.range_mm`. Ognuno fallisce in due
  confronti, kB = 0.145 e 0.165 contro 0.155, quindi 4 esiti negativi in tutto.
  - Differiscono al più di 1.3·10⁻¹² in relativo.
  - Tutte le grandezze intere sono identiche, e così `proton_exits.csv` e la colonna `Edep`,
    byte per byte.
  - La causa è l'ordine di somma nella verifica, non la simulazione (§1, run di kB). Previsione
    fallita (d), soglia invariata.

**Esito (`c2`, `checks_version` 9, D23):** rilette in ordine di `event_id`, passano tutte le
soglie dei sei run, comprese le quattro identità dei run di kB. La previsione di `c2` sulla
distanza da `c1` è fallita in parte: lo scarto dei float arriva a 1.8·10⁻⁸ relativo invece di
10⁻¹¹ (§7).

## Run

Build `release` al commit `782d9c2` (`git describe` v0.3.0-21-g782d9c2), Geant4 11.4.1,
`QGSP_BIC_HP`, G4NDL 4.7.1, 8 thread, seme 20260907, `scripts/run_capped.sh -m 5G`, un run alla
volta. Dati sull'SSD in `sim/<run>` e `validation/<run>/c1`. Le config di kB = 0.145 / 0.155 /
0.165 sono generate fuori dal repo da `config/mono_14MeV.json`, cambiando solo
`scintillator.birks_kB_mm_per_MeV`. La copia completa è in `run.json`.

| Run | Eventi | Tempo | RSS massimo | Su disco | Byte/evento |
|---|---|---|---|---|---|
| `2_mono_14MeV_1e7` | 10⁷ | 264 s | 169 MiB | 1138 MB | 113.8 |
| `2_flat_10_20MeV_1e7` | 10⁷ | 271 s | 169 MiB | 1209 MB | 120.9 |
| `2_flat_19p5_20MeV_1e7` | 10⁷ | 269 s | 167 MiB | 1129 MB | 112.9 |
| `2_mono_14MeV_kB0p145_1e6` | 10⁶ | 27.9 s | 158 MiB | 114.3 MB | 114.3 |
| `2_mono_14MeV_kB0p155_1e6` | 10⁶ | 29.1 s | 159 MiB | 114.3 MB | 114.3 |
| `2_mono_14MeV_kB0p165_1e6` | 10⁶ | 29.0 s | 162 MiB | 114.3 MB | 114.3 |

Tempo da `run.json` (inizializzazione compresa), RSS da `/usr/bin/time -v`, byte dagli
`events_t*.root` (c).

**Prestazioni contro le stime** (§7):
- **velocità:** circa 3.7·10⁴ eventi/s, cioè 4.6·10³ per thread, poco sotto i 5.5·10³ stimati.
  Il tempo (4.4–4.5 minuti) cade comunque nella banda stimata di 4–8 minuti.
- **memoria:** RSS sotto 0.3 GB, come stimato.
- **byte per evento:** 112.9–120.9 contro circa 118 dello scan. `kBytesPerEvent` = 256 B
  resta, con un margine di più di 2 volte.
- **disco:** 3.82 GB in tutto, contro i circa 4 GB stimati.

`riptide_sim_check` (`--bins-u 50 --bins-z 50`):
- run da 10⁷: 52–66 s, RSS al più 172 MiB;
- run di kB: 6.7 s e 63 MiB (c).

## Ordine dei fatti

1. **Previsioni v1** (`83b959d`), scritte dal codice di `7de4fec` prima di ogni run.
2. **Pilota v2** (schema 2, `checks_version` 7). Nel run a 8 thread, un protone contenuto ha
   `truth_p_Edep` − `truth_Ep` = −2.027 MeV (evento 246185). La causa è uno scattering
   Coulombiano singolo con un rinculo di C12 da 2.027 MeV, trovata ripetendo il run con la
   traccia degli step. L'oracolo v1 del bilancio (min ≥ −10⁻⁹ MeV) sarebbe fallito.
3. **D22** (`60c37bc`): colonna `truth_p_n_nuclear_recoil` (schema 3) e `checks_version` 8.
   I protoni con rinculi nucleari escono dal bilancio e da PSTAR, e vengono contati.
4. **Previsioni v2** (`85283ba`): cambia solo l'esclusione dei rinculi più la banda del loro
   conteggio. Nessun'altra soglia cambia.
5. **Pilota v3**: stessi eventi del pilota v2 evento per evento, `protons_nuclear_recoil` 0 e 1,
   tutto pulito. Il pilota non entra nei risultati.
6. **Produzione `c1`**: i sei run sopra, codice al commit `782d9c2`.

Dopo il pilota v3 sono stati cancellati, col via:
- dall'SSD: il run di traccia `scratch_trace_246185`;
- i piloti v2 (schema 2) con le loro verifiche.

## 1. Oracoli deterministici — (a)

| Verifica | Soglia | mono | piatto | bordo | kB 0.145 / 0.155 / 0.165 |
|---|---|---|---|---|---|
| `must_be_zero.total` (prima interazione) | = 0 | 0 | 0 | 0 | 0 / 0 / 0 |
| `deposits.must_be_zero.total` | = 0 | 0 | 0 | 0 | 0 / 0 / 0 |
| `deposits.skipped`, `pstar.skipped` | `null` | sì | sì | sì | sì |
| Σ classi = eventi | uguale | 10⁷ | 10⁷ | 10⁷ | 10⁶ |
| bilancio `max` [MeV] | ≤ 0.01 | 9.21·10⁻⁴ | 9.17·10⁻⁴ | 8.71·10⁻⁴ | 8.24·10⁻⁴ |
| bilancio `min` [MeV] | ≥ −10⁻⁹ | −5.3·10⁻¹⁵ | −7.1·10⁻¹⁵ | −1.1·10⁻¹⁴ | −5.3·10⁻¹⁵ |
| `protons_nuclear_recoil` | 1–150 (10⁷) | 40 | 51 | 57 | 3 (senza soglia) |
| `pstar.samples` + `below_table` = campioni del bilancio | uguale | 2 749 918 + 213 | 2 647 313 + 174 | 2 025 377 + 105 | 275 426 + 19 |
| `pstar.sha256` | `eec195fe…081b` | sì | sì | sì | sì |
| `setup.birks_kB_mm_per_MeV` | kB della config | 0.155 | 0.155 | 0.155 | 0.145 / 0.155 / 0.165 |
| `Ep_residuals.max_abs` [MeV] | < 10⁻⁹·E_n,max | 1.95·10⁻¹⁴ | 2.84·10⁻¹⁴ | 3.20·10⁻¹⁴ | 1.60·10⁻¹⁴ |
| `model_mix.Ep_residuals.max_abs` [MeV] | < 10⁻¹¹·E_n,max | 0 | 4.36·10⁻¹³ | 4.86·10⁻¹³ | 0 |
| eventi non HP a 14 MeV | = 0 | 0 | — | — | 0 |
| `awr` di H-1 | 0.999167 | sì | sì | sì | sì |

I tre numeri della prima interazione coincidono con quelli della fase 1b:
- `Ep_residuals`, `model_mix` e l'offset sono gli stessi a tutte le cifre stampate;
- `np_first_without_proton` vale 260 / 242 / 170, come nella 1b;
- `np_model_ambiguous` è 0.

Il rinculo nucleare ha un tasso di circa 1.5·10⁻⁵ per protone contenuto a 14 MeV, e di circa
2.8·10⁻⁵ al bordo, dentro la banda (b) di D22. Il tasso cresce con E_p come previsto, perché il
cammino del protone è più lungo.

**Offset relativo dalla cinematica fisica** (D14), al netto dell'arrotondamento, contro il range
della §1 di `1b_predictions.md`, con margine 10⁻⁹:

| Run | Range previsto (·10⁻⁴) | `rounding_free` min / max (·10⁻⁴) | Esito |
|---|---|---|---|
| mono | 5.36581425 – 5.48613868 | 5.36581436 / 5.48613626 | PASS |
| piatto | 5.33221807 – 5.50321193 | 5.33281530 / 5.50303403 | PASS |
| bordo | 5.33221807 – 5.50321193 | 5.33229844 / 5.50313941 | PASS |
| kB (tutti e tre) | come mono | 5.36581641 / 5.48613653 | PASS |

### Run di sensibilità a kB: stessi eventi

| Verifica | Soglia | 0.145 contro 0.155 | 0.165 contro 0.155 | Esito |
|---|---|---|---|---|
| colonna `Edep` di `deposit_spectra.csv` | identica byte per byte | identica | identica | PASS |
| `proton_exits.csv` | identico byte per byte | identico | identico | PASS |
| `event_classes`, `protons_nuclear_recoil` (3), `Edep_overflow` (28 256) | identici | identici | identici | PASS |
| `contained_balance_MeV` | identico | conteggi, `min`, `max` identici; `mean` diverso di 1 ULP | come a sinistra | **FAIL** |
| `pstar.range_mm` | identico | `samples` identico; medie diverse fino a 6·10⁻¹⁴ in relativo | come a sinistra | **FAIL** |
| Δ `pstar.light_MeV.relative_residual` | \|Δ\| < 0.01% | +0.0024% | −0.0022% | PASS |

**Causa dei due fallimenti** (§8 di `2_predictions.md`), cercata sugli stessi dati e senza nuovi
run. Le prove:
- **Ripartizione fra i thread:** gli eventi sono gli stessi, ma Geant4 MT li distribuisce fra i
  thread in modo dinamico. Le dimensioni dei file sono un esempio: con kB = 0.155,
  `events_t3.root` pesa 9.68 MB, mentre con 0.145 pesa 13.75 MB e con 0.165 pesa 15.94 MB.
  `riptide_sim_check` legge i file per thread in ordine, quindi somma gli stessi double in un
  ordine diverso.
- **Scala dei residui:** in tutto `checks.json`, tolti i campi che dipendono da kB per
  costruzione (luce, `Evis`, `kB_zero_max_abs_difference_MeV`), i confronti danno:
  - per ogni grandezza intera, nessuna differenza;
  - per le grandezze in virgola mobile, al più 1.3·10⁻¹² in relativo (pull di Legendre, a₁);
  - nelle colonne di range di `pstar_comparison.csv`, al più 10⁻¹⁴, con i conteggi dei protoni
    per bin identici.

  È la scala dell'arrotondamento di somme su 10⁵–10⁶ termini.
- **Il caso di controllo:** i conteggi di prima interazione (§2) mostrano lo stesso effetto
  rispetto alla fase 1b.

La previsione (d) «identici» valeva per i depositi di ogni evento, verificati bit per bit dal test
d'integrazione. Non vale per le medie, che dipendono dall'ordine di somma. Il fallimento resta
scritto così e la soglia non cambia. Il fatto fisico previsto, cioè stessi eventi, stessi
depositi e stessa traccia, è verificato da ogni campo intero e dai due file confrontati byte per
byte (a). Con la lettura in ordine di `event_id` (`c2`, §7) le due verifiche passano.

## 2. Test statistici contro `jeff33` — (c)

Soglie della §2 di `1b_predictions.md`:
- mono: α 1.11·10⁻³, |pull| < 3.261;
- piatto e bordo: α 9.09·10⁻⁴, |pull| < 3.317.

| Verifica | mono | piatto | bordo |
|---|---|---|---|
| spettro u, χ²/ndf (p) | 54.98/49 (0.258) | 39.60/49 (0.829) | 58.35/49 (0.169) |
| profondità z, χ²/ndf (p) | 59.83/49 (0.138) | 50.47/49 (0.415) | 64.78/49 (0.065) |
| frazioni, χ²/ndf (p) | 6.33/4 (0.176) | 1.82/4 (0.769) | 0.19/4 (0.996) |
| pull a₁…a₆ | +0.43 −0.25 −1.83 +1.41 −2.79 +0.06 | −0.13 +0.13 +0.24 +0.15 +0.08 −0.17 | −1.37 −0.28 −0.56 +1.60 −0.38 +1.14 |
| `model_mix` pull | — | +0.71 (40 837 contro 40 753.5) | −0.53 (819 885 contro 820 163.7) |
| `model_mix` χ²/ndf (p) | — | 7.41/10 (0.686) | 5.51/10 (0.855) |

Tutto PASS. Gli stessi numeri valgono per `c2` / `c3` della fase 1b.

La tabella di `1b_results.md` riportava a₁ e a₂ del bordo come −1.38 e −0.29, un errore di
arrotondamento: i valori nei CSV, in `c3` come qui, sono −1.37498 e −0.28474. Corretto con questo
documento.

I run di kB da 10⁶ eventi non hanno soglie statistiche. Danno spettro u 51.81/49 (p 0.365),
profondità 46.17/49 (p 0.589) e frazioni 14.74/4 (p 0.0053, n–C non elastico −3.69σ). Sono i
primi 10⁶ eventi del mono, sotto la soglia di Bonferroni di un run da 10⁷ ma vicini, e
riportati per completezza.

**Confronto con la prima interazione della fase 1b** (previsione (d), non soglia: «CSV
identici»). Byte per byte coincidono:
- nel mono, 4 dei 5 CSV (tutti tranne `legendre.csv`);
- nel piatto e nel bordo, nessuno.

Il confronto cella per cella mostra però che **tutti i conteggi osservati sono identici** nei
tre run, per ogni bin di u, di z, delle frazioni e di `model_mix`. Le differenze stanno solo
nelle grandezze in virgola mobile sommate evento per evento: le attese sommate con l'E_n di ogni
evento (piatto e bordo) e i momenti di Legendre (tutti). Sono al più 2.4·10⁻¹³ in relativo, e
3.6·10⁻¹¹ per un momento di Legendre vicino a zero. Nel mono le attese non dipendono
dall'evento, e per questo coincidono byte per byte. La causa è la stessa dei run di kB: ordine
di somma (a). La prima interazione di ogni evento coincide con quella del run 1b, come previsto.

## 3. Confronto con PSTAR — (c)

**Luce, `pstar.light_MeV.relative_residual`:**

| Run | Banda | Osservato | Esito |
|---|---|---|---|
| mono 14 MeV | [−0.070%, −0.055%] | −0.0620% ± 0.00003% | PASS |
| piatto 10–20 MeV | [−0.070%, −0.052%] | −0.0600% ± 0.00003% | PASS |
| bordo 19.5–20 MeV | [−0.070%, −0.050%] | −0.0565% ± 0.00003% | PASS |

Il bordo è meno negativo del piatto, come previsto (b): ha più protoni sopra 14 MeV.

**Range, `pstar.range_mm.relative_residual`:** banda [−0.20%, 0.00%] per tutti e tre.

| Run | Osservato | Esito |
|---|---|---|
| mono | −0.1058% ± 0.0009% | PASS |
| piatto | −0.1030% ± 0.0010% | PASS |
| bordo | −0.1012% ± 0.0011% | PASS |

**Luce per bin di E_p.** La soglia è: entro 0.005% assoluto dal valore dello scan per lo stesso
bin, per i bin con `Ep_low` ≥ 2 MeV e almeno 100 protoni.

| Bin [MeV] | Scan (tabella §3) | mono | piatto | bordo |
|---|---|---|---|---|
| 2.0–2.4 | +0.08% | — | +0.0779% | +0.0779% |
| 3.2–3.6 (mono 3.08–3.36) | ≈ 0 | +0.0047% | −0.0038% | −0.0038% |
| 4.8–5.2 (mono 4.76–5.04) | −0.055% | −0.0546% | −0.0564% | −0.0564% |
| 7.6–8.0 (mono 7.56–7.84) | −0.0735% | −0.0735% | −0.0735% | −0.0735% |
| 10.0–10.4 (mono 9.80–10.08) | −0.070% | −0.0704% | −0.0698% | −0.0698% |
| 14.0–14.4 | −0.060% | — | −0.0591% | −0.0591% |
| 19.2–19.6 | −0.048% | — | −0.0479% | −0.0479% |
| continuità 19.6–20.0 | entro 0.005% dal bin prima | — | −0.0473% | −0.0473% |

Tutti PASS. Sotto 2 MeV il residuo è positivo e decresce con E_p in tutti e tre i run, come
previsto (d) dallo stopping nucleare escluso dalla luce di Geant4 (D21):
- mono: da +1.12% nel primo bin a +0.10% a 1.7 MeV;
- piatto e bordo: da +0.78% a +0.10%.

**Difetto delle previsioni.** La soglia per bin chiede il confronto con il valore dello scan per
**ogni** bin ≥ 2 MeV: 42 nel mono, 45 nel piatto e nel bordo. I valori per bin dello scan però
erano conservati fuori dal repo, in una cartella temporanea che è andata persa prima dei run. Si
possono confrontare solo i bin che cadono sui 7 punti ≥ 2 MeV della tabella della §3, sopra.
Gli altri bin non sono valutabili. Non è un esito positivo né negativo, e la soglia resta
com'è. Gli altri bin vanno da −0.0735% a +0.0643% (mono) e a +0.0779% (piatto e bordo).
Formano una funzione liscia di E_p che coincide fra piatto e bordo a 4 cifre in ogni bin
confrontato (c).

**Range per bin** (attesa, non soglia):
- sotto 3 MeV è nettamente positivo: +7.7% nel primo bin del mono e +5.5% nel piatto e nel
  bordo, poi +0.4% a 2 MeV, come nello scan;
- sopra 3 MeV sta fra −0.16% e +0.11% in tutti i run: da circa +0.1% a 3 MeV scende a un
  minimo di circa −0.15% verso 10 MeV.

L'errore per bin sopra 3 MeV è di qualche 10⁻³ %, quindi la forma è una misura (c). La causa
dell'eccesso sotto 3 MeV resta non studiata (d).

## 4. Sensibilità a kB — attese, non soglie

| Grandezza | kB = 0.145 | kB = 0.165 | Etichetta |
|---|---|---|---|
| `pstar.light_MeV.observed_mean`, variazione contro 0.155 | +3.758% | −3.483% | (c) |
| media della tabella §4 sui protoni contenuti, pesata con la luce | +3.784% | −3.500% | (b) |
| media della stessa tabella senza peso | +4.34% | −3.97% | (b) |
| attesa PSTAR, L(E_p; kB) sugli stessi protoni | +3.756% | −3.481% | (b) |
| media di `det_Evis` (da `deposit_spectra.csv`, centri dei bin di 0.28 MeV) | +3.28% | −3.06% | (c) |

Per costruzione, la variazione di una media di luce è la media delle variazioni pesata con la
luce. Con quel peso la previsione (b) della §4 cade entro 0.03 punti dall'osservato; la media
senza peso, scritta in modo ambiguo nelle previsioni, no. La differenza residua viene
dall'interpolazione lineare di una tabella a 6 punti. `det_Evis` si sposta meno della luce dei
protoni, come previsto: elettroni e ioni contribuiscono in parte con una luce meno sensibile a
kB. Questo è il valore da riportare come sistematica.

**Sistematica di kB** sull'intervallo [0.145, 0.165] di D18, per una sorgente mono a 14 MeV
(c):
- luce media dei protoni: +3.8% / −3.5%;
- `det_Evis` medio: +3.3% / −3.1%.

## 5. Attese che non sono soglie — (c)

**Classi dell'evento (D19)**, contro lo scan a 10⁴ eventi (errore dello scan 0.3–0.5 punti):

| Classe | mono (scan) | piatto (scan) | bordo |
|---|---|---|---|
| `no_interaction` | 37.58% (37.3%) | 36.36% (35.7%) | 37.95% |
| `nc_inelastic_no_np` | 17.49% (17.6%) | 18.18% (18.6%) | 19.55% |
| `nc_elastic_only` | 17.14% (17.0%) | 18.68% (18.5%) | 21.84% |
| `mixed` | 13.92% (14.1%) | 13.66% (14.2%) | 10.47% |
| `single_np` | 9.64% (9.8%) | 9.14% (9.1%) | 7.57% |
| `multiple_np` | 4.02% (3.9%) | 3.79% (3.8%) | 2.51% |
| `other` | 0.20% (0.22%) | 0.19% (0.17%) | 0.12% |

Lo scan è entro il suo errore tranne `mixed` nel piatto (−0.5 punti, circa 1.6σ dello scan). Il
bordo ha meno `single_np` del piatto, come previsto (b).

Altre attese, per 10⁴ eventi, contro lo scan:

| Grandezza | mono (scan) | piatto (scan) | bordo |
|---|---|---|---|
| `np_events` | 2778 (2808) | 2677 (2723) | 2065 |
| protoni contenuti confrontati con PSTAR | 2750 (2768) | 2647 (2685) | 2025 |
| `protons_exiting` | 13.1 (16) | 14.0 (17) | 19.1 |
| `protons_hadronic` | 14.6 (24) | 15.6 (21) | 20.5 |
| `events_secondary_n_hadronic` | 979 (999) | 976 (—) | 1060 |
| `Edep_overflow` | 281 (267) | 0.31 (2) | 3.1 |
| `Evis_overflow` | 0.002 (0) | 0 (0) | 0 |
| `expected_rare.p_Edep_above_Ep` | 5.6 (4) | 5.5 (—) | 4.3 |

Gli scan da 10⁴ eventi hanno errori di Poisson del 20–25% sui conteggi piccoli. L'unico scarto
oltre 2σ dello scan è `protons_hadronic` nel mono: 24 contro 14.6 attesi, cioè +2.5σ, ed è una
fluttuazione dello scan. `Evis_overflow` vale 2 eventi su 10⁷ nel mono.

## 6. Confronto con `endfb80` — attese, non soglie (c)

| Verifica | mono | piatto | bordo |
|---|---|---|---|
| frazioni, χ²/4 | 305.5 | 304.0 | 292.9 |
| n–C non el. / el. | +14.91σ / −10.86σ | +16.17σ / −10.26σ | +14.60σ / −12.02σ |
| spettro u, χ²/49 | 64.13 | 43.40 | 69.48 |
| pull di Legendre, massimo \|·\| | 3.19 (a₃) | 1.20 (a₃) | 3.00 (a₄) |

Sono gli stessi valori di `c2` della fase 1b: la prima interazione è la stessa. Resta il limite
di D13: a 10⁷ eventi la simulazione non distingue le due revisioni di H-1.

## 7. Verifiche `c2`, lettura in ordine di `event_id` (D23)

Il checker con `checks_version` 9 rilegge i sei run come `c2`, con gli stessi argomenti di `c1`
(`--bins-u 50 --bins-z 50`). Non ci sono nuovi run di simulazione e le soglie non cambiano.

**Previsione, scritta prima delle verifiche:**

- run di kB: le quattro verifiche di identità fallite in `c1` passano.
  `contained_balance_MeV` e `pstar.range_mm` sono identici fra i tre run. Restano identici anche
  le colonne non `Evis` di `deposit_spectra.csv` e `proton_exits.csv`. Ogni altra verifica dà lo
  stesso esito di `c1`;
- tutti i sei run, rispetto a `c1`:
  - ogni intero (contatori, classi, conteggi dei CSV) è identico;
  - ogni float differisce di al più 10⁻¹¹ in relativo, perché cambia solo l'ordine delle somme;
  - ogni soglia delle §1–§3 ha lo stesso esito.

**Risultati** (sei verifiche, ognuna in circa 1 minuto per 10⁷ eventi, RSS ≤ 200 MB):

| Verifica | Soglia o previsione | Osservato | Esito |
|---|---|---|---|
| §1–§3, tutte le soglie, sei run | come in `c1` | stesso esito di `c1` per ogni soglia | PASS |
| run di kB, `contained_balance_MeV` | identico | identico, `mean` compreso | PASS |
| run di kB, `pstar.range_mm` | identico | identico | PASS |
| run di kB, `Edep` e `proton_exits.csv` | identici byte per byte | identici | PASS |
| interi di `checks.json` e dei CSV, `c2` contro `c1` | identici | identici, sei run | come previsto |
| float, `c2` contro `c1`, run mono e di kB | ≤ 10⁻¹¹ relativo | ≤ 8.7·10⁻¹³ (Legendre) | come previsto |
| float, `c2` contro `c1`, run piatto e bordo | ≤ 10⁻¹¹ relativo | 267 valori oltre, fino a 1.8·10⁻⁸ | **previsione fallita** (d) |

Le righe di `residuals_Ep.csv` (una per evento n–p) ora escono in ordine di `event_id`. Come
insieme sono identiche a quelle di `c1`, in tutti i run.

**Perché la previsione dei float è fallita**, dagli stessi file:
- I valori oltre 10⁻¹¹ sono tutti nei due run a spettro piatto: pull, χ² e p-value di
  `spectrum_u`, `depth_z`, `fractions` e `model_mix`, più un coefficiente di Legendre.
- In quei run i conteggi attesi sono somme di pesi che dipendono da E_n, su 10⁷ eventi.
  L'arrotondamento di una somma cresce col numero di termini: lo scarto relativo arriva a
  10⁻¹⁰ sui χ².
- Un pull è (osservato − atteso)/σ. Quando osservato e atteso si avvicinano, lo scarto relativo
  si amplifica: 1.8·10⁻⁸ corrisponde a 3.7·10⁻¹¹ in assoluto. Lo scarto assoluto più grande
  su un pull è 3.1·10⁻¹⁰; su un χ² è 3.4·10⁻⁹.
- La previsione aveva preso la scala di `c1` dei run mono di kB (10⁻¹²). Non aveva tenuto conto
  dei pesi in E_n né della cancellazione nei pull.

L'esito di ogni soglia non cambia. D'ora in poi `c2` è il riferimento di questi run, perché non
dipende dalla divisione fra i thread.

## Punti aperti

- **Verifica indipendente dall'ordine:** risolta da D23 (`c2`, §7).
- **Warning `-Warray-bounds` nella build `release`** di `sim/tests/integration_reader.cpp`, da
  GCC 13 dentro `nlohmann::json`. C'erano già prima di D23. Da eliminare, perché la regola è
  zero warning.
- **Riferimenti per bin dello scan perduti.** Vanno conservati nel repo, per esempio in un CSV
  sotto `docs/validation/predictions/`, d'ora in poi per ogni previsione per bin.
- **Range sotto 3 MeV** (+5–8% nel primo bin): causa non studiata.
- **kB oltre 500 MeV/cm** (D18): la fine del range resta un'estrapolazione; i run di
  sensibilità coprono solo l'intervallo di Pöschl.
- **Offset HP di H-1** (D14), **bordo a 20 MeV** (D16), **eccesso di `ionIoni` sui rinculi C12**
  (D14) e **C-13** (D15): invariati, da decidere prima dell'analisi.
