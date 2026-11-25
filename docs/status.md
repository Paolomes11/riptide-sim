# Stato

**v0.8.0 — 2026-11-25 — Fase 6 (Caso A su dati Geant4) completata: con l'accettanza di D40
la famiglia Geant4 = toy passa 117/117 e la calibrazione assoluta 32/32
(`validation/6_results.md`). Prossima: fase 7 (Caso B), da progettare.**

`master` è a v0.8.0 (tag `v0.8.0`): la fase 6 fino a D40 è unita con la PR #7 (merge
`209707a`), il passo D e la chiusura con la PR #8 (merge `252d7be`); i due branch sono
cancellati. La fase 5
è unita con la PR #3 (tag `v0.7.0`), il report di progetto in `docs/riptide_report.pdf` con la
PR #4. La fase 5 comprende taglio `has_proton` (D34),
cinematica inversa del rinculo, motore del Caso A in `analysis/stats` (D32, D33), generatore di
pseudo-esperienze e finestra su θ_p (D35), app `riptide_stats_bench` e `riptide_stats_calibrate`.
La calibrazione C++ passa tutte le 71 soglie contro il toy (`validation/5_results.md`).
Il pull n–C non elastico lungo +X della fase 3 è chiuso come fluttuazione: con un seme
indipendente scende da +3.73 a +1.79 (`validation/3_nC_results.md`).

La fase 6 comprende decisioni D36–D40, run, 15 record del Caso A senza accettanza, valutazione,
diagnosi a 14 MeV e 15 record con l'accettanza (passo D). Senza accettanza Geant4 dava un
posterior circa il 10% più largo del toy e un bias di −0.023 MeV a 14 MeV, N = 1000. La causa è
l'accettanza in θ_p, che la likelihood non conteneva. Con q(θ_p | E_n, sel) di D40 il bias a
14 MeV scende a −0.0014 ± 0.0047 e tutte le voci gemelle e assolute passano.

## Fatto

- **Fase 0 (v0.1.0):** superbuild CMake con preset `debug`, `release`, `ci`; C++17; warning
  severi; `-ffast-math` vietato; `common/` (unità, schema, radice dei dati); `analysis/stats`
  scheletro.
- **Dati:** `scripts/data_root.sh` (`check` / `prepare`) per controllare e preparare l'SSD.
- **Dipendenze:** nlohmann/json 3.11.3 e Lyra 1.6.1 in `external/` (D12).
- **Config e CLI:** `RunConfig` da JSON (`config/mono_14MeV.json`, `config/flat_10_20MeV.json`),
  riga di comando di `riptide_sim` con Lyra.
- **Geant4:** mondo nel vuoto (D10), blocco di `G4_PLASTIC_SC_VINYLTOLUENE`, `QGSP_BIC_HP`,
  sorgente mono o piatta, run MT.
- **Semi per evento** da (event_id, seme del run), con invalidazione della cache delle sezioni
  d'urto del neutrone: output indipendente dal numero di thread (D11).
- **Verità MC:** prima interazione nello scintillatore (classificata per sottotipo di processo e
  nucleo bersaglio) e protone di rinculo del primo n–p.
- **Output:** ntupla `events` con schema v1 (19 colonne), un file per thread; `run.json` con
  config, semi, versioni e provenienza, scritto per ultimo come marcatore di completamento.
- **Test:** `debug` 56/56, `ci` 54/54, zero warning, clang-format pulito. `sim_integration`
  (2×1000 eventi a 14 MeV, circa 2 s) verifica `run.json`, unicità degli `event_id`, coerenza
  delle colonne di verità e uguaglianza colonna per colonna fra 1 e 2 thread; senza SSD viene
  saltato.

### Fase 1b (v0.3.0, branch `feat/sim-checks`, `checks_version` 5)

Verifiche fisiche non circolari e prestazioni, vedi `roadmap.md`.

- `scripts/mount_ssd.sh`: monta l'SSD per UUID (UUID fuori dal repo), senza toccare `/etc/fstab`.
- **D13:** valutazioni di riferimento JEFF-3.3 (da cui deriva G4NDL: verifica
  dell'implementazione) ed ENDF/B-VIII.0, sull'SSD con URL e sha256 in
  `architecture/data_layout.md`. Per H-1 le due librerie sono la stessa analisi LANL: il
  confronto n–p indipendente è rimandato. Per il carbonio, JEFF-3.3 C-0 ed ENDF/B-VIII.0
  C-12 + C-13.
- `common/np_kinematics.hpp`: E_p(E_n, θ_p) relativistica, mappa E_p ↔ cos θ*, polinomi di
  Legendre.
- `sim/validation` (senza Geant4, anche in `ci`):
  - lettore ENDF-6 (sezioni d'urto MF3, coefficienti di Legendre MF4, leggi d'interpolazione
    1–5);
  - statistiche delle verifiche (χ² con probabilità di sopravvivenza, frazioni multinomiali,
    momenti di Legendre, fusione dei bin con attesa < 5);
  - distribuzioni attese (spettro E_p da Legendre, profondità del primo vertice, frazioni della
    prima interazione);
  - sha256, corda del fascio nel blocco, `CheckAccumulator` (attese sommate evento per
    evento), report JSON/CSV, riga di comando di `riptide_sim_check`;
  - test sui file reali dell'SSD, saltati se il disco manca.
- `riptide::sim_reader`: lettura di `run.json` e delle ntuple, anche in streaming
  (`forEachEventRow`).
- App `riptide_sim_check` (layout dell'output in `architecture/data_layout.md`), più lo smoke
  test `sim_check_smoke` sul run a 1 thread di `sim_integration`. Sul run da 1000 eventi a
  14 MeV:
  - RSS 45 MB in 0.45 s; a 10⁷ eventi si stima circa 0.1 GB (16 B per evento n–p);
  - n_H1 = 5.26·10²² cm⁻³, n_C = 4.73·10²² cm⁻³, n_H2/n_H1 = 1.15·10⁻⁴, L = 100 mm;
  - contatori di coerenza a zero.

  Questi numeri servono solo come controllo del funzionamento, non sono un risultato.
- **Residuo di E_p spiegato (D14).** Il modello HP di Geant4 usa per H-1 la massa dell'atomo
  (AWR = 0.999167 in G4NDL 4.7.1) e crea poi il rinculo con la massa del protone. Le verifiche
  modellano esattamente questa convenzione (`sim/validation/hp_recoil.hpp`), con la massa letta
  da G4NDL, e `checks.json` è passato alla versione 2. Nello smoke test:
  - il residuo di E_p è al più 1.2·10⁻¹⁴ MeV, cioè arrotondamento;
  - lo scarto dalla cinematica fisica è (5.37–5.49)·10⁻⁴ in relativo, con media 5.43·10⁻⁴;
  - i contatori u > 1 e θ_p > π/2 ora devono essere zero, e lo sono.
- **Carbonio HP (D15).** Verificato sui file: l'inelastico HP del carbonio esclude la cattura
  (coincide con JEFF-3.3 C-0 MT3 − MT102; integrale su 10–20 MeV entro 2.5·10⁻⁷). G4NDL non ha
  C-12 e compone il carbonio come 0.9893 · C-0 + 0.0107 · C-13, contando il C-13 due volte:
  l'attesa `jeff33` ora è composta allo stesso modo (`checks_version` 3). Lettore dei file di
  sezione d'urto G4NDL (`sim/validation/hp_data.hpp`) e oracolo permanente G4NDL ↔ JEFF-3.3
  per H-1, C-0 e C-13 (`sim/tests/hp_data_reference_test.cpp`). Residuo noto e non modellato:
  circa 3·10⁻⁴ sul non elastico del carbonio vicino a 14 MeV, dall'inelastico di C-13.
  Per H-1, G4NDL non ha un file inelastico: sotto 20 MeV le prime interazioni possibili sono
  solo l'elastico e la cattura.
- **Bordo a 20 MeV (D16):** in 19.5–20 MeV l'elastico n–p è HP o CHIPS secondo una rampa
  lineare; le verifiche contano la miscela attesa (`model_mix`).
- **Previsioni prima dei run** (`validation/1b_predictions.md`, commit `6de0124`): soglie
  Bonferroni con α_family = 0.01 a run, range dell'offset di D14, residui noti.
- **Run:** tre run da 10⁷ eventi (mono 14 MeV, piatto 10–20 MeV, bordo 19.5–20 MeV), 8 thread,
  circa 4.5 minuti e 145 MiB ciascuno; 54–62 B/evento, da cui `kBytesPerEvent` = 128 B.
- **Esito (D17, `validation/1b_results.md`):** al primo esito il mono falliva su due verifiche,
  per due difetti delle attese: il break-up n–d su H-2 non era previsto, e l'oracolo dell'offset
  risentiva dell'arrotondamento a E_p di pochi eV. Corrette le attese (H-2 da G4NDL,
  `checks_version` 5) e rilanciate le verifiche sugli stessi eventi, con soglie invariate e
  previsioni committate prima: tutte le verifiche passano nei tre run.
- **Test:** `debug` 172/172, `ci` 157/157, zero warning, clang-format pulito.

### Fase 2 (v0.4.0, branch `feat/sim-deposits`, `checks_version` 9)

Deposito, Birks, traccia del protone e classificazione dell'evento, vedi `roadmap.md`.

- **D18:** kB = 0.155 mm/MeV per il plastico vinil-toluenico, da letteratura citata.
- **Config v2:** blocco `scintillator` (kB e fonte) e funzione di passo del protone.
- **Libreria:** registrazione dei depositi per evento (totale e del protone di rinculo, con e
  senza Birks), traccia del protone (lunghezza, punto finale, uscita, interazioni adroniche),
  classificazione dell'evento in 7 classi dai contatori d'interazione.
- **Schema v2 (D19):** colonne di deposito, energia visibile, traccia del protone e classe
  dell'evento. Il protone può depositare fino a circa 10⁻⁴ MeV più di `truth_Ep` per
  G4NuclearStopping a fine percorso; tolleranza 0.01 MeV.
- **Verifiche (`checks_version` 6):** sezione `deposits` di `checks.json` (contatori che devono
  essere zero, rari attesi, classi, bilancio dei protoni contenuti), `deposit_spectra.csv`,
  `proton_exits.csv`; per i run con schema v1 la sezione è `skipped`.
- **Oracolo kB = 0:** terzo run `ctest_it_kB0` in `sim_integration`, con E_vis = E_dep entro
  10⁻⁹ MeV e depositi identici al run a 1 thread.
- **Smoke test** (1000 eventi a 14 MeV): contatori che devono essere zero tutti a zero; 292
  protoni contenuti con bilancio massimo 1.9·10⁻⁴ MeV e nessun bilancio negativo. Solo
  controllo del funzionamento, non un risultato.
- **PSTAR (D20, 2026-11-20):** tabella NIST del plastico vinil-toluenico sull'SSD con sha256,
  lettore in `sim/validation` (S e R_CSDA interpolati log-log, L(E_p; kB) con Simpson su S
  totale) e sezione `deposits.pstar` di `checks.json` (`checks_version` 7) con
  `pstar_comparison.csv`: residui di `truth_p_length` contro R_CSDA e di `truth_p_Evis` contro
  L, per bin di E_p, per i protoni contenuti senza interazioni adroniche.
- **Osservazione (d) dello smoke test** (292 protoni, solo funzionamento, non un risultato):
  range relativo −1.1·10⁻⁴ ± 9·10⁻⁴; luce relativa +0.57% ± 0.04% sopra L di PSTAR. Spiegata e
  rimossa da D21.
- **Regressione v7 sui run 1b (`c3`, 2026-11-20):** CSV identici byte per byte a `c2`, in
  `checks.json` cambiano solo le versioni e `check_id`, `deposits` saltata (schema 1). Dettagli
  in `validation/1b_results.md`.
- **Banda di overflow spiegata (D14, 2026-11-20):** i 27 eventi a circa +7.5 keV sopra E_n sono
  l'energia creata dall'offset HP in ogni urto n–p. Σε_np vale in media 7.4 keV; il resto della
  dispersione è il neutrone termalizzato che esce dal blocco. Il bilancio per evento chiude
  entro 4·10⁻¹⁵ MeV con un diagnostico fuori dal repo, che riproduce il run bit per bit (c).
  La cattura n–p conserva l'energia: la previsione di +1.3 keV per cattura era sbagliata.
- **Luce continua dei protoni (D21, 2026-11-20):** Birks per step dipendeva dalla funzione di
  passo, da +0.51% a −0.39% sulla griglia, senza limite (c). Ora per i protoni la luce di uno step
  è L(K_pre) − L(K_pre − E_dep) dallo S elettronico di Geant4. Sullo stesso scan il residuo contro
  PSTAR vale −0.062% ± 0.001% per tutte le config (c): è la differenza fra lo S di Geant4 e lo S
  totale di PSTAR. La funzione di passo resta (0.2, 0.1 mm).
- **Intervallo di kB (D18, 2026-11-20):** [0.145, 0.165] mm/MeV, ±2σ di Pöschl, solo da quella
  fonte. Laplace 2020 letto: non dà un kB.
- **Pilota v2 (2026-11-20, non un run di validazione):** 10⁵ eventi a 1 thread in 13.0 s e
  4·10⁵ a 8 thread in 9.6 s, circa 113 B/evento. Nel secondo un protone contenuto ha
  `truth_p_Edep` − `truth_Ep` = −2.027 MeV: uno scattering Coulombiano singolo ha creato un
  rinculo di C12 da 2.027 MeV (c, run ripetuto con la traccia degli step).
- **Rinculi nucleari (D22, 2026-11-20):** colonna `truth_p_n_nuclear_recoil` (schema 3); le
  verifiche dei depositi (`checks_version` 8) escludono quei protoni dal bilancio e da PSTAR, li
  contano in `protons_nuclear_recoil` e saltano i run con schema < 3.
  `validation/2_predictions.md` è alla versione 2, corretta dopo il pilota e detta tale.
- **Pilota v3 (2026-11-20, non un run di validazione):** `pilot2_v3_mono_14MeV_1e5_1t` in
  13.3 s, 89 MiB, 114 B/evento; `pilot2_v3_mono_14MeV_4e5_8t` in 9.5 s, 145 MiB, 115 B/evento.
  Verifica c1: contatori che devono essere zero a 0 in tutti e due; `protons_nuclear_recoil` 0 e
  1 (lo stesso evento del pilota v2, che riproduce evento per evento); minimo del bilancio
  −5.3·10⁻¹⁵ MeV; luce contro PSTAR −0.062%, range −0.107% e −0.103%. Stima per la produzione
  (b): circa 4 minuti e 1.15 GB per ogni run da 10⁷ a 8 thread.
- **Run di validazione (2026-11-21, `validation/2_results.md`):** tre da 10⁷ eventi (mono
  14 MeV, piatto 10–20 MeV, bordo 19.5–20 MeV) e tre da 10⁶ a kB = 0.145 / 0.155 / 0.165, a 8
  thread, contro le previsioni v2 (`85283ba`) con soglie invariate.
  - **Prestazioni (c):** circa 4.5 minuti e meno di 170 MiB di RSS per 10⁷ eventi;
    112.9–120.9 B/evento, quindi `kBytesPerEvent` = 256 B resta; 3.82 GB in tutto.
  - **Oracoli (a):** contatori che devono essere zero a 0 in tutti i run; bilancio dei protoni
    contenuti fra −1.1·10⁻¹⁴ e 9.2·10⁻⁴ MeV; rinculi nucleari (D22) 40–57 per 10⁷ eventi, nella
    banda 1–150.
  - **Statistiche contro `jeff33` (c):** tutte sotto le soglie Bonferroni. I CSV di prima
    interazione corrispondono a quelli della fase 1b.
  - **PSTAR (c):** luce da −0.062% (mono) a −0.057% (bordo) rispetto a L di PSTAR, range
    −0.10%. Nel confronto per bin con lo scan, entro 0.005%.
  - **Sensibilità a kB (b):** ±0.01 mm/MeV su kB sposta la luce media dei protoni di +3.8% /
    −3.5%. Depositi e traccia identici, byte per byte.
  - **Verifiche `c1`:** quattro identità fallite nei run di kB, a 10⁻¹² relativo, per l'ordine
    di somma nel checker; la previsione era sbagliata (d).
- **Lettura in ordine di `event_id` (D23, `checks_version` 9):** il reader fonde i file per
  thread in ordine di `event_id`, quindi le somme non dipendono più dalla divisione degli eventi
  fra i thread. Lo smoke test confronta byte per byte i run a 1 e a 2 thread. Le verifiche `c2`
  sugli stessi run passano tutte, comprese le quattro identità. Lo scarto da `c1` è fino a
  1.8·10⁻⁸ relativo sui pull dei run piatti, oltre i 10⁻¹¹ previsti (d).
- **Decisioni D24–D26 (2026-11-21):** la cinematica HP passa da `sim/validation` a
  `common/hp_kinematics.hpp`; `physicalProtonEnergyFromHp` toglie l'offset di E_p di D14, da
  applicare in `analysis/io` prima di σ_Ep (D24). Sorgente ≤ 19.5 MeV nelle fasi 5–6 (D25);
  eccesso di `ionIoni` come sistematica del fondo da carbonio (D26).
- **Test:** `debug` 241/241, `ci` 222/222, zero warning nel codice del progetto, clang-format
  pulito. In `release` restano due `-Wnull-dereference` da un header di Geant4
  (`G4NavigationHistory.icc`), presenti anche prima di D21, e quattro `-Warray-bounds` da
  nlohmann/GCC 13 in `sim/tests/integration_reader.cpp`, presenti anche con il file di
  `83b959d`.

### Fase 3 (v0.5.0, branch `feat/sim-geometry`, `checks_version` 10)

Geometria reale e sorgente a campo lontano, vedi `roadmap.md`.

- **D27:**
  - blocco di BC-408 da 60 × 20 × 20 mm centrato nell'origine, lato lungo su X, mondo nel vuoto
    senza supporti;
  - sorgente a campo lontano: direzione fissa o isotropa, punto di partenza uniforme su un disco
    ⟂ Ω_n fuori dal blocco;
  - il fascio a matita resta per la regressione.
- **Config v3:** `source.beam` `"pencil"` o `"far_field"` (con `direction` vettore o
  `"isotropic"` e `margin_mm`). I config v2 si leggono come pencil. Nuovi config
  `3_ff_mono_14MeV_x`, `3_ff_mono_14MeV_isotropic`, `3_ff_flat_10_19p5MeV_isotropic`.
- **Generatore:** posizione e direzione campionate per evento con l'RNG seminato per evento
  (D11), quindi l'output non dipende dal numero di thread. Funzioni pure di campionamento con
  test di direzione, disco, accettanza e corda media di Cauchy.
- **Schema 4:** colonne `truth_n_start_{x,y,z}` (mm), il punto di partenza del neutrone.
- **Verifiche (`checks_version` 10):**
  - corda L_i per evento, con frazioni di prima interazione e profondità s = z/L_i calcolate su
    quella corda;
  - nuova sezione `geometric_acceptance`: pull binomiale dell'accettanza e pull della corda
    media;
  - nuovi contatori che devono essere zero, `start_not_on_source` e `interaction_off_chord`;
  - la previsione media le attese sulla distribuzione delle corde (10⁸ raggi);
  - **oracolo di regressione (c):** su `2_mono_14MeV_kB0p155_1e6` (`c4`, 50 bin) e
    `2_flat_10_20MeV_1e7` (`c3`), v10 ridà i numeri di v9 a parte le chiavi nuove e i bordi
    della profondità divisi per L.
- **Previsioni (`validation/3_predictions.md`):**
  - v1 a `8429cc3`, v2 a `ff6f837` con i riferimenti PSTAR per bin nel repo
    (`predictions/3_pstar_light_reference_*.csv`);
  - attese `jeff33` in `predictions/3_ff_*_1e7.json`;
  - soglie Bonferroni |pull| < 3.29 (+X) e < 3.32 (isotropi).
- **Piloti da 10⁵ eventi a 1 thread (non run di validazione):** `pilot3_*_1e5_1t/c1` sull'SSD.
  - Contatori che devono essere zero a 0.
  - Frazione di protoni che escono 3.35% / 3.98% / 4.85% (+X / isotropa mono / isotropa
    piatta), contro lo 0.47% del cubo della fase 2.
  - 7–8 s per 10⁵ eventi, 91 MB, 110–119 B/evento. Stima (b) per un run da 10⁷ a 8 thread:
    2–6 minuti, circa 1.2 GB.
- **Run di validazione (2026-11-22, `validation/3_results.md`):** tre da 10⁷ eventi, build di
  `ded6c88`, verifiche `c1`.
  - Tutte le soglie passano: oracoli, χ² e pull contro `jeff33`, bande e bin di PSTAR.
  - Frazione di protoni che escono 3.17% / 4.37% / 4.76%, che cresce con E_p fino al 17% sopra
    15.6 MeV nel piatto isotropo.
  - 2.4–2.8 minuti e meno di 0.2 GB per run, 110–119 B/evento: `kBytesPerEvent` = 256 resta.
  - n–C non elastico lungo +X a +3.73σ da `jeff33` (il χ² delle frazioni passa con p = 0.0027).
    Chiuso il 2026-11-23 come fluttuazione: con il seme 20260908 il pull è +1.79, e il
    diagnostico per profondità e per facce non mostra struttura (`validation/3_nC_results.md`).
  - La previsione (d) sul verso del range contro PSTAR non è confermata.
- **Test:** `debug` 266/266, `ci` 247/247, clang-format pulito.

### Fase 4 (v0.6.0, branch `feat/analysis-io`, `analysis_record_version` 2)

I/O dell'analisi, vedi `roadmap.md`.

- **D28:** `run.json` passa a `run_record_version` 2, con il blocco `hp_elastic_H1` (file,
  sha256, AWR, `target_mass_MeV`). Per i run v1 la massa si prende da `--hp-from-check <id>` di
  una verifica dello stesso run, mai da un default.
- **D29:**
  - `riptide::io_core` è senza ROOT e compila anche in `ci`;
  - `riptide::io_root` contiene solo il lettore `RDataFrame`. Rifiuta le colonne mancanti,
    riordina per `event_id` e verifica che gli id siano unici;
  - `common/event_sample.hpp` (struct-of-arrays) e `common/random.hpp` (SplitMix64, motore per
    evento, Box–Muller propria) servono anche a `analysis/stats`;
  - opzione `RIPTIDE_BUILD_ROOT_IO`, ON per default e OFF nel preset `ci`.
- **D30:** catena di tagli `single_np`, `contained`, `evis_min=<MeV>` con cut-flow; il default
  per la fase 5 è `single_np,contained`.
- **D31:**
  - `Ep_phys` dalla cinematica HP inversa, poi `reco_Ep = Ep_phys + N(0, σ_Ep)`;
  - `reco_theta_p = truth_theta_p + N(0, σ_θ)` per il Caso A;
  - `reco_p_dir` spostato nel piano tangente per i Casi B e C;
  - motore per evento da (seme, `event_id`);
  - addendum: un evento senza traccia del protone ha `reco_p_dir` = (0, 0, 0) e nessuna
    estrazione per la direzione.
- **`riptide_analysis_prepare`:** scrive `analysis.json` (impostazioni, fonte della massa,
  cut-flow, riassunti dei residui) e `cutflow.csv` sotto `analysis/<run_id>/<analysis_id>/`. Il
  campione resta in memoria.
- **Previsioni (`validation/4_predictions.md`):** v1 a `d19fd67`, prima dei run.
- **Run di validazione (2026-11-22, `validation/4_results.md`):** i tre run da 10⁷ eventi della
  fase 3, letti in sola lettura.
  - Primo giro `a1`: `p_dir.pull_mean` era `null` per un NaN negli eventi senza protone. È un
    difetto del codice, corretto in `6f8817e`.
  - Secondo giro `a2`: tutte le soglie passano, senza cambiarle. I conteggi di `single_np` sono
    esatti, offset HP nella banda analitica, 15 pull sotto 3.3, riproducibilità bit per bit.
  - 9–10 s e 1.37 GB al massimo per un run da 10⁷ eventi, circa 110 B/evento sopra una base di
    0.3 GB.
- **Test:** `debug` 321/321, `ci` 299/299, clang-format pulito.

### Fase 5 (v0.7.0, branch `feat/stats-case-a`, `calibration_version` 1)

Port del motore statistico del Caso A, vedi `roadmap.md`.

- **D32:** cinematica esatta `riptide::protonEnergy` con masse n ≠ p e l'inversa
  `recoilCos2FromEnergy` / `protonAngleFromEnergy` in `common/np_kinematics.hpp`; griglie in
  double, E_n su un dominio parametrico (default [10, 20] MeV, 500 punti) e θ_p su [0, π/2]
  (500 punti). Niente float32.
- **D33:** `CaseAModel` con la tabella Ep(E_n, θ_p) precalcolata e la marginalizzazione di θ_p con
  logsumexp; posterior per evento e combinato con il prior contato una volta; prior su θ_p come
  argomento (sin 2θ_p con celle esatte, oppure piatto), prior su E_n piatto e proprio.
- **D34:** taglio `has_proton` (`truth_valid_np == 1`) nella catena D30; il default del Caso A
  diventa `single_np,has_proton,contained`. Le estrazioni degli altri eventi non cambiano.
- **D35:** generatore delle pseudo-esperienze del modello del toy, motore per evento; finestra su
  θ_p con limite d'errore C = 50, identica bit per bit alla somma completa sul benchmark.
- **`analysis/stats`:** griglie, prior, logsumexp, modello, combinazione, validazione (momenti,
  intervalli a code uguali, coverage, pull, curva del bias, checklist), pseudo-esperienze e
  calibrazione. Nessun ROOT né Geant4; compila in `ci`.
- **Oracoli (c):** R1, R2, R3 e i test di limite del Caso A in Catch2, verdi.
- **Previsioni (`validation/5_predictions.md`):** v1 a `c9e7669`, prima dei run.
- **Calibrazione (2026-11-22, `validation/5_results.md`):** `riptide_stats_calibrate` `c1`.
  - Tutte le 71 soglie passano, scarto massimo 1.75σ.
  - Checklist dal prior: pull −0.002 / 1.001, coverage 0.680 / 0.900 / 0.953; bias per bin
    entro la soglia dal toy.
  - `std_sqrt_N` 1.36 / 0.945 / 0.925 a N = 10 / 100 / 1000.
  - Prior piatto su θ_p: bias −0.097 ± 0.006 MeV a 15 MeV (toy −0.101 ± 0.009).
  - Spostamento dal prior su [10, 30]: +0.93 / +0.35 / +0.11 / +0.03 / +0.01 std a
    N = 1 / 3 / 10 / 30 / 100.
  - 39.6 s e 19 MB per il run completo.
- **Prestazioni (c):** 0.174 s per 1000 eventi a thread singolo, contro 0.71 s del toy. Il target
  D35 (≤ 0.10 s) è mancato; il tempo residuo è dominato da `exp` scalare.
- **Test:** `debug` 367/367, `ci` 345/345, clang-format pulito.

### Fase 6 (v0.8.0, branch `feat/case-a-geant4` e `feat/case-a-acceptance-runs`, `case_a_version` 2)

Caso A sugli eventi di Geant4, vedi `roadmap.md`. Fonti dei numeri: `validation/6_predictions.md`
e, sulla radice dei dati, `validation/6_eval/`, `validation/6_diag_14MeV/`,
`validation/6_acceptance/`, `validation/6a_eval/`; sintesi in `validation/6_results.md`.

- **D36–D39 (af81a89):** campioni (piatto 10–19.5 MeV +X da 10⁷, mono 11/17/19 MeV da 10⁶, più
  il 14 MeV della fase 3), esperienze disgiunte da una permutazione, toy gemello `--toy` e effetti
  attesi da Geant4 che il toy non ha.
- **Codice:** split in esperienze disgiunte (1738102), sorgente del run da `run.json` (00fa20d),
  funzioni pure (f080749), app `riptide_case_a_geant4` (7e676ca), config (6055ab1).
- **Previsioni (`validation/6_predictions.md`, 994ece0), prima dei run:** famiglia gemella
  K = 117 (soglia z = 3.93), famiglia assoluta K = 32 (z = 3.60), spostamento atteso del bias
  globale sul piatto +0.19 MeV.
- **Run (2026-11-23, build 994ece0):** i quattro nuovi run passano `c1`, con l'oracolo HP su
  `rounding_free_relative_*` (D17). 15 record in `stats/6_<run_id>{,_pilot,_toy}`.
  Selezionati: piatto 799076, 11 MeV 98693, 14 MeV 826770, 17 MeV 68736, 19 MeV 61069.
  - **Correzione:** a 14 MeV le previsioni riportano 826859, che è la catena della fase 4 senza
    `has_proton` (89 eventi in più). Il valore giusto è 826770.
- **Valutazione (`validation/6_eval/`):**
  - famiglia gemella 107/117: 9 FAIL su `std_sqrt_N` a N = 100 e 1000 (G4 circa il 10% più
    largo del toy a tutte le energie, z fino a +8.7) e `bias_mean` a 14 MeV, N = 1000
    (G4 −0.0227, toy +0.0072, z −4.25);
  - famiglia assoluta: toy 32/32, G4 28/32 (`pull_mean` a N = 1000 a 11/14/17 MeV e a N = 100 a
    11 MeV): lì il toy è circolare (D38);
  - bias globale sul piatto: osservato +0.166, atteso +0.180, residuo −0.014 con tolleranza
    0.101, PASS; coverage G4 0.703 / 0.919 / 0.964, toy 0.707 / 0.922 / 0.965;
  - le previsioni non attendevano nessuna voce di contrazione oltre soglia: **anomalia**.
  - Minore: `std_sqrt_N` del toy a 19 MeV cambia del 12% fra N = 100 e 1000 (previsto ≤ 10%).
- **Diagnosi a 14 MeV (`validation/6_diag_14MeV/`):**
  - D4 riproduce `case_a.json` esattamente;
  - cinematica respinta: scarto in E_p ~10⁻⁵ MeV;
  - larghezze per evento uguali fra G4 e toy;
  - accettanza confermata: rapporto θ_p G4/toy 0.58 a 0–5°, circa 1.06 a θ intermedi, 0.92
    sopra 80°; il toy ripesato alla distribuzione di G4 riproduce G4 (`std_sqrt_N` 1.010 / 0.942
    contro 1.018 / 0.932, bias −0.021 contro −0.023).
- **Misura dell'accettanza (`validation/6_acceptance/`):** il deficit in avanti viene da
  `single_np` (n–p multiplo e misto 55% a 0–5°), `contained` perde al più il 5.6%,
  l'anisotropia del primo urto vale ±6%. La forma dipende da E_n in modo non liscio
  (χ² lineare 420/306).
- **D40 (0385bc6):** q(θ_p | E_n, sel) al posto del prior sin 2θ_p, dal rapporto r(E_n, θ_p) su
  bin 1 MeV × 5° di un run piatto indipendente (seme 20260908), interpolato bilinearmente;
  toy gemello per accept-reject; percorso di default invariato. Previsioni nella §10 di
  `6_predictions.md`.
- **Codice di D40:** prior per riga in `CaseAModel` (d16022f), `AcceptanceRatio` (077d297), toy
  con accept-reject (9271b03), `acceptance.json` v1 (30b01b3), app `riptide_acceptance` e
  `--acceptance` in `riptide_case_a_geant4`, `case_a.json` v2 (6729a54). Percorso di default
  verificato identico bit per bit sui record G4 e toy a 14 MeV; 41.6 s contro 40 s.
- **Passo D (2026-11-25, build 7ec81bb):**
  - run piatto `6_ff_flat_10_19p5MeV_x_1e7_s20260908`: `c1` PASS, 798388 selezionati (−688
    dal piatto della fase 6, entro ±2700); tabella `6acc_s20260908`;
  - 15 record `6a_<run_id>{,_pilot,_toy}`: famiglia gemella 117/117 (max|z| 3.16), assoluta
    G4 32/32 e toy 32/32, bias globale del piatto PASS (residuo −0.019, tolleranza 0.102);
  - bias G4 a N = 1000: −0.002 / −0.001 / −0.002 / +0.014 MeV a 11 / 14 / 17 / 19 MeV;
  - il toy converge verso G4 (`std_sqrt_N` a 14 MeV, N = 100, da 0.920 a 1.026); G4 cresce del
    +1.3…+3.6%, con il ±3% previsto **mancato** a 19 MeV, N = 1000;
  - prova a metà (q dai soli `event_id` pari): spostamento 0.246·σ_c contro 0.3, PASS;
  - tempi ≤ 42 s per `case_a`, RSS ≤ 1.42 GB.
- **Risultati (`validation/6_results.md`, f54b787).**
- **Test:** `debug` 406/406, `ci` 382/382, clang-format pulito.

**Prossimo passo:** progettare la fase 7 (Caso B) in una sessione dedicata: decisioni D41+ e
previsioni prima del codice.

## Punti aperti

- **Fase 3, deciso (D27):** blocco 60 × 20 × 20 mm, mondo nel vuoto senza supporti né
  materiali attorno. Se una geometria più realistica li aggiungerà, andranno rivisti fondo e
  accettanza.
- SSD esterno montato su `/mnt/external_ssd`, ma non in modo permanente: prima di ogni run
  `scripts/data_root.sh check` (procedura in `architecture/data_layout.md`).
- **Fase 2, deciso (D25):** sopra 19.5 MeV lo stato finale elastico è in parte CHIPS (D16):
  sorgente ≤ 19.5 MeV nelle fasi 5–6, bordo come sistematica nella fase 9.
- **Fase 2, deciso (D24, D28):** la verità di Geant4 ha un E_p più alto della cinematica fisica
  di circa 5.4·10⁻⁴ in relativo (D14). `analysis/io` lo corregge con `physicalProtonEnergyFromHp`,
  con la massa HP da `run.json` v2 o da una verifica dello stesso run (D28).
- **Fase 4, deciso (D34):** 89 / 17 / 19 eventi `single_np` contenuti non hanno un protone
  tracciato (`truth_Ep_phys` = 0). Escono con il taglio `has_proton`; la frazione persa (~10⁻⁴)
  è una sistematica per la fase 9.
- **Fase 3, sotto la sensibilità (d):** un eccesso residuo di n–C non elastico lungo +X dello
  0.1–0.2% non è escluso. Non tocca il Caso A, perché il taglio `single_np` toglie gli eventi
  n–C; va rivisto nella fase 9, quando gli eventi n–C entreranno nel modello.
- **Fase 6 (d):** a 19 MeV, con q, il bias di G4 a N ≥ 100 è positivo (+0.014/+0.021 MeV,
  entro soglia) e `std_sqrt_N` cresce del 3.6%. Primo candidato da verificare: il bordo della
  tabella q a 19.5 MeV. Da riprendere con le sorgenti fino a 20 MeV (fase 9).
- **Fase 6 (b):** `physicalRecoilAngleLogPrior` usa la massa non relativistica di θ_p, mentre il
  toy campiona con la cinematica relativistica (scarto ±0.7% a 14 MeV, ±1% a 20 MeV). Con q non
  pesa; senza q va corretto o dichiarato.
- **Fase 6 (a):** la forma di q in θ_p dipende da E_n in modo non liscio (χ² lineare 420/306);
  l'origine non è verificata.
- **Fase 6 (d):** `std_sqrt_N` del toy a 19 MeV cambia dell'11–12% fra N = 100 e 1000, con q e
  senza (criterio ≤ 10%).
- **Fase 5:** il target D35 sulla verosimiglianza (≤ 0.10 s per 1000 eventi) è mancato, con
  0.174 s. Un `exp` vettoriale o il parallelismo sugli eventi richiedono una decisione separata
  e vanno rivisti solo se la fase 6 o i Casi B/C ne avranno bisogno.
- **Fase 4:** nella build `release` restano `-Wnull-dereference` e `-Warray-bounds` da header
  esterni (libstdc++, nlohmann/json, Geant4); nessuno viene dal codice del progetto.
- **Fase 2:** il doppio conteggio del C-13 nel modello HP (circa +0.5% sul non elastico del
  carbonio) va tenuto fra le sistematiche del fondo da carbonio (D15).
- Confronto n–p indipendente da JEFF-3.3/ENDF-B-VIII.0 (D13) e verifica di H-2 contro una
  valutazione indipendente (D17): non fatti in 1b.
- **Fase 2, deciso (D26):** a fine percorso `ionIoni` deposita più del K iniziale dello step. Riguarda 56
  rinculi C12 su 486 e 4 Be9, con K fra 0.11 e 0.40 MeV; l'eccesso arriva al 15%, per un totale
  di 0.74 MeV su 492 MeV di rinculi C12 (D14). La causa non è verificata; non dipende dalla
  funzione di passo del protone (D21, P4). Pesa su `det_Edep` più che sulla luce. Resta una
  sistematica del fondo da carbonio, da studiare solo se nella fase 9 conterà.
- **Fase 2:** kB è misurato fino a 500 MeV/cm; la fine del range dei protoni va oltre, quindi lì
  la luce è un'estrapolazione non coperta dall'intervallo di D18. Un fit di Birks sui dati EJ-200
  di Laplace 2020 la sonderebbe.
- **Fase 2:** le verifiche dei depositi devono tollerare `det_Edep` − E_n fino a circa
  6·10⁻⁴ · E_n senza un γ di cattura (D14).
- **Fase 2:** il range sotto 3 MeV supera R_CSDA di PSTAR del 5–8% nel primo bin; la causa non
  è studiata.
- I riferimenti per bin delle previsioni (scan) vanno conservati nel repo: quelli della fase 2
  sono andati persi, e si sono potuti confrontare solo 7 punti della tabella.
