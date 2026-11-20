# Stato

**v0.3.0 — 2026-11-19 — Fase 1b (verifiche fisiche e prestazioni) completata.**

Fase 2 in corso (branch `feat/sim-deposits`, aggiornato 2026-11-20): passi 1–5 del piano fatti,
più i confronti con PSTAR (`checks_version` 7), la spiegazione della banda di overflow (D14), la
luce continua dei protoni (D21), l'intervallo sistematico di kB (D18) e le previsioni dei run
(`validation/2_predictions.md`, scritte prima dei run). Il pilota ha trovato i rinculi nucleari
del protone (D22): schema 3, `checks_version` 8, previsioni alla versione 2.

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

### Fase 2 (in corso, branch `feat/sim-deposits`)

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
- **Test:** `debug` 235/235, `ci` 219/219, zero warning nel codice del progetto, clang-format
  pulito. In `release` restano due `-Wnull-dereference` da un header di Geant4
  (`G4NavigationHistory.icc`), presenti anche prima di D21, e quattro `-Warray-bounds` da
  nlohmann/GCC 13 in `sim/tests/integration_reader.cpp`, presenti anche con il file di
  `83b959d`.

**Prossimi passi, nell'ordine:**

1. Col via: tre run da 10⁷ e tre di sensibilità a kB, come in `validation/2_predictions.md` v2,
   uno alla volta con `scripts/run_capped.sh -m 5G`, ciascuno con la verifica c1.
2. Risultati `validation/2_results.md` (con i piloti, l'oracolo della v1 fallito e D22) e
   v0.4.0.

## Punti aperti

- Dimensioni dello scintillatore: cubo di 10 cm provvisorio nei config; definitive prima della
  fase 3.
- SSD esterno montato su `/mnt/external_ssd`, ma non in modo permanente: prima di ogni run
  `scripts/data_root.sh check` (procedura in `architecture/data_layout.md`).
- **Fase 2:** sopra 19.5 MeV lo stato finale elastico è in parte CHIPS (D16): il bordo della
  sorgente va evitato o trattato come sistematica.
- **Fase 2:** la verità di Geant4 ha un E_p più alto della cinematica fisica di circa 5.4·10⁻⁴
  in relativo, circa 7.6 keV a E_p = 14 MeV. Ne seguirebbe un bias relativo dello stesso ordine
  su E_n ricostruito. Va deciso prima dell'analisi se correggerlo, se modellarlo o se riportarlo
  come sistematica (D14).
- **Fase 2:** il doppio conteggio del C-13 nel modello HP (circa +0.5% sul non elastico del
  carbonio) va tenuto fra le sistematiche del fondo da carbonio (D15).
- Confronto n–p indipendente da JEFF-3.3/ENDF-B-VIII.0 (D13) e verifica di H-2 contro una
  valutazione indipendente (D17): non fatti in 1b.
- **Fase 2:** a fine percorso `ionIoni` deposita più del K iniziale dello step. Riguarda 56
  rinculi C12 su 486 e 4 Be9, con K fra 0.11 e 0.40 MeV; l'eccesso arriva al 15%, per un totale
  di 0.74 MeV su 492 MeV di rinculi C12 (D14). La causa non è verificata; non dipende dalla
  funzione di passo del protone (D21, P4). Pesa su `det_Edep` più che sulla luce.
- **Fase 2:** kB è misurato fino a 500 MeV/cm; la fine del range dei protoni va oltre, quindi lì
  la luce è un'estrapolazione non coperta dall'intervallo di D18. Un fit di Birks sui dati EJ-200
  di Laplace 2020 la sonderebbe.
- **Fase 2:** le verifiche dei depositi devono tollerare `det_Edep` − E_n fino a circa
  6·10⁻⁴ · E_n senza un γ di cattura (D14).
- **Fase 2:** `kBytesPerEvent` (256 B) resta: i piloti v3 misurano 114–115 B/evento; il valore
  va registrato in `2_results.md`.
- Sull'SSD aspettano l'OK per la cancellazione: `sim/scratch_trace_246185` e i piloti v2
  `sim/pilot2_mono_14MeV_{1e5_1t,4e5_8t}` con le loro cartelle in `validation/` (schema 2, le
  verifiche dei depositi non li leggono più).
