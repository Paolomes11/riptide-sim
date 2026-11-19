# Stato

**v0.3.0 — 2026-11-19 — Fase 1b (verifiche fisiche e prestazioni) completata.**

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
- kB di Birks per BC-408: da prendere da letteratura citata (fase 2).
