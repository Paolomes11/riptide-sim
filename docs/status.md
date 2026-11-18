# Stato

**v0.2.0 — 2026-11-17 — Fase 1a (simulazione minima con verità MC) completata.**

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

## In corso: fase 1b (branch `feat/sim-checks`, aggiornato al 2026-11-18)

Verifiche fisiche non circolari e prestazioni, vedi `roadmap.md`. Fatti i passi 0–5 di 9:

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
- **Test:** `debug` 134/134, `ci` 124/124, zero warning, clang-format pulito.

Prossimi passi:

6. previsioni scritte prima dei run (`docs/validation/1b_predictions.md`), comprese le soglie;
7. run, solo col via esplicito, sotto `scripts/run_capped.sh`;
8. byte/evento misurato;
9. risultati e chiusura con v0.3.0.

## Punti aperti

- Dimensioni dello scintillatore: cubo di 10 cm provvisorio nei config; definitive prima della
  fase 3.
- Byte per evento: stima provvisoria (`kProvisionalBytesPerEvent` = 1 KiB) da misurare in 1b.
- SSD esterno montato su `/mnt/external_ssd`, ma non in modo permanente: prima di ogni run
  `scripts/data_root.sh check` (procedura in `architecture/data_layout.md`).
- Comportamento di `G4NeutronHP` a 19.5–20 MeV (bordo del dominio).
- Frazioni della prima interazione, cattura (verificato in parte il 2026-11-18 sui sorgenti di
  Geant4 11.4.1 e su G4NDL 4.7.1):
  - in `G4HadronPhysicsQGSP_BIC_HP::Neutron()` l'inelastico (`neutronInelastic`, modello HP fino
    a 20 MeV) e la cattura (`nCaptureHP`, dataset `G4NeutronHPCaptureData`) sono processi
    distinti;
  - `G4NDL4.7.1/Inelastic/CrossSection` non ha un file `1_1_Hydrogen`, che invece c'è in
    `Capture/CrossSection`: per H-1 sotto 20 MeV le prime interazioni possibili sono solo
    l'elastico e la cattura;
  - resta da verificare numericamente, prima di scrivere le previsioni (passo 6 della 1b), che
    la sezione d'urto inelastica HP del carbonio (`Inelastic/CrossSection/6_nat_Carbon`) escluda
    la cattura, confrontandola con JEFF-3.3 C-0.
- Residui di E_p nello smoke test (248 eventi n–p a 14 MeV), cioè E_p simulata meno
  E_p(E_n, θ_p): media +3.7 keV, rms 4.3 keV, massimo 7.4 keV, tutti positivi. Il residuo
  relativo è costante, ΔE_p/E_p = (5.37–5.49)·10⁻⁴ su tutti gli eventi: è un fattore di scala
  deterministico, non rumore né moto termico. Il valore è vicino a m_e/m_p ≈ 5.45·10⁻⁴.
  Un'ipotesi da verificare è che la massa del bersaglio sia quella atomica dell'idrogeno e non
  quella del protone. Va capito e previsto al passo 6, prima dei run.
- kB di Birks per BC-408: da prendere da letteratura citata (fase 2).
