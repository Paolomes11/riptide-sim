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

## Prossimo passo

Fase 1b (verifiche fisiche non circolari e prestazioni), vedi `roadmap.md`. Per ogni verifica la
previsione si scrive prima del run; i run girano sull'SSD con `scripts/run_capped.sh`.

## Punti aperti

- Dimensioni dello scintillatore: cubo di 10 cm provvisorio nei config; definitive prima della
  fase 3.
- Byte per evento: stima provvisoria (`kProvisionalBytesPerEvent` = 1 KiB) da misurare in 1b.
- SSD esterno montato su `/mnt/external_ssd`, ma non in modo permanente: prima di ogni run
  `scripts/data_root.sh check` (procedura in `architecture/data_layout.md`).
- Comportamento di `G4NeutronHP` a 19.5–20 MeV (bordo del dominio).
- kB di Birks per BC-408: da prendere da letteratura citata (fase 2).
