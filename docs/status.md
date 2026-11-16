# Stato

**v0.1.0 — 2026-11-16 — Fase 0 (infrastruttura) completata.**

## Fatto

- Superbuild CMake con preset `debug`, `release`, `ci`; C++17; warning severi; `-ffast-math`
  vietato (`#error`); ASan/UBSan e `-Werror` nel preset `ci`.
- `common/`: `units.hpp` (MeV, mm, ns, rad = 1, verificate contro CLHEP con test a 0–1 ULP),
  `event_schema.hpp` (bozza v0: 3 colonne meta, 16 di verità), `data_paths.hpp` (radice dei dati
  con controlli di disco di sistema e spazio).
- `sim/riptide_sim --version`: stampa versione del progetto, di Geant4 (11.4.1 MT) e dello schema.
- `analysis/stats`: scheletro con controllo della semantica IEEE.
- Test: `debug` 14/14, `release` 14/14, `ci` 13/13 (senza Geant4), zero warning; clang-format
  pulito.

## Prossimo passo

Fase 1 (simulazione minima con verità MC), vedi `roadmap.md`.

## Punti aperti

- Dimensioni definitive dello scintillatore (prima della fase 3).
- SSD esterno non montato su `/mnt/external_ssd`: i run veri richiedono il mount o
  `RIPTIDE_DATA_DIR` su un altro disco.
- Comportamento di `G4NeutronHP` a 19.5–20 MeV (bordo del dominio).
- kB di Birks per BC-408: da prendere da letteratura citata (fase 2).
