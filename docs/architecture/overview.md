# Architettura

## Componenti

```
riptide-sim/
├── common/      header-only: unità, schema degli eventi, radice dei dati, versione
├── sim/         applicazione Geant4: genera eventi e li scrive su file
│   ├── config/  configurazione del run e run.json, senza Geant4
│   ├── truth/   classificazione della verità per evento, senza Geant4
│   ├── reader/  lettura di run.json e delle ntuple di un run (con Geant4), anche in streaming
│   ├── validation/  lettura ENDF e G4NDL, cinematica del rinculo HP, attese, statistiche e
│   │                report delle verifiche fisiche, riga di comando di riptide_sim_check;
│   │                senza Geant4
│   └── apps/    riptide_sim (simulazione), riptide_sim_check (verifiche fisiche di un run)
├── analysis/
│   ├── stats/   motore statistico, C++ puro (niente ROOT, niente Geant4)
│   ├── io/      (fase 4) lettura ROOT via RDataFrame → struct-of-arrays
│   └── apps/    (fase 5) riptide_reco, riptide_validate
├── config/      (fase 1) configurazioni di run JSON versionate
├── scripts/     wrapper shell; Python solo per figure da risultati ridotti
├── docs/
└── cmake/       opzioni del compilatore per target
```

Dipendenze, solo verso il basso:

```
common  ←  sim/config  ←  sim/truth  ←  sim
                           sim/truth  ←  sim/reader  (Geant4)  ←  test d'integrazione, app di verifica
common  ←  sim/config  ←  sim/validation  ←  riptide_sim_check (con sim/reader e Geant4)
common  ←  analysis/stats  ←  analysis/apps
common  ←  analysis/io     ←  analysis/apps
```

Note:

- `riptide_sim_check` usa Geant4 solo per costruire il materiale NIST dello scintillatore (nome
  in `sim/include/riptide/sim/detector_materials.hpp`, condiviso con la simulazione) e per il
  lettore ROOT; non crea un run manager.
- La conversione riga dell'ntupla → `TruthRecord` (`truthRecord`) sta in `sim/truth`, accanto
  alla conversione inversa `eventRow`.
- L'oracolo G4NDL ↔ JEFF-3.3 (`sim/tests/hp_data_reference_test.cpp`, D15) sta in `sim_tests`,
  che linka Geant4 e `sim/validation`: trova G4NDL con `G4FindDataDir` e lo legge con
  `G4ParticleHPManager::GetDataStream`, come il modello HP, e lo analizza con il lettore puro
  `sim/validation/hp_data.hpp`.
- `CliError` e il parsing dei numeri sono in `sim/config/cli_parse.hpp`, condivisi dalle due
  righe di comando.

Simulazione e analisi comunicano **solo** attraverso i file su disco e lo schema in
`common/include/riptide/common/event_schema.hpp`. Lo schema ha una versione (`kSchemaVersion`):
ogni cambio di colonne la incrementa e l'analisi rifiuta file con versione diversa.

## Build

- CMake ≥ 3.25 con preset (`CMakePresets.json`), generatore Ninja, cartelle `build/<preset>`.
- `debug` e `release`: build completa (richiede Geant4 11.4). `ci`: senza Geant4, cioè `common`,
  le librerie pure di `sim/` (`config`, `truth`, `validation`) e `analysis/stats`, con ASan/UBSan
  e `-Werror`.
- Opzioni del compilatore applicate per target tramite `riptide::warnings` e
  `riptide::options`: il codice di terze parti (Geant4, ROOT, Catch2) non è toccato. Gli header
  di Geant4 arrivano come include di sistema.
- Test con Catch2 3, scoperti da `catch_discover_tests`; codice di ritorno 4 = test saltato (es.
  ambiente senza il disco atteso).

## Efficienza (dal primo giorno)

- Simulazione MT, un file per thread, nessuna unione obbligatoria.
- Lo schema contiene solo tipi scalari (niente vettori per evento finché non servono): lettura
  colonnare veloce.
- Analisi su buffer struct-of-arrays contigui; nessun ciclo sugli eventi dentro le funzioni di
  griglia; parallelismo solo dopo aver profilato.
- Niente `-ffast-math` (vedi `decisions.md`, D4).

## Riproducibilità

- Ogni evento porta `event_id` e i due semi del suo RNG.
- Ogni run scrive `run.json`: configurazione, seme principale, git sha, versioni di Geant4, ROOT e
  schema, numero di eventi.
- Ogni analisi scrive `analysis.json` con i propri parametri e semi.

Il layout dei dati su disco è in `data_layout.md`.
