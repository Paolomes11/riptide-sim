# riptide-sim

Simulazione Geant4 e ricostruzione statistica in C++ del neutrone incidente a partire dal
singolo scattering n–p in uno scintillatore plastico (RIPTIDE), dominio 10–20 MeV.

Successore del toy model Python
[RIPTIDE_toy_model](https://github.com/Paolomes11/RIPTIDE_toy_model): lì il generatore degli
eventi e il modello di ricostruzione coincidevano, quindi la validazione era in parte circolare.
Qui gli eventi vengono da Geant4 (`G4NeutronHP`, dati valutati) e il motore statistico è
riscritto in C++. L'ottica di lettura è fuori scope.

- **Simulazione** (`sim/`): genera eventi e li scrive su file (TTree, un file per thread).
- **Analisi** (`analysis/`): legge i file, ricostruisce (Casi A, B, C) e valida (bias,
  risoluzione, pull, coverage, contrazione 1/√N, robustezza al prior).
- Comunicano solo tramite i file e lo schema in `common/`.

## Requisiti

- CMake ≥ 3.25, Ninja, compilatore C++17 (g++ ≥ 11 o clang ≥ 14), Catch2 3.
- Per la simulazione: Geant4 11.4 (MT) con i dati G4NDL, trovabile da CMake (`geant4-config`
  nel `PATH` o `CMAKE_PREFIX_PATH`).
- Per il lettore dell'analisi (`analysis/io/root`): ROOT ≥ 6.36 compilato in C++17 con
  `dataframe`, trovabile da CMake allo stesso modo (`root-config` nel `PATH` o il prefisso
  d'installazione in `CMAKE_PREFIX_PATH`). Si spegne con `-DRIPTIDE_BUILD_ROOT_IO=OFF`.

## Build e test

```bash
cmake --preset release && cmake --build --preset release && ctest --preset release
./build/release/sim/riptide_sim --version

# senza Geant4/ROOT (come in CI): common, sim/{config,truth,validation}, analysis/stats e
# analysis/io/core, con sanitizer e -Werror
cmake --preset ci && cmake --build --preset ci && ctest --preset ci
```

## Uso

```bash
scripts/data_root.sh check --need 1
scripts/run_capped.sh -m 5G -- ./build/release/sim/riptide_sim \
    --config config/mono_14MeV.json --events 1000 --threads 2 --run-id mono14_001
```

| Opzione | Significato |
|---|---|
| `--config <file>` | configurazione del run (JSON, vedi `config/`), obbligatoria |
| `--events <n>` | numero di eventi, obbligatorio, > 0 |
| `--run-id <id>` | nome della cartella sotto `<radice>/sim/`, obbligatorio; non deve esistere |
| `--data-dir <dir>` | radice dei dati (default: `$RIPTIDE_DATA_DIR`, poi `/mnt/external_ssd/riptide-sim`) |
| `--threads <n>` | thread worker (default 1) |
| `--seed <s>` | seme del run, 0..2147483647 (default 20260907) |
| `--allow-system-disk` | permette di scrivere sul disco di sistema (solo test brevi) |
| `--version` | stampa le versioni ed esce |

Il run scrive `<radice>/sim/<run_id>/events_t<k>.root` (un file per thread, TTree `events`) e,
alla fine, `run.json` con config, seme e versioni. Lo stesso seme dà gli stessi eventi con
qualsiasi numero di thread. Formato in `docs/architecture/data_layout.md`.

Il test d'integrazione (`ctest --preset debug -R sim_integration`) scrive due run brevi
sull'SSD; se l'SSD non è disponibile viene saltato.

## Dati

I dati non stanno nel repository né sul disco di sistema. Radice: `--data-dir` >
`$RIPTIDE_DATA_DIR` > `/mnt/external_ssd/riptide-sim`; deve esistere, essere su un disco diverso
da `/` e avere spazio sufficiente. Dettagli in `docs/architecture/data_layout.md`.

Prima di ogni run, controllo e preparazione dell'SSD. `data_root.sh` non usa `sudo`: se l'SSD
non è montato, stampa i comandi per montarlo. `mount_ssd.sh` monta l'SSD per UUID (chiede la
password di `sudo`) e poi prepara la radice.

```bash
scripts/mount_ssd.sh                   # UUID da RIPTIDE_SSD_UUID o ~/.config/riptide-sim/ssd_uuid
scripts/data_root.sh check --need 10   # 10 GiB + 2 di margine
scripts/data_root.sh prepare           # crea sim/ e analysis/ nella radice
```

Run lunghi con tetto di memoria:

```bash
scripts/run_capped.sh -m 5G -- ./build/release/sim/riptide_sim ...
```

## Documentazione

- `docs/status.md` — stato attuale
- `docs/roadmap.md` — fasi 0–10
- `docs/toy_model_findings.md` — risultati ereditati dal toy model
- `docs/decisions.md` — registro delle decisioni
- `docs/architecture/` — architettura e layout dei dati

## Stato

v0.2.0: fase 1a (simulazione minima con verità MC) completata; fase 1b (verifiche fisiche e
prestazioni) in corso su `feat/sim-checks`, passi 0–3 di 9. Dettagli in `docs/status.md`.

## Licenza

Apache-2.0, vedi `LICENSE` e `NOTICE`.
