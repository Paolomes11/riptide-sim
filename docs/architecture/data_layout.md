# Layout dei dati

I dati non stanno mai nel repository né sul disco di sistema.

## Radice

Risolta da `riptide::data::resolveDataRoot`, in quest'ordine:

1. opzione `--data-dir` della riga di comando;
2. variabile d'ambiente `RIPTIDE_DATA_DIR`;
3. default `/mnt/external_ssd/riptide-sim`.

Un valore vuoto conta come non impostato; un percorso relativo è un errore. Prima di scrivere,
`requireWritableRoot` verifica che la radice:

- esista e sia una cartella (non viene **mai** creata: un SSD non montato lascerebbe il mount
  point vuoto sul disco di sistema);
- sia scrivibile;
- non stia sullo stesso dispositivo di `/`;
- abbia spazio libero ≥ stima del run + 2 GiB (`estimateRunBytes(n_eventi, byte_per_evento)`).

## Struttura

```
$RIPTIDE_DATA_DIR/
├── sim/<run_id>/
│   ├── events_t<k>.root      # un file per thread, TTree "events"
│   └── run.json              # config, seme, git sha, versioni, n. eventi
└── analysis/<run_id>/<analysis_id>/
    ├── results.root | results.csv
    └── analysis.json         # parametri, semi dello smearing, versione dello schema letto
```

Nel repository, ignorata da git: `outputs/figures/<analysis_id>/` (figure da risultati ridotti).

Gli identificativi `run_id` e `analysis_id` contengono solo `[A-Za-z0-9._-]`, non iniziano con
`.` o `-` e sono lunghi al più 128 caratteri (`isValidId`): niente attraversamento di cartelle.

## Schema degli eventi (bozza v0)

Albero `events`. Meta: `event_id`, `rng_seed_0`, `rng_seed_1`. Verità (`truth_`): energia e
direzione del neutrone, vertice e tempo della prima interazione, tipo della prima interazione,
numero di scattering n–p, energia e direzione del protone di rinculo, θ_p, `truth_valid_np`.
Unità: MeV, mm, ns, rad. Elenco completo con unità e descrizioni in `event_schema.hpp`. Le
colonne `det_` (fase 2) e `reco_` (analisi) arriveranno con le rispettive fasi.
