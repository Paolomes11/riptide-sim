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

## Collegare e preparare l'SSD

Prima di ogni run che scrive dati si lancia `scripts/data_root.sh`. Lo script ripete i controlli
di `requireWritableRoot`: se fallisce lui, si fermerebbe anche la simulazione. Non usa mai
`sudo`, non monta niente e non crea la radice; quando manca qualcosa stampa i comandi da lanciare
a mano.

1. **Controllo**, in sola lettura:

   ```bash
   scripts/data_root.sh check [--data-dir DIR] [--need GIB]
   ```

   La radice si risolve come sopra. Lo script verifica, nell'ordine:
   - che il mount point dell'SSD sia montato (`$RIPTIDE_SSD_MOUNT`, default `/mnt/external_ssd`),
     se la radice sta sotto di esso;
   - che la radice esista e sia scrivibile;
   - che non stia sul dispositivo di `/`;
   - che ci siano almeno `GIB` + 2 GiB liberi.
2. **Mount**, se il controllo lo chiede. È un'operazione manuale, con l'UUID preso da
   `lsblk -o NAME,SIZE,FSTYPE,UUID`:

   ```bash
   sudo mount UUID=<UUID> /mnt/external_ssd
   ls -la /mnt/external_ssd                       # verificare che sia il disco giusto
   sudo mkdir -p /mnt/external_ssd/riptide-sim    # solo la prima volta
   sudo chown "$USER":"$USER" /mnt/external_ssd/riptide-sim
   ```

   Se si definisce `RIPTIDE_SSD_UUID`, lo script lo inserisce già nel comando suggerito.
   Facoltativo, per un mount permanente: aggiungere a `/etc/fstab` la riga
   `UUID=<UUID> /mnt/external_ssd ext4 defaults,nofail 0 2`. Grazie a `nofail` il boot non si
   blocca se il disco manca.
3. **Preparazione**. Ripete il controllo e poi crea `sim/` e `analysis/` dentro la radice;
   lanciarlo più volte non cambia nulla:

   ```bash
   scripts/data_root.sh prepare
   ```

Codici di uscita:

| Codice | Significato |
|---|---|
| 0 | tutto a posto |
| 2 | uso errato, oppure percorso relativo |
| 3 | radice mancante o SSD non montato |
| 4 | radice non scrivibile |
| 5 | radice sul disco di sistema |
| 6 | spazio insufficiente |

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
