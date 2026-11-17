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
2. **Mount**, se il controllo lo chiede. Lo fa `scripts/mount_ssd.sh`, lanciato dall'utente
   perché chiede la password di `sudo`:

   ```bash
   scripts/mount_ssd.sh
   ```

   L'UUID non sta nel repository. Lo script lo legge da `RIPTIDE_SSD_UUID` oppure dal file non
   tracciato `~/.config/riptide-sim/ssd_uuid` (una riga; `$XDG_CONFIG_HOME` se definita), da
   creare una volta sola con l'UUID preso da `lsblk -o NAME,SIZE,FSTYPE,UUID`. Lo script:
   - se il mount point è già montato con quell'UUID non fa nulla, quindi si può lanciare a ogni
     sessione;
   - si ferma se sul mount point c'è un altro disco, se l'SSD è montato altrove, se il
     filesystem non è ext4 o se il mount point non montato contiene file (sarebbero dati finiti
     sul disco di sistema);
   - monta con `sudo mount UUID=<UUID>`, unico comando privilegiato, e poi lancia
     `data_root.sh prepare`;
   - non tocca `/etc/fstab` e non crea la radice dei dati.

   Codici di uscita: 0 a posto, 2 UUID mancante o non valido, 3 SSD non collegato, 7 stato non
   sicuro, 8 mount fallito, anche quando manca un terminale per la password di `sudo` (per
   esempio dal prompt `!` di una sessione con l'assistente: va lanciato da un terminale
   normale); gli altri vengono da `data_root.sh prepare`.

   A mano, i comandi equivalenti sono:

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
│   └── run.json              # provenienza del run e marcatore di completamento
└── analysis/<run_id>/<analysis_id>/
    ├── results.root | results.csv
    └── analysis.json         # parametri, semi dello smearing, versione dello schema letto
```

Nel repository, ignorata da git: `outputs/figures/<analysis_id>/` (figure da risultati ridotti).

Gli identificativi `run_id` e `analysis_id` contengono solo `[A-Za-z0-9._-]`, non iniziano con
`.` o `-` e sono lunghi al più 128 caratteri (`isValidId`): niente attraversamento di cartelle.
I `run_id` `ctest_it_1t` e `ctest_it_2t` sono riservati al test d'integrazione
(`sim/tests/run_integration.sh`), che cancella e ricrea le loro cartelle a ogni esecuzione.

## File degli eventi

Il merging delle ntuple di Geant4 è spento: ogni thread worker scrive il proprio file.
`riptide_sim` apre `events.root` e Geant4 aggiunge `_t<k>`, dove `k` è l'id del thread da 0 a
N−1, quindi un run a N thread produce `events_t0.root` … `events_t<N−1>.root`. Il master non
crea l'analysis manager e non scrive nessun file. Gli eventi sono distribuiti fra i file in modo
non deterministico, perciò un evento si identifica dal suo `event_id` e mai dalla posizione nel
file. Per leggere un run vanno letti i file elencati in `event_files` di `run.json`.

La ntupla di Geant4 non ha attributi, quindi non può portare la versione dello schema: questa sta
in `run.json`.

## run.json

`run.json` raccoglie la provenienza del run (`sim/config/.../run_record.hpp`) ed è anche il suo
**marcatore di completamento**. `riptide_sim` lo scrive per ultimo, solo dopo che `BeamOn` è
terminato: prima scrive `run.json.tmp` e poi lo rinomina, così un lettore non vede mai un file
scritto a metà. Se `run.json` esiste già, il run non lo sovrascrive.

Il lettore:
- rifiuta una cartella di run senza `run.json`, perché è un run interrotto;
- rifiuta un `schema_version` o un `run_record_version` che non conosce.

Chiavi (versione 1 del formato):

| chiave | tipo | contenuto |
|---|---|---|
| `run_record_version` | int | versione del formato di `run.json` |
| `run_id` | string | nome della cartella del run |
| `project_version` | string | versione del progetto (CMake) |
| `git_describe` | string | `git describe --tags --always --dirty` dei sorgenti alla build; `-dirty` indica modifiche non committate, `unknown` che git non è disponibile |
| `geant4_version` | string | `major.minor.patch` di Geant4 |
| `schema_version` | int | versione dello schema dei file degli eventi |
| `config_version` | int | versione del formato del config |
| `seed` | uint | seme del run, da cui derivano i semi per evento |
| `events` | uint | eventi richiesti e simulati |
| `threads` | int | thread worker |
| `started_utc` | string | inizio del run, ISO 8601 in UTC (`2026-09-07T12:34:56Z`) |
| `elapsed_s` | double | tempo reale da `Initialize` alla fine di `BeamOn` [s] |
| `event_files` | array di string | `events_t0.root` … `events_t<N−1>.root` |
| `config` | object | config effettivo; rileggibile con `parseRunConfig` |

`git_describe` viene rigenerato a ogni build (`cmake/RiptideGitDescribe.cmake`) e descrive i
sorgenti al momento della build, non al momento del run.

## Schema degli eventi (v1)

Albero `events`, 19 colonne in quest'ordine. Ci sono solo due tipi, `Int32` (`Int_t`) e `Double`
(`Double_t`), perché il writer di Geant4 non ha colonne intere a 64 bit.

| Colonne | Tipo | Contenuto |
|---|---|---|
| `event_id`, `rng_seed_0`, `rng_seed_1` | Int32 | numero dell'evento e i due semi che lo riproducono da solo (il seme del run è limitato a 31 bit) |
| `truth_En` | Double | energia cinetica del neutrone incidente [MeV] |
| `truth_omega_n_{x,y,z}` | Double | direzione del neutrone incidente |
| `truth_r0_{x,y,z}`, `truth_t0` | Double | vertice [mm] e tempo [ns] della prima interazione nello scintillatore |
| `truth_first_interaction` | Int32 | codice `FirstInteraction`: 0 nessuna, 1 n–p elastico, 2 n–C elastico, 3 n–C inelastico, 4 altro |
| `truth_n_np_scatters` | Int32 | numero di scattering elastici su H del neutrone primario |
| `truth_Ep`, `truth_p_dir_{x,y,z}`, `truth_theta_p` | Double | protone di rinculo del **primo** n–p: energia [MeV], direzione, angolo [rad] rispetto a Ω_n del fascio |
| `truth_valid_np` | Int32 | 1 se le colonne del protone sono riempite |

Unità: MeV, mm, ns, rad. Nomi, unità e descrizioni stanno in `event_schema.hpp`, che è l'unica
fonte: la riga scritta dalla sim (`sim/truth/.../event_row.hpp`) è verificata colonna per colonna
contro lo schema da un test.

Il primo n–p può venire dopo un'altra interazione, per esempio un n–C elastico. In quel caso il
neutrone non viaggia più lungo Ω_n, e θ_p può superare π/2. Per il campione pulito del singolo
scattering si seleziona quindi `truth_first_interaction == 1`. Le colonne `det_` (fase 2) e `reco_`
(analisi) arriveranno con le rispettive fasi.
