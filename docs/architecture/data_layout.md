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
├── analysis/<run_id>/<analysis_id>/
│   ├── results.root | results.csv
│   └── analysis.json         # parametri, semi dello smearing, versione dello schema letto
├── validation/<run_id>/<check_id>/
│   ├── spectrum_u.csv | depth_z.csv | fractions.csv | legendre.csv | residuals_Ep.csv
│   ├── model_mix.csv
│   └── checks.json           # esiti numerici delle verifiche e marcatore di completamento
└── reference/                # valutazioni nucleari di riferimento (D13), sola lettura
    ├── endfb80/              # file ENDF-6 estratti; download/ con gli originali
    └── jeff33/
```

Nel repository, ignorata da git: `outputs/figures/<analysis_id>/` (figure da risultati ridotti).

Gli identificativi `run_id`, `analysis_id` e `check_id` contengono solo `[A-Za-z0-9._-]`, non iniziano con
`.` o `-` e sono lunghi al più 128 caratteri (`isValidId`): niente attraversamento di cartelle.
I `run_id` `ctest_it_1t` e `ctest_it_2t` sono riservati al test d'integrazione
(`sim/tests/run_integration.sh`), che cancella e ricrea le loro cartelle a ogni esecuzione.
Il `check_id` `ctest_smoke` è riservato allo smoke test di `riptide_sim_check`
(`sim/tests/run_check_smoke.sh`), che ricrea solo `validation/ctest_it_1t/ctest_smoke`. Lo
smoke test controlla `checks_version` = 5, i file JEFF-3.3 letti (H-1, C-0, C-13, D15) e lo
sha256 del primo, l'AWR di H-1 (0.999167, D14), il residuo di E_p all'arrotondamento, lo scarto
dalla cinematica fisica in [5.3, 5.5]·10⁻⁴, i contatori di coerenza a zero e, a 14 MeV, nessun
evento n–p con la cinematica non HP (D16).

## Valutazioni di riferimento

File in formato ENDF-6 usati dalle verifiche della fase 1b (D13). Si scaricano una volta, si
rendono di sola lettura, e il programma di verifica ne riporta lo sha256 in `checks.json`.

ENDF/B-VIII.0, da NNDC:

| File | Origine | sha256 |
|---|---|---|
| `download/ENDF-B-VIII.0_neutrons.zip` | https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_neutrons.zip (310158456 byte) | `a1cced11710eb36690c5ad7b64d414465c5849f5428904049c864e9beca28894` |
| `download/endf-manual-viii.0.pdf` | https://www.nndc.bnl.gov/endf-b8.0/endf-manual-viii.0.pdf | `77a0fee413c3b1d5d74a161ed9fe7f77bbcbc58a654304851b7b2b400183d022` |
| `n-001_H_001.endf` | dallo zip (MAT 125) | `e7b76c20c1d61f2cb2621de716d324047fa1f9583bef568f643f8d98094260fb` |
| `n-006_C_012.endf` | dallo zip (MAT 625) | `daa7901e0067f9ac6d9b2392514ba00e71e37d00c8a8ecbbde0620237815f981` |
| `n-006_C_013.endf` | dallo zip (MAT 628) | `abaa8cc34398e6ada0a8c451efd5eef1019a37c7a0a708a4280af71ab77621de` |

JEFF-3.3, da IAEA-NDS (https://www-nds.iaea.org/public/download-endf/JEFF-3.3/n/), scaricato a
mano dal browser perché il sito rifiuta gli accessi automatici. JEFF-3.3 non ha un file per C-12,
ma solo quello del carbonio naturale (C-0), che è anche quello usato da G4NDL (D13). Le
verifiche leggono H-1, C-0 e C-13 (D15).

| File | Origine | sha256 |
|---|---|---|
| `download/n_0125_1-H-1.zip` | IAEA-NDS (12665 byte) | `3ce3c497c56676d14affc680d787db276e3804ba820033a25312382925ce3b9c` |
| `download/n_0600_6-C-0.zip` | IAEA-NDS (171589 byte) | `e2f727f80492b6e59dc618ab5edcdde3483afb2e93fc5fbf6683259ce315ff1c` |
| `download/n_0628_6-C-13.zip` | IAEA-NDS (469364 byte) | `cb20c0afadee4bc7667f015d8e807e8e6b45c5039099663bdee121158d3e6e93` |
| `n_0125_1-H-1.dat` | dallo zip (MAT 125) | `d7dc8efb7bb1284013e7ae233d3d9c7fd1c61dad6fd97bcd29542fdb531ceca8` |
| `n_0600_6-C-0.dat` | dallo zip (MAT 600) | `30ca87a83148abe20653c892bc78ffe762b0cf6757e4ad597197140cfb0c5b0c` |
| `n_0628_6-C-13.dat` | dallo zip (MAT 628) | `c51b1895594c2ee354c2c0c353bf90b7e229fac1154464c88971f200a642d244` |

In ogni file l'intestazione MF1 riporta "JEFF-3.3" e ZA/MAT corrispondono al nome.

H-2 non ha un file valutato fra questi: entrambe le librerie lo prendono dalle sezioni d'urto
di G4NDL 4.7.1 lette da HP (D17), sotto `$G4NEUTRONHPDATA` e non nella radice dei dati.
Elastico e cattura hanno una discontinuità a 14.25 MeV (energia ripetuta), che il totale
riproduce.

| File | sha256 |
|---|---|
| `Elastic/CrossSection/1_2_Hydrogen.z` | `ad1e901e2545bd67792d78a2bfb013daca24d47d006d95690e882ef692fde004` |
| `Inelastic/CrossSection/1_2_Hydrogen.z` | `ebe306167a4dca0a8c5aca39f7b17deeb79c396b286bd0625fb134175290bd1d` |
| `Capture/CrossSection/1_2_Hydrogen.z` | `e40c5a72fedd83643107c3de866be1538b5048c0d84ddca366e5bdd48422f786` |

## Verifiche fisiche (`riptide_sim_check`)

`riptide_sim_check --run-id <run> --check-id <id>` legge un run completo una sola volta e scrive
in `validation/<run_id>/<check_id>/`. Se la cartella esiste già si ferma: nessuna sovrascrittura.
Non dà verdetti: le soglie di esito vengono dalle previsioni scritte prima dei run.

Le attese si sommano evento per evento alla E_n di ciascun evento e sono calcolate per due
librerie:

- `jeff33`: H-1, C-0 con la densità di C-12 e C-13 con la sua, come il modello HP compone il
  carbonio (D15);
- `endfb80`: H-1, C-12 e C-13.

Le densità atomiche sono quelle del materiale NIST `G4_PLASTIC_SC_VINYLTOLUENE` costruito da
Geant4, quindi del mezzo effettivamente simulato. Il deuterio non ha un file di riferimento:
nella verità è classificato n–p (Z = 1), ma le attese usano solo H-1, e `checks.json` riporta
n_H2/n_H1 (circa 1.15·10⁻⁴) come sistematica nota.

I CSV, a precisione piena (`max_digits10`):

- `spectrum_u.csv`, `depth_z.csv`: `low,high,observed,expected_<lib>…`, una riga per bin
  originale (u = (1 − μ)/2, con μ ricavato da E_p con la cinematica HP di D14; profondità del primo vertice lungo il fascio, in mm);
- `fractions.csv`: `code,category,observed,expected_<lib>…` per i 5 codici di
  `FirstInteraction`;
- `legendre.csv`: `l,value,error,expected_<lib>,pull_<lib>…` per l ≥ 1 (solo intestazione se la
  verifica è saltata);
- `residuals_Ep.csv`: `mu,Ep_residual_MeV,Ep_physical_offset_MeV`, una riga per evento n–p
  classificato HP (D16).
  Le colonne sono:
  - `mu`, il coseno CM da θ_p;
  - `Ep_residual_MeV`, l'E_p simulata meno la cinematica HP (arrotondamento);
  - `Ep_physical_offset_MeV`, l'E_p simulata meno la cinematica fisica n–p (convenzione HP, D14).
- `model_mix.csv`: `En_low_MeV,En_high_MeV,classified,observed_non_hp,expected_non_hp,variance`,
  una riga per bin di E_n nella finestra di sovrapposizione HP/CHIPS (D16) intersecata con la
  sorgente; solo intestazione se l'intersezione è vuota.

`checks.json` (versione `checks_version` = 5) si scrive per ultimo, con tmp + rename. Chiavi:

| Chiave | Contenuto |
|---|---|
| `checks_version`, `run_id`, `check_id` | identificazione |
| `versions` | `riptide`, `git`, `geant4`, `schema`, `run_record` (del run letto) |
| `reference_files` | per file: `library`, `nuclide`, `path`, `sha256` |
| `hp_elastic_H1` | `file`, `sha256`, `rep_flag`, `awr`, `frame_flag`, `target_mass_MeV` dell'intestazione G4NDL di H-1 (D14) |
| `hp_deuterium` | per file: `file`, `sha256` delle sezioni d'urto G4NDL di H-2 (elastica, anelastica, cattura; D17) |
| `medium` | `material`, `density_g_cm3`, `n_H1_per_mm3`, `n_H2_per_mm3`, `n_C12_per_mm3`, `n_C13_per_mm3`, `H2_over_H1` |
| `beam` | `L_mm` (corda del fascio nel blocco), `En_min_MeV`, `En_max_MeV`, `bins_u`, `bins_z` |
| `results` | `events`, `np_events`, `np_hp_events`, `spectrum_u`, `depth_z`, `fractions`, `legendre`, `Ep_residuals`, `Ep_physical_offset`, `model_mix`, `counters` |

In `results`, ogni verifica riporta per libreria `skipped` (`null` oppure il motivo). Se la
verifica è stata fatta, riporta anche χ², `ndf`, `p_value` e i pull; per gli istogrammi pure i
bin fusi (attesa < 5) da cui viene il χ². `counters.must_be_zero` raccoglie i contatori che
solo un bug può rendere diversi da zero, con il loro `total`. Fra questi ci sono anche
u > 1 + 10⁻⁹ e θ_p > π/2, impossibili con il bersaglio a riposo (D14). `counters.expected_rare`
raccoglie la fisica nota, che viene riportata e non tagliata: n–d senza protone ed eventi n–p
compatibili con entrambe le cinematiche (`np_model_ambiguous`).
`Ep_residuals` riporta `samples`, `mean`, `rms` e `max_abs` in MeV. `Ep_physical_offset` ha le
stesse chiavi e in più `relative_mean`, `relative_min` e `relative_max` dello scarto relativo,
con `rounding_free_relative_min` e `rounding_free_relative_max`: gli stessi estremi al netto
dell'arrotondamento, cioè min(rel + R/E_p,fis) e max(rel − R/E_p,fis) con R =
`Ep_residuals.max_abs`. Sono questi che si confrontano con il range previsto, perché a E_p di pochi
eV un residuo di arrotondamento di 10⁻¹⁵ MeV sposta da solo lo scarto relativo di ~10⁻¹⁰.

Ogni evento n–p con protone valido (`np_events`) è classificato dalla cinematica (D16), con
tolleranza 10⁻¹¹·E_n su E_p:
- HP, se E_p coincide con la cinematica HP di D14: entra in `np_hp_events`, spettro u, Legendre,
  `Ep_residuals` e `Ep_physical_offset`;
- non HP (CHIPS), se coincide con la cinematica fisica n–p: entra solo in `model_mix`;
- ambiguo, se coincide con entrambe: trattato come HP e contato in `np_model_ambiguous`;
- nessuna delle due: `must_be_zero.np_kinematics_unmatched`.
Un evento non HP con E_n fuori dalla finestra 19.5–20 MeV va in
`must_be_zero.np_non_hp_outside_window`.

`model_mix` riporta `window_MeV`, `classified`, `observed_non_hp`, `expected_non_hp` (Σ p(E_n),
con p la rampa lineare di D16), `variance` (Σ p(1 − p)), `skipped` e `pull` (`null` se saltato),
`bins`, `chi2_skipped` e, se fatto, `chi2`, `ndf`, `p_value` sui bin di E_n, più `Ep_residuals`
degli eventi non HP contro la cinematica fisica.

Dalla versione 5 (2026-11-19):
- le attese contengono H-2 dalle sezioni d'urto di G4NDL (D17): l'elastico va in n–p senza
  protone di rinculo, il break-up e la cattura in "other"; c'è la chiave `hp_deuterium`.
- `Ep_physical_offset` ha `rounding_free_relative_min` e `rounding_free_relative_max` (D17).

Dalla versione 4 (2026-11-19):
- classificazione HP/non HP degli eventi n–p e chiavi `np_hp_events` e `model_mix`;
- `must_be_zero` ha `np_kinematics_unmatched` e `np_non_hp_outside_window`;
- `expected_rare` ha `np_model_ambiguous`;
- c'è `model_mix.csv`.

Dalla versione 2 (2026-11-18):
- c'è la chiave `hp_elastic_H1`;
- in `results` c'è `Ep_physical_offset`;
- i contatori u > 1 e θ_p > π/2 sono passati in `must_be_zero`;
- `residuals_Ep.csv` ha la terza colonna.

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
