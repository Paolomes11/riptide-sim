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
│   ├── cutflow.csv           # eventi rimasti dopo ogni taglio (D30)
│   └── analysis.json         # parametri, semi dello smearing, massa HP e sua fonte, schema letto (D31)
├── validation/<run_id>/<check_id>/
│   ├── spectrum_u.csv | depth_z.csv | fractions.csv | legendre.csv | residuals_Ep.csv
│   ├── model_mix.csv
│   ├── deposit_spectra.csv | proton_exits.csv   # solo per run con schema ≥ 3
│   ├── pstar_comparison.csv  # solo per run con schema ≥ 3
│   └── checks.json           # esiti numerici delle verifiche e marcatore di completamento
└── reference/                # dati di riferimento, sola lettura
    ├── endfb80/              # file ENDF-6 estratti (D13); download/ con gli originali
    ├── jeff33/
    └── pstar/                # stopping power e range dei protoni nel PVT (D20)
```

Nel repository, ignorata da git: `outputs/figures/<analysis_id>/` (figure da risultati ridotti).

Gli identificativi `run_id`, `analysis_id` e `check_id` contengono solo `[A-Za-z0-9._-]`, non iniziano con
`.` o `-` e sono lunghi al più 128 caratteri (`isValidId`): niente attraversamento di cartelle.
I `run_id` `ctest_it_1t`, `ctest_it_2t` e `ctest_it_kB0` sono riservati al test d'integrazione
(`sim/tests/run_integration.sh`), che cancella e ricrea le loro cartelle a ogni esecuzione.
`ctest_it_kB0` ripete il run a 1 thread con kB = 0 (config generato al volo, fuori dal repo):
il reader controlla `det_Evis` = `det_Edep` evento per evento entro 10⁻⁹ MeV e `det_Edep`
identico, bit per bit, a quello di `ctest_it_1t`, perché kB entra solo nella luce visibile.
Il `check_id` `ctest_smoke` è riservato allo smoke test di `riptide_sim_check`
(`sim/tests/run_check_smoke.sh`), che ricrea solo `validation/ctest_it_1t/ctest_smoke` e
`validation/ctest_it_2t/ctest_smoke`. I due run contengono gli stessi eventi, divisi in uno e in
due file: i CSV devono essere identici byte per byte e `checks.json` identico a meno di `run_id`
(D23). Lo smoke test controlla `checks_version` = 10, la modalità del fascio, `geometric_acceptance` (raggi = eventi, colpi ≤ raggi, colpi = raggi nel pencil), i file JEFF-3.3 letti (H-1, C-0, C-13, D15) e lo
sha256 del primo, l'AWR di H-1 (0.999167, D14), il residuo di E_p all'arrotondamento, lo scarto
dalla cinematica fisica in [5.3, 5.5]·10⁻⁴, i contatori di coerenza a zero e, a 14 MeV, nessun
evento n–p con la cinematica non HP (D16). Controlla anche la sezione `deposits`: fatta, con kB
nominale, con le classi che sommano agli eventi e con `must_be_zero.total` = 0, e la sezione
`deposits.pstar`: fatta, con lo sha256 di `pstar/pstar_216.txt` uguale a quello di `sha256sum`,
con `samples` + `below_table` uguale ai protoni di `contained_balance_MeV` e con un bin per bin
dei depositi. Nessuna soglia fisica: quelle stanno nelle previsioni.

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

## Stopping power dei protoni (PSTAR)

Tabella PSTAR del NIST per il materiale 216, "PLASTIC SCINTILLATOR (vinyltoluene-based)" (D20),
scaricata il 2026-11-20 in `reference/pstar/`. La tabella in testo si ottiene con un POST a
https://physics.nist.gov/cgi-bin/Star/apdata.pl con `prog=PSTAR`, `matno=216`,
`ShowDefault=on`, `NumofEnergies=0`, `character=space` e tutte le colonne attive (`electronic`,
`nuclear`, `total`, `csda`, `project`, `detour`).

| File | Origine | sha256 |
|---|---|---|
| `pstar_216.txt` | apdata.pl, tabella in testo (9308 byte) | `eec195fe1257334039598db5d24e6040c91689200e1731c05f24301a0f3c081b` |
| `download/ap_table_216.html` | https://physics.nist.gov/cgi-bin/Star/ap_table.pl, stessa tabella in HTML (36078 byte) | `a01d863fc2d1cba1fffa72e2ccd010ad9e8302222d2cb06936b54488c47aa4e8` |
| `download/compos_ap216.html` | https://physics.nist.gov/cgi-bin/Star/compos.pl?ap216, composizione, densità e I (7014 byte) | `49c5ca730f42e3d77c8db3d6622432d763b961e91d030634d62f5a33b8321416` |

`pstar_216.txt` ha un'intestazione di 8 righe e 132 righe di dati con 7 colonne: E [MeV],
stopping elettronico, nucleare e totale [MeV cm²/g], range CSDA e proiettato [g/cm²], fattore
di detour. Il lettore converte in MeV/mm e mm con la densità 1.032 g/cm³.

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
  originale (u = (1 − μ)/2, con μ ricavato da E_p con la cinematica HP di D14; profondità s = z/L_i
  del primo vertice lungo la corda L_i dell'evento, adimensionale su [0, 1]);
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

- `deposit_spectra.csv`: `low,high,Edep,Evis`, eventi con `det_Edep` > 0 e `det_Evis` > 0 in
  bin uguali su [0, E_n,max]; i valori sopra E_n,max (catture, con il loro Q) vanno nei
  contatori di overflow di `checks.json`, non nell'ultimo bin.
- `proton_exits.csv`: `Ep_low,Ep_high,protons,exiting`, protoni di rinculo validi e quelli che
  escono dal blocco, per bin di `truth_Ep` con gli stessi bordi.
- `pstar_comparison.csv`: `Ep_low,Ep_high,protons,length_mean,R_csda_mean,length_residual_mean,
  length_residual_sd,Evis_mean,L_mean,Evis_residual_mean,Evis_residual_sd`, per bin di
  `truth_Ep` con gli stessi bordi: medie di `truth_p_length` [mm] e di R_CSDA(E_p), di
  `truth_p_Evis` [MeV] e di L(E_p; kB), e media e deviazione standard campionaria dei residui
  (osservato − atteso). Campo vuoto per le medie di un bin senza protoni e per le deviazioni
  con meno di due protoni.

Il checker legge le righe dei file per thread fuse in ordine crescente di `event_id`
(`mergeByEventId`, D23): ogni somma, e quindi ogni output, dipende solo dagli eventi e non da come
Geant4 li ha divisi fra i thread.

`checks.json` (versione `checks_version` = 10) si scrive per ultimo, con tmp + rename. Chiavi:

| Chiave | Contenuto |
|---|---|
| `checks_version`, `run_id`, `check_id` | identificazione |
| `versions` | `riptide`, `git`, `geant4`, `schema` del programma di verifica (`kSchemaVersion` compilato, non lo schema del run: quello è in `run.json`), `run_record` del run letto |
| `reference_files` | per file: `library`, `nuclide`, `path`, `sha256` |
| `hp_elastic_H1` | `file`, `sha256`, `rep_flag`, `awr`, `frame_flag`, `target_mass_MeV` dell'intestazione G4NDL di H-1 (D14) |
| `hp_deuterium` | per file: `file`, `sha256` delle sezioni d'urto G4NDL di H-2 (elastica, anelastica, cattura; D17) |
| `medium` | `material`, `density_g_cm3`, `n_H1_per_mm3`, `n_H2_per_mm3`, `n_C12_per_mm3`, `n_C13_per_mm3`, `H2_over_H1` |
| `beam` | `mode` (`pencil` o `far_field`), `En_min_MeV`, `En_max_MeV`, `bins_u`, `bins_z`; pencil: `L_mm` (corda del fascio nel blocco); far-field: `direction` (vettore oppure `"isotropic"`), `disk_radius_mm`, `disk_distance_mm`, `expected_acceptance`, `expected_mean_chord_mm` (D27) |
| `results` | `events`, `np_events`, `np_hp_events`, `geometric_acceptance`, `spectrum_u`, `depth_z`, `fractions`, `legendre`, `Ep_residuals`, `Ep_physical_offset`, `model_mix`, `counters`, `deposits` |

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

Ogni evento ha la sua corda L_i = `rayBoxChord(truth_n_start, truth_omega_n)` nel blocco (D27).
Nel fascio pencil è la corda del config, uguale per tutti, e una direzione diversa da quella del
fascio va in `must_be_zero.omega_not_beam` (l'evento resta verificato lungo il fascio). Nella
sorgente a campo lontano, invece:
- una direzione non unitaria, o diversa da quella fissa del config, va in `omega_not_beam`;
- un punto di partenza fuori dal disco ⟂ Ω_n (piano e raggio entro 10⁻⁶ mm) va in
  `must_be_zero.start_not_on_source`.
In entrambi i casi l'evento è escluso da ogni altra verifica. Un raggio che manca il blocco ha
prima interazione Unknown con certezza: entra in `fractions` con attesa 1 in Unknown, e una prima
interazione diversa va in `must_be_zero.interaction_off_chord`. Le attese delle frazioni e della
profondità usano L_i evento per evento: la profondità è binnata come frazione s = z/L_i, e l'attesa
di ogni evento è l'esponenziale troncato su [0, L_i] diviso negli stessi bin. Nel pencil s è z/L
e i numeri coincidono con la versione 9.

`geometric_acceptance` riporta `rays` (eventi con raggio verificato), `hits` (raggi con corda > 0),
`expected_fraction` (accettanza in forma chiusa: area proiettata media su πR²), `skipped` e
`pull` binomiale dei colpi, `chord_mm` (`mean`, `sd`, `min`, `max` delle corde dei colpi, media e
varianza di Welford, e `expected_mean`: la corda fissa, oppure 4V/S per la sorgente isotropa,
formula di Cauchy) e `mean_chord_skipped` / `mean_chord_pull` = (media − attesa)/(sd/√colpi). Nel
pencil entrambi i pull sono saltati (`null`): ogni raggio colpisce e la corda è una sola.

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

`deposits` riguarda le colonne v2 e v3 (fase 2). Per un run con schema 1 o 2 vale solo
`{"skipped": "schema version N is older than 3, the first these checks read"}`, e le altre
verifiche restano invariate. Lo schema 2 è escluso perché il lettore darebbe 0 a
`truth_p_n_nuclear_recoil` e metterebbe nel bilancio i protoni con un rinculo nucleare (D22). Altrimenti ha `skipped` = `null` e:
- `setup`: `half_size_mm`, `En_max_MeV`, `birks_kB_mm_per_MeV`, `bins` e le tolleranze
  `deposit_sum_tolerance_MeV` (10⁻⁹), `nuclear_stopping_slack_MeV` (0.01, D19),
  `surface_tolerance_mm` (10⁻⁶);
- `events`, `np_events`, `event_classes` (conteggi per nome della classe, D19),
  `events_secondary_n_hadronic`, `protons_hadronic`, `protons_nuclear_recoil` (protoni con
  `truth_p_n_nuclear_recoil` > 0, D22), `protons_exiting`;
- `contained_balance_MeV`: `samples`, `mean`, `min`, `max`, `positive`, `negative` di
  `truth_p_Edep` − `truth_Ep` per i protoni contenuti senza interazioni adroniche e senza
  rinculi nucleari (D22);
- `kB_zero_max_abs_difference_MeV`: max |`det_Evis` − `det_Edep`| su tutti gli eventi (con kB = 0
  è l'oracolo, altrimenti misura il quenching);
- `Edep_overflow`, `Evis_overflow`;
- `counters.must_be_zero`: depositi negativi, Evis > Edep (evento e protone), deposito del
  protone sopra il totale o sopra E_p oltre lo slack, protone valido senza traccia, deposito senza
  interazione del primario, classe scritta diversa da quella ricalcolata, prima interazione non
  contata, uscita con punto finale fuori dalla superficie, punto finale fuori dal blocco,
  Evis ≠ Edep con kB = 0, e `total`;
- `counters.expected_rare`: `p_Edep_above_Ep` (entro lo slack, G4NuclearStopping a fine range)
  e `contained_end_on_surface`;
- `pstar`: confronto con la tabella PSTAR (D20) degli stessi protoni di
  `contained_balance_MeV`. Per ciascuno l'attesa è calcolata alla sua `truth_Ep`, non al centro
  del bin: R_CSDA(E_p) interpolato log-log per `truth_p_length`, e L(E_p; kB) = ∫ dE / (1 + kB·S),
  con S totale della tabella e kB del run, per `truth_p_Evis`. Chiavi: `skipped` (`null`),
  `file`, `sha256`, `density_g_cm3` (quella della tabella; se la densità del mezzo Geant4 ne
  differisce oltre 10⁻¹² in relativo, il programma si ferma), `samples` (protoni confrontati),
  `below_table` (contenuti con E_p sotto la prima riga, 1 keV, esclusi), `range_mm` e
  `light_MeV` su tutti i bin e `bins` (per bin: `Ep_low_MeV`, `Ep_high_MeV`, `range_mm`,
  `light_MeV`). Ogni riassunto ha `samples`, `observed_mean`, `expected_mean`, `residual_mean`,
  `residual_sd`, `residual_sem`, `relative_residual` (Σd / Σatteso) e `relative_sem`; `null`
  quando non è definito (medie senza protoni, deviazioni ed errori con meno di due).

`riptide_sim_predict` scrive su stdout la previsione delle stesse verifiche prima del run (le copie
committate stanno in `docs/validation/predictions/`). In `prediction.chords`:
- `mode`, `acceptance`, `expected_mean_chord_mm`; nel pencil anche `L_mm`;
- nel campo lontano, la distribuzione delle corde di `--chord-samples` raggi (default 10⁸) estratti
  con il seme comune: `samples`, `seed`, `bins`, `bin_width_mm` (10⁻² mm), `sample_hits`,
  `sample_acceptance_pull`, `sample_mean_chord_mm`, `sample_sd_chord_mm`, `min_chord_mm`,
  `max_chord_mm`. I pesi dei colpi sommano all'accettanza in forma chiusa, quindi la frazione dei
  raggi che mancano non ha errore di campionamento.
`prediction.chord_sample_shift` ripete la previsione con un secondo campione (seme + 1) e ne
riporta `seed`, `fractions_in_sigma`, `legendre_in_sigma` e `u_noncentrality`: quanto pesa la
dimensione finita del campione; nel pencil è la stringa `"pencil beam: one chord, no sampling"`.

Dalla versione 10 (2026-11-21):
- la corda è calcolata per evento da `truth_n_start_{x,y,z}` (schema 4) e la profondità è la
  frazione s = z/L_i, con `depth_z.high` = 1 invece della corda in mm;
- c'è `results.geometric_acceptance`; `beam` ha `mode` e, per il campo lontano, la sorgente;
- `must_be_zero` ha `start_not_on_source` e `interaction_off_chord`.
Oracolo: sulle corse `c2` della fase 2 (pencil, 50 bin) i CSV e `checks.json` sono identici bit per
bit a meno delle chiavi nuove, delle versioni e dei bordi di `depth_z` divisi per L (scarto
≤ 1.3·10⁻¹⁶ in relativo).

Dalla versione 8 (2026-11-20): c'è `results.deposits.protons_nuclear_recoil`; i protoni con un
rinculo nucleare sono esclusi da `contained_balance_MeV` e da `pstar` (D22); le verifiche dei
depositi richiedono lo schema 3.

Dalla versione 7 (2026-11-20): c'è `results.deposits.pstar` e c'è `pstar_comparison.csv`; il
programma richiede `reference/pstar/pstar_216.txt` per i run con schema ≥ 2.

Dalla versione 6 (2026-11-19): c'è `results.deposits` e ci sono `deposit_spectra.csv` e
`proton_exits.csv`.

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

## Preparazione dell'analisi (`riptide_analysis_prepare`)

```bash
riptide_analysis_prepare --run-id <run> --analysis-id <id> --seed 20260907 \
    --sigma-Ep 0.10 --sigma-theta 0.08 --cuts single_np,has_proton,contained \
    [--hp-from-check <check_id>] [--max-events <n>] [--data-dir <radice>]
```

Il programma legge i file degli eventi del run e applica la catena di tagli (D30). I tagli
disponibili sono `single_np` (classe vera), `has_proton` (`truth_valid_np` = 1, D34),
`contained` (`truth_p_exits` = 0) ed `evis_min=<MeV>`; la catena del Caso A è
`single_np,has_proton,contained`. Poi converte
`truth_Ep` nell'E_p fisico con la massa HP (D28) ed estrae la risoluzione con il seme dato (D31).
Tutti i parametri dell'analisi sono obbligatori e nessuno ha un valore di default.
`--hp-from-check` serve per i run con `run_record_version` 1. Se la cartella
`analysis/<run_id>/<analysis_id>/` esiste già, il programma esce con codice 1 senza toccarla.
Il campione resta in memoria e non va su disco. Il programma scrive solo:

- `cutflow.csv`, con le colonne `index,cut,threshold_MeV,passed,fraction`. La riga 0 è `all`, cioè
  gli eventi letti. Seguono i tagli nell'ordine in cui sono applicati; `fraction` è relativa alla
  riga 0.
- `analysis.json`, scritto per ultimo: la sua presenza segna un'analisi completa. Le chiavi
  (`analysis_record_version` = 2) sono:
  - `analysis_id`, `project_version`, `git_describe`.
  - `input`: `run_id`, `run_record_version`, `schema_version`, `sim_git_describe`,
    `events_in_run`, `max_events` (0 = tutti).
  - `hp_target_mass`: `source` (`run.json` o `checks.json`, nel secondo caso con `check_id`),
    `file`, `sha256`, `awr`, `target_mass_MeV`.
  - `response`: `seed`, `sigma_Ep_MeV`, `sigma_theta_rad`, `per_event_engine`.
  - `selection`: `text`, `events_read`, `events_selected`, e `cuts`, con `name`,
    `threshold_MeV` (solo per `evis_min`) e `passed` per ciascun taglio.
  - `summaries.hp_conversion`, calcolati sugli eventi selezionati:
    - `events`;
    - `events_without_offset`;
    - `relative_offset` {`min`, `max`, `mean`}, dove l'offset è (`truth_Ep` − `truth_Ep_phys`) /
      `truth_Ep_phys`;
    - `lab_angle_relative_max`, lo scarto relativo massimo di `truth_Ep_phys` dalla cinematica n–p
      fisica in θ_p.
  - `summaries.resolution`:
    - `Ep_MeV` e `theta_p_rad`, ciascuno con {`mean`, `std`, `max_abs`, `pull_mean`,
      `pull_variance`} dei residui reco − truth. Il pull della media è m/(σ/√N), quello della
      varianza è (s²/σ² − 1)/√(2/(N−1)).
    - `p_dir`, con {`events_without_track`, `angle_max`, `norm_deviation_max`,
      `scaled_angle_squared_mean`, `pull_mean`}. δ è l'angolo fra `reco_p_dir` e la traccia vera.
      δ²/(2σ_θ²) segue Exp(1) e il suo pull è (media − 1)·√N, con N il numero di eventi con
      traccia.
    - Gli eventi senza protone (`truth_p_dir` = (0, 0, 0), cioè `truth_valid_np` = 0) hanno
      `reco_p_dir` = (0, 0, 0). Entrano nei residui di Ep e θ_p, non in `p_dir`, e sono contati
      in `events_without_track`.

    I pull sono `null` quando σ = 0, oppure, per `p_dir`, quando nessun evento ha una traccia. Un
    residuo non finito interrompe il programma, perché nei massimi sparirebbe.

  La versione 1 (fino a `d19fd67`) non aveva `events_without_track`. Normalizzava anche le
  direzioni nulle, per cui `p_dir.pull_mean` usciva `null` nei run con eventi senza protone.

Gli `analysis_id` `ctest_prepare_a`, `ctest_prepare_b` e `ctest_prepare_c` sono riservati al test
d'integrazione (`analysis/apps/tests/run_prepare_it.sh`). Il test ricrea le loro cartelle sotto
`analysis/ctest_it_1t/` e controlla tre cose:

- con lo stesso seme, `analysis.json` (a meno di `analysis_id`) e `cutflow.csv` sono identici;
- con un altro seme cambia `analysis.json` ma non il cut-flow;
- una cartella esistente è rifiutata.

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

Chiavi (versione 2 del formato; la versione 1 non ha `hp_elastic_H1`, D28):

| chiave | tipo | contenuto |
|---|---|---|
| `run_record_version` | int | versione del formato di `run.json` |
| `run_id` | string | nome della cartella del run |
| `project_version` | string | versione del progetto (CMake) |
| `git_describe` | string | `git describe --tags --always --dirty` dei sorgenti alla build; `-dirty` indica modifiche non committate, `unknown` che git non è disponibile |
| `geant4_version` | string | `major.minor.patch` di Geant4 |
| `schema_version` | int | versione dello schema dei file degli eventi |
| `config_version` | int | versione del formato del config: 1 per i run della fase 1b, 2 da kB di Birks in poi (`config/README.md`) |
| `seed` | uint | seme del run, da cui derivano i semi per evento |
| `events` | uint | eventi richiesti e simulati |
| `threads` | int | thread worker |
| `started_utc` | string | inizio del run, ISO 8601 in UTC (`2026-09-07T12:34:56Z`) |
| `elapsed_s` | double | tempo reale da `Initialize` alla fine di `BeamOn` [s] |
| `event_files` | array di string | `events_t0.root` … `events_t<N−1>.root` |
| `config` | object | config effettivo; rileggibile con `parseRunConfig`. Dalla versione 2 contiene kB di Birks con la sua fonte e la funzione di passo applicata |
| `hp_elastic_H1` | object | dalla versione 2: `file`, `sha256`, `rep_flag`, `awr`, `frame_flag`, `target_mass_MeV` dell'intestazione G4NDL di H-1 letta all'avvio del run, lo stesso blocco di `checks.json` (D14, D28). L'analisi ne usa `target_mass_MeV` per la conversione di D24 |

`git_describe` viene rigenerato a ogni build (`cmake/RiptideGitDescribe.cmake`) e descrive i
sorgenti al momento della build, non al momento del run.

## Schema degli eventi (v4)

Albero `events`: le 19 colonne della v1, le 18 della v2, quella della v3, poi le 3 della v4, in quest'ordine. Ci sono solo due tipi, `Int32` (`Int_t`) e `Double`
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
| `det_Edep`, `det_Evis` | Double | energia depositata nello scintillatore da tutte le particelle e luce visibile con Birks [MeV] (v2); per i protoni la luce continua L(K_pre) − L(K_pre − E_dep), per le altre particelle `G4EmSaturation` (D21, dal commit 3af4c3a) |
| `truth_np_vertex_{x,y,z}` | Double | vertice del primo n–p [mm] (v2) |
| `truth_p_end_{x,y,z}` | Double | punto finale del protone di rinculo nello scintillatore [mm] (v2) |
| `truth_p_length` | Double | cammino del protone nello scintillatore [mm] (v2) |
| `truth_p_exits` | Int32 | 1 se il protone esce dallo scintillatore (v2) |
| `truth_p_Edep`, `truth_p_Evis` | Double | deposito e luce visibile della sola traccia del protone [MeV] (v2); la luce è continua, L(E_p) per un protone che si ferma (D21) |
| `truth_p_n_hadronic` | Int32 | interazioni adroniche del protone nello scintillatore (v2) |
| `truth_n_nC_elastic`, `truth_n_nC_inelastic`, `truth_n_other` | Int32 | interazioni n–C elastiche, n–C anelastiche e altre del neutrone primario (v2) |
| `truth_n_secondary_n_hadronic` | Int32 | interazioni adroniche dei neutroni non primari nello scintillatore (v2) |
| `truth_event_class` | Int32 | codice `EventClass` (D19): 0 nessuna, 1 singolo n–p, 2 n–p multiplo, 3 solo n–C elastico, 4 n–C anelastico senza n–p, 5 misto, 6 altro (v2) |
| `truth_p_n_nuclear_recoil` | Int32 | step non adronici del protone di rinculo che hanno creato un rinculo nucleare, cioè un secondario con numero barionico > 0 (scattering Coulombiano singolo, D22) (v3) |
| `truth_n_start_{x,y,z}` | Double | punto di partenza del neutrone primario [mm] (v4): fisso nel fascio pencil, estratto per evento sul disco della sorgente a campo lontano (D27); con `truth_omega_n` dà la corda dell'evento nello scintillatore |

Unità: MeV, mm, ns, rad. Nomi, unità e descrizioni stanno in `event_schema.hpp`, che è l'unica
fonte: la riga scritta dalla sim (`sim/truth/.../event_row.hpp`) è verificata colonna per colonna
contro lo schema da un test.

Il primo n–p può venire dopo un'altra interazione, per esempio un n–C elastico. In quel caso il
neutrone non viaggia più lungo Ω_n, e θ_p può superare π/2. Per il campione pulito del singolo
scattering si seleziona quindi `truth_first_interaction == 1`. Le colonne `reco_` arriveranno con
l'analisi.

Le colonne del vertice n–p e del protone (`truth_np_vertex_*`, `truth_p_*`) valgono 0 se
`truth_valid_np` = 0. Il protone è la traccia nata nel primo n–p, senza i suoi secondari; i
dettagli e il margine di `G4NuclearStopping` su `truth_p_Edep` stanno in D19.

Le versioni successive aggiungono colonne solo in coda e non cambiano quelle esistenti. Il
lettore (`sim/reader`) accetta le versioni da 1 a quella della build e collega solo le colonne
presenti nel file. Le altre restano a zero, quindi in un file v1 le colonne v2 valgono 0 e in un file v2 `truth_p_n_nuclear_recoil` vale 0. In un file anteriore alla v4 `truth_n_start_*` vale 0: quei run sono tutti pencil, e il punto di partenza si prende da `source.position_mm` del config.
