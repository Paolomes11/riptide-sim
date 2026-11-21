# Configurazioni di run

File JSON letti da `riptide_sim --config`. Le chiavi portano l'unità nel nome (`_mm`, `_MeV`);
chiavi mancanti, sconosciute o non valide fanno fallire il run (`sim/config/run_config.hpp`).

| Chiave | Significato |
|---|---|
| `config_version` | versione del formato: un nuovo run vuole `2` |
| `geometry.scint_size_mm` | lati del blocco di scintillatore, centrato nell'origine |
| `source.mode` | `mono` (con `En_MeV`) o `flat` (con `En_min_MeV`, `En_max_MeV`) |
| `source.direction` | versore del fascio (tolleranza 10⁻⁶, poi normalizzato) |
| `source.position_mm` | posizione del cannone, fuori dal blocco |
| `physics_list` | nome per `G4PhysListFactory` |
| `scintillator.birks_kB_mm_per_MeV` | costante di Birks del materiale [mm/MeV], ≥ 0; `0` spegne il quenching |
| `scintillator.birks_source` | provenienza del valore (id della decisione, oggi `D18`) |
| `proton_step_function.dRoverRange` | funzione di passo di Geant4 per muoni e adroni, protoni compresi: frazione del range, in (0, 1] |
| `proton_step_function.finalRange_mm` | lunghezza sotto cui il passo non si riduce più, > 0 |

Birks si applica passo per passo, quindi la luce di un protone che si ferma dipende dalla funzione
di passo. Per questo è obbligatoria e scritta per esteso. I file usano i default di Geant4 11.4.1
(0.2 e 0.1 mm), che vengono impostati in modo esplicito: Geant4 non permette di rileggere il valore
in uso. Nel log di inizializzazione ogni processo di ionizzazione stampa `StepFunction=(...)`.

La versione `1`, senza `scintillator` né `proton_step_function`, è quella dei run della fase 1b.
Si legge ancora dal loro `run.json`, ma non si può usare per un nuovo run.

L'energia deve stare in (0, 20] MeV: oltre 20 MeV finiscono i dati valutati di `G4NeutronHP`.

- `mono_14MeV.json`: fascio monoenergetico a 14 MeV.
- `flat_10_20MeV.json`: spettro piatto 10–20 MeV.
- `flat_19p5_20MeV.json`: spettro piatto 19.5–20 MeV, la finestra in cui l'elastico HP e CHIPS
  si sovrappongono (D16); serve a verificarne la rampa.

**Provvisorio:** il blocco è un cubo di 10 cm in tutti i file. Le dimensioni definitive
verranno dal progetto di ottimizzazione.
