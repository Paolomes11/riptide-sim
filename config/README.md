# Configurazioni di run

File JSON letti da `riptide_sim --config`. Le chiavi portano l'unità nel nome (`_mm`, `_MeV`);
chiavi mancanti, sconosciute o non valide fanno fallire il run (`sim/config/run_config.hpp`).

| Chiave | Significato |
|---|---|
| `config_version` | versione del formato: un nuovo run vuole `3` |
| `geometry.scint_size_mm` | lati del blocco di scintillatore, centrato nell'origine |
| `source.beam` | `pencil` (fascio a matita) o `far_field` (sorgente a campo lontano, D27) |
| `source.mode` | `mono` (con `En_MeV`) o `flat` (con `En_min_MeV`, `En_max_MeV`) |
| `source.direction` | versore del fascio (tolleranza 10⁻⁶, poi normalizzato); in `far_field` anche la stringa `"isotropic"` |
| `source.position_mm` | solo `pencil`: posizione del cannone, fuori dal blocco |
| `source.margin_mm` | solo `far_field`: distanza in più del disco di partenza oltre la mezza diagonale del blocco, > 0 |
| `physics_list` | nome per `G4PhysListFactory` |
| `scintillator.birks_kB_mm_per_MeV` | costante di Birks del materiale [mm/MeV], ≥ 0; `0` spegne il quenching |
| `scintillator.birks_source` | provenienza del valore (id della decisione, oggi `D18`) |
| `proton_step_function.dRoverRange` | funzione di passo di Geant4 per muoni e adroni, protoni compresi: frazione del range, in (0, 1] |
| `proton_step_function.finalRange_mm` | lunghezza sotto cui il passo non si riduce più, > 0 |

Birks si applica passo per passo, quindi la luce di un protone che si ferma dipende dalla funzione
di passo. Per questo è obbligatoria e scritta per esteso. I file usano i default di Geant4 11.4.1
(0.2 e 0.1 mm), che vengono impostati in modo esplicito: Geant4 non permette di rileggere il valore
in uso. Nel log di inizializzazione ogni processo di ionizzazione stampa `StepFunction=(...)`.

Con `far_field` ogni evento ha una direzione Ω_n (fissa, o uniforme sulla sfera con
`"isotropic"`) e un punto di partenza uniforme su un disco perpendicolare a Ω_n, centrato sulla
retta per il centro del blocco, a distanza mezza diagonale + `margin_mm` a monte. Il raggio del
disco è il minimo che contiene l'ombra del blocco per una direzione fissa, la mezza diagonale per
la sorgente isotropa (D27). Il mondo si allarga fino a contenere tutto il disco.

Le versioni precedenti si leggono ancora dal `run.json` dei loro run, ma non si possono usare per
un nuovo run:

- `2`: run della fase 2, senza `source.beam`, letti come `pencil`;
- `1`: run della fase 1b, senza `scintillator` né `proton_step_function`.

L'energia deve stare in (0, 20] MeV: oltre 20 MeV finiscono i dati valutati di `G4NeutronHP`.

Fase 7b, Caso B su eventi Geant4 (D45): copie di `6_ff_flat_10_19p5MeV_x.json` in cui cambia solo
la direzione fissa; la direzione +X riusa il run della fase 6. 10⁷ eventi ciascuno.

- `7_ff_flat_10_19p5MeV_z.json`: direzione +Z, lungo un lato corto.
- `7_ff_flat_10_19p5MeV_obl.json`: direzione obliqua (3, 2, 1)/√14, scritta con 16 cifre.

Fase 6, Caso A su eventi Geant4 (D36): stesso blocco, fisica e sorgente di
`3_ff_mono_14MeV_x.json`, direzione fissa +X; cambia solo lo spettro.

- `6_ff_flat_10_19p5MeV_x.json`: spettro piatto 10–19.5 MeV (D25), 10⁷ eventi; checklist per
  evento con la verità distribuita sul dominio.
- `6_ff_mono_11MeV_x.json`, `6_ff_mono_17MeV_x.json`, `6_ff_mono_19MeV_x.json`: 11, 17 e 19 MeV,
  10⁶ eventi ciascuno; con il run `3_ff_mono_14MeV_x_1e7` danno quattro energie per la contrazione
  1/√N e la robustezza al prior.

Fase 3, blocco 60 × 20 × 20 mm con il lato lungo su X (D27), sorgente `far_field`:

- `3_ff_mono_14MeV_isotropic.json`: 14 MeV, direzioni isotrope.
- `3_ff_flat_10_19p5MeV_isotropic.json`: spettro piatto 10–19.5 MeV (D25), direzioni isotrope.
- `3_ff_mono_14MeV_x.json`: 14 MeV, direzione fissa +X lungo il lato lungo.

Fase 2, cubo di 10 cm con fascio `pencil`, tenuti per la regressione:

- `mono_14MeV.json`: fascio monoenergetico a 14 MeV.
- `flat_10_20MeV.json`: spettro piatto 10–20 MeV.
- `flat_19p5_20MeV.json`: spettro piatto 19.5–20 MeV, la finestra in cui l'elastico HP e CHIPS
  si sovrappongono (D16); serve a verificarne la rampa.
