# Configurazioni di run

File JSON letti da `riptide_sim --config`. Le chiavi portano l'unità nel nome (`_mm`, `_MeV`);
chiavi mancanti, sconosciute o non valide fanno fallire il run (`sim/config/run_config.hpp`).

| Chiave | Significato |
|---|---|
| `config_version` | versione del formato, oggi `1` |
| `geometry.scint_size_mm` | lati del blocco di scintillatore, centrato nell'origine |
| `source.mode` | `mono` (con `En_MeV`) o `flat` (con `En_min_MeV`, `En_max_MeV`) |
| `source.direction` | versore del fascio (tolleranza 10⁻⁶, poi normalizzato) |
| `source.position_mm` | posizione del cannone, fuori dal blocco |
| `physics_list` | nome per `G4PhysListFactory` |

L'energia deve stare in (0, 20] MeV: oltre 20 MeV finiscono i dati valutati di `G4NeutronHP`.

- `mono_14MeV.json`: fascio monoenergetico a 14 MeV.
- `flat_10_20MeV.json`: spettro piatto 10–20 MeV.

**Provvisorio:** il blocco è un cubo di 10 cm in entrambi i file. Le dimensioni definitive
verranno dal progetto di ottimizzazione.
