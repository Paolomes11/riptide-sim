# Risultati della fase 6 — Caso A su dati Geant4

Data: 2026-11-25. Previsioni e soglie: `docs/validation/6_predictions.md`, versione 1 in
`994ece0` (§§1–9) e aggiunta della §10 (D40) in `0385bc6`, entrambe committate **prima** dei
run a cui si riferiscono. Dopo i run le soglie non sono state toccate.

**Esito.**
- Senza accettanza, la prima valutazione ha dato un'**anomalia**: la famiglia gemella passa
  107/117, la calibrazione assoluta di Geant4 28/32.
- La diagnosi a 14 MeV ne ha trovato la causa: la likelihood non contiene l'accettanza in θ_p
  della selezione.
- D40 aggiunge al modello il prior q(θ_p | E_n, sel), misurato su un run indipendente. Con q
  (passo D) passano:
  - la famiglia gemella, 117/117;
  - la calibrazione assoluta, 32/32 per Geant4 e per il toy;
  - il bias globale del piatto;
  - la prova a metà.
- Una previsione non è rispettata: a 19 MeV, N = 1000, `std_sqrt_N` di Geant4 cresce del 3.6%
  contro il ±3% previsto (§8).

Etichette: **(a)** fatto stabilito; **(b)** stima con ipotesi esplicite; **(c)** verificato con
Monte Carlo; **(d)** ipotesi da mettere alla prova.

Tutti i comandi sono passati da `scripts/run_capped.sh -m 5G`, uno alla volta, dopo
`scripts/data_root.sh check`. Output sull'SSD:
- `stats/<stats_id>/case_a.json`;
- le valutazioni, fatte con script fuori dal repo che applicano le formule di
  `6_predictions.md` così come sono scritte, in `validation/6_eval/`, `validation/6_diag_14MeV/`,
  `validation/6_acceptance/` e `validation/6a_eval/`.

## 1. Run e verifiche `c1`

Build `release` a `994ece0` (`git describe` v0.7.0-16-g994ece0) per i quattro run della §1 di
`6_predictions.md`, e a `7ec81bb` per il run piatto del passo D. Geant4 11.4.1, `QGSP_BIC_HP`,
BC-408 60 × 20 × 20 mm, direzione +X.

| Run | p χ² (u / z / frazioni) | max\|pull\| | a₁ | Protoni che escono | Esito |
|---|---|---|---|---|---|
| `6_ff_flat_10_19p5MeV_x_1e7` | 0.299 / 0.790 / 0.975 | 1.44 | −0.01435 | 36727 | PASS **(a)** |
| `6_ff_mono_11MeV_x_1e6` | 0.176 / 0.520 / 0.935 | 1.65 | −0.01045 | 2755 | PASS **(a)** |
| `6_ff_mono_17MeV_x_1e6` | 0.709 / 0.240 / 0.393 | 1.91 | −0.01407 | 4194 | PASS **(a)** |
| `6_ff_mono_19MeV_x_1e6` | 0.411 / 0.040 / 0.523 | 1.63 | −0.01568 | 4638 | PASS **(a)** |
| `6_ff_flat_10_19p5MeV_x_1e7_s20260908` | 0.530 / 0.869 / 0.821 | 1.87 | −0.0151 | 36902 | PASS **(a)** |

- Gli oracoli deterministici sono a zero o entro le soglie esatte su tutti e cinque i run:
  contatori, corda, bilancio dei depositi, kB, `observed_non_hp`.
- **Oracolo HP, precisazione (a).** La §1 non dice quali estremi dell'offset HP − fisico
  confrontare. Si usano `rounding_free_relative_{min,max}`, come stabilito in D17 e come nella
  fase 3. Per esempio, a 11 MeV il minimo osservato è 5.3827725·10⁻⁴, contro il limite
  5.3827717·10⁻⁴.
- Gli a₁ osservati sono tutti entro 2σ dalle attese della §1 **(c)**.
- Il confronto con PSTAR non è stato rifatto: le bande valgono solo per il 14 MeV e il piatto
  isotropo della fase 3.

## 2. Selezione

| Run | Selezionati | Attesi (§2) | Rapporto |
|---|---|---|---|
| piatto | 799076 | 8.01·10⁵ | 0.998 |
| 11 MeV | 98693 | 9.95·10⁴ | 0.992 |
| 14 MeV | 826770 | 826770 | 1.000 |
| 17 MeV | 68736 | 6.88·10⁴ | 0.999 |
| 19 MeV | 61069 | 6.13·10⁴ | 0.996 |

- Tutti i run sono ben dentro la banda di ±10% **(b)**, e nessuna riga è ridotta **(a)**.
- **Correzione (a).** La §2 dà 826 859 selezionati a 14 MeV: è la catena della fase 4 senza
  `has_proton` (89 eventi in più, D34). Il valore giusto per la catena
  `single_np,has_proton,contained` è 826770. Le attese degli altri run sono scalate da quel
  numero, quindi si spostano di 10⁻⁴ in relativo, molto meno della banda.
- **Run del passo D (b).** Il piatto con seme 20260908 seleziona 798388 eventi: 688 meno del
  piatto della fase 6, entro i ±2700 previsti dalla §10.

## 3. Valutazione senza accettanza (prima valutazione)

Record `6_<run_id>{,_pilot,_toy}`. Valutazione in `validation/6_eval/eval_6.txt`.

**Famiglia Geant4 = toy (K = 117, z = 3.93): 107/117 (a).** Le 10 voci che falliscono:

| Voce | G4 | toy | z |
|---|---|---|---|
| 11 MeV N=1000 `std_sqrt_N` | 0.809 | 0.730 | +8.73 |
| 11 MeV N=100 `std_sqrt_N` | 0.860 | 0.761 | +7.23 |
| 17 MeV N=1000 `std_sqrt_N` | 1.084 | 0.987 | +6.53 |
| 19 MeV N=1000 `std_sqrt_N` | 1.162 | 1.050 | +6.13 |
| 14 MeV N=100 `std_sqrt_N` | 1.018 | 0.920 | +5.22 |
| 17 MeV N=100 `std_sqrt_N` | 1.187 | 1.077 | +4.63 |
| 14 MeV N=1000 `std_sqrt_N` | 0.931 | 0.871 | +4.55 |
| 14 MeV N=1000 `bias_mean` | −0.0227 | +0.0072 | −4.25 |

Le altre due sono `std_sqrt_N` sopra soglia. In sintesi, Geant4 è il 6–11% più largo del toy a
tutte le energie, e a 14 MeV ha un bias di −0.023 MeV.

**Calibrazione assoluta (K = 32, z = 3.60) (a).**
- Il toy passa 32/32.
- Geant4 passa 28/32. Falliscono i `pull_mean` a N = 1000 a 11, 14 e 17 MeV (fra −0.81 e
  −0.68, z fino a −5.2) e quello a N = 100 a 11 MeV (−0.33).
- Su quelle voci il toy è circolare (D38).

**Bias globale del piatto (b): PASS.**
- Δ osservato +0.166, Δ atteso +0.180, residuo −0.014, tolleranza 0.101.
- La stima prima del run era ≈ +0.19.
- Coverage: G4 0.703 / 0.919 / 0.964, toy 0.707 / 0.922 / 0.965.

**Lettura.** La §5 non prevedeva nessuna voce di contrazione oltre soglia. Il verso del bias non
era previsto e la grandezza di `std_sqrt_N` è fuori dal fattore 2 della §3: secondo le regole
della §3 è un'**anomalia**, da indagare prima della chiusura.

## 4. Diagnosi a 14 MeV

In `validation/6_diag_14MeV/diag_14.txt`. Tutto su `3_ff_mono_14MeV_x_1e7`, N = 100 e 1000.

1. **Riproduzione (a).** Lo script ricalcola le righe di contrazione di `case_a.json` con
   differenza 0.
2. **Cinematica, respinta (a).** Il residuo fra la verità di Geant4 corretta per l'HP e la
   cinematica del modello ha rms ~10⁻⁵ in E_p e ~10⁻⁶ in θ, trascurabile rispetto a σ_Ep e
   σ_θ.
3. **Larghezze per evento uguali (c).** Per bin di θ_p ricostruito, `mean_std` di G4 e del toy
   coincidono entro l'1%.
4. **Accettanza in θ_p, confermata (c).**
   - Il rapporto delle distribuzioni di θ_p fra G4 e toy vale 0.58 a 0–5°, 0.74 a 5–10°,
     circa 1.06 a θ intermedi e 0.92 sopra 80°.
   - Il toy ripesato alla distribuzione di Geant4 (H3) riproduce Geant4:
     - `std_sqrt_N` 1.010 / 0.942, contro 1.018 / 0.932 di G4;
     - bias a N = 1000 −0.021, contro −0.023.
   - Le varianti H1 e H2 della diagnosi danno gli stessi numeri di Geant4 entro 10⁻⁴, quindi
     non spiegano la differenza.

**Conclusione (c).** Gli eventi con θ_p piccolo hanno E_p ≈ E_n e sono i più informativi.
Geant4 li perde per `single_np`, perché il neutrone diffuso riurta. La likelihood con
sin 2θ_p li dà per presenti, quindi il posterior di Geant4 è più largo di quello del toy e
leggermente spostato.

## 5. Misura dell'accettanza

In `validation/6_acceptance/README.md`, con le previsioni scritte prima del run.

**Chiusura (a).** I conteggi per bin di θ_p coincidono con quelli della diagnosi, e q/iso
riproduce il rapporto G4/toy.

**Cut-flow a 14 MeV (a), ε_tot = 0.792:**

| θ_p | anisotropia | ε_single | ε_cont | q/iso |
|---|---|---|---|---|
| 0–5° | 1.063 | 0.453 | 0.965 | 0.586 |
| 5–10° | 1.068 | 0.581 | 0.957 | 0.751 |
| 40–45° | 1.003 | 0.869 | 0.965 | 1.061 |
| 85–90° | 0.957 | 0.751 | 1.000 | 0.908 |

- **Deficit in avanti (a).** Viene da ε_single: a 0–5°, il 35% degli eventi è n–p multiplo e
  il 20% è misto. Previsione confermata.
- **Contenimento (a).** ε_cont perde al più il 5.6%, con il minimo a 20–25°. Previsione
  confermata in parte.
- **Anisotropia (a).** Il primo urto ha un'anisotropia del ±7% rispetto alla massa non
  relativistica, contro il ±5% previsto: previsione mancata di poco.
- **Dipendenza da E_n (a).** La forma dipende da E_n in modo non liscio: χ²/dof di 519/324
  (costante), 420/306 (lineare) e 393/288 (quadratico). Da qui la tabella in bin di E_n di D40,
  invece di un modello parametrico.
- **Prior fisico non relativistico (b).** `physicalRecoilAngleLogPrior` usa la massa non
  relativistica, mentre il toy campiona con la cinematica relativistica. Lo scarto va
  da ±0.7% a 14 MeV a ±1% a 20 MeV. La tabella q di D40 lo assorbe.

## 6. Passo D di D40: risultati con l'accettanza

**Catena (a).**
1. Run piatto `6_ff_flat_10_19p5MeV_x_1e7_s20260908`, con `c1` come nella §1.
2. `riptide_acceptance` con `--stats-id 6acc_s20260908`: 798388 eventi, tabella da 10 bin di
   E_n × 18 bin di θ_p, `ratio_max` 1.080.
3. 15 record `6a_<run_id>{,_pilot,_toy}` con `--acceptance 6acc_s20260908`.

In ogni `case_a.json` lo script controlla `settings.acceptance.{stats_id, run_id, half}`, e la
provenienza è OK su tutti i record. Valutazione in `validation/6a_eval/eval_6a.txt`, con lo
script della fase 6 puntato agli id `6a_*`.

**Famiglia Geant4 = toy: 117/117 (c), max|z| 3.16.** Le voci più grandi:

| Voce | G4 | toy | z |
|---|---|---|---|
| 17 MeV N=10 `pull_width` | 1.114 | 0.950 | +3.16 |
| 11 MeV N=100 `std_sqrt_N` | 0.872 | 0.828 | +3.14 |
| 17 MeV N=10 coverage 0.90 | 0.860 | 0.918 | −2.71 |
| 14 MeV N=1000 `std_sqrt_N` | 0.951 | 0.979 | −2.03 |

Nessuna voce supera la soglia. Le due vicine a 3σ_c sono isolate e di verso opposto.
Previsione confermata.

**Calibrazione assoluta: G4 32/32 e toy 32/32 (c).**
- La voce più grande di G4 è il `pull_mean` a 19 MeV, N = 100: +0.17, z +2.83.
- I `pull_mean` di G4 a N = 1000, che a 11, 14 e 17 MeV fallivano fra −0.81 e −0.68, ora
  valgono −0.07, −0.05 e −0.06. A 19 MeV si passa da −0.31 a +0.38 (z +2.24).

**Bias a N = 1000 (b), contro le formule fissate della §5:**

| E_n | G4 senza q | G4 con q | toy con q | Δ G4 − toy | Δ atteso |
|---|---|---|---|---|---|
| 11 | −0.0209 | −0.0020 ± 0.0041 | +0.0007 ± 0.0037 | −0.0027 ± 0.0055 | +0.0019 |
| 14 | −0.0226 | **−0.0014 ± 0.0047** | +0.0075 ± 0.0050 | −0.0089 ± 0.0068 | +0.0011 |
| 17 | −0.0230 | −0.0019 ± 0.0051 | +0.0003 ± 0.0040 | −0.0022 ± 0.0065 | −0.0011 |
| 19 | −0.0119 | +0.0140 ± 0.0066 | +0.0047 ± 0.0057 | +0.0093 ± 0.0087 | −0.0026 |

- A 14 MeV la previsione |`bias_mean`| ≲ 0.01 MeV è **confermata**.
- Il Δ da accettanza in E_n (g·`std_sqrt_N`²) resta sotto la sensibilità, come previsto.
- **A 19 MeV** il bias di G4 a N = 100 e 1000 è positivo, +0.021 ± 0.008 e +0.014 ± 0.007, con
  il `pull_mean` più alto della famiglia. Sono tutte voci entro soglia. Si annota come punto da
  rivedere **(d)**: 19 MeV è vicino al bordo alto della tabella q, che arriva a 19.5 MeV. La
  causa non è verificata.

**Bias globale del piatto (b): PASS.**
- Δ osservato +0.1645, Δ atteso +0.1832, residuo −0.0187, tolleranza 0.1015.
- Coverage: G4 0.703 / 0.918 / 0.964, toy 0.707 / 0.921 / 0.965, |Δ| ≤ 0.004 (atteso ≲ 0.01).

**Prior su θ_p (c).** `physical` qui vuol dire q. Bias di G4 / toy:

| E_n | `physical` | `flat` |
|---|---|---|
| 11 | +0.003 / −0.002 | −0.090 / −0.094 |
| 14 | +0.009 / −0.008 | −0.104 / −0.118 |
| 17 | +0.001 / −0.005 | −0.124 / −0.123 |
| 19 | +0.008 / +0.000 | −0.127 / −0.129 |

Tutto Geant4 = toy, come previsto.

**Prior su E_n (c).** `shift_mean` decresce con N per ogni mono e per entrambe le sorgenti. A
N = 1 lo spostamento è massimo a 19 MeV (1.71) e minimo a 11 MeV (0.29). Previsioni
confermate.

## 7. Prova a metà

`riptide_acceptance --half even --stats-id 6acc_s20260908_even` usa solo gli eventi con
`event_id` pari: 398902 eventi. Con questa tabella si è rifatto il `case_a` Geant4 completo a
14 MeV, `6ah_even_3_ff_mono_14MeV_x_1e7` (`validation/6a_eval/eval_half.txt`).

| | `bias_mean` N=1000 | `pull_mean` | coverage |
|---|---|---|---|
| q da tutto il run | −0.00140 | −0.046 | 0.65 / 0.975 / 1.0 |
| q dalla metà pari | +0.00028 | +0.009 | 0.65 / 0.975 / 1.0 |

- Lo spostamento vale +0.00168 MeV. Con σ_c = 0.00683 è 0.246·σ_c, contro la soglia di 0.3:
  **PASS (c)**.
- Il margine è piccolo. Per costruzione la metà della statistica raddoppia la varianza della
  tabella. A questa statistica il rumore di q pesa quindi meno di un terzo dell'errore del bias,
  ma non è trascurabile.
- La metà dispari non è stata usata: la §10 chiede solo la metà pari.

## 8. Larghezza: confronto con la fase 6

`std_sqrt_N` a N = 100 / 1000, prima e dopo q:

| E_n | G4 N=100 | G4 N=1000 | toy N=100 | toy N=1000 |
|---|---|---|---|---|
| 11 | 0.860 → 0.872 (+1.3%) | 0.809 → 0.822 (+1.6%) | 0.760 → 0.828 (+8.9%) | 0.730 → 0.810 (+10.9%) |
| 14 | 1.018 → 1.034 (+1.6%) | 0.931 → 0.951 (+2.1%) | 0.920 → 1.026 (+11.5%) | 0.871 → 0.979 (+12.3%) |
| 17 | 1.187 → 1.209 (+1.9%) | 1.084 → 1.111 (+2.5%) | 1.077 → 1.196 (+11.1%) | 0.987 → 1.093 (+10.7%) |
| 19 | 1.271 → 1.305 (+2.7%) | 1.162 → **1.204 (+3.6%)** | 1.193 → 1.356 (+13.6%) | 1.050 → 1.204 (+14.6%) |

- **Il toy converge verso Geant4 (c).** A 14 MeV, N = 100, il toy va da 0.920 a 1.026, contro
  ≈ 1.01 previsto. Il meccanismo previsto, H3 della diagnosi, è confermato.
- **G4 cambia poco: previsione (b) mancata in un punto.**
  - Sette valori su otto stanno entro ±3%.
  - A 19 MeV, N = 1000, l'aumento è del +3.6%.
  - L'aumento cresce con E_n a entrambi gli N, quindi sembra un effetto sistematico di q
    piuttosto che una fluttuazione **(d)**: la likelihood con q dà meno peso alle regioni di
    θ_p soppresse.
  - I due record usano gli stessi eventi, quindi l'errore della differenza è più piccolo di
    `std_sqrt_N_error` di ciascuno, ma non è stato stimato. La soglia non si cambia: il ±3%
    resta una previsione mancata.
- **Contrazione 1/√N (c).**
  - `std_sqrt_N` decresce da N = 10 a 100 per tutti.
  - Fra 100 e 1000 cambia meno del 10%, tranne il toy a 19 MeV: 1.356 → 1.204, cioè −11%.
  - Lo stesso caso era già fuori senza q (−12%). Non è legato a D40; è un punto minore della
    contrazione a 19 MeV con 200 / 40 esperimenti **(d)**.

## 9. Tempi e memoria (a, senza soglia)

Tempo trascorso e RSS da `/usr/bin/time -v`.

| Passo | Tempo | RSS |
|---|---|---|
| sim piatto 10⁷, 8 thread (seme 20260908) | 2:45 | 178 MB |
| `c1` sul piatto 10⁷ | 52 s | 122 MB |
| `riptide_acceptance` (tutto / metà pari) | 8.5 s / 8.6 s | 1.41 / 1.40 GB |
| `case_a` piatto: pilota / completo / toy | 11.1 / 12.7 / 4.4 s | 1.41 GB / 1.41 GB / 80 MB |
| `case_a` 11 MeV: pilota / completo / toy | 4.9 / 35.0 / 34.1 s | 0.41 GB / 0.41 GB / 48 MB |
| `case_a` 14 MeV: pilota / completo / toy | 13.5 / 42.2 / 34.2 s | 1.41 GB / 1.41 GB / 82 MB |
| `case_a` 17 MeV: pilota / completo / toy | 4.8 / 34.9 / 41.0 s | 0.41 GB / 0.41 GB / 47 MB |
| `case_a` 19 MeV: pilota / completo / toy | 4.7 / 34.1 / 33.3 s | 0.41 GB / 0.41 GB / 46 MB |
| `case_a` 14 MeV con q della metà pari | 42.3 s | 1.42 GB |

- Tutti i piloti soddisfano il criterio (< 3 GB, < 60 s).
- Il termine angolare per riga di E_n costa al più il +5% rispetto ai ~40 s della fase 6,
  entro il ≤ +30% previsto **(b)**.
- `riptide_acceptance` impiega 8.5 s contro i ~10 s stimati, con la stessa memoria.
- L'RSS dei mono 10⁶ (0.41 GB) supera la stima < 0.3 GB della §8; non ha una soglia, e la
  ripartizione fra lettura e tabella q non è stata misurata.

## 10. Checklist

Nell'ordine bias → risoluzione → pull → coverage → contrazione → prior. I numeri sono quelli con
q, passo D.

| Voce | Esito | Etichetta |
|---|---|---|
| Bias | 14 MeV N = 1000 −0.0014 ± 0.0047; piatto entro la formula fissata; 19 MeV +0.014 ± 0.007, entro soglia ma da rivedere | (c), (d) per il 19 MeV |
| Risoluzione | Geant4 = toy su tutte le `std_sqrt_N`; G4 cambia del +1.3…+3.6% rispetto a senza q, con il ±3% previsto mancato a 19 MeV, N = 1000 | (c), (b) mancata |
| Pull | `pull_mean` / `pull_width` entro le due famiglie; massimo +0.38 a 19 MeV, N = 1000 (z 2.24) | (c) |
| Coverage | entro le due famiglie; piatto G4 − toy ≤ 0.004 | (c) |
| Contrazione 1/√N | decresce ovunque; toy a 19 MeV −11% fra 100 e 1000 (già senza q) | (c), (d) |
| Robustezza al prior | θ_p: `physical` (q) e `flat` Geant4 = toy; E_n: `shift_mean` Geant4 = toy, decrescente con N; prova a metà 0.246·σ_c | (c) |

## 11. Punti aperti lasciati dalla fase 6

- **19 MeV (d).** Con q, il bias di G4 a N ≥ 100 è positivo (+0.014/+0.021) e la larghezza
  cresce del 3.6%. Il bordo della tabella q (19.5 MeV) è il primo candidato da verificare. Va
  ripreso quando si useranno sorgenti fino a 20 MeV (fase 9).
- **Prior fisico non relativistico (b).** `physicalRecoilAngleLogPrior` va reso relativistico,
  oppure va dichiarato che vale solo dove q non si usa.
- **Dipendenza non liscia di q da E_n (a).** L'origine non è verificata; la struttura in E_n
  della sezione d'urto n–C è un candidato.
- **Contrazione del toy a 19 MeV (d).** −11/−12% fra N = 100 e 1000, sopra il 10%. Si
  ripresenta sia con q sia senza; una statistica di esperimenti più alta lo chiarirebbe.
