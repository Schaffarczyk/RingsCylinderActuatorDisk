# VortexCode V4 — Strukturwechsel: Knoten-Kollokation + Wirbelkern
**Claude, 2026-07-27 — nach van Kuik & Lignarolo, Wind Energy 19 (2016), App. A+B**

## Idee

V1–V3 (Mittelpunkts-Kollokation, kernfrei) haben einen strukturellen
Zielkonflikt: Stabilität erzwingt orts­abhängige Dämpfung, die Kante bleibt
unterbestimmt, cP behält ein Fluktuationsband (~10⁻³·|cP| bei moderatem N).
V4 wechselt die Struktur nach dem Vorbild von van Kuik & Lignarolo:

1. **Ringe mit Zirkulationen Γᵢ sitzen an den Knoten** (xᵢ, rᵢ);
   **Ring 1 ist an der Scheibenkante (0,1) fixiert** und von den
   Randbedingungen ausgenommen (vgl. Paper Abschn. 4.3).
2. **Wirbelkern-Regularisierung:** Für Auswertungsabstände
   ρ < ½·β·dsᵢ werden geglättete (Marshall-)Ausdrücke benutzt.
   **Abweichung vom Paper:** statt festem δ = 0.002 ein
   **zellweitenproportionaler Kern** β·dsᵢ mit β = e^(−1/2) ≈ 0.607
   (konsistent mit dem analytischen Streifenintegral). Der feste Kern
   verfälscht cP um ~10⁻², sobald ds ≫ δ variiert — im Prototyp gemessen.
3. **Halbunendlicher Zylinder** bei x = 30·R_wake (Radius R_wake,
   Stärke γ_wake); v_r-Vorzeichen korrekt (Befund C, SubA-V3).
4. **Einphasiges Schema mit EINHEITLICHER Dämpfung 0.05:**
   Γᵢ ← Γᵢ + 0.05·(Γ*ᵢ − Γᵢ) mit Γ*ᵢ = −cT/(2vᵢ)·dsᵢ, und
   rᵢ ← rᵢ + 0.05·(ψ_wake − ψᵢ)/(rᵢ(1+vzᵢ)). Der Kern beseitigt die
   Sägezahn-Mode — keine z-abhängige Dämpfung mehr nötig.
   Das v_n-„fine-tuning" des Papers (Phase 2, d = 0.0025) wurde im
   Prototyp getestet: in beiden Vorzeichen instabil und **überflüssig**,
   da die Ψ-Phase bereits auf Maschinenniveau konvergiert; v_n dient
   als Diagnostik (Spalte 4 in sls-V4.DAT).

## Validierung (Python-Prototyp adcode_v4.py, Betz cT = −8/9)

| N | max|Δψ| (Ende) | cP | cP − (−16/27) |
|---|---|---|---|
| 400 | 3·10⁻¹³ | −0.5920118 | +5.8·10⁻⁴ |
| 800 | 2·10⁻¹² | −0.5923623 | +2.3·10⁻⁴ |
| 1600 | 8·10⁻¹² | −0.5925434 | **+4.9·10⁻⁵** |

Eigenschaften: **eindeutiger Fixpunkt** (kein cP-Band, keine Drift,
Residuum → Maschinengenauigkeit), monotone Konvergenz in ~1500
Iterationen, Fehlerskalierung beschleunigt sich Richtung 1/N².
v_n nach Konvergenz: Kante ~3·10⁻², Blattmitte ~2·10⁻³; nur die
Zylinder-Naht (letzte ~10 Ringe) behält v_n ~ 0.5 (bekanntes
Modellartefakt des abrupt beginnenden Zylinders).
Vergleich: V3 (Mittelpunkte) bei N = 1000–4000 → Band ±(0.5…4)·10⁻⁴
mit Restdrift; vK&L-Paper: „a few ‰", cP auf 1 ‰.

## Dateien und Benutzung

- `VortexCode-V4.f` — Hauptprogramm (Modul memv4 + initv4 + ringsv4 +
  cppv4); verlinkt `mem.f` und `SubA-V3.f` (elliptische Integrale,
  Zylinderfunktionen; Kopien liegen bei)
- `inpa-V4.dat` — Eingabe (cT, N, LtubeFac, beta, dampA/G, maxiter,
  epsPSI, c1, c2, ncpout)
- `compVortexCode-V4.cmd` — Kompilierskript
- `adcode_v4.py` — Python-Prototyp (Validierung; braucht numpy/scipy)
- Ausgaben: `conv-V4.DAT` (niter, max|Δψ|, rms, cP, errcP),
  `sls-V4.DAT` (x, r, γ, v_n, Δψ)

    compVortexCode-V4.cmd
    VortexCode-V4.exe

Erwartung (N = 1600): max|Δψ| fällt monoton; nach ~1500 Iterationen
~10⁻¹⁰ (epsPSI), cP ≈ −0.59254, errcP ≈ 8·10⁻⁵. Laufzeit: O(N²) pro
Iteration wie V3; cP-Auswertung (3000·N) nur alle ncpout Iterationen.

**Hinweis:** Wie bei V3 stand in der Sandbox kein gfortran zur
Verfügung — Quellen sind statisch geprüft (Zeilenlängen, Balance) und
der Algorithmus ist 1:1 gegen den validierten Python-Prototyp
geschrieben. Compilermeldungen bitte zurückspielen.

## Einordnung / Grenzen

- V4 löst die **Iterations- und Eindeutigkeitsprobleme** (Band, Drift,
  Dämpfungs-Tuning) — der verbleibende cP-Fehler ist reiner, glatter
  Diskretisierungsfehler und mit Richardson sauber extrapolierbar.
- Die **Kantensingularität selbst** bleibt kerngeglättet (γ an Ring 1
  beschränkt); wer die 1/√s-Zone auflösen will, braucht weiterhin die
  analytische Kantenzelle oder die mehrwertige Parametrisierung
  (README-Kante-und-FineTuning.md, Optionen 2–3).
- Erwartbare nächste Schritte: N-Studie mit V4 (Richardson dann mit
  sauberem p), Propeller-Fall cT = 1, weiche γ-Rampe an der
  Zylinder-Naht.






