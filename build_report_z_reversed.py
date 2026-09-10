"""Assemble the reviewed derivation with tables from actual execution records."""
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parent
validation = json.loads((ROOT/"validation_z_reversed.json").read_text())
comparison = json.loads((ROOT/"comparison_z_reversed.json").read_text())

convergence = ["| Nq | epsilon_old | epsilon_z | Delta_epsilon |",
               "|----|-------------|-----------|---------------|"]
for row in comparison["convergence"]:
    convergence.append(f'| {row["Nq"]} | {row["epsilon_old"]:.12f} | {row["epsilon_z"]:.12f} | {row["Delta_epsilon"]:+.3e} |')

special = ["| q | OLD (numérico) | Z-REVERSED (numérico) |",
           "|---|----------------|----------------------|"]
for point in ("-pi/2","-pi/4","0","pi/4","pi/2"):
    values = [", ".join(f"{x:.6f}" for x in comparison["special_points"][name][point])
              for name in ("OLD","Z-REVERSED")]
    special.append(f"| {point} | ({values[0]}) | ({values[1]}) |")

labels = {
    "unitarity":"Unitariedade de U",
    "raw_formula_vs_geometry":"RAW: fórmula versus ligações",
    "eq27":"Chamada matricial Eq.(27)",
    "bloch_vs_closed":"Bloch da geometria versus fórmulas fechadas",
    "q2_cancellation":"Cancelamento de q2",
    "hermiticity":"Hermiticidade",
    "majorana":"Relação Majorana",
    "matrix_period_2pi":"Periodicidade matricial 2pi",
    "matrix_covariance_pi":"Covariância matricial pi",
    "spectrum_period_pi":"Periodicidade espectral pi",
    "raw_bloch_spectrum":"Espectros RAW versus Bloch",
    "same_q_pairing":"Pairing no mesmo q",
    "polynomial_relative":"Polinômio característico (normalizado)",
    "only_z_matrix_difference":"Diferença OLD versus Z somente nos termos z",
    "finite_chain_bloch":"Cadeia finita versus Bloch",
    "finite_chain_raw":"Cadeia finita versus RAW",
    "finite_chain_full_spectrum":"Espectro total finito versus união dos blocos",
    "finite_chain_energy_per_site":"Energia por sítio finita em ambas representações",
}
errors = ["| Verificação | Erro máximo |", "|-------------|-------------|"]
for key,value in validation["max_errors"].items():
    errors.append(f"| {labels[key]} | {value:.3e} |")

text = (ROOT/"report_sections_1_8.md").read_text()+"\n"+(ROOT/"report_sections_9_17.md").read_text()
for marker,value in (("{{CONVERGENCE_TABLE}}","\n".join(convergence)),
                     ("{{SPECIAL_POINTS_TABLE}}","\n".join(special)),
                     ("{{VALIDATION_TABLE}}","\n".join(errors))):
    assert text.count(marker) == 1
    text = text.replace(marker,value)
assert "{{" not in text
assert chr(956) not in text
headers = re.findall(r"^## (\d+)\.",text,flags=re.M)
assert headers == [str(x) for x in range(1,18)], headers
for target in re.findall(r"\]\(<([^>]+)>\)",text):
    assert Path(re.sub(r":\d+$","",target)).exists(), target
destination = ROOT/"Quadrupolar_Spin_Liquid_Z_Reversed.md"
destination.write_text(text,encoding="utf-8")
print(f"Wrote {destination.name}: {len(text)} characters, all 17 sections, live result tables and existing links.")
