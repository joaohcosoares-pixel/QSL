"""Bundle results with portable links and unmodified source snapshots."""
from pathlib import Path
import hashlib
import json
import re
import zipfile

root = Path(__file__).resolve().parent
old_source = Path("/home/joaohenrique/Estudo/Pesquisa/Pesquisa/QSL/Quadrupolar_Spin_Liquid_OLD_reference.md")
pdf_source = Path("/home/joaohenrique/Estudo/Pesquisa/Vanuildo (Mat. Condensada)/Quadrupolar Spin Liquid in 1D (cópia).pdf")
filenames = [
    "quadrupolar_spin_liquid_z_reversed.py", "compare_old_vs_z_reversed.py", "draw_z_lattice.py",
    "qsl_lattice_z_reversed.png", "dispersion_old_vs_z_reversed.png", "Fq_old_vs_z_reversed.png",
    "ground_state_energy_old_vs_z.csv", "energy_convergence_old_vs_z.csv",
    "validation_z_reversed.json", "comparison_z_reversed.json", "geometry_flux_audit.json",
    "validation_z_reversed.log", "comparison_z_reversed.log", "README_Z_REVERSED.md",
]
entries = {name:(root/name).read_bytes() for name in filenames}
for path in sorted((root/"auditoria_fontes").glob("original-*.png")):
    entries[f"auditoria_fontes/{path.name}"] = path.read_bytes()
entries["fontes/Quadrupolar_Spin_Liquid_OLD_reference.md"] = old_source.read_bytes()
entries["fontes/PDF_original.pdf"] = pdf_source.read_bytes()
report = (root/"Quadrupolar_Spin_Liquid_Z_Reversed.md").read_text()
report = report.replace(str(old_source),"fontes/Quadrupolar_Spin_Liquid_OLD_reference.md")
report = report.replace(str(pdf_source),"fontes/PDF_original.pdf")
report = report.replace(str(root)+"/","")
for link in re.findall(r"\]\(<([^>]+)>\)",report):
    assert link in entries, link
entries["Quadrupolar_Spin_Liquid_Z_Reversed.md"] = report.encode("utf-8")
manifest = {name:{"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest()} for name,data in entries.items()}
entries["MANIFEST_SHA256.json"] = (json.dumps(manifest,indent=2,ensure_ascii=False)+"\n").encode()
destination = root/"Z_Reversed_entregaveis.zip"
with zipfile.ZipFile(destination,"w",compression=zipfile.ZIP_DEFLATED) as archive:
    for name,data in entries.items():
        archive.writestr(name,data)
with zipfile.ZipFile(destination) as archive:
    assert archive.testzip() is None
    assert set(archive.namelist()) == set(entries)
print(f"Verified {destination.name}: {len(entries)} entries, {destination.stat().st_size} bytes; all report links resolve within archive.")
