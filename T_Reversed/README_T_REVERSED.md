# T-REVERSED: entregáveis e reprodução

Leia `Quadrupolar_Spin_Liquid_T_Reversed.md` para a derivação e os resultados.
O setor parte exclusivamente de OLD e muda u13t e u24t de −1 para +1.

Execute nesta pasta:

```bash
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/qsl_t_mpl python3 -W error run_t_reversed.py
```

Dependências: Python 3.10+, NumPy, SciPy, SymPy, Matplotlib e Pillow.
As versões exatas da execução estão em todos os metadados. O executor compila
os scripts, roda a suíte de 1000 casos também sob `-O`, a auditoria simbólica,
cadeia finita, operador explícito de 64 estados, convergência, benchmark com
aquecimento e cinco repetições, figuras 300 dpi, CSVs, JSONs e o relatório.

O construtor `build_report_t_reversed.py` recusa hashes desatualizados em JSON
e CSV. `report_manifest_t_reversed.json` registra os hashes dos artefatos.
`sources/` preserva o PDF e referências sem edição. Os relatórios legados
nessa pasta são fontes de regressão; não constituem a base física de T.

No ponto isotrópico, a igualdade OLD/T é da energia de matéria por sítio no
limite termodinâmico. A projeção física de spins e a classificação topológica
não fazem parte da validação realizada. A divergência da Eq. (26) literal é
preservada e explicada no relatório, conforme a hierarquia de fontes solicitada.
