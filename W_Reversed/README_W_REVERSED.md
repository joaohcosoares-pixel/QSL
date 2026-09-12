# Auditoria W-REVERSED

O relatório principal é `Quadrupolar_Spin_Liquid_W_Reversed.md`, com 29 seções,
dedução completa e tabelas produzidas da execução. Todos os gráficos, integrais
e buscas principais usam **q em [-pi,pi]**.

No ponto Kx=Ky=Kz=Kt=Kw=1, OLD é favorecido sobre W:
epsilon_old = -1.063544409973365; epsilon_w = -0.838804985931099;
Delta = +0.224739424042266. A tolerância efetiva de comparação é 1e-10.
Veja no relatório a divergência documental da Eq.(26), sem alteração do original.

## Reproduzir

Python 3.10 ou superior; NumPy, SciPy, SymPy e Matplotlib. As versões efetivamente
utilizadas estão nos JSON. Na pasta destes entregáveis, execute em sequência:

```bash
OPENBLAS_NUM_THREADS=1 python3 -W error quadrupolar_spin_liquid_w_reversed.py
OPENBLAS_NUM_THREADS=1 python3 -W error draw_w_lattice.py
OPENBLAS_NUM_THREADS=1 python3 -W error compare_old_vs_w_reversed.py
OPENBLAS_NUM_THREADS=1 python3 -W error build_report_w_reversed.py
```

O primeiro comando realiza os testes simbólicos, 1000 casos aleatórios e cadeias
finitas. O terceiro gera dispersões, energia, convergência, comparação finita,
consolidação Z e benchmark com cinco repetições. A execução não suprime warnings.
Os scripts de desenho definem uma pasta temporária gravável para o cache do Matplotlib.

Qualquer edição de um script exige reexecutar os geradores antes do builder:
os resultados registram SHA-256 de **todos** os scripts da pasta. O builder
rejeita hashes inconsistentes, diferenças CSV/JSON e fontes locais alteradas.
`report_manifest_w_reversed.json` registra também os hashes dos artefatos finais.

As validações críticas permanecem ativas com `python3 -O -W error`.
Para outros acoplamentos, `compare_old_vs_w_reversed.py --couplings Kx Ky Kz Kt Kw`
gera nova comparação numérica; o relatório analítico deste pacote exige acoplamentos
unitários e recusa resultados anisotrópicos. A busca genérica de gap é numérica;
no ponto isotrópico há também uma prova analítica do mínimo global.

## Arquivos e fontes

- `quadrupolar_spin_liquid_w_reversed.py`: geometria → RAW → Eq.(27) → Bloch; testes.
- `symbolic_w_audit.py`: transcrição independente do PDF, determinante e raízes.
- `draw_w_lattice.py`: geometria e fluxos, comparações de arestas e desenho.
- `compare_old_vs_w_reversed.py`: eigvalsh, quadratura, tracking visual e benchmark.
- `build_report_w_reversed.py` e `report_template_w.md`: relatório autocontido e proteção de procedência.
- `sources/`: snapshots inalterados da referência OLD e dos alvos de regressão Z.
- `geometry_audit_notes.md`, `algebra_audit_notes.md`: notas das verificações independentes.
- CSV/JSON/PNG: resultados novos, metadados, figuras a 300 dpi e tabelas de execução.

O PDF original e arquivos na pasta de pesquisa foram apenas lidos. Seus hashes
estão nos JSON. Os snapshots locais permitem reproduzir a regressão sem depender
dos caminhos originais; o PDF original deve acompanhar uma nova conferência visual
da fonte. A versão Z nomeada com `(1)` não existia; a versão Obsidian encontrada
está identificada nos metadados. Upper-reversed foi omitido da consolidação por
ausência de conjunto numérico auditado e compatível nesta execução.

Os resultados descrevem matéria no gauge fixado. Não demonstram o mínimo entre
todos os gauges, classificação topológica ou projeção física de paridade finita.
