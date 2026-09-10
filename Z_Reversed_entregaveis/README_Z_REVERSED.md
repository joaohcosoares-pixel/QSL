# Entregáveis Z-REVERSED

Abra `Quadrupolar_Spin_Liquid_Z_Reversed.md` para a derivação completa de 17 seções. A referência OLD é obrigatória e foi preservada: a discrepância da Eq.(26) está transcrita e demonstrada no relatório. A implementação segue a geometria original e inverte somente as duas famílias z.

No ponto isotrópico:

| quantidade | OLD | Z-REVERSED |
|------------|-----|------------|
| energia por sítio | -1.063544 | -1.063544 |
| zero-energy band gap completo | 2.000000 | 2.000000 |
| distância mínima ao zero | 1.000000 | 1.000000 |
| fluxos x, ordem W124,W234,W123,W134 | +,+,+,+ | -,-,+,+ |
| fluxos w, ordem W213,W143,W214,W243 | +,+,+,+ | -,-,+,+ |

A diferença de energia é exatamente zero, demonstrada analiticamente. Os setores não são gauge-equivalentes. A comparação não determina um estado fundamental global entre todos os gauges.

Para reproduzir os testes e os resultados, mantenha estes três scripts na mesma pasta:

```text
quadrupolar_spin_liquid_z_reversed.py
compare_old_vs_z_reversed.py
draw_z_lattice.py
```

Dependências: Python, NumPy, SymPy, SciPy e Matplotlib. Versões efetivamente utilizadas: NumPy 1.26.4, SymPy 1.12, SciPy 1.11.4 e Matplotlib 3.6.3.

```bash
python3 -W error quadrupolar_spin_liquid_z_reversed.py
python3 -W error compare_old_vs_z_reversed.py
```

O primeiro executa verificações simbólicas, 1000 casos assimétricos em ambos os setores e as cadeias finitas. O segundo calcula convergência, quadratura, gaps, bandas e gera os PNGs a 300 dpi e os CSVs. As últimas execuções de ambos terminaram com saída 0 e sem warnings. Mensagens ocorridas durante a preparação estão discriminadas no relatório.

Os JSONs e logs preservam os resultados brutos; `energy_convergence_old_vs_z.csv` registra as cinco malhas. `ground_state_energy_old_vs_z.csv` usa o gap completo e `TIE` para o empate isotrópico.

O pacote ZIP contém o relatório final com links relativos, figuras e cópias intactas da referência OLD e do PDF original para facilitar compartilhamento. Basta abrir o Markdown em um leitor com suporte a equações LaTeX. O manifesto inclui hashes dos arquivos entregues.

Veredito físico: APROVADO COM RESSALVAS, pela errata documental explicitada. Veredito do código: A, no escopo analisado. É seguro discutir com o orientador apresentando essas ressalvas.
