# Validação da revisão final

Arquivo entregue: `compare_gauge_ground_state.py`.
Referência: arquivo `compare_gauge_ground_state.py` fornecido em Downloads.
Ambiente de execução: Python 3, NumPy 1.26.4 e Matplotlib.

## Resultado isotrópico calculado

Para Kx = Ky = Kz = Kt = Kw = 1:

| Quantidade | Resultado |
| --- | ---: |
| epsilon_old | -1.063544 |
| epsilon_new | -0.957196 |
| Delta_epsilon = epsilon_new - epsilon_old | +0.106349 |
| OLD: zero-energy band gap | 2.000000 |
| NEW: zero-energy band gap | 0.000000 |

OLD é gapped; NEW é gapless. Setor energeticamente favorecido: OLD.
Os valores foram calculados; o arredondamento é apenas de apresentação.
O gap informado é o intervalo espectral total em torno de E=0.

## Modificações realizadas

1. Adição de `tracked_band_structure`, com `np.linalg.eigh` e avaliação das
   24 permutações de quatro bandas por `itertools.permutations`, maximizando
   a soma dos overlaps quadráticos entre passos consecutivos.
2. Adição de `_align_degenerate_eigenvectors`: alinhamento unitário da base
   em subespaços degenerados com a base do passo anterior, por SVD do NumPy.
   A tolerância identifica degenerescências na escala do arredondamento de
   ponto flutuante e é independente de GAP_TOL. Os autovalores não são alterados.
3. Adição de asserts em cada q: igualdade do espectro rastreado ordenado com
   `eigvalsh(H(q))` e verificação dos autovetores após o alinhamento.
   TRACKING_TOL = 1e-12; as tolerâncias originais foram preservadas.
4. A figura de dispersão recebe exclusivamente as bandas rastreadas. Preserva
   OLD à esquerda, NEW à direita, BZ física, escala comum de energia, linha
   E=0 e 300 dpi. As cores identificam trajetórias de autoestados. Não há
   legenda de ordenação energética; o rótulo vertical passou de E_n(q) a E(q).
5. Adição de recálculo independente de F(q) e das energias após o tracking,
   sempre pelo caminho original de `eigvalsh`, com asserts explícitos para
   F_OLD, F_NEW, epsilon_old, epsilon_new e Delta_epsilon antes/depois.
6. Saída de acoplamentos, energias e gaps com seis casas decimais, sinal
   explícito em Delta_epsilon e remoção da tabela e das classificações
   duplicadas. Diagnósticos de erro mantêm notação científica. O CSV mantém
   sua precisão original.
7. Substituição do rótulo impresso `gap` por `zero-energy band gap`, sem
   modificar o cálculo. Inclusão de mensagens de aprovação dos novos testes.

Não foram adicionados SciPy, suavização, interpolação ou fórmulas analíticas
em substituição à diagonalização.

## Verificações executadas

- Script original e script final executados com sucesso. Todos os testes
  antigos passaram, incluindo hermiticidade, espectro real/finito,
  normalização por quatro sítios e comparação das BZ física e estendida.
- Comparação textual confirmou 18 funções originais integralmente
  preservadas, incluindo coeficientes, Hamiltonianos, F(q), integração,
  normalização, gap, testes antigos, figura de F(q) e escrita do CSV.
  Acoplamentos, malhas, limites das BZ e tolerâncias originais são iguais.
- Os espectros ordenados e os arrays de F(q) do original e da revisão
  são exatamente iguais em toda a malha física.
- Tracking validado nos 4001 pontos de cada setor, totalizando 8002
  verificações. Maior diferença entre os espectros: 4,441e-15.
- epsilon_old, epsilon_new e Delta_epsilon são exatamente iguais aos do
  script original, antes do arredondamento. Diferença de cada valor: zero.
- Sanity check isotrópico aprovado, com os resultados acima.
- Testes com fases arbitrárias e bases degeneradas rotacionadas preservaram
  as trajetórias NEW. Percorrer a malha no sentido inverso produziu as mesmas
  trajetórias, salvo uma permutação global das colunas.
- Um Hamiltoniano auxiliar de teste confirmou a continuação em um cruzamento
  exato. Outro confirmou a preservação de um pequeno anticruzamento resolvido.
  Esses Hamiltonianos são exclusivos dos testes externos à implementação física.
- As figuras foram inspecionadas visualmente e têm metadados de 300 dpi.
  Os cruzamentos reais e as cúspides de F_NEW(q) estão preservados.
- O CSV e o PNG de F(q) regenerados são idênticos byte a byte aos produzidos
  pelo script original.

## Execução

Com NumPy e Matplotlib instalados, execute no diretório desejado para as saídas:

```bash
python3 compare_gauge_ground_state.py
```

O script gera `dispersion_old_vs_new.png`, `Fq_old_vs_new.png` e
`ground_state_energy_comparison.csv` no diretório corrente.
