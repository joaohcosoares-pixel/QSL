# Auditoria independente da geometria ORIGINAL e W-REVERSED

Esta nota separa transcrição documental, evidência geométrica e dedução. As fontes foram apenas lidas, sem alterações. A zona principal solicitada é sempre `PROJECT_BZ = (-pi, pi)` para a variável `q = k.a1`.

## Fontes e inspeção visual

- PDF original: `/home/joaohenrique/Pesquisa/Vanuildo (Mat. Condensada)/Quadrupolar Spin Liquid in 1D (cópia).pdf`; SHA256 `e8d763f2b40c0df4cbd764b1d467d278e348076f589922752da54a5df5761e1e`.
- Referência documental OLD: `/home/joaohenrique/Pesquisa/Pesquisa/QSL/Quadrupolar_Spin_Liquid_OLD_reference_obsidian.md`; SHA256 `3c8dec124881e58a66e5f86cb95b07124bd82130a6aea79159bddae8d38b6156`.
- O PDF possui oito páginas. As páginas 3 a 8 foram renderizadas diretamente com Poppler; as figuras e equações relevantes foram inspecionadas visualmente. A extração textual do manuscrito é pouco confiável e não foi usada para decidir sinais.
- Página 3: Eq.(10), posições das quatro sub-redes, vetores a1/a2 e começo do Hamiltoniano real.
- Página 4: continuação da Eq.(11), definição de fluxo orientado Eqs.(13)-(14), e figura do gauge OLD.
- Página 5: convenção das setas, lista de dez valores de gauge Eq.(15), Fourier Eq.(16), começo da expansão Eq.(17).
- Página 6: termos individuais e agrupamento de Fourier da Eq.(17).
- Página 7: final da Eq.(17), matriz de Bloch e Eqs.(21)-(26).
- Página 8: Eq.(27) e energia fundamental por sítio Eq.(28).

Os scripts anteriores servem somente para organização de software. Os dez bonds abaixo foram conferidos nos pontos terminais da Eq.(11) e nas setas da figura da página 4.

## Posições e tradução de célula

Os sítios de uma célula de referência estão em

\[
r_1=0,\quad r_2=a_1,\quad r_3=a_1+a_2,\quad r_4=a_2.
\]

A tradução que conserva o rótulo das quatro sub-redes é `T = 2 a1`. Denotamos por `m_n` o sítio de sub-rede `m` na célula `n`. Uma representação fiel da conectividade usa `a1=(2,0)`, `a2=(1,2)` e portanto

\[
r(1_n)=(4n,0),\quad r(2_n)=(4n+2,0),\quad
r(3_n)=(4n+3,2),\quad r(4_n)=(4n+1,2).
\]

A escala do desenho é convencional; nenhuma seta é inferida da escala. Nos dois desenhos as coordenadas são idênticas. A repetição das famílias horizontais x/w é alternada, com ambas apontando para a direita embaixo e para a esquerda em cima no OLD.

## Gauge original congelado e alteração permitida

A figura da página 4 e a Eq.(15) na página 5 concordam com a lista:

| Chave | Bond armazenado | OLD | W-REVERSED |
|---|---|---:|---:|
| u12x | 1 para 2, x | +1 | +1 |
| u34x | 3 para 4, x | +1 | +1 |
| u23y | 2 para 3, y | +1 | +1 |
| u41y | 4 para 1, y | -1 | -1 |
| u42z | 4 para 2, z | +1 | +1 |
| u31z | 3 para 1, z | +1 | +1 |
| u43w | 4 para 3, w | +1 | -1 |
| u12w | 1 para 2, w | -1 | +1 |
| u13t | 1 para 3, t | -1 | -1 |
| u24t | 2 para 4, t | -1 | -1 |

Uma seta de `i` para `j` significa `u_ij,gamma=+1`; a antissimetria é `u_ji,gamma=-u_ij,gamma`. Os índices do bond armazenado não são necessariamente o sentido da seta: no OLD, `u12w=-1` representa a seta de `2_(n-1)` para `1_n`.

## Dez bonds reais, incluindo os deslocamentos de Fourier

`da1,da2` representam o deslocamento físico do segundo extremo menos o primeiro, e `cell_shift` representa a célula final menos a inicial. São quantidades diferentes. Para todo bond:

\[
(da1,da2)=2\,cell\_shift\,(1,0)+r_{end}-r_{start},
\]

com as posições internas expressas na base `(a1,a2)`.

| start | end | gamma | stored u OLD | da1 | da2 | cell_shift | Ligação em células |
|---:|---:|:---:|---:|---:|---:|---:|---|
| 1 | 2 | x | +1 | +1 | 0 | 0 | 1_n, 2_n |
| 3 | 4 | x | +1 | -1 | 0 | 0 | 3_n, 4_n |
| 2 | 3 | y | +1 | 0 | +1 | 0 | 2_n, 3_n |
| 4 | 1 | y | -1 | 0 | -1 | 0 | 4_n, 1_n |
| 4 | 2 | z | +1 | +1 | -1 | 0 | 4_n, 2_n |
| 3 | 1 | z | +1 | +1 | -1 | +1 | 3_n, 1_(n+1) |
| 4 | 3 | w | +1 | -1 | 0 | -1 | 4_n, 3_(n-1) |
| 1 | 2 | w | -1 | -1 | 0 | -1 | 1_n, 2_(n-1) |
| 1 | 3 | t | -1 | +1 | +1 | 0 | 1_n, 3_n |
| 2 | 4 | t | -1 | +1 | +1 | +1 | 2_n, 4_(n+1) |

Na Eq.(11), o bond z31 aparece originalmente como `3(R+a2-a1),1(R)`; transladar os dois extremos por `2a1` produz `3_n,1_(n+1)`. Analogamente o bond t24 aparece como `2(R-a1),4(R+a2)` e é equivalente a `2_n,4_(n+1)`. Essas translações comuns dos dois extremos não alteram os deslocamentos físicos.

## Oito fluxos derivados de caminhos anti-horários

Os caminhos abaixo têm área orientada positiva no desenho, calculada pela fórmula do polígono. Os quatro caminhos do bloco x e os quatro do bloco w possuem área `2` nas unidades gráficas. O nome da plaqueta registra apenas as sub-redes; as células explicitadas no percurso evitam confundir triângulos diferentes.

| Plaqueta | Percurso anti-horário | Produto de variáveis armazenadas | Fatores OLD | OLD | Fatores W | W-REVERSED | Mudou? |
|---|---|---|---|---:|---|---:|---|
| W124 | 1_0 → 2_0 → 4_0 → 1_0 | u12x (-u42z) u41y | (+1)(-1)(-1) | +1 | (+1)(-1)(-1) | +1 | Não |
| W234 | 2_0 → 3_0 → 4_0 → 2_0 | u23y u34x u42z | (+1)(+1)(+1) | +1 | (+1)(+1)(+1) | +1 | Não |
| W123 | 1_0 → 2_0 → 3_0 → 1_0 | u12x u23y (-u13t) | (+1)(+1)(+1) | +1 | (+1)(+1)(+1) | +1 | Não |
| W134 | 1_0 → 3_0 → 4_0 → 1_0 | u13t u34x u41y | (-1)(+1)(-1) | +1 | (-1)(+1)(-1) | +1 | Não |
| W213 | 2_-1 → 1_0 → 3_-1 → 2_-1 | (-u12w) (-u31z) (-u23y) | (+1)(-1)(-1) | +1 | (-1)(-1)(-1) | -1 | Sim |
| W143 | 1_0 → 4_0 → 3_-1 → 1_0 | (-u41y) u43w u31z | (+1)(+1)(+1) | +1 | (+1)(-1)(+1) | -1 | Sim |
| W214 | 2_-1 → 1_0 → 4_0 → 2_-1 | (-u12w) (-u41y) (-u24t) | (+1)(+1)(+1) | +1 | (-1)(+1)(+1) | -1 | Sim |
| W243 | 2_-1 → 4_0 → 3_-1 → 2_-1 | u24t u43w (-u23y) | (-1)(+1)(-1) | +1 | (-1)(-1)(-1) | -1 | Sim |

O algoritmo usa a lista de bonds transladada, adiciona explicitamente a aresta reversa com sinal oposto e procura cada passo do caminho nesse grafo. Os resultados de fluxo não são armazenados como entrada do cálculo. A expectativa de que somente ciclos que contêm w mudem é verificada depois do produto.

Os setores não são gauge-equivalentes porque quatro fluxos locais diferem. Isso exclui até transformações locais dependentes de célula: os fatores de sítio sempre se cancelam em um ciclo fechado. A busca adicional das 16 matrizes `diag(s1,s2,s3,s4)` não encontra solução. É também possível ver a contradição diretamente: preservar `u12x` exige `s1*s2=+1`, enquanto inverter `u12w` exige `s1*s2=-1` numa transformação periódica de quatro sub-redes.

## Hamiltoniano real e diferença exclusiva em w

Suprimindo o sobrescrito y nos operadores e somando sobre todas as células:

\[
\begin{aligned}
H_{OLD}=i\sum_n\{&K_x(\theta_{1n}\theta_{2n}+\theta_{3n}\theta_{4n})
+K_y(\theta_{2n}\theta_{3n}-\theta_{4n}\theta_{1n})\\
&+K_z(\theta_{4n}\theta_{2n}+\theta_{3n}\theta_{1,n+1})
+K_w(\theta_{4n}\theta_{3,n-1}-\theta_{1n}\theta_{2,n-1})\\
&-K_t(\theta_{1n}\theta_{3n}+\theta_{2n}\theta_{4,n+1})\}.
\end{aligned}
\]

O Hamiltoniano `H_W` mantém os termos x,y,z,t exatamente iguais e substitui o termo da família w por

\[
iK_w\sum_n(-\theta_{4n}\theta_{3,n-1}+\theta_{1n}\theta_{2,n-1}).
\]

Logo,

\[
\Delta H_W=2iK_w\sum_n(\theta_{1n}\theta_{2,n-1}-\theta_{4n}\theta_{3,n-1}).
\]

## Fourier RAW, sem usar as fórmulas finais como fonte

Para cada bond orientado `(i,j)`, sua contribuição à matriz RAW é

\[
(H_{raw})_{ij}\mathrel{+}=iK_\gamma u_{ij,\gamma}
e^{i(q_1 da1+q_2 da2)},\qquad
(H_{raw})_{ji}\mathrel{+}=[(H_{raw})_{ij}]^*.
\]

É necessário somar ambas as contribuições quando um par tem mais de um bond. Para o OLD, a lista acima produz:

\[
\begin{aligned}
f_{12}^{raw}&=i(K_xe^{iq_1}-K_we^{-iq_1}),\\
f_{13}^{raw}&=-ie^{iq_2}(K_te^{iq_1}+K_ze^{-iq_1}),\\
f_{14}^{raw}&=iK_ye^{iq_2},\\
f_{23}^{raw}&=iK_ye^{iq_2},\\
f_{24}^{raw}&=-ie^{iq_2}(K_te^{iq_1}+K_ze^{-iq_1}),\\
f_{34}^{raw}&=i(K_xe^{-iq_1}-K_we^{iq_1}).
\end{aligned}
\]

Os dois termos z armazenados apontam da cadeia superior para a inferior; inverter sua orientação para preencher os elementos 13 e 24 fornece `-iKz exp[-iq1+iq2]`. O fator transversal nesses elementos é, portanto, **exp(+iq2)**. A mesma fase aparece nos termos y e t que cruzam as cadeias.

Invertendo apenas os dois valores w da bond list, obtém-se

\[
f_{12,W}^{raw}=i(K_xe^{iq_1}+K_we^{-iq_1}),\qquad
f_{34,W}^{raw}=i(K_xe^{-iq_1}+K_we^{iq_1}),
\]

e os quatro demais elementos RAW permanecem idênticos.

Com `Theta_raw=U Theta_Bloch` e `U=diag(1,1,exp(-iq2),exp(-iq2))`, o Hamiltoniano transformado é `U† H_raw U`. Isso elimina o fator `exp(+iq2)` dos elementos que cruzam as cadeias. Os elementos 12 e 34 não recebem alteração de fase. A combinação inversa de fases no RAW e transformação `U H_raw U†` poderia ocultar um erro; por isso os elementos RAW devem ser verificados antes da Eq.(27).

## Divergência documental da Eq.(26)

Na página 7, a Eq.(26) literal escreve

\[
f_{34}^{(26),literal}=i(K_xe^{iq_1}-K_we^{-iq_1}).
\]

A referência OLD preserva corretamente essa transcrição. Entretanto, o bond x34 da Eq.(11) vai de `R+a1+a2` para `R+a2`; seu deslocamento é **-a1**. Logo sua contribuição em 34 é `iKx exp(-iq1)`. O bond w43 fornece `-iKw exp(+iq1)` ao elemento 34 por conjugação. Portanto a geometria e o Hamiltoniano real impõem

\[
f_{34}^{geometria}=i(K_xe^{-iq_1}-K_we^{iq_1}).
\]

A Eq.(17) agrupada da página 6, na linha `theta3† [...] theta4`, escreve justamente esta última expressão. O último termo da Eq.(17), no topo da página 7, também tem o conjugado correto em 43. A transformação da Eq.(27) não pode sanar a discrepância, pois seus fatores diagonais em 3 e 4 são iguais e se cancelam.

Assim, a divergência foi localizada entre Eq.(17)/geometria e a Eq.(26) final. O novo cálculo usa a primeira por hierarquia física e mantém a segunda documentada como literal. A diferença é

\[
f_{34}^{geometria}-f_{34}^{(26),literal}=2(K_x+K_w)\sin q_1.
\]

Ela desaparece em pontos especiais como `q1=0` e `q1=pi`, que seriam testes insuficientes. Para o OLD auditado e para W vale `f34=-f12*`; a Eq.(26) literal, em contraste, copia `f12`. Não se alterou o arquivo OLD de referência para esconder isso.

## Outro sinal local no manuscrito e seu limite

Na página 6, segunda linha da expansão de termos individuais, o termo `theta4† theta1` aparece visualmente com o fator `-iKy exp(+iq2)`, ao lado de `+iKy exp(+iq2) theta1† theta4`. A forma agrupada da **mesma** Eq.(17), mais abaixo, escreve `-iKy exp(-iq2) theta4† theta1`, que é o conjugado correto e concorda com a Eq.(11). Trata-se de um erro local da linha intermediária, sem ambiguidade após o agrupamento. Essa constatação reforça a necessidade de derivar o RAW das posições, verificá-lo separadamente e não copiar uma única linha manuscrita sem conferir Hermiticidade.

## Zona principal e o que a geometria implica

A célula com quatro sub-redes tem tradução `2a1`. Para o momento discreto de célula `p`, o parâmetro do manuscrito é `q1=p/2`, antes das fases internas. Isso não autoriza trocar a zona principal pedida: todos os resultados principais devem usar `q1 in [-pi,pi]`. A eventual covariância por deslocamento de pi deve ser tratada como redundância espectral, com contagem e normalização explicitadas no relatório científico e conferidas na cadeia finita. Na página 8 a Eq.(28) é apresentada textualmente como energia fundamental **por sítio**, e não como energia total.
