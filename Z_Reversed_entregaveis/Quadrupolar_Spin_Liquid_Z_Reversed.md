# Quadrupolar Spin Liquid in 1D — Z-Reversed Gauge

Este estudo parte do gauge ORIGINAL e inverte somente as duas famílias periódicas de ligações z. A referência documental obrigatória é [Quadrupolar_Spin_Liquid_OLD_reference.md](<fontes/Quadrupolar_Spin_Liquid_OLD_reference.md>), lida em conjunto com a geometria, o Hamiltoniano real e a Fourier do [PDF original disponibilizado](<fontes/PDF_original.pdf>). A referência OLD permanece intacta. Seus enunciados documentais são distinguidos das deduções abaixo.

**Ressalva de fonte:** a referência reproduz deliberadamente a Eq.(26) impressa, sem correção. Essa equação diverge da geometria e da Eq.(11), e também da Eq.(17) na mesma cópia do PDF. Este relatório registra a transcrição literal e a errata demonstrada; o Hamiltoniano físico implementado segue a geometria original. Assim, “OLD” abaixo designa o gauge original reconstruído da geometria, com a correção documentada de fases em \(f_{34}\), e não a matriz obtida pela aplicação literal da Eq.(26) inconsistente.

As marcações identificam a origem de cada afirmação: **[DERIVADO DO PDF ORIGINAL]**, **[DEDUÇÃO MATEMÁTICA]**, **[RESULTADO NUMÉRICO]** e **[TESTE INDEPENDENTE]**. As imagens de fonte preservam a escrita original; toda a notação transcrita e o código novo usam \(u\).

## 1. Original geometry and conventions

[DERIVADO DO PDF ORIGINAL]

Os sítios da célula original são

\[
1:\mathbf R,\qquad 2:\mathbf R+\mathbf a_1,\qquad
4:\mathbf R+\mathbf a_2,\qquad
3:\mathbf R+\mathbf a_1+\mathbf a_2.
\]

A repetição de uma célula de quatro sub-redes ocorre por \(\mathbf T=2\mathbf a_1\). Escrevemos \(\mathbf R_n=n\mathbf T\) e \(m_n\) para o sítio \(m\) nessa célula. Esse índice explicita a periodicidade existente: não introduz outra célula nem modifica \(\mathbf a_1\) ou \(\mathbf a_2\).

![Geometria original, página 3 do PDF](<auditoria_fontes/original-3.png>)

![Setas do gauge ORIGINAL, página 4 do PDF](<auditoria_fontes/original-4.png>)

A convenção orientada é

\[
\boxed{u_{ij,\gamma}=-u_{ji,\gamma},\qquad
 i\longrightarrow j\iff u_{ij,\gamma}=+1.}
\]

A orientação escolhida para **armazenar uma variável** pode ser oposta à seta. Por exemplo, \(u_{41,y}=-1\) representa a seta \(1\to4\), enquanto \(u_{34,x}=+1\) representa a seta superior \(3\to4\), para a esquerda.

| bond | type | visual arrow | stored variable | original value |
|------|------|--------------|-----------------|----------------|
| \(1_n,2_n\) | x | \(1_n\to2_n\), inferior para a direita | `u12x` | +1 |
| \(3_n,4_n\) | x | \(3_n\to4_n\), superior para a esquerda | `u34x` | +1 |
| \(2_n,3_n\) | y | \(2_n\to3_n\), inferior para superior | `u23y` | +1 |
| \(4_n,1_n\) | y | \(1_n\to4_n\), oposta ao armazenamento | `u41y` | −1 |
| \(4_n,2_n\) | z | \(4_n\to2_n\), diagonal descendente | `u42z` | +1 |
| \(3_{n-1},1_n\) | z | \(3_{n-1}\to1_n\), diagonal descendente | `u31z` | +1 |
| \(4_n,3_{n-1}\) | w | \(4_n\to3_{n-1}\), superior para a esquerda | `u43w` | +1 |
| \(1_n,2_{n-1}\) | w | \(2_{n-1}\to1_n\), inferior para a direita | `u12w` | −1 |
| \(1_n,3_n\) | t | \(3_n\to1_n\), oposta ao armazenamento | `u13t` | −1 |
| \(2_{n-1},4_n\) | t | \(4_n\to2_{n-1}\), oposta ao armazenamento | `u24t` | −1 |

Esses dez sinais concordam simultaneamente com as setas e com a Eq.(15), cuja imagem também registra a Fourier usada nas notas:

![Eq.(15), Fourier e início da Eq.(17), página 5](<auditoria_fontes/original-5.png>)

[TESTE INDEPENDENTE]

O desenho anterior localizado em `plot_qsl_graph.py` orientava as horizontais superiores para a direita e, portanto, não era uma reprodução do ORIGINAL. O novo desenho foi reconstruído do PDF: primeiro se estabeleceu a base original, com superiores para a esquerda; a diferença entre essa base e Z contém apenas a reversão de z. Nenhum sinal do setor anterior de horizontais invertidas entra na definição física abaixo.

## 2. Definition of the Z-reversed configuration

[DEDUÇÃO MATEMÁTICA]

“Inverter z” significa inverter as duas famílias de setas reais, incluindo sua repetição em todas as células:

\[
4_n\to2_n\quad\longmapsto\quad2_n\to4_n,
\qquad
3_{n-1}\to1_n\quad\longmapsto\quad1_n\to3_{n-1}.
\]

Mantendo as variáveis armazenadas e os deslocamentos originais, isso equivale exatamente a

\[
\boxed{u_{42,z}:+1\to-1,\qquad u_{31,z}:+1\to-1.}
\]

| Variável | OLD | Z-REVERSED | Alterada? |
|----------|-----|------------|-----------|
| \(u_{12,x}\) | +1 | +1 | não |
| \(u_{34,x}\) | +1 | +1 | não |
| \(u_{23,y}\) | +1 | +1 | não |
| \(u_{41,y}\) | −1 | −1 | não |
| \(u_{42,z}\) | +1 | −1 | sim, z |
| \(u_{31,z}\) | +1 | −1 | sim, z |
| \(u_{43,w}\) | +1 | +1 | não |
| \(u_{12,w}\) | −1 | −1 | não |
| \(u_{13,t}\) | −1 | −1 | não |
| \(u_{24,t}\) | −1 | −1 | não |

Na figura reconstruída, \(b_i=(2i,0)\), \(t_i=(2i+1,2)\), \(\mathbf a_1=(2,0)\), \(\mathbf a_2=(1,2)\), e os rótulos centrais são \(1=b_2\), \(2=b_3\), \(4=t_2\), \(3=t_3\). Somente as setas verdes z apontam no sentido oposto ao ORIGINAL.

![Rede com apenas z invertido, 300 dpi](<qsl_lattice_z_reversed.png>)

[TESTE INDEPENDENTE]

A rotina `assert_only_z_reversed()` verifica explicitamente `u34x == +1` e `u43w == +1` nos dois setores, a inversão de `u42z` e `u31z`, e a igualdade de todas as demais chaves. Qualquer mudança não-z produz `AssertionError`. A rotina de desenho também compara suas setas com o grafo real transladado de ambos os gauges; esse teste evita um desenho visualmente plausível com sinais errados.

## 3. Plaquette fluxes

[DERIVADO DO PDF ORIGINAL]

A definição usada nas notas e neste relatório é

\[
W_p=\prod_{\langle ij\rangle_\gamma\in p}u_{ij,\gamma},
\qquad W_p=e^{i\Phi_p},
\]

com circulação anti-horária. Portanto \(W_p=+1\) corresponde a \(\Phi_p=0\) e \(W_p=-1\) a \(\Phi_p=\pi\), módulo \(2\pi\), **nessa convenção de produto orientado**. Em um triângulo, inverter toda a circulação troca o sinal do produto porque cada um dos três fatores troca de sinal; fixar a circulação é indispensável.

[DEDUÇÃO MATEMÁTICA]

Os quatro triângulos do bloco x estão na mesma célula. Cada sinal abaixo é obtido da ligação efetivamente percorrida, usando \(u_{ji}=-u_{ij}\) quando necessário.

**Caminho \(1_n\to2_n\to4_n\to1_n\):**

\[
\begin{aligned}
W_{124}&=u_{12,x}u_{24,z}u_{41,y}
=u_{12,x}(-u_{42,z})u_{41,y},\\
W_{124}^{\rm OLD}&=(+1)(-1)(-1)=+1,\\
W_{124}^{Z}&=(+1)(+1)(-1)=-1.
\end{aligned}
\]

**Caminho \(2_n\to3_n\to4_n\to2_n\):**

\[
\begin{aligned}
W_{234}&=u_{23,y}u_{34,x}u_{42,z},\\
W_{234}^{\rm OLD}&=(+1)(+1)(+1)=+1,\\
W_{234}^{Z}&=(+1)(+1)(-1)=-1.
\end{aligned}
\]

**Caminho \(1_n\to2_n\to3_n\to1_n\):**

\[
\begin{aligned}
W_{123}&=u_{12,x}u_{23,y}u_{31,t}
=u_{12,x}u_{23,y}(-u_{13,t}),\\
W_{123}^{\rm OLD}&=(+1)(+1)(+1)=+1,\\
W_{123}^{Z}&=(+1)(+1)(+1)=+1.
\end{aligned}
\]

**Caminho \(1_n\to3_n\to4_n\to1_n\):**

\[
\begin{aligned}
W_{134}&=u_{13,t}u_{34,x}u_{41,y},\\
W_{134}^{\rm OLD}&=(-1)(+1)(-1)=+1,\\
W_{134}^{Z}&=(-1)(+1)(-1)=+1.
\end{aligned}
\]

O bloco w imediatamente à esquerda tem vértices \(2_{n-1},1_n,4_n,3_{n-1}\). Os índices de célula não podem ser descartados ao construir suas ligações.

**Caminho \(2_{n-1}\to1_n\to3_{n-1}\to2_{n-1}\):**

\[
\begin{aligned}
W_{213}&=u_{21,w}u_{13,z}u_{32,y}
=(-u_{12,w})(-u_{31,z})(-u_{23,y}),\\
W_{213}^{\rm OLD}&=(+1)(-1)(-1)=+1,\\
W_{213}^{Z}&=(+1)(+1)(-1)=-1.
\end{aligned}
\]

**Caminho \(1_n\to4_n\to3_{n-1}\to1_n\):**

\[
\begin{aligned}
W_{143}&=u_{14,y}u_{43,w}u_{31,z}
=(-u_{41,y})u_{43,w}u_{31,z},\\
W_{143}^{\rm OLD}&=(+1)(+1)(+1)=+1,\\
W_{143}^{Z}&=(+1)(+1)(-1)=-1.
\end{aligned}
\]

**Caminho \(2_{n-1}\to1_n\to4_n\to2_{n-1}\):**

\[
\begin{aligned}
W_{214}&=u_{21,w}u_{14,y}u_{42,t}
=(-u_{12,w})(-u_{41,y})(-u_{24,t}),\\
W_{214}^{\rm OLD}&=(+1)(+1)(+1)=+1,\\
W_{214}^{Z}&=(+1)(+1)(+1)=+1.
\end{aligned}
\]

**Caminho \(2_{n-1}\to4_n\to3_{n-1}\to2_{n-1}\):**

\[
\begin{aligned}
W_{243}&=u_{24,t}u_{43,w}u_{32,y}
=u_{24,t}u_{43,w}(-u_{23,y}),\\
W_{243}^{\rm OLD}&=(-1)(+1)(-1)=+1,\\
W_{243}^{Z}&=(-1)(+1)(-1)=+1.
\end{aligned}
\]

| Plaqueta | Tipos percorridos | OLD | Z-REVERSED | Mudou? |
|----------|-------------------|-----|------------|--------|
| \(W_{124}\) | x,z,y | +1 | −1 | sim |
| \(W_{234}\) | y,x,z | +1 | −1 | sim |
| \(W_{123}\) | x,y,t | +1 | +1 | não |
| \(W_{134}\) | t,x,y | +1 | +1 | não |
| \(W_{213}\) | w,z,y | +1 | −1 | sim |
| \(W_{143}\) | y,w,z | +1 | −1 | sim |
| \(W_{214}\) | w,y,t | +1 | +1 | não |
| \(W_{243}\) | t,w,y | +1 | +1 | não |

Assim, em cada bloco x ou w, os dois triângulos que contêm z mudam de fluxo; os dois que contêm t conservam o fluxo. Os triângulos se sobrepõem como ciclos do grafo; não se pressupõe que todos os oito sejam invariantes independentes.

[TESTE INDEPENDENTE] [RESULTADO NUMÉRICO]

O código localiza cada aresta do caminho na lista real, aplica a orientação de percurso e multiplica os fatores. Os valores de fluxo não são usados para gerar os próprios produtos. Os oito caminhos fecham; na escala \(\mathbf a_1=(2,0),\mathbf a_2=(1,2)\), todos têm área orientada \(+2\), confirmando circulação anti-horária. O arquivo `geometry_flux_audit.json` preserva caminhos, fatores, produtos e resultados.

## 4. Gauge equivalence

[DEDUÇÃO MATEMÁTICA]

Sob a transformação local

\[
\theta_i\mapsto s_i\theta_i,\qquad
u_{ij}\mapsto s_i u_{ij}s_j,\qquad s_i\in\{-1,+1\},
\]

cada sinal local aparece duas vezes em um ciclo fechado. Escrevendo diretamente em termos de \(u\),

\[
W_p'=(s_i u_{ij}s_j)(s_j u_{jk}s_k)(s_k u_{ki}s_i)
=s_i^2s_j^2s_k^2W_p=W_p.
\]

Como \(W_{124}^{\rm OLD}=+1\) e \(W_{124}^{Z}=-1\),

\[
\boxed{\text{OLD e Z-REVERSED não são gauge-equivalentes.}}
\]

Esse argumento vale inclusive para sinais \(s_{m,n}\) dependentes da célula e para células maiores. Não se limita à busca de quatro sinais.

Para a transformação periódica de quatro sub-redes,

\[
G=\operatorname{diag}(s_1,s_2,s_3,s_4),
\]

preservar as ligações x exige \(s_1s_2=+1\) e \(s_3s_4=+1\); preservar y exige \(s_2s_3=+1\) e \(s_4s_1=+1\). Logo,

\[
s_1=s_2=s_3=s_4,\qquad
s_4s_2=s_3s_1=+1.
\]

Entretanto, inverter z exigiria \(s_4s_2=s_3s_1=-1\). As exigências são incompatíveis. As ligações t e w conservadas apenas acrescentam restrições já satisfeitas pelos sinais todos iguais.

[TESTE INDEPENDENTE] [RESULTADO NUMÉRICO]

A busca exaustiva percorreu os \(2^4=16\) elementos de \(\{-1,+1\}^4\) e verificou todas as dez ligações. Número de soluções OLD \(\to\) Z: **zero**. Setas diferentes poderiam representar o mesmo setor em outro problema; aqui os fluxos demonstram setores diferentes. Uma eventual igualdade energética não desfaz essa distinção.

## 5. Real-space Hamiltonian

[DERIVADO DO PDF ORIGINAL]

Partimos de

\[
H=i\sum_{\gamma\in\{x,y,z,t,w\}}
\sum_{\langle ij\rangle_\gamma}K_\gamma u_{ij,\gamma}\theta_i^y\theta_j^y,
\qquad K_\gamma\in\mathbb R.
\]

Para tornar as expressões legíveis, nesta seção \(\theta_m(\mathbf r)\equiv\theta_m^y(\mathbf r)\). A soma em \(\mathbf R\) percorre as origens das células originais. Inserindo os sinais confirmados na Eq.(11), os **dez termos OLD** são

\[
\begin{aligned}
H_{\rm OLD}=i\sum_{\mathbf R}\Big[&
 K_x\theta_1(\mathbf R)\theta_2(\mathbf R+\mathbf a_1)\\
&+K_y\theta_2(\mathbf R+\mathbf a_1)\theta_3(\mathbf R+\mathbf a_1+\mathbf a_2)\\
&+K_x\theta_3(\mathbf R+\mathbf a_1+\mathbf a_2)\theta_4(\mathbf R+\mathbf a_2)\\
&-K_y\theta_4(\mathbf R+\mathbf a_2)\theta_1(\mathbf R)\\
&+\underbrace{K_z\theta_4(\mathbf R+\mathbf a_2)\theta_2(\mathbf R+\mathbf a_1)}_{z:\,u_{42,z}=+1}\\
&+K_w\theta_4(\mathbf R+\mathbf a_2)\theta_3(\mathbf R+\mathbf a_2-\mathbf a_1)\\
&+\underbrace{K_z\theta_3(\mathbf R+\mathbf a_2-\mathbf a_1)\theta_1(\mathbf R)}_{z:\,u_{31,z}=+1}\\
&-K_w\theta_1(\mathbf R)\theta_2(\mathbf R-\mathbf a_1)\\
&-K_t\theta_1(\mathbf R)\theta_3(\mathbf R+\mathbf a_1+\mathbf a_2)\\
&-K_t\theta_2(\mathbf R-\mathbf a_1)\theta_4(\mathbf R+\mathbf a_2)
\Big].
\end{aligned}
\]

[DEDUÇÃO MATEMÁTICA]

Mantendo a ordem dos operadores e todos os seus argumentos, os **dez termos Z-REVERSED** são

\[
\begin{aligned}
H_Z=i\sum_{\mathbf R}\Big[&
 K_x\theta_1(\mathbf R)\theta_2(\mathbf R+\mathbf a_1)\\
&+K_y\theta_2(\mathbf R+\mathbf a_1)\theta_3(\mathbf R+\mathbf a_1+\mathbf a_2)\\
&+K_x\theta_3(\mathbf R+\mathbf a_1+\mathbf a_2)\theta_4(\mathbf R+\mathbf a_2)\\
&-K_y\theta_4(\mathbf R+\mathbf a_2)\theta_1(\mathbf R)\\
&-\underbrace{K_z\theta_4(\mathbf R+\mathbf a_2)\theta_2(\mathbf R+\mathbf a_1)}_{z:\,u_{42,z}=-1}\\
&+K_w\theta_4(\mathbf R+\mathbf a_2)\theta_3(\mathbf R+\mathbf a_2-\mathbf a_1)\\
&-\underbrace{K_z\theta_3(\mathbf R+\mathbf a_2-\mathbf a_1)\theta_1(\mathbf R)}_{z:\,u_{31,z}=-1}\\
&-K_w\theta_1(\mathbf R)\theta_2(\mathbf R-\mathbf a_1)\\
&-K_t\theta_1(\mathbf R)\theta_3(\mathbf R+\mathbf a_1+\mathbf a_2)\\
&-K_t\theta_2(\mathbf R-\mathbf a_1)\theta_4(\mathbf R+\mathbf a_2)
\Big].
\end{aligned}
\]

A diferença é, portanto,

\[
\boxed{H_Z-H_{\rm OLD}=-2iK_z\sum_{\mathbf R}
\left[
\theta_4(\mathbf R+\mathbf a_2)\theta_2(\mathbf R+\mathbf a_1)
+\theta_3(\mathbf R+\mathbf a_2-\mathbf a_1)\theta_1(\mathbf R)
\right].}
\]

Confirma-se a hipótese: mudam diretamente apenas os termos associados a \(u_{42,z}\) e \(u_{31,z}\). Reordenar operadores em etapas seguintes pode acrescentar sinais de anticomutação; esses sinais não são novas mudanças físicas de gauge.

A lista independente que alimenta a reconstrução computacional contém os deslocamentos **fim menos início**, \(\mathbf d=da_1\mathbf a_1+da_2\mathbf a_2\):

| início | fim | tipo | variável | \(u^{\rm OLD}\) | \(u^Z\) | \(da_1\) | \(da_2\) |
|--------|-----|------|----------|----------------|----------|------------|------------|
| 1 | 2 | x | `u12x` | +1 | +1 | +1 | 0 |
| 2 | 3 | y | `u23y` | +1 | +1 | 0 | +1 |
| 3 | 4 | x | `u34x` | +1 | +1 | −1 | 0 |
| 4 | 1 | y | `u41y` | −1 | −1 | 0 | −1 |
| 4 | 2 | z | `u42z` | +1 | −1 | +1 | −1 |
| 4 | 3 | w | `u43w` | +1 | +1 | −1 | 0 |
| 3 | 1 | z | `u31z` | +1 | −1 | +1 | −1 |
| 1 | 2 | w | `u12w` | −1 | −1 | −1 | 0 |
| 1 | 3 | t | `u13t` | −1 | −1 | +1 | +1 |
| 2 | 4 | t | `u24t` | −1 | −1 | +1 | +1 |

## 6. Majorana Fourier transform

[DERIVADO DO PDF ORIGINAL]

A Fourier é mantida na forma original, sem mudar o sinal do expoente nem seu prefator:

\[
\boxed{\theta_m^y(\mathbf r)=\sqrt{\frac4N}
\sum_{\mathbf k\in\mathrm{BZ}/2}
\left[e^{i\mathbf k\cdot\mathbf r}\theta_m^y(\mathbf k)
+e^{-i\mathbf k\cdot\mathbf r}\theta_m^{y\dagger}(\mathbf k)\right].}
\]

Definimos \(q_1=\mathbf k\cdot\mathbf a_1\), \(q_2=\mathbf k\cdot\mathbf a_2\). As notas restringem \(\mathbf k=k\hat{\mathbf x}\). Conservamos inicialmente \(q_1,q_2\) na álgebra para auditar as fases antes da Eq.(27).

A Eq.(16) das notas registra, para o termo orientado escrito com deslocamento no primeiro operador,

\[
iK u\sum_{\mathbf R}\theta_m(\mathbf R+\boldsymbol\delta)\theta_n(\mathbf R)
=2\sum_{\mathbf k\in\mathrm{BZ}/2}
\left[iK u e^{-i\mathbf k\cdot\boldsymbol\delta}\theta_m^\dagger\theta_n
-iK u e^{+i\mathbf k\cdot\boldsymbol\delta}\theta_n^\dagger\theta_m\right].
\]

O fator 2 acima pertence à escrita em meia BZ. Na Eq.(17), o documento passa à BZ completa e escreve os elementos matriciais com prefator \(iK\), como nas Eqs.(21)–(26). Adotamos essa normalização documental de bandas e a normalização energética especificada no pedido; não multiplicamos as matrizes finais por 2.

[DEDUÇÃO MATEMÁTICA]

A origem das fases pode ser auditada termo a termo. Para uma ligação \(i\) em \(\mathbf r_i\) e \(j\) em \(\mathbf r_j\), o produto que acompanha \(\theta_i^\dagger(\mathbf k)\theta_j(\mathbf k)\) contém

\[
e^{-i\mathbf k\cdot\mathbf r_i}e^{+i\mathbf k\cdot\mathbf r_j}
=e^{+i\mathbf k\cdot(\mathbf r_j-\mathbf r_i)}
=e^{+i\mathbf k\cdot\mathbf d}.
\]

O produto de operadores na ordem oposta satisfaz, para sub-redes diferentes,

\[
\theta_i\theta_j^\dagger=-\theta_j^\dagger\theta_i,
\]

e produz o conjugado hermitiano. Logo, na convenção da matriz em BZ completa das notas,

\[
\boxed{(H_{\rm raw})_{ij}\mathrel{+}=iK_\gamma u_{ij,\gamma}
 e^{i(q_1da_1+q_2da_2)},\quad
(H_{\rm raw})_{ji}\mathrel{+}=-iK_\gamma u_{ij,\gamma}
 e^{-i(q_1da_1+q_2da_2)}.}
\]

Para compreender a extensão da meia BZ, chame o bilinear entre colchetes de \(B(\mathbf k)\). A relação Majorana \(\Theta(-\mathbf k)=\Theta^\dagger(\mathbf k)\), com os índices tratados componente a componente, e a reordenação acima fornecem \(B(-\mathbf k)=B(\mathbf k)\). Assim, a escrita documental \(2\sum_{\mathrm{BZ}/2}B(\mathbf k)\) corresponde a \(\sum_{\mathrm{BZ}}B(\mathbf k)\). Momentos autoconjugados de uma cadeia finita são tratados sem duplicação pela matriz real; são de medida nula na integral contínua. Esse argumento explica a contagem de momentos, sem redefinir a Fourier ou a escala de energia pedidas.

Para acompanhar ambos os setores, introduzimos apenas o sinal auxiliar

\[
\sigma=u_{31,z}=u_{42,z},\qquad
\sigma=+1\ (\mathrm{OLD}),\qquad\sigma=-1\ (Z).
\]

Essa abreviação não substitui as variáveis \(u_{ij,\gamma}\) do gauge. As duas ligações z produzem:

\[
\begin{aligned}
\mathbf d_{31,z}&=\mathbf a_1-\mathbf a_2,
&(H_{\rm raw})_{31}^{(z)}&=i\sigma K_z e^{i(q_1-q_2)},\\
\mathbf d_{42,z}&=\mathbf a_1-\mathbf a_2,
&(H_{\rm raw})_{42}^{(z)}&=i\sigma K_z e^{i(q_1-q_2)}.
\end{aligned}
\]

Ao colocá-las nas posições superiores \(13\) e \(24\),

\[
\boxed{(H_{\rm raw})_{13}^{(z)}=(H_{\rm raw})_{24}^{(z)}
=-i\sigma K_z e^{-iq_1+iq_2}.}
\]

Portanto,

\[
f_{13,\rm raw}^{(z),\,\mathrm{OLD}}
=f_{24,\rm raw}^{(z),\,\mathrm{OLD}}
=-iK_z e^{-iq_1+iq_2},
\qquad
f_{13,\rm raw}^{(z),\,Z}
=f_{24,\rm raw}^{(z),\,Z}
=+iK_z e^{-iq_1+iq_2}.
\]

As ligações t têm \(u_{13,t}=u_{24,t}=-1\) e deslocamento \(\mathbf a_1+\mathbf a_2\), de modo que

\[
f_{13,\rm raw}^{(t)}=f_{24,\rm raw}^{(t)}
=-iK_t e^{+iq_1+iq_2}.
\]

A tabela seguinte mostra a soma de **todas** as contribuições superiores antes da Eq.(27):

| Elemento | Contribuição real 1 | Contribuição real 2 | Soma RAW para sinal \(\sigma\) |
|----------|--------------------|--------------------|--------------------------------|
| \(12\) | \(u_{12,x}:\ iK_xe^{iq_1}\) | \(u_{12,w}:\ -iK_we^{-iq_1}\) | \(i(K_xe^{iq_1}-K_we^{-iq_1})\) |
| \(13\) | \(u_{13,t}:\ -iK_te^{i(q_1+q_2)}\) | \(u_{31,z}\), invertida para \(13\): \(-i\sigma K_ze^{-iq_1+iq_2}\) | \(-ie^{iq_2}(K_te^{iq_1}+\sigma K_ze^{-iq_1})\) |
| \(14\) | \(u_{41,y}\), invertida para \(14\): \(+iK_ye^{iq_2}\) | — | \(+iK_ye^{iq_2}\) |
| \(23\) | \(u_{23,y}:\ +iK_ye^{iq_2}\) | — | \(+iK_ye^{iq_2}\) |
| \(24\) | \(u_{24,t}:\ -iK_te^{i(q_1+q_2)}\) | \(u_{42,z}\), invertida para \(24\): \(-i\sigma K_ze^{-iq_1+iq_2}\) | \(-ie^{iq_2}(K_te^{iq_1}+\sigma K_ze^{-iq_1})\) |
| \(34\) | \(u_{34,x}:\ +iK_xe^{-iq_1}\) | \(u_{43,w}\), invertida para \(34\): \(-iK_we^{iq_1}\) | \(i(K_xe^{-iq_1}-K_we^{iq_1})\) |

Em particular, os **seis coeficientes RAW Z-REVERSED** são

\[
\boxed{\begin{aligned}
f_{12,\rm raw}^{Z}&=i(K_xe^{iq_1}-K_we^{-iq_1}),\\
f_{13,\rm raw}^{Z}&=-ie^{iq_2}(K_te^{iq_1}-K_ze^{-iq_1}),\\
f_{14,\rm raw}^{Z}&=iK_ye^{iq_2},\\
f_{23,\rm raw}^{Z}&=iK_ye^{iq_2},\\
f_{24,\rm raw}^{Z}&=-ie^{iq_2}(K_te^{iq_1}-K_ze^{-iq_1}),\\
f_{34,\rm raw}^{Z}&=i(K_xe^{-iq_1}-K_we^{iq_1}).
\end{aligned}}
\]

Todos os quatro elementos que conectam cadeias têm o fator **\(e^{+iq_2}\)**. Seu sinal decorre dos deslocamentos reais, inclusive quando uma ligação z armazenada aponta da cadeia superior para a inferior. O sinal negativo do deslocamento transversal nessa orientação é conjugado ao preencher a entrada superior da matriz.

## 7. Eq.(27)

[DERIVADO DO PDF ORIGINAL]

A transformação impressa é a substituição

\[
\Theta\longmapsto U\Theta,
\qquad
\Theta=(\theta_1^y,\theta_2^y,\theta_3^y,\theta_4^y)^T,
\qquad
\boxed{U=\operatorname{diag}(1,1,e^{-iq_2},e^{-iq_2}).}
\]

[DEDUÇÃO MATEMÁTICA]

Para \(q_2\in\mathbb R\),

\[
U^\dagger=\operatorname{diag}(1,1,e^{iq_2},e^{iq_2}),
\qquad
U^\dagger U=\operatorname{diag}(1,1,1,1)=I_4.
\]

A direção correta da transformação matricial segue da **substituição ativa dentro da forma quadrática**:

\[
\mathcal H_{\rm raw}[\Theta]=\Theta^\dagger H_{\rm raw}\Theta,
\qquad
\mathcal H_{\rm raw}[U\Theta]
=(U\Theta)^\dagger H_{\rm raw}(U\Theta)
=\Theta^\dagger\underbrace{(U^\dagger H_{\rm raw}U)}_{H_{\rm Bloch}}\Theta.
\]

Portanto,

\[
\boxed{H_{\rm Bloch}=U^\dagger H_{\rm raw}U.}
\]

Se o objetivo for escrever o **mesmo operador em coordenadas novas**, a identificação equivalente é

\[
\boxed{\Theta_{\rm raw}=U\Theta_{\rm Bloch},\qquad
\Theta_{\rm Bloch}=U^\dagger\Theta_{\rm raw}.}
\]

Definir, em vez disso, \(\Theta'=U\Theta_{\rm raw}\) e reescrever o mesmo operador em \(\Theta'\) daria \(UH_{\rm raw}U^\dagger\). Por isso, essa segunda identificação não deve ser confundida com a substituição ativa acima. A Eq.(27), sua matriz U e o cancelamento exigido são preservados.

Para os quatro elementos cross-chain, \(a\in\{1,2\}\), \(b\in\{3,4\}\), escrevemos \((H_{\rm raw})_{ab}=e^{iq_2}c_{ab}(q_1)\). Então

\[
(H_{\rm Bloch})_{ab}
=(U^\dagger)_{aa}(H_{\rm raw})_{ab}U_{bb}
=1\cdot e^{iq_2}c_{ab}\cdot e^{-iq_2}=c_{ab}.
\]

Explicitamente,

\[
\begin{aligned}
f_{13}^{Z}&=e^{iq_2}\big[-i(K_te^{iq_1}-K_ze^{-iq_1})\big]e^{-iq_2}
=-i(K_te^{iq_1}-K_ze^{-iq_1}),\\
f_{14}^{Z}&=e^{iq_2}(iK_y)e^{-iq_2}=iK_y,\\
f_{23}^{Z}&=e^{iq_2}(iK_y)e^{-iq_2}=iK_y,\\
f_{24}^{Z}&=e^{iq_2}\big[-i(K_te^{iq_1}-K_ze^{-iq_1})\big]e^{-iq_2}
=-i(K_te^{iq_1}-K_ze^{-iq_1}).
\end{aligned}
\]

As entradas inferiores são os conjugados. Para as duas entradas da mesma cadeia,

\[
f_{12}^{Z}=1\cdot f_{12,\rm raw}^{Z}\cdot1,
\qquad
f_{34}^{Z}=e^{iq_2}f_{34,\rm raw}^{Z}e^{-iq_2}=f_{34,\rm raw}^{Z}.
\]

Assim \(q_2\) desaparece de **toda** a matriz. Daqui em diante \(q=q_1=\mathbf k\cdot\mathbf a_1\).

[TESTE INDEPENDENTE]

A unitariedade e o cancelamento de \(q_2\) foram verificados simbolicamente. A matriz RAW construída pela lista de ligações reais foi também transformada por U numericamente e comparada às fórmulas fechadas, para ambos os setores. A fórmula fechada não é a fonte de geometria da rotina física `h_bloch`.

## 8. Final Bloch Hamiltonian

[DEDUÇÃO MATEMÁTICA]

A matriz final Z é

\[
\boxed{H_Z(q)=
\begin{pmatrix}
0&f_{12}^{Z}&f_{13}^{Z}&f_{14}^{Z}\\
(f_{12}^{Z})^*&0&f_{23}^{Z}&f_{24}^{Z}\\
(f_{13}^{Z})^*&(f_{23}^{Z})^*&0&f_{34}^{Z}\\
(f_{14}^{Z})^*&(f_{24}^{Z})^*&(f_{34}^{Z})^*&0
\end{pmatrix}.}
\]

Os coeficientes foram obtidos das dez ligações reais, e não por substituição presumida em uma matriz final:

| coefficient | OLD | Z-INVERTED | changed? | reason |
|-------------|-----|------------|----------|--------|
| \(f_{12}\) | \(i(K_xe^{iq}-K_we^{-iq})\) | \(i(K_xe^{iq}-K_we^{-iq})\) | não | x: \(1_n\to2_n\); w armazenado: \(1_n\to2_{n-1}\); nenhum z |
| \(f_{13}\) | \(-i(K_te^{iq}+K_ze^{-iq})\) | \(-i(K_te^{iq}-K_ze^{-iq})\) | sim | t: \(1_n,3_n\) conservado; z: \(3_{n-1},1_n\) invertido e reordenado |
| \(f_{14}\) | \(iK_y\) | \(iK_y\) | não | y: \(4_n,1_n\), entrada superior obtida por conjugação |
| \(f_{23}\) | \(iK_y\) | \(iK_y\) | não | y: \(2_n,3_n\) conservado |
| \(f_{24}\) | \(-i(K_te^{iq}+K_ze^{-iq})\) | \(-i(K_te^{iq}-K_ze^{-iq})\) | sim | t: \(2_{n-1},4_n\) conservado; z: \(4_n,2_n\) invertido e reordenado |
| \(f_{34}\) | \(i(K_xe^{-iq}-K_we^{iq})\) | \(i(K_xe^{-iq}-K_we^{iq})\) | não | x: \(3_n,4_n\); w: \(4_n,3_{n-1}\), reordenado para entrada 34 |

Em particular,

\[
\boxed{f_{13}^{Z}-f_{13}^{\rm OLD}
=f_{24}^{Z}-f_{24}^{\rm OLD}=2iK_ze^{-iq},}
\]

com todas as outras entradas superiores idênticas. As entradas inferiores variam pelo conjugado dessas diferenças.

Para explicitar todos os elementos, defina

\[
a=i(K_xe^{iq}-K_we^{-iq}),\qquad
b_Z=-i(K_te^{iq}-K_ze^{-iq}),\qquad y=K_y.
\]

A expressão anterior equivale a

\[
\boxed{H_Z(q)=
\begin{pmatrix}
0&a&b_Z&iy\\
a^*&0&iy&b_Z\\
b_Z^*&-iy&0&-a^*\\
-iy&b_Z^*&-a&0
\end{pmatrix}.}
\]

Aqui a relação \(f_{34}=-f_{12}^*\) foi calculada das expressões: \(-a^*=i(K_xe^{-iq}-K_we^{iq})\). O estudo das consequências espectrais aparece nas seções seguintes.

[DERIVADO DO PDF ORIGINAL — TRANSCRIÇÃO LITERAL]

A referência OLD obrigatória preserva a Eq.(26) como impressa:

\[
\boxed{f_{34}^{\mathrm{Eq.(26)\ literal}}
=i(K_xe^{+iq}-K_we^{-iq}).}
\]

Essa transcrição não foi alterada no documento-fonte. Seu status é documental, não de identidade matemática com a reconstrução geométrica.

[DEDUÇÃO MATEMÁTICA — ERRATA DOCUMENTADA]

A ligação superior x armazenada \(3\to4\) tem deslocamento \(-\mathbf a_1\), e por isso contribui \(+iK_xe^{-iq}\). A ligação superior w armazenada \(4\to3\) tem deslocamento \(-\mathbf a_1\); ao preencher a entrada 34, sua contribuição é conjugada: \(-iK_we^{+iq}\). Portanto,

\[
f_{34}^{\rm geometria}=i(K_xe^{-iq}-K_we^{+iq}),
\qquad
f_{34}^{\rm geometria}-f_{34}^{\mathrm{Eq.(26)\ literal}}
=2(K_x+K_w)\sin q.
\]

A Eq.(17) concorda com a expressão geométrica. A Eq.(27) não pode resolver a discrepância, pois suas fases nas sub-redes 3 e 4 são iguais e deixam \(f_{34}\) inalterado. A igualdade com a expressão literal ocorre apenas em pontos ou acoplamentos especiais, e não serve como identidade geral.

[TESTE INDEPENDENTE]

A auditoria preliminar já encontrou erro máximo \(2{,}48\times10^{-16}\) entre a projeção de uma cadeia ORIGINAL de 28 sítios e a expressão geométrica para \(f_{34}\), enquanto a expressão literal diferiu em aproximadamente \(3{,}51\) para o mesmo conjunto assimétrico de acoplamentos. A validação completa posterior reconstrói todas as entradas para OLD e Z, conforme a seção 15. Desse modo, a correção empregada é sustentada pela geometria e por teste independente; não é escolhida para obter um resultado energético desejado.

## 9. Symmetries

[DEDUÇÃO MATEMÁTICA] Todos os acoplamentos desta análise são reais. Os seis coeficientes acima determinam a parte inferior da matriz por conjugação, e portanto

\[
H_\sigma(q)^\dagger=H_\sigma(q),\qquad \sigma=+1\ (\mathrm{OLD}),\quad \sigma=-1\ (Z).
\]

Cada ligação real contribui com \(iK_\gamma u_{ij,\gamma}e^{i\mathbf k\cdot\mathbf d}\), com \(K_\gamma u_{ij,\gamma}\) real. Logo,

\[
f_{mn}(-q)=-f_{mn}(q)^*,\qquad
\boxed{H_\sigma(-q)=-H_\sigma(q)^*.}
\]

Esta é a relação Majorana, verificada nas expressões novas. Ela, isoladamente, relaciona o espectro em \(q\) ao de \(-q\); o pareamento no mesmo \(q\) exige uma propriedade adicional, demonstrada abaixo.

[DEDUÇÃO MATEMÁTICA] As relações especiais que sobrevivem em ambos os setores são

\[
f_{23}=f_{14}=iK_y,\qquad f_{24}=f_{13},\qquad f_{34}=-f_{12}^*.
\]

Escrevendo \(a=f_{12}\), \(b=f_{13}\), a matriz tem a forma

\[
H_\sigma=\begin{pmatrix}
0&a&b&iK_y\\
a^*&0&iK_y&b\\
b^*&-iK_y&0&-a^*\\
-iK_y&b^*&-a&0
\end{pmatrix}.
\]

Para

\[
J=\begin{pmatrix}
0&0&-1&0\\0&0&0&-1\\1&0&0&0\\0&1&0&0
\end{pmatrix},\qquad J^\dagger J=I,\quad J^2=-I,
\]

a multiplicação explícita fornece

\[
\boxed{JH_\sigma(q)^*J^\dagger=-H_\sigma(q).}
\]

Assim, se \(H_\sigma v=Ev\), então \(H_\sigma(Jv^*)=-E(Jv^*)\) no mesmo \(q\). Essa relação foi demonstrada para a matriz dos dois setores; não foi importada do caso com horizontais superiores invertidas. Ela não implica, por si só, degenerescência entre estados de mesmo sinal de energia.

[DEDUÇÃO MATEMÁTICA] A célula física tem translação \(\mathbf T=2\mathbf a_1\). As coordenadas internas em unidades de \(\mathbf a_1\) são \(r=(0,1,1,0)\). Definindo

\[
D_\pi=\operatorname{diag}(1,-1,-1,1),
\]

os coeficientes que contêm \(e^{\pm iq}\) mudam de sinal quando \(q\to q+\pi\), enquanto \(f_{14},f_{23}\) permanecem iguais. Portanto,

\[
\boxed{H_\sigma(q+\pi)=D_\pi H_\sigma(q)D_\pi^\dagger,}
\qquad H_\sigma(q+2\pi)=H_\sigma(q).
\]

A matriz é periódica em \(2\pi\) nessa base, e o espectro é periódico em \(\pi\). A BZ física \([-\pi/2,\pi/2]\) decorre de \(\mathbf T\), sem mudança de célula. Um período menor pode aparecer em parâmetros especiais, sem alterar essa BZ de referência.

[TESTE INDEPENDENTE] Hermiticidade, Majorana, periodicidade matricial, covariância em \(\pi\), relações especiais e antiunitária foram verificadas simbolicamente. As relações matriciais e espectrais pertinentes também foram testadas nos 1000 casos aleatórios.

## 10. Characteristic polynomial

[DEDUÇÃO MATEMÁTICA] Para tornar a expansão legível, escreva

\[
a=a_R+ia_I,\qquad b=b_R+ib_I,
\]

\[
\begin{aligned}
a_R&=-(K_x+K_w)\sin q,& a_I&=(K_x-K_w)\cos q,\\
b_R&=(K_t-\sigma K_z)\sin q,& b_I&=-(K_t+\sigma K_z)\cos q.
\end{aligned}
\]

Essas quatro quantidades são reais. Expandindo \(\det(\lambda I-H_\sigma)\) diretamente com SymPy, obtém-se

\[
\begin{aligned}
P_\sigma={}&\lambda^4
-2(a_R^2+a_I^2+b_R^2+b_I^2+K_y^2)\lambda^2\\
&+a_R^4+a_I^4+b_R^4+b_I^4+K_y^4\\
&+2a_R^2a_I^2+2a_R^2b_R^2+2a_R^2b_I^2+2a_R^2K_y^2\\
&-2a_I^2b_R^2-2a_I^2b_I^2+2a_I^2K_y^2\\
&+2b_R^2b_I^2+2b_R^2K_y^2-2b_I^2K_y^2.
\end{aligned}
\]

Defina

\[
S_\sigma=a_R^2+a_I^2+b_R^2+b_I^2+K_y^2,
\qquad R_\sigma=a_I^2(b_R^2+b_I^2)+b_I^2K_y^2.
\]

O mesmo resultado se simplifica exatamente para

\[
\boxed{P_\sigma(\lambda,q)=\lambda^4-2S_\sigma\lambda^2+S_\sigma^2-4R_\sigma.}
\]

Em acoplamentos originais, sem componentes auxiliares indefinidas,

\[
S_\sigma=K_x^2+K_y^2+K_z^2+K_t^2+K_w^2
+2(\sigma K_tK_z-K_xK_w)\cos2q,
\]

\[
R_\sigma=\cos^2q\left[(K_x-K_w)^2
\big(K_t^2+K_z^2+2\sigma K_tK_z\cos2q\big)
+K_y^2(K_t+\sigma K_z)^2\right].
\]

Em particular, o polinômio **Z** completo é

\[
\boxed{\begin{aligned}
P_Z={}&\lambda^4
-2\left[\sum_{\gamma=x,y,z,t,w}K_\gamma^2
-2(K_tK_z+K_xK_w)\cos2q\right]\lambda^2\\
&+\left[\sum_{\gamma=x,y,z,t,w}K_\gamma^2
-2(K_tK_z+K_xK_w)\cos2q\right]^2\\
&-4\cos^2q\left[(K_x-K_w)^2(K_t^2+K_z^2-2K_tK_z\cos2q)
+K_y^2(K_t-K_z)^2\right].
\end{aligned}}
\]

Não há termo cúbico nem **termo linear em \(\lambda\)**. O pareamento \(\pm E\) no mesmo \(q\) resulta do determinante e da relação antiunitária, inclusive para acoplamentos distintos. Resolvendo a quadrática em \(\lambda^2\),

\[
\boxed{E=\pm\sqrt{S_\sigma+2\sqrt{R_\sigma}},\qquad
E=\pm\sqrt{S_\sigma-2\sqrt{R_\sigma}}.}
\]

A forma \(R_\sigma=a_I^2(b_R^2+b_I^2)+b_I^2K_y^2\) demonstra sua não negatividade. A hermiticidade assegura a não negatividade das duas soluções para \(E^2\). Para o cálculo de energia, o código diagonaliza a matriz com `eigvalsh`, sem impor sinais ou arredondar radicandos destas fórmulas fechadas.

[TESTE INDEPENDENTE] O polinômio foi calculado por `Matrix.charpoly` e sua diferença para a expressão compacta foi expandida a zero, simbolicamente. O programa também verifica sua compatibilidade com os coeficientes numéricos obtidos diretamente da matriz em ambos os setores. O log registra a expansão realmente produzida por SymPy.

## 11. Isotropic spectrum

[DEDUÇÃO MATEMÁTICA] Para \(K_x=K_y=K_z=K_t=K_w=1\), ponha \(s=\sin q\), \(c=\cos q\). As matrizes efetivamente obtidas são

\[
H_{\rm OLD}(q)=\begin{pmatrix}
0&-2s&-2ic&i\\
-2s&0&i&-2ic\\
2ic&-i&0&2s\\
-i&2ic&2s&0
\end{pmatrix},
\]

\[
\boxed{H_Z(q)=\begin{pmatrix}
0&-2s&2s&i\\
-2s&0&i&2s\\
2s&-i&0&2s\\
-i&2s&2s&0
\end{pmatrix}.}
\]

Em OLD,

\[
S_O=5,\quad R_O=4\cos^2q,\qquad
P_O=\lambda^4-10\lambda^2+25-16\cos^2q.
\]

Como \(\cos q\ge0\) na BZ física,

\[
E_O(q)=\left(-\sqrt{5+4\cos q},-\sqrt{5-4\cos q},
\sqrt{5-4\cos q},\sqrt{5+4\cos q}\right).
\]

Em Z,

\[
S_Z=1+8\sin^2q,\qquad R_Z=0,\qquad
P_Z=[\lambda^2-(1+8\sin^2q)]^2,
\]

\[
\boxed{E_Z(q)=\left(-\sqrt{1+8\sin^2q},-\sqrt{1+8\sin^2q},
\sqrt{1+8\sin^2q},\sqrt{1+8\sin^2q}\right).}
\]

Cada energia de Z tem multiplicidade dois em todo \(q\). A verificação adicional

\[
H_Z(q)^2=(1+8\sin^2q)I
\]

foi feita simbolicamente no script principal. Em OLD, as bandas de mesmo sinal se encontram nas extremidades da BZ; não são degeneradas no interior.

| q | OLD: autovalores crescentes exatos | Z: autovalores crescentes exatos |
|---|----------------------------------|--------------------------------|
| \(0\) | \((-3,-1,1,3)\) | \((-1,-1,1,1)\) |
| \(\pm\pi/4\) | \((-\sqrt{5+2\sqrt2},-\sqrt{5-2\sqrt2},\sqrt{5-2\sqrt2},\sqrt{5+2\sqrt2})\) | \((-\sqrt5,-\sqrt5,\sqrt5,\sqrt5)\) |
| \(\pm\pi/2\) | \((-\sqrt5,-\sqrt5,\sqrt5,\sqrt5)\) | \((-3,-3,3,3)\) |

[RESULTADO NUMÉRICO] A diagonalização direta fornece:

| q | OLD (numérico) | Z-REVERSED (numérico) |
|---|----------------|----------------------|
| -pi/2 | (-2.236068, -2.236068, 2.236068, 2.236068) | (-3.000000, -3.000000, 3.000000, 3.000000) |
| -pi/4 | (-2.797933, -1.473626, 1.473626, 2.797933) | (-2.236068, -2.236068, 2.236068, 2.236068) |
| 0 | (-3.000000, -1.000000, 1.000000, 3.000000) | (-1.000000, -1.000000, 1.000000, 1.000000) |
| pi/4 | (-2.797933, -1.473626, 1.473626, 2.797933) | (-2.236068, -2.236068, 2.236068, 2.236068) |
| pi/2 | (-2.236068, -2.236068, 2.236068, 2.236068) | (-3.000000, -3.000000, 3.000000, 3.000000) |

As pequenas diferenças entre \(+q\) e \(-q\) nos arquivos brutos são da ordem do arredondamento. As fórmulas isotrópicas foram comparadas com a matriz geométrica nos 2001 pontos usados para visualização, com erro máximo de \(5.33\times10^{-15}\).

## 12. Band dispersions

[RESULTADO NUMÉRICO]

![Dispersões OLD e Z, mesma escala vertical e BZ física](<dispersion_old_vs_z_reversed.png>)

O painel esquerdo é **Original gauge (OLD)**; o direito é **z bonds reversed (Z-REVERSED)**. A linha horizontal marca \(E=0\). A figura usa 300 dpi e os mesmos limites verticais. Em Z, as quatro bandas aparecem sobrepostas aos pares devido à degenerescência exata; nenhuma banda foi removida.

[TESTE INDEPENDENTE] As cores acompanham autovetores pelo algoritmo

\[
\pi_* =\underset{\pi\in S_4}{\operatorname{argmax}}
\sum_{m=1}^4\left|\langle v_m(q_i),v_{\pi(m)}(q_{i+1})\rangle\right|^2.
\]

As 24 permutações são testadas. Dentro de subespaços degenerados, usa-se transporte unitário por SVD para estabilizar a base sem alterar os autovalores. O espectro rastreado, reordenado por energia para comparação, coincide com `eigvalsh` até \(4.22\times10^{-15}\). Tracking entra exclusivamente no gráfico. As integrais e a ocupação usam uma diagonalização independente com `eigvalsh`.

[DEDUÇÃO MATEMÁTICA] Para evitar ambiguidade, neste relatório e no CSV,

\[
\boxed{g_{\rm zero}=\min_{q,n:E_n>0}E_n(q)-\max_{q,n:E_n<0}E_n(q).}
\]

Esse **zero-energy band gap** é a largura total do intervalo sem bandas em torno de zero. Devido ao pareamento demonstrado,

\[
g_{\rm zero}=2\delta_0,\qquad \delta_0=\min_{q,n}|E_n(q)|.
\]

[RESULTADO NUMÉRICO] Uma malha de 4001 pontos localiza candidatos a mínimos de \(\delta_0\); minimização escalar contínua refina cada candidato. Os valores são

| setor | \(\delta_0\) | zero-energy band gap \(g_{\rm zero}\) | classificação |
|-------|--------------|-------------------------------------------------|---------------|
| OLD | 1.000000 | 2.000000 | gapped |
| Z-REVERSED | 1.000000 | 2.000000 | gapped |

[DEDUÇÃO MATEMÁTICA] As fórmulas isotrópicas certificam que esses mínimos são globais: \(5-4\cos q\ge1\) e \(1+8\sin^2q\ge1\), com igualdade em \(q=0\). O algoritmo numérico não contém esses valores como resultados fixos. Para acoplamentos arbitrários, a busca de mínimos por malha e refinamento não constitui uma prova global; a certificação analítica aqui se refere ao ponto isotrópico.

O gap das bandas de \(h=iA/2\) é a grandeza definida acima. Os autovalores positivos de \(iA\) são o dobro das bandas positivas de h; não se deve misturar essas convenções ao comparar excitações de Majorana.

## 13. Occupied-band sum F(q)

[DEDUÇÃO MATEMÁTICA] Define-se, separadamente em cada setor,

\[
F_\sigma(q)=\sum_{n:E_{\sigma,n}(q)<0}E_{\sigma,n}(q).
\]

No ponto isotrópico há exatamente duas bandas negativas em toda a BZ. Substituindo as bandas derivadas,

\[
\boxed{F_O(q)=-\sqrt{5+4\cos q}-\sqrt{5-4\cos q},}
\]

\[
\boxed{F_Z(q)=-2\sqrt{1+8\sin^2q}.}
\]

[RESULTADO NUMÉRICO]

![Soma das bandas ocupadas, sem suavização](<Fq_old_vs_z_reversed.png>)

O gráfico foi calculado de autovalores novos de `eigvalsh`, mantendo a multiplicidade das bandas. Não houve suavização, ajuste de curvas ou uso do tracking no cálculo de F.

[DEDUÇÃO MATEMÁTICA] As funções não coincidem ponto a ponto. Por exemplo,

\[
F_O(0)=-4,\quad F_Z(0)=-2;
\qquad F_O(\pi/2)=-2\sqrt5,\quad F_Z(\pi/2)=-6.
\]

Não surgem cúspides: os radicandos são estritamente positivos e nenhuma banda cruza \(E=0\). O encontro de bandas OLD em \(\pm\pi/2\) ocorre em \(\pm\sqrt5\), distante do nível zero, e não produz uma singularidade na soma ocupada. Os cruzamentos entre as curvas \(F_O\) e \(F_Z\) comparam setores distintos e também não são cruzamentos de bandas em zero.

## 14. Ground-state energy per site

[DERIVADO DO PDF ORIGINAL] A Eq.(28) da referência obrigatória soma todos os autovalores negativos, dividindo pelo número de sítios. Um autovalor isolado não é a energia fundamental do sistema.

[DEDUÇÃO MATEMÁTICA] Há \(N_{\rm sites}=4N_{\rm cell}\), e a translação \(2a_1\) fornece densidade de momentos \(N_{\rm cell}/\pi\). No limite termodinâmico,

\[
\frac{E_0}{N_{\rm cell}}=\frac1\pi\int_{-\pi/2}^{\pi/2}F(q)\,dq,
\qquad
\boxed{\epsilon=\frac{E_0}{N_{\rm sites}}=
\frac1{4\pi}\int_{-\pi/2}^{\pi/2}F(q)\,dq.}
\]

Para \(\{\gamma_i,\gamma_j\}=2\delta_{ij}\), a convenção finita é

\[
H=\frac i4\gamma^TA\gamma,\qquad h=\frac{iA}{2},\qquad
E_0=\sum_{\lambda(h)<0}\lambda(h)=-\frac14\operatorname{Tr}|iA|.
\]

Um dímero \(H=iK\gamma_1\gamma_2\) tem energias muitos-corpos \(\pm K\), e h tem autovalores \(\pm K\); esse teste fixa a ausência de um fator extra \(1/2\) na soma dos autovalores negativos de h.

[RESULTADO NUMÉRICO] A tabela vem das execuções efetivas, com regra dos trapézios e pontos incluindo as duas extremidades da BZ:

| Nq | epsilon_old | epsilon_z | Delta_epsilon |
|----|-------------|-----------|---------------|
| 1001 | -1.063544409973 | -1.063544409973 | +2.220e-16 |
| 2001 | -1.063544409973 | -1.063544409973 | +0.000e+00 |
| 4001 | -1.063544409973 | -1.063544409973 | -2.220e-16 |
| 8001 | -1.063544409973 | -1.063544409973 | +0.000e+00 |
| 16001 | -1.063544409973 | -1.063544409973 | +0.000e+00 |

As três últimas malhas são estáveis em 12 casas decimais, limite usado pelo diagnóstico de convergência. As cinco malhas também concordam nessa precisão. Para apresentação física usamos seis casas:

\[
\boxed{\epsilon_{\rm OLD}=-1.063544,\qquad
\epsilon_Z=-1.063544,\qquad\Delta\epsilon_Z=0.000000.}
\]

As diferenças de \(\pm2.22\times10^{-16}\) em algumas malhas são arredondamento e não seleção energética. A quadratura adaptativa fornece \(-1.0635444099733649\) (OLD) e \(-1.063544409973365\) (Z), com estimativas de erro por sítio de \(2.66\times10^{-14}\) e \(2.51\times10^{-14}\). Os dígitos adicionais são sustentados pela convergência e pela verificação analítica a seguir.

### Demonstração do empate exato

[DEDUÇÃO MATEMÁTICA] Defina \(g(t)=\sqrt{1+8\sin^2t}\). Como

\[
5-4\cos q=1+8\sin^2(q/2),\qquad
5+4\cos q=1+8\cos^2(q/2),
\]

e as funções são pares,

\[
\begin{aligned}
-\int_{-\pi/2}^{\pi/2}F_O(q)\,dq
&=2\int_0^{\pi/2}\left[g(q/2)+g(\pi/2-q/2)\right]dq\\
&=4\int_0^{\pi/4}\left[g(t)+g(\pi/2-t)\right]dt\\
&=4\int_0^{\pi/2}g(t)\,dt.
\end{aligned}
\]

Por outro lado,

\[
-\int_{-\pi/2}^{\pi/2}F_Z(q)\,dq=4\int_0^{\pi/2}g(q)\,dq.
\]

Portanto as duas integrais são **exatamente iguais**, apesar de dispersões e fluxos diferentes. Usando a integral elíptica completa de segunda espécie na convenção de parâmetro m,

\[
\mathbb E(m):=\int_0^{\pi/2}\sqrt{1-m\sin^2t}\,dt,
\qquad\int_0^{\pi/2}g(t)\,dt=3\mathbb E(8/9),
\]

obtemos

\[
\boxed{\epsilon_O=\epsilon_Z=-\frac3\pi\mathbb E(8/9),\qquad
\Delta\epsilon_Z=0.}
\]

[RESULTADO NUMÉRICO] A avaliação da expressão elíptica é \(-1.0635444099733649\), compatível com as duas integrais calculadas diretamente das matrizes. A energia por célula é quatro vezes a energia por sítio, aproximadamente \(-4.254178\).

## 15. Numerical validation

[TESTE INDEPENDENTE] A implementação separa três fontes de verificação:

1. A lista de ligações reais transcrita da Eq.(11), usada para construir a matriz RAW e a cadeia finita.
2. As seis fórmulas fechadas, usadas como comparação independente da montagem matricial.
3. A derivação simbólica do polinômio, as fórmulas isotrópicas e a igualdade analítica das integrais.

Essa independência é de construção algébrica e numérica, não de experimentos físicos distintos. O teste finito detecta fases e deslocamentos errados, mas não substituiria uma auditoria visual da lista de ligações; por isso a figura e a Eq.(11) foram conferidas primeiro.

### Cadeia periódica finita

[TESTE INDEPENDENTE] Para cada ligação dirigida da célula n, soma-se

\[
A_{(n,i),(n+\Delta n,j)}\mathrel{+}=2K_\gamma u_{ij,\gamma},
\qquad A_{(n+\Delta n,j),(n,i)}\mathrel{-}=2K_\gamma u_{ij,\gamma},
\]

com índices periódicos e

\[
\Delta n=\frac{da_1-(r_j-r_i)}2,\qquad r=(0,1,1,0).
\]

Os \(\Delta n\) são inteiros conferidos a partir da geometria. A base de Fourier completa é

\[
V_{(n,m),a}(q,q_2)=\frac{\delta_{ma}}{\sqrt{N_{\rm cell}}}
e^{i[q(2n+r_m)+q_2\eta_m]},\qquad \eta=(0,0,1,1),
\]

\[
q=\frac{\pi\ell}{N_{\rm cell}}.
\]

Com \(q_2=0\), compara-se \(V^\dagger(iA/2)V\) com \(H_{\rm Bloch}(q)\). Com \(q_2=0.371\), compara-se com \(H_{\rm raw}(q,q_2)\), verificando também os deslocamentos entre cadeias. A união dos espectros dos blocos deve coincidir com todos os autovalores da matriz finita; a soma negativa por sítio deve coincidir nas duas representações.

Foram executadas cadeias de **5, 8 e 17 células**, isto é, **20, 32 e 68 sítios**, em **ambos** os setores, incluindo todos os momentos permitidos. Os acoplamentos foram \((1.13,0.79,1.31,0.92,0.67)\).

### Casos aleatórios e erros máximos

[TESTE INDEPENDENTE] Foram executados **1000 conjuntos aleatórios**, cada um nos dois setores: 2000 avaliações de setor, semente `20260908`. Os módulos dos cinco acoplamentos são distintos, amostrados entre 0.15 e 2.0 com separação mínima 0.025. Metade dos conjuntos usa acoplamentos positivos; na outra metade os sinais também são sorteados. \(q,q_2\) são sorteados em \([-3\pi,3\pi]\).

[RESULTADO NUMÉRICO] Erros máximos efetivamente registrados:

| Verificação | Erro máximo |
|-------------|-------------|
| Unitariedade de U | 2.220e-16 |
| RAW: fórmula versus ligações | 3.587e-15 |
| Chamada matricial Eq.(27) | 0.000e+00 |
| Bloch da geometria versus fórmulas fechadas | 3.928e-15 |
| Cancelamento de q2 | 3.928e-15 |
| Hermiticidade | 8.951e-16 |
| Relação Majorana | 3.928e-15 |
| Periodicidade matricial 2pi | 4.021e-15 |
| Covariância matricial pi | 3.928e-15 |
| Periodicidade espectral pi | 5.329e-15 |
| Espectros RAW versus Bloch | 6.217e-15 |
| Pairing no mesmo q | 5.773e-15 |
| Polinômio característico (normalizado) | 6.678e-15 |
| Diferença OLD versus Z somente nos termos z | 4.965e-16 |
| Cadeia finita versus Bloch | 1.102e-15 |
| Cadeia finita versus RAW | 1.510e-15 |
| Espectro total finito versus união dos blocos | 3.109e-15 |
| Energia por sítio finita em ambas representações | 2.220e-16 |

O erro do polinômio é normalizado pela maior magnitude de seus coeficientes ou 1; os demais são absolutos. O erro zero no teste direto da chamada Eq.(27) confirma a implementação da chamada, mas não é evidência independente por si só. As comparações RAW versus geometria, fórmulas fechadas, cancelamento de \(q_2\) e projeção finita fornecem essa verificação adicional.

### Cobertura dos requisitos

| teste requerido | resultado |
|-----------------|-----------|
| todos os u inteiros em \(\{-1,+1\}\) | passou |
| anticontaminação x,w superiores e todos os não-z | passou com assertions explícitos |
| quatro fluxos x e quatro fluxos w | calculados de percursos reais; passou |
| fechamento e orientação anti-horária dos percursos | passou por área orientada positiva |
| unitariedade de U | identidade simbólica e teste numérico |
| RAW fechado versus RAW das ligações | passou nos dois setores |
| Eq.(27) e cancelamento de \(q_2\) | passou simbólica e numericamente |
| hermiticidade e Majorana | passou simbólica e numericamente |
| período matricial \(2\pi\), covariância e espectro em \(\pi\) | passou |
| polinômio, ausência de termo linear e pairing | passou |
| cadeia finita, espectro total e energia por sítio | passou nos seis casos de setor/tamanho |
| diferença OLD versus Z somente em f13 e f24 | passou |
| busca das 16 transformações locais | nenhuma solução; invariância dos fluxos verificada para as 16 |
| 1000 casos assimétricos | passou |
| fórmulas isotrópicas versus `eigvalsh` | passou |
| tracking preserva o espectro | passou; não usado na energia |
| convergência, BZ deslocada e BZ estendida | passou |

### Invariância da BZ e normalização

[DEDUÇÃO MATEMÁTICA] Pela periodicidade espectral de \(\pi\), a mesma energia resulta de qualquer intervalo de largura \(\pi\). Usando a BZ estendida, deve-se dividir sua integral por \(8\pi\), pois a largura foi duplicada:

\[
\epsilon=\frac1{4\pi}\int_{a}^{a+\pi}F(q)dq
=\frac1{8\pi}\int_{-\pi}^{\pi}F(q)dq.
\]

[RESULTADO NUMÉRICO] A BZ deslocada por 0.371 e a BZ estendida foram testadas com o mesmo passo de integração. O maior erro em energia por sítio entre essas escolhas foi \(2.22\times10^{-16}\). Isso verifica a normalização sem redefinir a BZ física.

### Reprodução, arquivos e warnings

[TESTE INDEPENDENTE] Os dois scripts finais foram efetivamente executados, na versão entregue, com warnings tratados como erros:

```bash
python3 -W error quadrupolar_spin_liquid_z_reversed.py
python3 -W error compare_old_vs_z_reversed.py
```

Ambos terminaram com **código de saída 0 e sem warnings na execução final**. O primeiro grava `validation_z_reversed.json`; o segundo gera os três PNGs, o CSV requerido, um CSV adicional de convergência e `comparison_z_reversed.json`. `draw_z_lattice.py` é o helper de desenho importado pelo segundo script e deve acompanhá-lo. Os scripts anteriores não foram sobrescritos.

Durante a preparação ocorreram duas mensagens de ambiente gráfico: criação de cache temporário do Matplotlib porque o diretório padrão não era gravável; e um `UserWarning` de incompatibilidade Matplotlib/pyparsing ao inicializar mathtext, acompanhado de mensagens no destrutor. As versões finais configuram um diretório gravável e usam rótulos Unicode sem mathtext. Nenhum warning foi suprimido. Uma assertion de desenvolvimento também detectou que o gerador de `charpoly` do SymPy não herdava as assumptions do símbolo externo; o teste foi corrigido substituindo explicitamente `characteristic.gen`, e todas as verificações simbólicas passaram depois disso.

Ambiente usado: NumPy 1.26.4, SymPy 1.12, SciPy 1.11.4 e Matplotlib 3.6.3. A API `eigvalsh` fornece autovalores hermitianos ordenados com multiplicidade, como requerido para a ocupação ([documentação NumPy](https://numpy.org/doc/stable/reference/generated/numpy.linalg.eigvalsh.html)). O determinante foi verificado pelo polinômio característico de matriz ([documentação SymPy](https://docs.sympy.org/latest/modules/matrices/matrices.html)), e `ellipe` usa a integral elíptica por parâmetro m definida acima ([documentação SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html)). Essas referências documentam somente as ferramentas de cálculo; a fonte física obrigatória é a referência OLD e sua geometria original.

Os arquivos de saída são:

- [Código físico e testes](<quadrupolar_spin_liquid_z_reversed.py>).
- [Comparação, gráficos e CSV](<compare_old_vs_z_reversed.py>).
- [Helper da geometria](<draw_z_lattice.py>).
- [CSV energético](<ground_state_energy_old_vs_z.csv>).
- [CSV de convergência](<energy_convergence_old_vs_z.csv>).
- [Resultados dos testes](<validation_z_reversed.json>).
- [Resultados da comparação](<comparison_z_reversed.json>).

O CSV requerido contém exatamente `Kx,Ky,Kz,Kt,Kw,epsilon_old,epsilon_z,Delta_epsilon,gap_old,gap_z,favored_sector`. No ponto isotrópico, `favored_sector=TIE`. `gap_old` e `gap_z` usam o gap completo de largura 2, não \(\delta_0=1\). Os dados brutos mantêm o arredondamento numérico real, sem forçar gaps inteiros. Para execuções opcionais com outros acoplamentos, a decisão de preferência também considera uma estimativa do erro de integração e pode retornar `UNRESOLVED`; não declara empate exato fora do caso demonstrado.

Os três PNGs foram confirmados como arquivos não vazios, com metadados equivalentes a 300 dpi, e inspecionados visualmente. O helper confere que as setas desenhadas coincidem com o grafo real transladado e que somente z muda em relação ao ORIGINAL reconstruído.

## 16. Physical conclusion

[DEDUÇÃO MATEMÁTICA + RESULTADO NUMÉRICO]

**“Entre OLD e Z-INVERTED, qual configuração possui menor energia no ponto isotrópico?”**

**Nenhuma possui energia menor: as duas configurações empatam exatamente na energia por sítio no limite termodinâmico, sob a normalização prescrita.**

Esse resultado foi obtido pela integração das quatro bandas, respeitando a ocupação negativa, e demonstrado por mudança de variável. O gap não foi usado para decidir a preferência energética.

Os setores não são gauge-equivalentes: os fluxos das plaquetas contendo z mudam de sinal. O empate de energia não restaura essa equivalência. As estruturas de bandas são diferentes; em Z, as bandas positivas e negativas são duplamente degeneradas em todo q no ponto isotrópico.

Esta comparação tampouco demonstra que alguma dessas configurações minimize a energia entre todos os setores de gauge possíveis, inclusive padrões de células maiores ou não periódicos.

## 17. What is proven and what is not

**VALIDADO**

- Gauge ORIGINAL conferido contra a figura e a Eq.(15), usando obrigatoriamente a referência OLD fornecida.
- Nova configuração altera somente `u42z` e `u31z`; x,y,t,w preservados.
- Oito produtos de fluxo derivados de caminhos anti-horários reais, e inequivalência de gauge local.
- Hamiltonianos reais OLD e Z, Fourier RAW, Eq.(27), cancelamento de \(q_2\), matriz de Bloch e polinômio.
- Hermiticidade, relação Majorana, simetrias especiais e periodicidades.
- Testes simbólicos, 1000 conjuntos assimétricos nos dois setores e projeções independentes de cadeias finitas.
- Dispersões isotrópicas, gap completo 2, \(F(q)\), convergência e empate energético exato.
- Arquivos gráficos e CSV efetivamente gerados e verificados.

**AINDA NÃO CONCLUÍDO**

- Estado fundamental global entre todos os setores de gauge.
- Classificação topológica, invariantes, estados de borda ou correlações quadrupolares deste novo setor.
- Diagrama de fases energético completo com cinco acoplamentos variáveis.
- Equivalência de energias para acoplamentos arbitrários: o empate demonstrado aqui é isotrópico.
- Correções de paridade/projeção sobre o espaço físico para cada cadeia finita. A energia apresentada é a densidade termodinâmica do Hamiltoniano de matéria no setor fixo, na convenção solicitada. Os testes finitos verificam a matriz quadrática e sua normalização; não implementam essa projeção adicional.

**Ressalva documental:** a Eq.(26) transcrita literalmente na referência obrigatória foi preservada e identificada como discrepante. Os resultados deste relatório correspondem ao gauge e à geometria ORIGINAL, com o coeficiente f34 derivado da Eq.(11)/Eq.(17). Eles não representam a matriz obtida pela adoção simultânea de todas as fórmulas transcritas sem resolver essa discrepância. A interpretação de Eq.(27) foi explicitada sem alterar U.

### Physical/mathematical verdict

**APROVADO COM RESSALVAS.** As derivações e validações são consistentes com a geometria ORIGINAL; a ressalva é a errata de fases na Eq.(26) transcrita e a distinção entre substituição ativa e coordenadas na Eq.(27).

### Code verdict

**A) correto e cientificamente robusto**, para os dois setores e a análise isotrópica definidos neste relatório, com testes independentes e limitações de escopo declaradas.

### Main physical result at isotropic point

**OLD:**

- flux pattern: x \((+1,+1,+1,+1)\); w \((+1,+1,+1,+1)\).
- gap: zero-energy band gap completo \(2.000000\); \(\min|E|=1.000000\).
- epsilon_old: \(-1.063544\) por sítio.

**Z-REVERSED:**

- flux pattern: x \((-1,-1,+1,+1)\); w \((-1,-1,+1,+1)\), nas ordens explicitadas na seção 3.
- gap: zero-energy band gap completo \(2.000000\); \(\min|E|=1.000000\).
- epsilon_z: \(-1.063544\) por sítio.

\[
\boxed{\Delta\epsilon=\epsilon_Z-\epsilon_{\rm OLD}=0\quad\text{exatamente}.}
\]

**Favored sector between these two: EMPATE — nenhum é favorecido energeticamente no ponto isotrópico.**

“É seguro utilizar estes resultados na próxima discussão com o orientador?”

**SIM, COM RESSALVAS.** Apresente a correção documentada da Eq.(26), a identificação dos espinores na Eq.(27) e o alcance restrito da comparação entre os dois setores. Dentro desse escopo, as contas e a execução foram verificadas.
