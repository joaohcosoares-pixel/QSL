# Auditoria algébrica independente — OLD e W-REVERSED

[DERIVADO DO PDF] Esta reconstrução usa os dez termos da Eq.(11), nas páginas 3–4 do PDF, e os dez sinais da Eq.(15), na página 5. As páginas 4 e 5 foram inspecionadas visualmente nos renders preservados em `auditoria_fontes/original-4.png` e `original-5.png`. A lista independente em `symbolic_w_audit.py` permite refutar a implementação principal sem importá-la. A única duplicação de constantes é essa pequena transcrição de fonte usada como fixture de auditoria.

O intervalo principal em todas as contas e comparações é

\[
\boxed{\mathrm{PROJECT\_BZ}=[-\pi,\pi].}
\]

## 1. Reconstrução de RAW e mudança de base

[DERIVADO DO PDF] Para um termo armazenado orientado \(i\to j\), com deslocamento final menos inicial \((d_1,d_2)\), a contribuição é

\[
(H_{\rm raw})_{ij}\mathrel{+}=iK_\gamma u_{ij,\gamma}e^{i(qd_1+q_2d_2)},\qquad
(H_{\rm raw})_{ji}\mathrel{+}=\left[(H_{\rm raw})_{ij}\right]^*.
\]

O sinal \(u=-1\) representa orientação física oposta à orientação de armazenamento. A conjugação do termo armazenado com índices \(43\), por exemplo, deve ser realizada antes de extraí-lo como \(f_{34}\).

| Termo | OLD \(u\) | W \(u\) | \(d_1,d_2\) |
|---|---:|---:|---|
| \(12,x\) | +1 | +1 | \((1,0)\) |
| \(23,y\) | +1 | +1 | \((0,1)\) |
| \(34,x\) | +1 | +1 | \((-1,0)\) |
| \(41,y\) | −1 | −1 | \((0,-1)\) |
| \(42,z\) | +1 | +1 | \((1,-1)\) |
| \(43,w\) | +1 | −1 | \((-1,0)\) |
| \(31,z\) | +1 | +1 | \((1,-1)\) |
| \(12,w\) | −1 | +1 | \((-1,0)\) |
| \(13,t\) | −1 | −1 | \((1,1)\) |
| \(24,t\) | −1 | −1 | \((1,1)\) |

[DEDUÇÃO MATEMÁTICA] Somando esses termos, sem começar pelas fórmulas finais:

| Elemento superior RAW | OLD | W-REVERSED |
|---|---|---|
| \(12\) | \(i(K_xe^{iq}-K_we^{-iq})\) | \(i(K_xe^{iq}+K_we^{-iq})\) |
| \(13\) | \(-ie^{iq_2}(K_te^{iq}+K_ze^{-iq})\) | igual |
| \(14\) | \(iK_ye^{iq_2}\) | igual |
| \(23\) | \(iK_ye^{iq_2}\) | igual |
| \(24\) | \(-ie^{iq_2}(K_te^{iq}+K_ze^{-iq})\) | igual |
| \(34\) | \(i(K_xe^{-iq}-K_we^{iq})\) | \(i(K_xe^{-iq}+K_we^{iq})\) |

Para os quatro elementos entre as duas cadeias, o fator RAW é **\(e^{+iq_2}\)**. Com

\[
U=\operatorname{diag}(1,1,e^{-iq_2},e^{-iq_2}),\qquad
\Theta_{\rm raw}=U\Theta_{\rm Bloch},\qquad H=U^\dagger H_{\rm raw}U,
\]

esse fator cancela. O script testa RAW antes da transformação e verifica também que a transformação incorreta \(UH_{\rm raw}U^\dagger\) preserva dependência espúria em \(q_2\). A validação não permite que dois erros de sinal se cancelem.

## 2. Discrepância documental da Eq.(26)

[DERIVADO DO PDF] A Eq.(26) literal preservada na referência OLD é

\[
f_{34}^{\rm literal}=i(K_xe^{iq}-K_we^{-iq}).
\]

[DEDUÇÃO MATEMÁTICA] A geometria da Eq.(11), em concordância com a Eq.(17), fornece

\[
f_{34}^{\rm OLD}=i(K_xe^{-iq}-K_we^{iq}),\qquad
f_{34}^{\rm OLD}-f_{34}^{\rm literal}=2(K_x+K_w)\sin q.
\]

Essa diferença não pode ser removida pela Eq.(27), cujas fases nas sub-redes 3 e 4 são iguais. A referência original permanece intocada. Aqui, OLD designa os sinais originais com a matriz reconstruída da geometria, não a aplicação literal da expressão documental discrepante.

## 3. Diferença dos Hamiltonianos

[DEDUÇÃO MATEMÁTICA] Escrevendo \(\theta_m=\theta_m^y\), os dois termos que mudam são

\[
H_{\rm OLD}^{(w)}=iK_w\sum_R\left[
\theta_4(R+a_2)\theta_3(R+a_2-a_1)
-\theta_1(R)\theta_2(R-a_1)\right],
\]

\[
H_{\rm W}^{(w)}=-H_{\rm OLD}^{(w)},\qquad
\boxed{\Delta H_W=2iK_w\sum_R\left[
\theta_1(R)\theta_2(R-a_1)
-\theta_4(R+a_2)\theta_3(R+a_2-a_1)\right].}
\]

Todos os outros oito termos ficam invariantes. A matriz extraída exclusivamente desses dois bonds é

\[
\Delta H(q)=2iK_w
\begin{pmatrix}
0&e^{-iq}&0&0\\
-e^{iq}&0&0&0\\
0&0&0&e^{iq}\\
0&0&-e^{-iq}&0
\end{pmatrix}.
\]

O script compara essa construção com a diferença das duas matrizes completas, sem usar as fórmulas fechadas de \(f_{12}\) e \(f_{34}\) na construção do lado independente.

## 4. Forma reduzida e determinante

[RESULTADO SIMBÓLICO] Para ambos os setores, a extração da matriz demonstra

\[
f_{23}=f_{14}=iK_y,\quad f_{24}=f_{13}=b,\quad
f_{12}=a,\quad f_{34}=-a^*.
\]

Com \(a=A+iB\), \(b=C+iD\) e \(y=K_y\),

\[
H=\begin{pmatrix}
0&a&b&iy\\
a^*&0&iy&b\\
b^*&-iy&0&-a^*\\
-iy&b^*&-a&0
\end{pmatrix}.
\]

| Variável | OLD | W-REVERSED |
|---|---|---|
| \(A\) | \(-(K_x+K_w)\sin q\) | \((K_w-K_x)\sin q\) |
| \(B\) | \((K_x-K_w)\cos q\) | \((K_x+K_w)\cos q\) |
| \(C\) | \((K_t-K_z)\sin q\) | igual |
| \(D\) | \(-(K_t+K_z)\cos q\) | igual |

Definindo

\[
S=A^2+B^2+C^2+D^2+y^2,\qquad
R=B^2(C^2+D^2)+y^2D^2,
\]

o determinante simbólico independente fornece

\[
\boxed{P(\lambda,q)=\lambda^4-2S\lambda^2+S^2-4R.}
\]

Os cinco coeficientes, em ordem decrescente, são \(1,0,-2S,0,S^2-4R\). Em particular:

\[
S_W=K_x^2+K_y^2+K_z^2+K_t^2+K_w^2
+2(K_tK_z+K_xK_w)\cos2q,
\]

\[
R_W=\cos^2q\left[
(K_x+K_w)^2(K_t^2+K_z^2+2K_tK_z\cos2q)
+K_y^2(K_t+K_z)^2\right].
\]

Para OLD, substitui-se nesta forma já derivada o sinal dos termos de \(K_w\):

\[
S_O=\sum_\gamma K_\gamma^2+2(K_tK_z-K_xK_w)\cos2q,
\]

\[
R_O=\cos^2q\left[
(K_x-K_w)^2(K_t^2+K_z^2+2K_tK_z\cos2q)
+K_y^2(K_t+K_z)^2\right].
\]

O espectro é

\[
\boxed{E=\left\{-\sqrt{S+2\sqrt R},-\sqrt{S-2\sqrt R},
\sqrt{S-2\sqrt R},\sqrt{S+2\sqrt R}\right\}.}
\]

Esta é uma consequência do determinante, e não um pareamento imposto ao resultado. Para acoplamentos reais, a Hermiticidade garante \(S\ge2\sqrt R\). Degenerescências de mesmo sinal ocorrem quando \(R=0\); a condição de fechamento em zero é \(S^2=4R\). Essas condições gerais incluem casos especiais de acoplamentos nulos.

## 5. Simetrias

[RESULTADO SIMBÓLICO] O script demonstra separadamente

\[
H^\dagger(q)=H(q),\qquad H(-q)=-H(q)^*,\qquad H(q+2\pi)=H(q),
\]

e

\[
H(q+\pi)=GH(q)G^\dagger,\qquad G=\operatorname{diag}(1,-1,-1,1).
\]

Portanto, o espectro tem período \(\pi\), embora a matriz usualmente tenha período \(2\pi\). Essa redundância **não altera** o intervalo principal solicitado \([-\pi,\pi]\).

O pareamento no mesmo momento admite ainda a prova antiunitária

\[
JH(q)^*J^\dagger=-H(q),\qquad
J=\begin{pmatrix}
0&0&-1&0\\0&0&0&-1\\1&0&0&0\\0&1&0&0
\end{pmatrix},\qquad J^2=-I.
\]

Esta relação leva \(E\) a \(-E\), portanto não implica degenerescência entre duas bandas positivas em um momento genérico. Não deve ser confundida com uma simetria antiunitária que comuta com o Hamiltoniano.

## 6. Limite isotrópico

[RESULTADO SIMBÓLICO] A substituição \(K_x=K_y=K_z=K_t=K_w=1\) na matriz derivada dos bonds, com \(c=\cos q\), dá

\[
H_W(q)=i\begin{pmatrix}
0&2c&-2c&1\\-2c&0&1&-2c\\2c&-1&0&2c\\-1&2c&-2c&0
\end{pmatrix},
\]

\[
S_W=1+8c^2,\qquad R_W=4c^2(1+4c^2),\qquad
\boxed{P_W(\lambda,q)=\lambda^4-2(1+8c^2)\lambda^2+1.}
\]

Em particular, \(\det H_W=1\) para todo \(q\); não existe fechamento em energia zero. As duas energias positivas ordenadas são

\[
e_+(q)=\sqrt{1+4c^2}+2|c|,\qquad
e_-(q)=\sqrt{1+4c^2}-2|c|,
\]

e o espectro completo é \(\{-e_+,-e_-,e_-,e_+\}\). As raízes propostas são substituídas no determinante no teste simbólico, e comparadas com `eigvalsh` em uma malha independente de 1001 pontos.

Como

\[
e_-e_+=1,\qquad e_-=\frac{1}{\sqrt{1+4c^2}+2|c|},
\]

o mínimo global, em \(|c|=1\), é

\[
\boxed{\min_{q,n}|E_n(q)|=\sqrt5-2=0.2360679774997897\ldots}
\]

nos pontos \(q=-\pi,0,\pi\). Definindo a largura do intervalo sem bandas centrado em zero,

\[
\boxed{\mathrm{zero\_energy\_band\_gap}=2\sqrt5-4=0.4721359549995794\ldots.}
\]

A degenerescência de mesmo sinal ocorre somente em \(q=\pm\pi/2\), onde o espectro é \((-1,-1,1,1)\). Os valores de bandas ordenadas podem ter mudanças de derivada nesses cruzamentos, mas a soma ocupada abaixo permanece suave.

| \(q\) | Espectro W ordenado |
|---|---|
| \(-\pi\) | \(-2-\sqrt5,\ 2-\sqrt5,\ \sqrt5-2,\ 2+\sqrt5\) |
| \(-\pi/2\) | \(-1,-1,1,1\) |
| \(-\pi/4\) | \(-\sqrt3-\sqrt2,\ \sqrt2-\sqrt3,\ \sqrt3-\sqrt2,\ \sqrt3+\sqrt2\) |
| \(0\) | \(-2-\sqrt5,\ 2-\sqrt5,\ \sqrt5-2,\ 2+\sqrt5\) |
| \(\pi/4\) | \(-\sqrt3-\sqrt2,\ \sqrt2-\sqrt3,\ \sqrt3-\sqrt2,\ \sqrt3+\sqrt2\) |
| \(\pi/2\) | \(-1,-1,1,1\) |
| \(\pi\) | \(-2-\sqrt5,\ 2-\sqrt5,\ \sqrt5-2,\ 2+\sqrt5\) |

Em \(\pi/4\), as energias positivas são \(0.3178372451957820\ldots\) e \(3.146264369941972\ldots\).

Para comparação, a mesma derivação em OLD fornece

\[
P_O=\lambda^4-10\lambda^2+25-16\cos^2q,
\quad E_O=\left\{\pm\sqrt{5+4|\cos q|},\pm\sqrt{5-4|\cos q|}\right\}.
\]

Assim, \(\min|E_O|=1\), e a largura de seu intervalo sem bandas centrado em zero é \(2\).

## 7. Normalização a partir de Majoranas e momentos discretos

[DEDUÇÃO MATEMÁTICA] Seja \(L=N_{\rm cell}\), com quatro Majoranas de matéria por célula e \(N=N_{\rm site}=4L\). Use \(\{\gamma_i,\gamma_j\}=2\delta_{ij}\),

\[
\widehat H=\frac{i}{4}\gamma^T A\gamma,\qquad
A_{ij}=2K_\gamma u_{ij,\gamma},\qquad A_{ji}=-A_{ij}.
\]

Cada par canônico de \(A\), com \(iA\) tendo autovalores \(\pm\omega_\alpha\), produz

\[
\widehat H_\alpha=\omega_\alpha(f_\alpha^\dagger f_\alpha-1/2),\qquad
E_0=-\frac12\sum_{\alpha=1}^{2L}\omega_\alpha.
\]

As bandas da convenção do PDF são autovalores de \(h=iA/2\), logo \(\pm E_\alpha=\pm\omega_\alpha/2\). Portanto,

\[
\boxed{E_0=\sum_{E_j(h)<0}E_j(h)=-\frac14\operatorname{Tr}|iA|.}
\]

O fator Majorana já está incorporado na relação \(h=iA/2\). Acrescentar outro fator \(1/2\) à soma negativa de \(h\) produziria um erro. Como verificação mínima, um bond isolado \(iK\gamma_1\gamma_2\) tem \(A_{12}=2K\), \(h\) com autovalores \(\pm|K|\) e energia fundamental \(-|K|\), exatamente a soma negativa de \(h\).

Os \(L\) momentos de célula são \(p_m=2\pi m/L\). A geometria original usa translação \(2a_1\), de modo que \(p=2q\), com fases internas das sub-redes mantidas. Assim há \(L\) valores independentes de \(q\) por período espectral \(\pi\). Não são \(2L\) momentos independentes quando se representa a zona solicitada de largura \(2\pi\).

Com \(F(q)=\sum_{E_n(q)<0}E_n(q)\), a Eq.(28) e a contagem discreta dão

\[
\frac{E_0(L)}{4L}=\frac1{4L}\sum_{m=0}^{L-1}F(q_m)
\longrightarrow\frac1{4\pi}\int_{\text{um período espectral}}F(q)\,dq.
\]

Na zona **principal solicitada**, com dois períodos espectrais e sem duplicar estados,

\[
\boxed{\epsilon=\frac{E_0}{N_{\rm site}}
=\frac1{8\pi}\int_{-\pi}^{\pi}F(q)\,dq.}
\]

A representação de Fourier original em BZ/2 tem \(L/2\) pares genéricos \((q,-q)\), quatro modos complexos por par e tratamento separado dos momentos auto-conjugados quando presentes. A representação na BZ completa inclui os parceiros Majorana e usa os elementos \(iK\), em lugar dos \(2iK\) da escrita em meia BZ. A cadeia finita acima elimina a ambiguidade de dupla contagem sem depender de como os pontos especiais são alocados à meia zona.

Esta normalização corresponde à energia de matéria no setor de gauge fixado, conforme Eq.(28). Uma projeção adicional no espaço físico de spins pode impor uma condição global de paridade em tamanhos finitos; não foi solicitada nem provada uma energia finita projetada. Tal ressalva não modifica a densidade de energia no limite termodinâmico.

## 8. Soma ocupada e energia exata isotrópica

[DEDUÇÃO MATEMÁTICA] Para W,

\[
\boxed{F_W(q)=-e_+(q)-e_-(q)=-2\sqrt{1+4\cos^2q}.}
\]

Para OLD,

\[
F_O(q)=-\sqrt{5+4\cos q}-\sqrt{5-4\cos q}.
\]

Os radicandos são estritamente positivos; nenhuma dessas somas ocupadas isotrópicas contém cúspide em energia zero. Não é necessário suavizar numericamente qualquer uma delas.

Defina a integral elíptica completa de segunda espécie pelo parâmetro \(m\), e não pelo módulo:

\[
\mathbb E(m)=\int_0^{\pi/2}\sqrt{1-m\sin^2t}\,dt.
\]

Então

\[
\boxed{\epsilon_W=-\frac{\sqrt5}{\pi}\mathbb E(4/5)
=-0.8388049859310991\ldots,}
\]

pois \(1+4\cos^2q=5[1-(4/5)\sin^2q]\). Para OLD, escrevendo \(5+4\cos q=1+8\cos^2(q/2)\) e \(5-4\cos q=1+8\sin^2(q/2)\), uma mudança de variável mostra

\[
\boxed{\epsilon_O=-\frac3\pi\mathbb E(8/9)
=-1.0635444099733649\ldots.}
\]

[RESULTADO NUMÉRICO] A avaliação independente com `scipy.special.ellipe` dá

\[
\boxed{\Delta\epsilon=\epsilon_W-\epsilon_O
=0.22473942404226577\ldots>0.}
\]

OLD é energeticamente favorecido **entre os dois setores comparados no ponto isotrópico**. O menor gap de W não decide a comparação de energia. Não se demonstrou que OLD seja o mínimo entre todas as configurações de gauge, tampouco uma ordem energética universal para acoplamentos anisotrópicos.

## 9. Escopo dos testes independentes

[TESTE INDEPENDENTE] `symbolic_w_audit.py` reconstrói RAW da transcrição da fonte, cancela \(q_2\), verifica os seis coeficientes estruturais e as simetrias, calcula o determinante em variáveis reais independentes, substitui a matriz efetivamente derivada nessa forma reduzida e calcula a diferença usando somente os bonds \(w\). O limite isotrópico é obtido por substituição após a derivação; suas raízes são testadas no determinante e contra `eigvalsh`. As validações críticas usam `raise AssertionError` explicitamente e permanecem ativas sob `python -O`.

O script também preserva a discrepância da Eq.(26) como identidade simbólica não nula. Ele não reutiliza JSON, CSV ou coeficientes numéricos previamente gerados. A auditoria principal deve completar esses testes com cadeias finitas, 1000 amostras anisotrópicas, convergência, hashes e benchmarks.

[TESTE INDEPENDENTE] A execução com `python3 -W error symbolic_w_audit.py` terminou em **PASS**. O erro máximo entre espectro isotrópico fechado e `eigvalsh` foi \(3.553\times10^{-15}\). Em 64 amostras anisotrópicas por setor, com sinais mistos, a comparação da transcrição independente com o módulo principal teve erro máximo RAW \(0\), Bloch \(1.832\times10^{-15}\) e diferença dos bonds \(w\) \(0\). O teste de gauge também exigiu igualdade das dez entradas da produção com a transcrição independente da fonte.
