r"""
Script de visualização científica em Python/Matplotlib para reprodução exata
do grafo de rede 1D definido no código TikZ fornecido.

Topologia:
- Cadeia inferior: nós b0 a b5 em (0, 0), (2, 0), ..., (10, 0)
- Cadeia superior: nós t0 a t5 em (1, 2), (3, 2), ..., (11, 2)
- Ligações horizontais: b_i -> b_{i+1} e t_i -> t_{i+1} (alternando bondx e bondw)
- Ligações y (verticais inclinadas): b_i -> t_i
- Ligações z (diagonais descendentes): t_i -> b_{i+1}
- Ligações t (diagonais ascendentes): b_i -> t_{i+1}
- Rótulos matemáticos: 1 abaixo de b2, 2 abaixo de b3, 4 acima de t2, 3 acima de t3
- Reticências horizontais (\cdots) em (-0.55, 1.0) e (11.55, 1.0)
"""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import matplotlib as mpl

# ---------------------------------------------------------------------------
# 1. Definição da paleta de cores para as classes de ligações (bonds)
# ---------------------------------------------------------------------------
# Paleta científica de alto contraste e acessível:
DEFAULT_BOND_COLORS: dict[str, str] = {
    "bondx": "#2B6CB0",  # Azul royal (ligação horizontal tipo x)
    "bondw": "#DD6B20",  # Laranja quente (ligação horizontal tipo w)
    "bondy": "#2F855A",  # Verde floresta (ligação vertical b_i -> t_i)
    "bondz": "#C53030",  # Vermelho carmim (diagonal t_i -> b_{i+1})
    "bondt": "#6B46C1",  # Roxo ametista (diagonal b_i -> t_{i+1})
}


def draw_qsl_lattice_graph(
    output_path: str = "qsl_lattice_graph.png",
    dpi: int = 300,
    usetex: bool = False,
    bond_colors: dict[str, str] | None = None,
    arrow_style: str = "-|>",
    node_radius_pt: float = 3.0,
    line_width: float = 1.5,
) -> tuple[plt.Figure, plt.Axes]:
    """
    Renderiza o grafo e salva a figura em PNG com fundo transparente.

    Parâmetros
    ----------
    output_path : str
        Caminho do arquivo PNG de saída.
    dpi : int
        Resolução da imagem em pontos por polegada (padrão: 300).
    usetex : bool
        Se True, usa o compilador LaTeX do sistema. Se False, usa o interpretador
        mathtext nativo do Matplotlib configurado com Computer Modern ('cm').
    bond_colors : dict[str, str], opcional
        Dicionário mapeando as classes ('bondx', 'bondw', 'bondy', 'bondz', 'bondt')
        para seus respectivos códigos hexadecimais de cor.
    arrow_style : str
        Estilo da ponta da seta ('-|>' para ponta preenchida LaTeX-like, ou '->').
    node_radius_pt : float
        Raio dos nós em pontos tipográficos (padrão: 3.0 pt, como no TikZ).
    line_width : float
        Espessura das linhas das arestas em pontos (padrão: 1.5 pt).

    Retorno
    -------
    tuple[plt.Figure, plt.Axes]
        Objeto da figura e dos eixos do Matplotlib.
    """
    colors = bond_colors or DEFAULT_BOND_COLORS

    # Configuração de fontes tipográficas para fidelidade matemática (estilo TeX)
    if usetex:
        mpl.rcParams["text.usetex"] = True
    else:
        mpl.rcParams["text.usetex"] = False
        mpl.rcParams["mathtext.fontset"] = "cm"
        mpl.rcParams["font.family"] = "serif"

    # Criar figura e eixo
    fig, ax = plt.subplots(figsize=(10.5, 2.8), dpi=dpi)

    # -----------------------------------------------------------------------
    # Coordenadas dos vértices:
    # b0 a b5: (0, 0) a (10, 0) com passo 2
    # t0 a t5: (1, 2) a (11, 2) com passo 2
    # -----------------------------------------------------------------------
    b = {i: (float(2 * i), 0.0) for i in range(6)}
    t = {i: (float(1 + 2 * i), 2.0) for i in range(6)}

    # Ajuste de encurtamento (shrink) para que as pontas das setas toquem
    # a borda dos nós circulares sem serem sobrepostas por eles.
    node_diameter_pt = 2.0 * node_radius_pt
    arrow_shrink_pt = node_radius_pt + 1.2

    # Formatação da seta no FancyArrowPatch
    if arrow_style == "-|>":
        fancy_style = "-|>,head_length=5.0,head_width=3.2"
    elif arrow_style == "->":
        fancy_style = "->,head_length=5.0,head_width=3.2"
    else:
        fancy_style = arrow_style

    # -----------------------------------------------------------------------
    # Construção da lista de arestas conforme a especificação TikZ
    # -----------------------------------------------------------------------
    edges: list[tuple[tuple[float, float], tuple[float, float], str]] = []

    # 1. Horizontais inferiores (para a direita: b_i -> b_{i+1})
    edges.extend([
        (b[0], b[1], "bondx"),
        (b[1], b[2], "bondw"),
        (b[2], b[3], "bondx"),
        (b[3], b[4], "bondw"),
        (b[4], b[5], "bondx"),
    ])

    # 2. Horizontais superiores (para a direita: t_i -> t_{i+1})
    edges.extend([
        (t[0], t[1], "bondx"),
        (t[1], t[2], "bondw"),
        (t[2], t[3], "bondx"),
        (t[3], t[4], "bondw"),
        (t[4], t[5], "bondx"),
    ])

    # 3. Ligações y: inferior -> superior (b_i -> t_i)
    for i in range(6):
        edges.append((b[i], t[i], "bondy"))

    # 4. Ligações z: superior -> inferior à direita (t_i -> b_{i+1})
    for i in range(5):
        edges.append((t[i], b[i + 1], "bondz"))

    # 5. Ligações t: inferior -> superior à direita (b_i -> t_{i+1})
    for i in range(5):
        edges.append((t[i + 1], b[i], "bondt"))

    # Desenho de todas as arestas direcionadas
    for p_start, p_end, bond_class in edges:
        arrow = FancyArrowPatch(
            posA=p_start,
            posB=p_end,
            arrowstyle=fancy_style,
            linewidth=line_width,
            color=colors[bond_class],
            shrinkA=arrow_shrink_pt,
            shrinkB=arrow_shrink_pt,
            capstyle="round",
            joinstyle="round",
            zorder=2,
        )
        ax.add_patch(arrow)

    # -----------------------------------------------------------------------
    # Desenho dos nós (círculos pretos preenchidos de tamanho uniforme)
    # -----------------------------------------------------------------------
    all_nodes = list(b.values()) + list(t.values())
    node_x, node_y = zip(*all_nodes)
    ax.plot(
        node_x,
        node_y,
        marker="o",
        color="black",
        linestyle="None",
        markersize=node_diameter_pt,
        zorder=3,
    )

    # -----------------------------------------------------------------------
    # Rótulos matemáticos:
    # - "1" abaixo do nó b2 (4, 0)
    # - "2" abaixo do nó b3 (6, 0)
    # - "4" acima do nó t2 (5, 2)
    # - "3" acima do nó t3 (7, 2)
    # - Reticências (\cdots) em (-0.55, 1.0) e (11.55, 1.0)
    # -----------------------------------------------------------------------
    label_fontsize = 13
    ellipsis_fontsize = 15

    ax.annotate(
        r"$1$",
        xy=b[2],
        xytext=(0, -7),
        textcoords="offset points",
        ha="center",
        va="top",
        fontsize=label_fontsize,
        color="black",
        zorder=4,
    )
    ax.annotate(
        r"$2$",
        xy=b[3],
        xytext=(0, -7),
        textcoords="offset points",
        ha="center",
        va="top",
        fontsize=label_fontsize,
        color="black",
        zorder=4,
    )
    ax.annotate(
        r"$4$",
        xy=t[2],
        xytext=(0, 7),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=label_fontsize,
        color="black",
        zorder=4,
    )
    ax.annotate(
        r"$3$",
        xy=t[3],
        xytext=(0, 7),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=label_fontsize,
        color="black",
        zorder=4,
    )

    ax.text(
        -0.55,
        1.0,
        r"$\cdots$",
        ha="center",
        va="center",
        fontsize=ellipsis_fontsize,
        color="black",
        zorder=4,
    )
    ax.text(
        11.55,
        1.0,
        r"$\cdots$",
        ha="center",
        va="center",
        fontsize=ellipsis_fontsize,
        color="black",
        zorder=4,
    )

    # -----------------------------------------------------------------------
    # Ajustes finais da imagem
    # -----------------------------------------------------------------------
    ax.set_aspect("equal")
    ax.axis("off")

    # Limites das coordenadas com margem adequada para reticências e rótulos
    ax.set_xlim(-1.2, 12.2)
    ax.set_ylim(-0.8, 2.8)

    # Salva com fundo estritamente transparente
    fig.savefig(
        output_path,
        transparent=True,
        bbox_inches="tight",
        pad_inches=0.08,
        dpi=dpi,
    )
    print(f"[OK] Grafo exportado com sucesso: {output_path}")

    return fig, ax


if __name__ == "__main__":
    draw_qsl_lattice_graph(output_path="qsl_lattice_graph.png")
