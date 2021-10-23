import argparse
import re
from dataclasses import dataclass
from typing import List, Optional

from ete3 import NodeStyle, TextFace, Tree, TreeStyle


@dataclass
class Args:
    """Arguments DataClass"""

    map_tree_infile: str
    plot_outfile: str
    plot_scale: int
    plot_margin: int
    plot_width: Optional[int]
    title: Optional[str]
    ladderize: bool
    add_brn_dup_trn: bool
    edit_mode: bool
    color_gene_num: str
    color_gain_num: str
    color_loss_num: str
    gain_symbol: str
    loss_symbol: str
    fsize_node_name: int
    fsize_leaf_name: int
    fsize_gene_num: int
    fsize_gain_num: int
    fsize_loss_num: int
    fsize_title: int
    fsize_legend: int


def main(args: Args = None):
    """Plot gain(brn/dup/trn)/loss map main function"""
    # Get command line arguments
    if args is None:
        args = get_args()

    # Convert to ultrametric tree for evenly spaced tree plot
    # (Adjust minimum branch length = 1.0)
    tree = Tree(args.map_tree_infile, format=1)
    max_depth = max([len(n.get_ancestors()) for n in tree.iter_leaves()])
    tree.convert_to_ultrametric(tree_length=max_depth, strategy="balanced")

    if args.ladderize:
        tree.ladderize()

    for node in tree.traverse():
        # Get node name, gene number, gain/loss number
        pattern = r"^(.+) \| (\d+) \[brn=(\d+) dup=(\d+) los=(\d+) trn=(\d+)\]"
        result = re.match(pattern, node.name)
        if result is None or None in result.groups():
            print(f"Input '{args.map_tree_infile}' is invalid dtl map file!!")
            exit(1)
        node_name, gene_num, brn_num, dup_num, los_num, trn_num = result.groups()
        node.name = node_name

        # Make gain information
        gain_num = int(brn_num) + int(dup_num) + int(trn_num)
        if args.add_brn_dup_trn:
            gain_info = f"{gain_num} ({brn_num}/{dup_num}/{trn_num})"
        else:
            gain_info = str(gain_num)

        # Define NodeStyle (no display style)
        node_style = NodeStyle()
        node_style["size"] = 0
        node_style["shape"] = "square"
        node.set_style(node_style)

        # Define node name TextFace (normal string style)
        node_name_face = TextFace(node_name, bold=True, fstyle="italic")
        # Define gene num TextFace (color rectangle style)
        gene_num_face = TextFace(" " + gene_num, fsize=args.fsize_gene_num, bold=True)
        gene_num_face.background.color = args.color_gene_num
        gene_num_face.border.width = 1
        gene_num_face.margin_left, gene_num_face.margin_right = 0, 2
        # Define gain num TextFace (color string with symbol)
        gain_num_face = TextFace(
            args.gain_symbol + gain_info,
            fgcolor=args.color_gain_num,
            fsize=args.fsize_gain_num,
        )
        gain_num_face.margin_bottom, gain_num_face.margin_right = 0, 2
        # Define gene num TextFace (color string with symbol)
        los_num_face = TextFace(
            args.loss_symbol + los_num,
            fgcolor=args.color_loss_num,
            fsize=args.fsize_loss_num,
        )
        los_num_face.margin_top, los_num_face.margin_right = 0, 2

        # Add defined TextFaces to Node(or Leaf)
        if node.is_leaf():
            node_name_face.fsize = args.fsize_leaf_name
            node_name_face.margin_left, node_name_face.margin_right = 3, 3
            node.add_face(gain_num_face, column=0, position="float")
            node.add_face(los_num_face, column=0, position="float")
            node.add_face(gene_num_face, column=1, position="float")
            node.add_face(node_name_face, column=0, position="branch-right")
        else:
            node_name_face.fsize = args.fsize_node_name
            node.add_face(gain_num_face, column=0, position="float")
            node.add_face(los_num_face, column=0, position="float")
            node.add_face(TextFace(""), column=1, position="float")
            node.add_face(gene_num_face, column=1, position="float")
            node.add_face(node_name_face, column=1, position="float")

    # Define TreeStyle
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_scale = False
    ts.branch_vertical_margin = args.plot_margin
    ts.scale = args.plot_scale

    # Define Title
    if args.title:
        ts.title.add_face(
            TextFace(args.title, fsize=args.fsize_title, bold=True), column=0
        )
        ts.title.add_face(TextFace(""), column=0)

    # Define Legend
    ts.legend_position = 3
    gene_legend_face = TextFace("Gene Number", fsize=args.fsize_legend, bold=True)
    gene_legend_face.border.width = 1
    gene_legend_face.background.color = args.color_gene_num
    gene_legend_face.margin_left, gene_legend_face.margin_right = 5, 0
    if args.add_brn_dup_trn:
        gain_legend_text = "Gain(Brn/Dup/Trn) Number"
    else:
        gain_legend_text = "Gain Number"
    gain_legend_face = TextFace(
        args.gain_symbol + gain_legend_text,
        fsize=args.fsize_legend,
        fgcolor=args.color_gain_num,
    )
    loss_legend_face = TextFace(
        args.loss_symbol + "Loss Number\n",
        fsize=args.fsize_legend,
        fgcolor=args.color_loss_num,
    )
    ts.legend.add_face(gene_legend_face, column=0)
    ts.legend.add_face(gain_legend_face, column=0)
    ts.legend.add_face(loss_legend_face, column=0)

    if args.edit_mode:
        # Launch interactive tree editor
        tree.show(tree_style=ts)

    # Plot gain/loss map tree
    tree.render(args.plot_outfile, tree_style=ts, w=args.plot_width)


def get_args(argv: Optional[List[str]] = None) -> Args:
    """Get arguments

    Returns:
        Args: Args Class
    """
    parser = argparse.ArgumentParser(description="Plot gain/loss map tool")

    ###########################################################################
    # Input/Output file path options
    ###########################################################################
    parser.add_argument(
        "-i",
        "--infile",
        required=True,
        type=str,
        help="Input gain/loss map result file",
        metavar="I",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        required=True,
        type=str,
        help="Output plot file (*.png|*.svg|*.pdf)",
        metavar="O",
    )
    ###########################################################################
    # Plot style options
    ###########################################################################
    parser.add_argument(
        "--plot_scale",
        type=int,
        help="Horizontal branch scale plot parameter (Default: Auto Defined)",
        default=None,
        metavar="",
    )
    default_plot_margin = 15
    parser.add_argument(
        "--plot_margin",
        type=int,
        help=f"Vertical margin plot parameter (Default: {default_plot_margin})",
        default=default_plot_margin,
        metavar="",
    )
    default_plot_width = None
    parser.add_argument(
        "--plot_width",
        type=int,
        help=f"Pixel width size plot parameter (Default: {default_plot_width})",
        default=default_plot_width,
        metavar="",
    )
    default_title = None
    parser.add_argument(
        "--title",
        type=str,
        help=f"Plot figure title (Default: {default_title})",
        default=default_title,
        metavar="",
    )
    parser.add_argument(
        "--ladderize",
        help="Plot ladderized style tree (Default: off)",
        action="store_true",
    )
    parser.add_argument(
        "--add_brn_dup_trn",
        help="Add Brn/Dup/Trn plot (Default: off)",
        action="store_true",
    )
    parser.add_argument(
        "--edit_mode",
        help="Launch interactive tree editor (Default: off)",
        action="store_true",
    )
    ###########################################################################
    # Plot color options
    ###########################################################################
    default_color_gene = "lightgreen"
    parser.add_argument(
        "--color_gene_num",
        type=str,
        help=f"Plot color of gene number (Default: '{default_color_gene}')",
        default=default_color_gene,
        metavar="",
    )
    default_color_gain = "red"
    parser.add_argument(
        "--color_gain_num",
        type=str,
        help=f"Plot color of gain number (Default: '{default_color_gain}')",
        default=default_color_gain,
        metavar="",
    )
    default_color_loss = "blue"
    parser.add_argument(
        "--color_loss_num",
        type=str,
        help=f"Plot color of loss number (Default: '{default_color_loss}')",
        default=default_color_loss,
        metavar="",
    )
    ###########################################################################
    # Plot gain/loss symbol options
    ###########################################################################
    default_gain_symbol = " ▲ "
    parser.add_argument(
        "--gain_symbol",
        type=str,
        help=f"Plot symbol of gain gene (Default: '{default_gain_symbol}')",
        default=default_gain_symbol,
        metavar="",
    )
    default_loss_symbol = " ▼ "
    parser.add_argument(
        "--loss_symbol",
        type=str,
        help=f"Plot symbol of loss gene (Default: '{default_loss_symbol}')",
        default=default_loss_symbol,
        metavar="",
    )
    ###########################################################################
    # Plot font size options
    ###########################################################################
    default_fsize_node_name = 8
    parser.add_argument(
        "--fsize_node_name",
        type=int,
        help=f"Plot font size of node name (Default: {default_fsize_node_name})",
        default=default_fsize_node_name,
        metavar="",
    )
    default_fsize_leaf_name = 10
    parser.add_argument(
        "--fsize_leaf_name",
        type=int,
        help=f"Plot font size of leaf name (Default: {default_fsize_leaf_name})",
        default=default_fsize_leaf_name,
        metavar="",
    )
    default_fsize_gene_num = 8
    parser.add_argument(
        "--fsize_gene_num",
        type=int,
        help=f"Plot font size of gene number (Default: {default_fsize_gene_num})",
        default=default_fsize_gene_num,
        metavar="",
    )
    default_fsize_gain_num = 8
    parser.add_argument(
        "--fsize_gain_num",
        type=int,
        help=f"Plot font size of gain number (Default: {default_fsize_gain_num})",
        default=default_fsize_gain_num,
        metavar="",
    )
    default_fsize_loss_num = 8
    parser.add_argument(
        "--fsize_loss_num",
        type=int,
        help=f"Plot font size of loss number (Default: {default_fsize_loss_num})",
        default=default_fsize_loss_num,
        metavar="",
    )
    default_fsize_title = 30
    parser.add_argument(
        "--fsize_title",
        type=int,
        help=f"Plot font size of title (Default: {default_fsize_title})",
        default=default_fsize_title,
        metavar="",
    )
    default_fsize_legend = 8
    parser.add_argument(
        "--fsize_legend",
        type=int,
        help=f"Plot font size of legend (Default: {default_fsize_legend})",
        default=default_fsize_legend,
        metavar="",
    )

    args = parser.parse_args(argv)

    # Define plot scale by whether brn/dup/trn plotting or not
    if args.plot_scale is None:
        if args.add_brn_dup_trn:
            args.plot_scale = 150
        else:
            args.plot_scale = 80

    # Outfile extension validation check
    if args.outfile[-4:] not in (".png", ".svg", ".pdf"):
        parser.error("-o/--outfile extension must be *.png|*.svg|*.pdf")

    return Args(
        args.infile,
        args.outfile,
        args.plot_scale,
        args.plot_margin,
        args.plot_width,
        args.title,
        args.ladderize,
        args.add_brn_dup_trn,
        args.edit_mode,
        args.color_gene_num,
        args.color_gain_num,
        args.color_loss_num,
        args.gain_symbol,
        args.loss_symbol,
        args.fsize_node_name,
        args.fsize_leaf_name,
        args.fsize_gene_num,
        args.fsize_gain_num,
        args.fsize_loss_num,
        args.fsize_title,
        args.fsize_legend,
    )


if __name__ == "__main__":
    main()
