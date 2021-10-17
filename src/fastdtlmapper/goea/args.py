import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


@dataclass
class Args:
    """Parse Arguments Class"""

    indir: Path
    plot_pvalue_thr: float
    plot_max_num: int
    plot_format: str
    plot_color: str
    use_adjusted_pvalue: bool


def get_args(argv: Optional[List[str]] = None) -> Args:
    """Get argument values

    Returns:
        argparse.Namespace: argument values
    """
    parser = argparse.ArgumentParser(
        description="Gain/Loss genes GOEA(GO Enrichment Analysis) tool"
    )
    parser.add_argument(
        "-i",
        "--indir",
        required=True,
        type=Path,
        help="FastDTLmapper result directory",
        metavar="IN",
    )
    default_plot_pvalue_thr = 0.05
    parser.add_argument(
        "--plot_pvalue_thr",
        type=float,
        default=default_plot_pvalue_thr,
        help=f"Plot GOterm pvalue threshold (Default: {default_plot_pvalue_thr})",
        metavar="",
    )
    default_max_plot_num = 10
    parser.add_argument(
        "--plot_max_num",
        type=int,
        default=default_max_plot_num,
        help=f"Plot GOterm max number (Default: {default_max_plot_num})",
        metavar="",
    )
    default_plot_format = "png"
    available_format = ["png", "svg", "jpg", "pdf"]
    parser.add_argument(
        "--plot_format",
        type=str,
        default=default_plot_format,
        choices=available_format,
        help=f"Plot file format [{'|'.join(available_format)}] "
        + f"(Default: '{default_plot_format}')",
        metavar="",
    )
    parser.add_argument(
        "--plot_color",
        type=str,
        default="",
        help="Plot specified hexcolor [e.g. '1affdb'] "
        + "(Default: yellow to red gradient color)",
        metavar="",
    )
    parser.add_argument(
        "--adjusted_pvalue",
        help="Use BH adjusted pvalue for plot threshold",
        action="store_true",
    )

    args = parser.parse_args(argv)

    # Plot hexcolor check
    def is_valid_hexcolor(hexcolor: str) -> bool:
        search_result = re.search(r"^(?:[0-9a-fA-F]{3}){1,2}$", hexcolor)
        return False if search_result is None else True

    if args.plot_color and not is_valid_hexcolor(args.plot_color):
        parser.error(
            f"argument --plot_color: invalid hexcolor code '{args.plot_color}'."
        )
    else:
        args.plot_color = "#" + args.plot_color if args.plot_color else ""

    return Args(
        args.indir,
        args.plot_pvalue_thr,
        args.plot_max_num,
        args.plot_format,
        args.plot_color,
        args.adjusted_pvalue,
    )
