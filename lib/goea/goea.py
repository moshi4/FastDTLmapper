import math
import os
import subprocess as sp
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
from urllib import request

import numpy as np
import pandas as pd
from goatools.godag_obosm import OboToGoDagSmall
from goatools.godag_plot import GODagPltVars, GODagSmallPlot
from goatools.obo_parser import GODag


@dataclass
class GOEA:
    """GO Enrichment Analysis Class"""

    target_list_file: Path
    all_list_file: Path
    association_file: Path
    obo_file: Path = Path("go-basic.obo")
    pvalue_thr: float = 0.05
    plot_max_num: int = 10
    plot_format: str = "png"
    use_adjusted_pvalue: bool = False

    def run(self, output_prefix: Path) -> List[Path]:
        """Run GOEA using goatools

        Args:
            output_prefix (Path): Output files prefix path

        Returns:
            List[Path]: Output file path list ({prefix}_[BP|MF|CC].tsv)
        """
        outdir = output_prefix.parent
        os.makedirs(outdir, exist_ok=True)

        goea_result_file_list = []
        for go_category in ("BP", "MF", "CC"):
            goea_result_file = Path(f"{output_prefix}_{go_category}.tsv")
            goea_cmd = (
                f"find_enrichment.py {self.target_list_file} {self.all_list_file} "
                + f"{self.association_file} --obo {self.obo_file} --ns {go_category} "
                + "--pval=-1 --method=fdr_bh --pval_field=fdr_bh "
                + f"--outfile={goea_result_file} "
            )
            print("=" * 100)
            sp.run(goea_cmd, shell=True)
            goea_result_file_list.append(goea_result_file)

        return goea_result_file_list

    def plot(
        self,
        goea_result_file: Path,
        output_prefix: Path,
    ) -> None:
        """Plot GOEA significant GOterms

        Args:
            goea_result_file (Path): GOEA result file path
            output_prefix (Path): Output files prefix path
        """
        for goea_type in ("enrichment", "purified"):
            # Extract goterm & pvalue
            goterm2pvalue = self._extract_goterm2pvalue(goea_result_file, goea_type)
            if len(goterm2pvalue) == 0:
                return
            # Get hexcolor from pvalue for color plot
            pvalue_abs_log10_list = [abs(math.log10(v)) for v in goterm2pvalue.values()]
            pvalue_hexcolor_list = self._convert_hexcolor_gradient(
                pvalue_abs_log10_list
            )
            goterm2hexcolor = {}
            for goterm, hexcolor in zip(goterm2pvalue.keys(), pvalue_hexcolor_list):
                goterm2hexcolor[goterm] = hexcolor
            # Plot GOterm with gradient color
            plot_outfile = Path(f"{output_prefix}_{goea_type}.{self.plot_format}")
            self._color_plot(plot_outfile, goterm2hexcolor, goterm2pvalue)

    def _extract_goterm2pvalue(
        self,
        goea_result_file: Path,
        extract_type: str,
    ) -> Dict[str, float]:
        """Extract GOterm & Pvalue from GOEA result file

        Args:
            goea_result_file (Path): GOEA result file
            extract_type (str): "enrichment" or "purified"

        Returns:
            Dict[str, float]: GOterm & Pvalue dict
        """
        if extract_type not in ("enrichment", "purified"):
            raise ValueError("extract_type must be 'enrichment' or 'purified'!!")

        pvalue_column = "p_fdr_bh" if self.use_adjusted_pvalue else "p_uncorrected"

        df = pd.read_table(goea_result_file)
        if extract_type == "enrichment":
            extract_df = df[
                (df["enrichment"] == "e") & (df[pvalue_column] < self.pvalue_thr)
            ]
        elif extract_type == "purified":
            extract_df = df[
                (df["enrichment"] == "p") & (df[pvalue_column] < self.pvalue_thr)
            ]
        extract_df = extract_df.head(self.plot_max_num)

        goterm2pvalue = {}
        for goterm, pvalue in zip(extract_df["# GO"], extract_df[pvalue_column]):
            goterm2pvalue[goterm] = pvalue

        return goterm2pvalue

    def _convert_hexcolor_gradient(
        self,
        value_list: List[float],
        default_min: Optional[float] = 0.05,
        default_max: Optional[float] = 10,
    ) -> List[str]:
        """Convert to gradient hexcolor (e.g. "#ffff00") list based on the size of the value

        Args:
            value_list (List[float]): List of float values
            default_min (float, optional): Default minimum value. Defaults to 0.05.
            default_max (float, optional): Default maximum value. Defaults to 10.

        Returns:
            List[str]: List of gradient hexcolor
        """
        # Make yellow(ffff00) -> red(ff0000) hexcolor gradient list
        step_num = 5
        hexcolor_gradient_list = [
            f"#ff{i:02x}00" for i in reversed(range(0, 256, step_num))
        ]

        # Define value range for gradient hexcolor
        min_value, max_value = min(value_list), max(value_list)
        if default_min:
            min_value = min(min_value, default_min)
        if default_max:
            max_value = max(max_value, default_max)
        value_range_list = list(
            np.linspace(min_value, max_value, int(255 / step_num) + 1)
        )[1:]

        # Get gradient hexcolor corresponding to input value
        convert_hexcolor_gradient_list = []
        for value in value_list:
            for idx, value_range in enumerate(value_range_list):
                if value <= value_range:
                    convert_hexcolor_gradient_list.append(hexcolor_gradient_list[idx])
                    break

        return convert_hexcolor_gradient_list

    def _color_plot(
        self,
        plot_outfile: str,
        goid2color: Dict[str, str],
        goid2pvalue: Dict[str, float] = {},
    ) -> None:
        """Plot GO DAG using self-defined GO color

        Args:
            plot_outfile (str): Output plot file path
            goid2color (Dict[str, str]): go id and hexcolor dict
            goid2pvalue (Dict[str, float], optional): go id and pvalue dict
        """
        # Get plot target GO DAG
        obodag = GODag(self.obo_file)
        godagsmall = OboToGoDagSmall(goids=goid2color.keys(), obodag=obodag).godag

        # Wrapping GO description line at appropriate location
        for v in godagsmall.go2obj.values():
            if len(v.name) < 20 or v.name.count("\n"):
                continue
            split_word = v.name.split(" ")
            split_cnt = math.ceil(len(split_word) / 2) - 1
            line_wrap_name = ""
            for cnt, word in enumerate(split_word):
                newline_or_space = "\n" if cnt == split_cnt else " "
                line_wrap_name += word + newline_or_space
            v.name = line_wrap_name

        # Add pvalue to the end of GO description
        for v in godagsmall.go2obj.values():
            if v.id in goid2pvalue.keys():
                v.name += f"\n{goid2pvalue[v.id]:.2e}"

        # Suppress useless header plot (e.g. L2 D2)
        godag_plg_vars = GODagPltVars()
        godag_plg_vars.fmthdr = "{GO}"

        # Create plot obj & add plot color
        godagplot = GODagSmallPlot(
            godagsmall, abodag=obodag, GODagPltVars=godag_plg_vars
        )
        godagplot.goid2color = goid2color

        # Plot color go dag
        godagplot.plt(plot_outfile, "pydot")

    @staticmethod
    def download_obo(obo_outfile: Path) -> None:
        """Download GO OBO file

        Args:
            obo_outfile (Path): Output OBO file path
        """
        if obo_outfile.is_file():
            return

        obo_download_url = "http://geneontology.org/ontology/go-basic.obo"
        with request.urlopen(obo_download_url) as u, open(obo_outfile, "wb") as f:
            f.write(u.read())

        if not obo_outfile.is_file():
            raise FileNotFoundError(f"'{obo_outfile}' not found. Download failure!!")
