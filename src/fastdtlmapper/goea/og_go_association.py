import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set


@dataclass
class OgGoAssociation:
    """OrthologGroup(OG) & GOterms Association Class"""

    og_gene_list_file: Path
    annotation_file_list: List[Path]
    go_define_ratio_thr: float = 0.5

    def load_og_go_association(self):
        """Load OG & GOterms association"""
        # Get OrthologGroup(OG) ID & Gene ID association
        self.og_id2gene_id_list = self._get_og_id2gene_id_association()
        # Get Gene ID & GO ID association
        self.gene_id2go_id_set = self._get_gene_id2go_id_association()
        # Get OG ID & GO ID association
        self.og_id2go_id_list = self._get_og_id2go_id_association(
            self.og_id2gene_id_list, self.gene_id2go_id_set
        )

    def _get_og_id2gene_id_association(self) -> Dict[str, List[str]]:
        """Get OG(OrthologGroup) ID & Gene ID association

        Returns:
            Dict[str, List[str]]: OG ID & Gene ID association dict
        """
        og_id2gene_id_list = {}
        with open(self.og_gene_list_file) as f:
            line_list = f.read().splitlines()
        for line in line_list:
            og_id = line[0:9]
            gene_id_list = line[11:].split(" ")
            og_id2gene_id_list[og_id] = gene_id_list
        return og_id2gene_id_list

    def _get_gene_id2go_id_association(self) -> Dict[str, Set[str]]:
        """Get Gene ID & GO ID association

        Returns:
            Dict[str, Set[str]]: Gene ID & GO ID association dict
        """
        gene_id2go_id_set = defaultdict(set)
        for annotation_file in self.annotation_file_list:
            with open(annotation_file) as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    if len(row) < 14 or row[13] == "-":
                        continue
                    gene_id, go_list = row[0], row[13].split("|")
                    gene_id2go_id_set[gene_id] |= set(go_list)
        return gene_id2go_id_set

    def _get_og_id2go_id_association(
        self,
        og_id2gene_id_list: Dict[str, List[str]],
        gene_id2go_id_set: Dict[str, Set[str]],
    ) -> Dict[str, List[str]]:
        """Get OG ID & GO ID association

        Args:
            og_id2gene_id_list (Dict[str, List[str]]): OG ID & Gene ID association dict
            gene_id2go_id_set (Dict[str, Set[str]]): Gene ID & GO ID association dict
            go_define_ratio_thr (float): OG's GO define ratio threshold (0.0 - 1.0)

        Returns:
            Dict[str, List[str]]: OG ID & GO ID association dict
        """
        og_id2go_id_list = defaultdict(list)
        for og_id, gene_id_list in og_id2gene_id_list.items():
            go_id_list = []
            for gene_id in gene_id_list:
                go_id_list.extend(list(gene_id2go_id_set[gene_id]))
            for go_id in set(go_id_list):
                go_id_count = go_id_list.count(go_id)
                if go_id_count >= len(gene_id_list) * self.go_define_ratio_thr:
                    og_id2go_id_list[og_id].append(go_id)
        return og_id2go_id_list

    def write_og2go_association(self, og2go_association_file: Path) -> None:
        """Write OG ID & GO ID association file

        Args:
            og2go_association_file (Path): OG & GO association output file path
        """
        self.load_og_go_association()
        association_info = ""
        for og_id, go_id_list in self.og_id2go_id_list.items():
            og_gene_num = len(self.og_id2gene_id_list[og_id])
            for cnt in range(og_gene_num):
                association_info += f"{og_id}_{cnt:05d}\t{';'.join(go_id_list)}\n"
        with open(og2go_association_file, "w") as f:
            f.write(association_info)
