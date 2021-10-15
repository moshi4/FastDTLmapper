from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Union

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree

from fastdtlmapper.util import UtilFasta, UtilGenbank, UtilTree


@dataclass
class InputCheck:
    """User input check Class"""

    fasta_or_genbank_dir: Union[str, Path]
    nwk_tree_file: Union[str, Path]
    valid_fasta_suffix_list: List[str] = field(
        default_factory=lambda: [".fa", ".faa", ".fasta"]
    )
    valid_genbank_suffix_list: List[str] = field(
        default_factory=lambda: [".gb", ".gbk", ".genbank"]
    )
    invalid_char_list: List[str] = field(
        default_factory=lambda: ["'", '"', "-", "_", "|"]
    )

    def __post_init__(self):
        self.fasta_or_genbank_dir = Path(self.fasta_or_genbank_dir)
        self.nwk_tree_file = Path(self.nwk_tree_file)

    def run(self):
        """Run all check process"""
        fasta_or_genbank_file_list = self._get_fasta_or_genbank_file_list()
        self._is_valid_fasta_or_genbank_format(fasta_or_genbank_file_list)
        self._is_valid_filename(fasta_or_genbank_file_list)
        self._check_rooted_tree()
        self._check_tree_seqfile_consistency(fasta_or_genbank_file_list)

    def _get_fasta_or_genbank_file_list(self) -> List[Path]:
        """Get fasta or genbank file list"""
        fasta_or_genbank_file_list = []
        for file in Path(self.fasta_or_genbank_dir).glob("*"):
            if (
                file.suffix in self.valid_fasta_suffix_list
                or file.suffix in self.valid_genbank_suffix_list
            ):
                fasta_or_genbank_file_list.append(file)
        return fasta_or_genbank_file_list

    def _is_valid_fasta_or_genbank_format(
        self, fasta_or_genbank_file_list: List[Path]
    ) -> None:
        """Check fasta|genbank format or not"""
        for file in fasta_or_genbank_file_list:
            # Check valid fasta file or not
            if file.suffix in self.valid_fasta_suffix_list:
                if not UtilFasta(file).is_valid_format:
                    print(f"ERROR: Input file '{file}' is invalid fasta format!!")
                    exit(1)
            # Check valid genbank file or not
            if file.suffix in self.valid_genbank_suffix_list:
                if not UtilGenbank(file).is_valid_format:
                    print(f"ERROR: Input file '{file}' is invalid genbank format!!")
                    exit(1)

    def _is_valid_filename(self, fasta_or_genbank_file_list: List[Path]) -> None:
        invalid_filename_list = []
        for file in fasta_or_genbank_file_list:
            for invalid_char in self.invalid_char_list:
                if file.name.count(invalid_char):
                    invalid_filename_list.append(file.name)
        for invalid_filename in invalid_filename_list:
            print(f"ERROR: Invalid filename '{invalid_filename}' detected!!")
        if len(invalid_filename_list) != 0:
            print(f"ERROR: {self.invalid_char_list} characters are invaild!!")
            exit(1)

    def _check_rooted_tree(self) -> None:
        """Check rooted tree or not"""
        if not UtilTree(self.nwk_tree_file).is_rooted_tree():
            print(f"ERROR: Input file '{self.nwk_tree_file}' is not rooted tree!!")
            exit(1)

    def _check_tree_seqfile_consistency(self, fasta_or_genbank_file_list: List[Path]):
        """Check tree vs seqfile consistency"""
        seq_filename_list = [f.stem for f in fasta_or_genbank_file_list]
        seq_file_num = len(seq_filename_list)

        tree: Tree = Phylo.read(self.nwk_tree_file, "newick")
        species_name_list = [leaf.name for leaf in tree.get_terminals()]
        species_num = len(species_name_list)
        # Check number of seqfile and tree species number
        if seq_file_num != species_num:
            print(f"ERROR: Number of input fasta or genbank file = {seq_file_num}")
            print(f"ERROR: Number of input newick tree species = {species_num}")
            print("ERROR: Both number must be equal!!")
            exit(1)
        # Check seqfile names & species names match
        if len(set(seq_filename_list) - set(species_name_list)) != 0:
            print("ERROR: Input filenames and newick tree species name must match!!")
            exit(1)
