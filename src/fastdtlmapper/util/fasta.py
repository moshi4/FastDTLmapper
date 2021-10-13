from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Union

from Bio import SeqIO


@dataclass
class UtilFasta:
    """Sequence Utility Class"""

    fasta_file: Union[str, Path]

    def __post_init__(self):
        self.fasta_file = Path(self.fasta_file)

    @property
    def seq_num(self) -> int:
        """Count number of fasta sequence"""
        return len(list(SeqIO.parse(self.fasta_file, "fasta")))

    @property
    def uniq_seq_num(self) -> int:
        """Count number of unique fasta sequence"""
        return len(set([r.seq for r in SeqIO.parse(self.fasta_file, "fasta")]))

    @property
    def is_valid_format(self) -> bool:
        """Check fasta format file or not"""
        if self.seq_num == 0:
            return False
        else:
            return True

    @property
    def species_name2seq_num(self) -> Dict[str, int]:
        """Get species_name & seq_num dict"""
        species_name2seq_num = defaultdict(int)
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            species_name = "_".join(record.id.split("_")[:-1])
            species_name2seq_num[species_name] += 1
        return species_name2seq_num

    @property
    def uniq_species_name_list(self) -> List[str]:
        """Get fasta unique species name list"""
        species_name_list = []
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            species_name = "_".join(record.id.split("_")[:-1])
            species_name_list.append(species_name)
        return sorted(list(set(species_name_list)))

    def add_serial_id(
        self,
        fasta_outfile: Union[str, Path],
        id_prefix: str = "GENE",
    ) -> None:
        """Add fasta id prefix tag

        Args:
            fasta_outfile (Union[str, Path]): Output serial id tagged fasta file path
            id_prefix (str, optional): Id prefix tag. Defaults to "GENE".
        """
        fix_records = []
        for cnt, record in enumerate(SeqIO.parse(self.fasta_file, format="fasta"), 1):
            serial_id_tag = f"{id_prefix}{cnt:06d}"
            record.id = f"{serial_id_tag} {record.id} "
            record.description = f"{serial_id_tag} {record.description} "
            fix_records.append(record)
        SeqIO.write(fix_records, fasta_outfile, "fasta-2line")
