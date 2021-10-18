from dataclasses import dataclass
from pathlib import Path
from typing import List, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@dataclass
class UtilGenbank:
    """Sequence Utility Class"""

    genbank_file: Union[str, Path]

    def __post_init__(self):
        self.genbank_file = Path(self.genbank_file)

    @property
    def is_valid_format(self) -> bool:
        """Check genbank format file or not"""
        if len(list(SeqIO.parse(self.genbank_file, "genbank"))) == 0:
            return False
        else:
            return True

    def convert_cds_fasta(
        self, fasta_outfile: Union[str, Path], seqtype: str, id_prefix: str
    ) -> None:
        """Convert genbank to CDS fasta file

        Args:
            fasta_outfile (str): Output CDS fasta file
            seqtype (str): Output sequence type ("nucleotide" or "protein")
            id_prefix (str, optional): sequence serial id prefix
        """
        gene_idx = 0
        cds_seq_record_list: List[SeqRecord] = []
        for record in SeqIO.parse(self.genbank_file, "genbank"):
            cds_feature_list = [f for f in record.features if f.type == "CDS"]
            for cds_feature in cds_feature_list:
                qualifiers = cds_feature.qualifiers
                protein_id = qualifiers.get("protein_id", ["NA"])[0]
                product = qualifiers.get("product", ["NA"])[0]
                # Delete invalid characters
                invalid_char_list = ["'", '"', "(", ")", "[", "]", ":", ";", "|", ","]
                transtable = str.maketrans({char: "_" for char in invalid_char_list})
                protein_id = protein_id.translate(transtable)

                if qualifiers.get("translation") is None:
                    continue

                gene_idx += 1
                seq_id = f"{id_prefix}_GENE{gene_idx:06d}"
                if protein_id != "NA":
                    seq_id += f"_{protein_id}"

                if seqtype == "nucleotide":
                    cds_seq = cds_feature.location.extract(record.seq)
                elif seqtype == "protein":
                    cds_seq = Seq(qualifiers.get("translation")[0])
                else:
                    raise KeyError(f"seqtype '{seqtype}' is not 'nucleotide|protein'")

                cds_seq_record = SeqRecord(seq=cds_seq, id=seq_id, description=product)
                cds_seq_record_list.append(cds_seq_record)

        SeqIO.write(cds_seq_record_list, fasta_outfile, "fasta-2line")
