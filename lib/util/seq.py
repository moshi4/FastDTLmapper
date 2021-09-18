from pathlib import Path
from typing import List, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def count_fasta_seq(fasta_file: str) -> int:
    """Count number of fasta sequence

    Args:
        fasta_file (str): Fasta file path

    Returns:
        int: Sequence number
    """
    return len(list(SeqIO.parse(fasta_file, "fasta")))


def add_fasta_annotation(
    fasta_infile: Union[str, Path],
    fasta_outfile: Union[str, Path],
    id_prefix: str = "GENE",
) -> None:
    """Add annotation to fasta

    Args:
        fasta_infile (Union[str, Path]): Input fasta file path
        fasta_outfile (Union[str, Path]): Output annotated fasta file path
        id_prefix (str, optional): Annotate id prefix. Defaults to "GENE".
    """
    fix_records = []
    for cnt, record in enumerate(SeqIO.parse(fasta_infile, format="fasta"), 1):
        unique_gene_id = f"{id_prefix}{cnt:05d}"
        record.id = f"{unique_gene_id} {record.id} "
        record.description = f"{unique_gene_id} {record.description} "
        fix_records.append(record)
    SeqIO.write(fix_records, fasta_outfile, "fasta-2line")


def gbk2cds_fasta(
    gbk_infile: str, fasta_outfile: str, seqtype: str, id_prefix: str = "GENE"
) -> None:
    """Convert genbank to CDS fasta file

    Args:
        gbk_infile (str): Input genbank file
        fasta_outfile (str): Output CDS fasta file
        seqtype (str): Output sequence type ("nucleotide" or "protein")
        id_prefix (str, optional): sequence serial id prefix
    """
    gene_idx = 0
    cds_seq_record_list: List[SeqRecord] = []
    for record in SeqIO.parse(gbk_infile, "genbank"):
        cds_feature_list = [f for f in record.features if f.type == "CDS"]
        for cds_feature in cds_feature_list:
            qualifiers = cds_feature.qualifiers
            protein_id = qualifiers.get("protein_id", ["NA"])[0]
            product = qualifiers.get("product", ["NA"])[0]

            if qualifiers.get("translation") is None:
                continue

            gene_idx += 1
            if seqtype == "nucleotide":
                seq_id = f"{id_prefix}{gene_idx:05d}"
                cds_seq = cds_feature.location.extract(record.seq)
            elif seqtype == "protein":
                seq_id = f"{id_prefix}{gene_idx:05d}|{protein_id}"
                cds_seq = Seq(qualifiers.get("translation")[0])
            else:
                raise KeyError(f"seqtype '{seqtype}' is not 'nucleotide|protein'")

            cds_seq_record = SeqRecord(seq=cds_seq, id=seq_id, description=product)
            cds_seq_record_list.append(cds_seq_record)

    SeqIO.write(cds_seq_record_list, fasta_outfile, "fasta-2line")
