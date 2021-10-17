import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Union


@dataclass
class OutPath:
    """Output Path Class"""

    rootdir: Union[str, Path]
    goea_mode: bool = False

    def __post_init__(self):
        self.rootdir = Path(self.rootdir)

        # 00. User fasta & tree directory
        self.user_data_dir = self.rootdir / "00_user_data"
        self.user_fasta_dir = self.user_data_dir / "fasta"
        self.user_tree_dir = self.user_data_dir / "tree"

        self.user_tree_file = self.user_tree_dir / "user_tree.nwk"
        self.um_tree_file = self.user_tree_dir / "ultrametric_tree.nwk"
        self.um_nodeid_tree_file = self.user_tree_dir / "ultrametric_nodeid_tree.nwk"

        # 01. OrthoFinder directory
        self.ortho_dir = self.rootdir / "01_orthofinder"
        self.ortho_group_fasta_dir = self.ortho_dir / "Orthogroup_Sequences"
        self.ortho_tmpwork_dir = self.user_fasta_dir / "OrthoFinder"

        self.og_gene_file = self.ortho_dir / "Orthogroups" / "Orthogroups.txt"

        # 02. DTL reconciliation directory
        self.dtl_rec_dir = self.rootdir / "02_dtl_reconciliation"

        # 03. All group node DTL event aggregate map directory
        self.aggregate_map_dir = self.rootdir / "03_aggregate_map_result"

        self.all_dtl_map_file = self.aggregate_map_dir / "all_dtl_map.nwk"
        self.all_gain_loss_map_file = self.aggregate_map_dir / "all_gain_loss_map.nwk"
        self.all_node_event_file = self.aggregate_map_dir / "all_og_node_event.tsv"
        self.all_trn_count_file = self.aggregate_map_dir / "all_transfer_gene_count.tsv"
        self.all_trn_genes_file = self.aggregate_map_dir / "all_transfer_gene_list.tsv"

        # 04. (Optional) GOEA analysis directory
        self.goea_dir = self.rootdir / "04_functional_analysis"

        self.go_annotation_dir = self.goea_dir / "go_annotation"
        self.go_annotation_workdir = self.go_annotation_dir / "work"
        self.go_enrichment_dir = self.goea_dir / "go_enrichment"
        self.result_summary_dir = self.goea_dir / "result_summary"
        self.result_summary_plot_dir = self.result_summary_dir / "significant_go_plot"

        self.obo_file = self.goea_dir / "go-basic.obo"
        self.og2go_association_file = self.go_enrichment_dir / "og2go_association.txt"
        self.significant_go_list_file = (
            self.result_summary_dir / "significant_go_list.tsv"
        )
        self.significant_go_count_file = (
            self.result_summary_dir / "significant_go_count.tsv"
        )

        # Log commands list and command stderr
        self.log_dir = self.rootdir / "log"
        self.parallel_cmds_dir = self.log_dir / "parallel_cmds"

        self.run_config_log_file = self.log_dir / "run_config.log"

        self.tmp_parallel_cmds_file = self.parallel_cmds_dir / "tmp_parallel_cmds.txt"

        self.mafft_output_log_file = self.parallel_cmds_dir / "mafft_output_log.csv"
        self.trimal_output_log_file = self.parallel_cmds_dir / "trimal_output_log.csv"
        self.iqtree_output_log_file = self.parallel_cmds_dir / "iqtree_output_log.csv"
        self.angst_output_log_file = self.parallel_cmds_dir / "angst_output_log.csv"

        self._makedirs()

    def _makedirs(self) -> None:
        """Create config output directories"""
        # Delete previous work dir
        shutil.rmtree(self.ortho_tmpwork_dir, ignore_errors=True)
        shutil.rmtree(self.go_enrichment_dir, ignore_errors=True)
        shutil.rmtree(self.result_summary_dir, ignore_errors=True)

        # Default output dir
        target_dir_list = [
            self.user_fasta_dir,
            self.user_tree_dir,
            self.ortho_dir,
            self.dtl_rec_dir,
            self.aggregate_map_dir,
            self.parallel_cmds_dir,
        ]
        if self.goea_mode:
            # Check rootdir is FastDTLmapper result dir or not
            if not self.user_fasta_dir.exists():
                err_msg = f"Input fasta directory '{self.user_fasta_dir}' not found!!\n"
                err_msg += "Please specify FastDTLmapper result directory as input."
                raise ValueError(err_msg)
            # Add GOEA output dir
            target_dir_list.extend(
                [
                    self.go_annotation_workdir,
                    self.go_enrichment_dir,
                    self.result_summary_plot_dir,
                ]
            )

        # Make configured output directory
        for target_dir in target_dir_list:
            os.makedirs(target_dir, exist_ok=True)
