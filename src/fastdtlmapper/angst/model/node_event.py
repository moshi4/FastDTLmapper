from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

from fastdtlmapper.angst.model.transfer import Transfer


@dataclass
class NodeEvent:
    """Node Event DataClass"""

    node_id: str
    brn_num: int = 0
    dup_num: int = 0
    los_num: int = 0
    trn_num: int = 0
    trn_detail_list: List[Transfer] = field(default_factory=list)
    gene_num: Optional[int] = None

    @property
    def gain_num(self) -> int:
        """Number of (Brn + Dup + Trn) event"""
        return self.brn_num + self.dup_num + self.trn_num

    @property
    def as_dtl_text(self) -> str:
        """DTL event text"""
        text = ""
        text += f"{self.gene_num} ["
        text += f"brn={self.brn_num} "
        text += f"dup={self.dup_num} "
        text += f"los={self.los_num} "
        text += f"trn={self.trn_num}]"
        return text

    @property
    def as_gain_loss_text(self) -> str:
        """Gain-Loss event text"""
        text = ""
        text += f"{self.gene_num} ["
        text += f"gain={self.gain_num} "
        text += f"los={self.los_num}]"
        return text

    @property
    def as_csv_format(self) -> str:
        """Node event csv format text"""
        trn_detail_text = ""
        for trn_detail in self.trn_detail_list:
            trn_detail_text += trn_detail.as_text + "|"
        trn_detail_text = trn_detail_text.rstrip("|")
        csv_text = (
            f"{self.node_id},{self.gene_num},{self.gain_num},{self.brn_num},"
            + f"{self.dup_num},{self.trn_num},{self.los_num},{trn_detail_text}"
        )
        return csv_text

    @property
    def as_tsv_format(self) -> str:
        """Node event tsv format text"""
        csv_text = self.as_csv_format
        return csv_text.replace(",", "\t")

    def __add__(self, other: NodeEvent):
        if not isinstance(other, self.__class__) and self.node_id != other.node_id:
            raise NotImplementedError()

        ret_obj = NodeEvent(self.node_id)
        ret_obj.brn_num = self.brn_num + other.brn_num
        ret_obj.dup_num = self.dup_num + other.dup_num
        ret_obj.los_num = self.los_num + other.los_num
        ret_obj.trn_num = self.trn_num + other.trn_num
        if self.gene_num is None or other.gene_num is None:
            ret_obj.gene_num = None
        else:
            ret_obj.gene_num = self.gene_num + other.gene_num
        ret_obj.trn_detail_list.extend(self.trn_detail_list)
        ret_obj.trn_detail_list.extend(other.trn_detail_list)

        return ret_obj
