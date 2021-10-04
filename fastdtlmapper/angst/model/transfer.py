from dataclasses import dataclass
from typing import Union


@dataclass
class Transfer:
    """Transfer DataClass"""

    donor_node_id: Union[int, str]
    recipient_node_id: Union[int, str]
    direction: bool = True

    @property
    def as_text(self) -> str:
        """Transfer as text format (e.g. 1 -> 2)"""
        if self.direction:
            return f"{self.donor_node_id} -> {self.recipient_node_id}"
        else:
            return f"{self.donor_node_id} <-> {self.recipient_node_id}"
