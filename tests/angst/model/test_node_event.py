import pytest
from fastdtlmapper.angst.model import NodeEvent, Transfer


@pytest.fixture()
def node_event() -> NodeEvent:
    """NodeEvent fixture"""
    return NodeEvent("N00X", 1, 3, 1, 0, [], 5)


def test_gain_num(node_event: NodeEvent):
    """gain_num test"""
    assert node_event.gain_num == 4


def test_as_dtl_text(node_event: NodeEvent):
    """as_dtl_text test"""
    assert node_event.as_dtl_text == "5 [brn=1 dup=3 los=1 trn=0]"


def test_as_gain_loss_text(node_event: NodeEvent):
    """as_gain_loss_text test"""
    assert node_event.as_gain_loss_text == "5 [gain=4 los=1]"


def test_as_csv_format(node_event: NodeEvent):
    """as_csv_format test"""
    assert node_event.as_csv_format == "N00X,5,4,1,3,0,1,"


def test_as_tsv_format(node_event: NodeEvent):
    """as_tsv_format test"""
    assert node_event.as_tsv_format == "N00X\t5\t4\t1\t3\t0\t1\t"


def test_add():
    """__add__ test"""
    node_event1 = NodeEvent("N00X", 1, 3, 1, 0, [Transfer(10, 20)], 5)
    node_event2 = NodeEvent("N00X", 0, 2, 1, 1, [Transfer(30, 40)], 5)

    assert node_event1 + node_event2 == NodeEvent(
        "N00X", 1, 5, 2, 1, [Transfer(10, 20), Transfer(30, 40)], 10
    )
