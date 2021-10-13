from fastdtlmapper.angst.model import NodeEvent, Transfer


def test_gain_num():
    """gain_num test"""
    node_event = NodeEvent(0, 1, 1, 1, 0, [], 0)
    assert node_event.gain_num == 2


def test_as_dtl_text():
    """as_dtl_text test"""
    node_event = NodeEvent(0, 1, 3, 1, 0, [], 5)
    assert node_event.as_dtl_text == "5 [brn=1 dup=3 los=1 trn=0]"


def test_as_gain_loss_text():
    """as_gain_loss_text test"""
    node_event = NodeEvent(0, 1, 3, 1, 0, [], 5)
    assert node_event.as_gain_loss_text == "5 [gain=4 los=1]"


def test_as_csv_format():
    """as_csv_format test"""
    node_event = NodeEvent(0, 1, 3, 1, 0, [], 5)
    assert node_event.as_csv_format == "0,5,4,1,3,0,1,"


def test_as_tsv_format():
    """as_tsv_format test"""
    node_event = NodeEvent(0, 1, 3, 1, 0, [], 5)
    assert node_event.as_tsv_format == "0\t5\t4\t1\t3\t0\t1\t"


def test_add():
    """__add__ test"""
    node_event1 = NodeEvent(0, 1, 3, 1, 0, [Transfer(10, 20)], 5)
    node_event2 = NodeEvent(0, 0, 2, 1, 1, [Transfer(30, 40)], 5)

    assert node_event1 + node_event2 == NodeEvent(
        0, 1, 5, 2, 1, [Transfer(10, 20), Transfer(30, 40)], 10
    )
