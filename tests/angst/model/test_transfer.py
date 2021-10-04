from fastdtlmapper.angst.model import Transfer


def test_as_text_with_direction():
    """as_text with direction test"""
    transfer = Transfer(1, 2)
    assert transfer.as_text == "1 -> 2"


def test_as_text_with_no_direction():
    """as_text with no direction test"""
    transfer = Transfer(1, 2, False)
    assert transfer.as_text == "1 <-> 2"
