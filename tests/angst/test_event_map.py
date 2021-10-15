from pathlib import Path

import pytest
from fastdtlmapper.angst import AngstEventMap, NodeEvent, Transfer


def test_event_map(species_tree_nodeid_file: Path, angst_result_dir: Path):
    """test event_map"""
    event_map = AngstEventMap(species_tree_nodeid_file, angst_result_dir)
    nodeid2node_event = event_map.nodeid2node_event
    assert nodeid2node_event["N001"] == NodeEvent("N001", 1, 0, 0, 0, [], 1)
    assert nodeid2node_event["N002"] == NodeEvent("N002", 0, 0, 0, 0, [], 1)
    assert nodeid2node_event["N003"] == NodeEvent("N003", 0, 0, 0, 0, [], 1)
    assert nodeid2node_event["N004"] == NodeEvent("N004", 0, 0, 0, 0, [], 1)
    assert nodeid2node_event["N005"] == NodeEvent("N005", 0, 0, 0, 0, [], 1)
    assert nodeid2node_event["N006"] == NodeEvent("N006", 0, 0, 0, 0, [], 1)
    assert nodeid2node_event["outsp1"] == NodeEvent("outsp1", 0, 0, 0, 0, [], 1)
    assert nodeid2node_event["outsp2"] == NodeEvent("outsp2", 0, 0, 0, 0, [], 1)
    assert nodeid2node_event["insp1"] == NodeEvent("insp1", 0, 0, 0, 0, [], 1)
    assert nodeid2node_event["insp2"] == NodeEvent("insp2", 0, 0, 0, 0, [], 1)
    assert nodeid2node_event["insp3"] == NodeEvent(
        "insp3", 0, 1, 0, 1, [Transfer("N002", "insp3")], 3
    )
    assert nodeid2node_event["insp4"] == NodeEvent("insp4", 0, 0, 0, 0, [], 1)
    assert nodeid2node_event["insp5"] == NodeEvent("insp5", 0, 0, 1, 0, [], 0)


def test_event_map_write_tree_gain_loss(
    species_tree_nodeid_file: Path, angst_result_dir: Path, tmp_path: Path
):
    """test event_map write_tree(gain-loss)"""
    event_map = AngstEventMap(species_tree_nodeid_file, angst_result_dir)
    map_tree_outfile = tmp_path / "gain_loss_map.nwk"
    event_map.write_tree(map_tree_outfile, "gain-loss")
    assert map_tree_outfile.is_file()


def test_event_map_write_tree_dtl(
    species_tree_nodeid_file: Path, angst_result_dir: Path, tmp_path: Path
):
    """test event_map write_tree(gain-loss)"""
    event_map = AngstEventMap(species_tree_nodeid_file, angst_result_dir)
    map_tree_outfile = tmp_path / "gain_loss_map.nwk"
    event_map.write_tree(map_tree_outfile, "dtl")
    assert map_tree_outfile.is_file()


def test_event_map_write_tree_ng(
    species_tree_nodeid_file: Path, angst_result_dir: Path, tmp_path: Path
):
    """test event_map write_tree NG case"""
    event_map = AngstEventMap(species_tree_nodeid_file, angst_result_dir)
    map_tree_outfile = tmp_path / "gain_loss_map.nwk"
    with pytest.raises(ValueError):
        event_map.write_tree(map_tree_outfile, "NG")
