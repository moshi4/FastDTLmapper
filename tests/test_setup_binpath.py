from pathlib import Path

import pytest
from fastdtlmapper.setup_binpath import SetupBinpath
from pytest import MonkeyPatch


def test_setup_binpath_ok(monkeypatch: MonkeyPatch):
    """test setup_binpath ok (no error)"""
    monkeypatch.setenv("PATH", "/usr/bin/")
    root_bin_path = Path(__file__).parent.parent / "src" / "fastdtlmapper" / "bin"
    SetupBinpath(root_bin_path)


def test_setup_binpath_ng(monkeypatch: MonkeyPatch):
    """test setup_binpath ng (error occur)"""
    monkeypatch.setenv("PATH", "")
    root_bin_path = Path("testdir")
    with pytest.raises(SystemExit) as e:
        SetupBinpath(root_bin_path)
    assert e.type == SystemExit
    assert e.value.code == 1
