"""Why there is almost no test code?

Due to the installation size of the 'interproscan' dependency program,
it is difficult to run the FastDTLgoea integration tests in CI environment.
Therefore, FastDTLgoea integration test code is not written.
"""
from pathlib import Path

import pytest
from fastdtlmapper.scripts.FastDTLgoea import run_interproscan
from pytest import MonkeyPatch


def test_run_interproscan_ng(monkeypatch: MonkeyPatch):
    """test run_interproscan NG case"""
    # Set PATH blank (Set no interproscan PATH env)
    monkeypatch.setenv("PATH", "")
    with pytest.raises(SystemExit) as e:
        run_interproscan(Path("dummy"), Path("dummy"), Path("dummy"))
    assert e.type == SystemExit
    assert e.value.code == 1
