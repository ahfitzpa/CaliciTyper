"""Unit tests for calicityper.config path resolution and environment overrides."""
from pathlib import Path
import pytest

from calicityper import config


def test_default_data_host_vs_container(monkeypatch):
    monkeypatch.delenv("APPTAINER_CONTAINER", raising=False)
    monkeypatch.delenv("SINGULARITY_CONTAINER", raising=False)
    assert config._default_data() == config.PKG / "data"

    monkeypatch.setenv("APPTAINER_CONTAINER", "/opt/calicityper.sif")
    assert config._default_data() == Path("/data")

    monkeypatch.delenv("APPTAINER_CONTAINER", raising=False)
    monkeypatch.setenv("SINGULARITY_CONTAINER", "/opt/calicityper.sif")
    assert config._default_data() == Path("/data")


def test_is_sapovirus():
    assert config.is_sapovirus("Caliciviridae_Sapovirus_GI")
    assert not config.is_sapovirus("GII")
    assert not config.is_sapovirus("")


def test_min_sites_for_region_override(monkeypatch):
    monkeypatch.setattr(config, "MIN_INFORMATIVE_SITES", 300)
    monkeypatch.setattr(config, "MIN_INFORMATIVE_SITES_ORF1", 170)
    monkeypatch.setattr(config, "MIN_INFORMATIVE_SITES_ORF2", None)
    assert config.min_sites_for("orf1") == 170
    assert config.min_sites_for("orf2") == 300


def test_exe_env_override_and_fallback(monkeypatch):
    monkeypatch.setenv("CALICITYPER_IQTREE", "/custom/iqtree3")
    assert config.exe("iqtree") == "/custom/iqtree3"

    monkeypatch.delenv("CALICITYPER_IQTREE", raising=False)
    monkeypatch.setattr(config.shutil, "which", lambda c: "/usr/bin/iqtree" if c == "iqtree" else None)
    assert config.exe("iqtree") == "/usr/bin/iqtree"

    monkeypatch.setattr(config.shutil, "which", lambda c: None)
    with pytest.raises(FileNotFoundError, match="CALICITYPER_IQTREE"):
        config.exe("iqtree")
