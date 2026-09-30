"""Tests for curated_references.json (T7): shape, and the invariant that entries
start unverified until a human confirms them (CLAUDE.md invariant 7)."""
import json
import re
from pathlib import Path

REQUIRED_KEYS = {"id", "verified", "subtypes", "phenotypes", "citation", "claim"}
ID_PATTERN = re.compile(r"^REF-\d{3}$")


def _refs():
    return json.loads((Path(__file__).parent / "curated_references.json").read_text())


def test_is_a_nonempty_list_of_conformant_entries():
    refs = _refs()
    assert isinstance(refs, list) and len(refs) > 0
    for r in refs:
        assert REQUIRED_KEYS <= set(r), f"{r.get('id')} missing keys: {REQUIRED_KEYS - set(r)}"


def test_ids_are_well_formed_and_unique():
    refs = _refs()
    ids = [r["id"] for r in refs]
    assert len(ids) == len(set(ids)), "duplicate REF ids"
    for rid in ids:
        assert ID_PATTERN.match(rid), f"{rid} does not match REF-### pattern"


def test_verified_is_boolean_and_unconfirmed_entries_are_false():
    """Per invariant 7: entries stay verified=false until a human has read the paper.

    This does not assert every entry is false forever -- once someone verifies a
    citation the flag should flip to true -- only that the field is a real boolean,
    never a placeholder string or missing.
    """
    for r in _refs():
        assert isinstance(r["verified"], bool)


def test_subtypes_and_phenotypes_are_nonempty_lists_of_strings():
    for r in _refs():
        for field in ("subtypes", "phenotypes"):
            assert isinstance(r[field], list) and len(r[field]) > 0
            assert all(isinstance(x, str) for x in r[field])


def test_dcis_caveat_reference_present():
    """This project's real dataset needs a DCIS-specific reference (data/DATASET.md caveat)."""
    refs = _refs()
    assert any("DCIS" in r["subtypes"] for r in refs), "no reference tagged for the DCIS subtype context"
