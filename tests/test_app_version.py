"""The deployed service must be able to say which release it is.

Ported from the sibling molAOP Builder, where `/health` carried its own
`"version"` literal and sat three releases behind for months because nothing
connected the literal to any release. The analyser started from a worse place:
it had no version string anywhere and no `/health` route at all, so a running
instance could not be identified except by comparing image digests — and an
archived copy under a DOI could not be matched to the results it produced.

`config.__version__` is the single source, pinned below to the newest released
heading in CHANGELOG.md, which is where a release is actually declared. A release
that forgets to bump it fails here.

`build` is separate and deliberately so: a semantic version cannot distinguish
two deployments of the same release, which is what "is the latest code live?"
asks. The Docker build bakes the commit SHA in as MOLAOP_IMAGE_SHA.

The citation metadata is pinned to the same constant, because `CITATION.cff` and
`.zenodo.json` are read by Zenodo from the tagged commit and a DOI cannot be
re-minted to correct a version that disagreed with the code.
"""
import json
import os
import re

import pytest

from config import __version__, get_build_ref

HERE = os.path.dirname(__file__)
ROOT = os.path.join(HERE, "..")
CHANGELOG = os.path.join(ROOT, "CHANGELOG.md")

# "## [5.0.0] - 2026-08-12" — deliberately excludes "## [Unreleased]".
# Three components required: the repo's historical headings are two-component
# ("## [4.0]") and are matched by the looser pattern below only for the tag/
# heading correspondence check, not for the version the code must equal.
RELEASE_HEADING = re.compile(r"^##\s*\[(\d+\.\d+\.\d+)\]", re.M)
ANY_RELEASE_HEADING = re.compile(r"^##\s*\[(\d+(?:\.\d+)+)\]", re.M)


def _changelog_text():
    with open(CHANGELOG, "r", encoding="utf-8") as fh:
        return fh.read()


def _newest_released_version():
    versions = RELEASE_HEADING.findall(_changelog_text())
    if not versions:
        pytest.skip(
            "no three-component release heading in CHANGELOG.md yet — the first "
            "such heading is written by the release-prep commit"
        )
    return versions[0]


def test_version_matches_newest_changelog_release():
    assert __version__ == _newest_released_version(), (
        f"config.__version__ is {__version__} but the newest release in "
        f"CHANGELOG.md is {_newest_released_version()}. Bump __version__ when "
        "cutting a release — /health reports it, and the Zenodo record is minted "
        "from the tagged commit."
    )


def test_health_reports_the_canonical_version_not_a_literal():
    """Guard against a literal being reintroduced into the route."""
    import inspect

    import app as app_module

    source = inspect.getsource(app_module)
    assert '"version": __version__' in source, (
        "/health must report config.__version__, not a hardcoded string"
    )


def test_build_ref_is_unknown_when_not_built_by_ci(monkeypatch):
    """No SHA in the environment must report "unknown", never a stand-in.

    A fabricated or blank build ref is worse than an honest "unknown": it makes
    an unidentifiable deployment look identified.
    """
    monkeypatch.delenv("MOLAOP_IMAGE_SHA", raising=False)
    assert get_build_ref() == "unknown"

    monkeypatch.setenv("MOLAOP_IMAGE_SHA", "")
    assert get_build_ref() == "unknown"


def test_build_ref_reports_the_injected_sha(monkeypatch):
    monkeypatch.setenv("MOLAOP_IMAGE_SHA", "44462be0000000000000000000000000000000aa")
    assert get_build_ref() == "44462be0000000000000000000000000000000aa"


def test_health_payload_carries_version_and_build(flask_client):
    response = flask_client.get("/health")
    assert response.status_code == 200
    body = response.get_json()

    assert body["version"] == __version__
    assert "build" in body, (
        "/health must identify the running build, not only the release version"
    )
    assert body["database"] is True


def test_health_reports_503_when_the_database_is_unreachable(flask_client, monkeypatch):
    """The status code is what a load balancer reads, so it has to move."""
    import app as app_module

    def _boom():
        raise RuntimeError("database is gone")

    monkeypatch.setattr(app_module.db_manager, "get_session", _boom)

    response = flask_client.get("/health")

    assert response.status_code == 503
    body = response.get_json()
    assert body["status"] == "unhealthy"
    assert body["database"] is False
    # The version must still be reported: identifying an unhealthy deployment is
    # more useful than identifying a healthy one.
    assert body["version"] == __version__


def test_dockerfile_and_ci_pass_the_commit_through():
    """The build arg is useless unless CI actually supplies it."""
    with open(os.path.join(ROOT, "Dockerfile"), encoding="utf-8") as fh:
        dockerfile = fh.read()
    assert "ARG GIT_SHA" in dockerfile
    assert "ENV MOLAOP_IMAGE_SHA=$GIT_SHA" in dockerfile

    workflow_path = os.path.join(ROOT, ".github", "workflows", "docker-publish.yml")
    with open(workflow_path, encoding="utf-8") as fh:
        workflow = fh.read()
    assert "GIT_SHA=${{ github.sha }}" in workflow, (
        "the image build step must pass the commit SHA as a build arg, or "
        "/health reports 'unknown' on every deployed image"
    )


def test_image_tag_matches_the_git_tag_verbatim():
    """`type=ref,event=tag` — without it no image tag equals the archived tag.

    The semver rules strip the leading v, so a v5.0.0 release produced 5.0.0 and
    5.0 but nothing a Zenodo record could be matched against by string equality.
    """
    workflow_path = os.path.join(ROOT, ".github", "workflows", "docker-publish.yml")
    with open(workflow_path, encoding="utf-8") as fh:
        workflow = fh.read()
    assert "type=ref,event=tag" in workflow


def test_citation_metadata_agrees_with_the_code_version():
    """CITATION.cff and .zenodo.json are frozen by the tag; they cannot drift.

    Skipped until the files exist, so this lands before them without failing.
    """
    cff_path = os.path.join(ROOT, "CITATION.cff")
    zenodo_path = os.path.join(ROOT, ".zenodo.json")
    if not os.path.exists(cff_path):
        pytest.skip("CITATION.cff not present yet")

    with open(cff_path, encoding="utf-8") as fh:
        cff = fh.read()
    match = re.search(r"^version:\s*['\"]?([0-9][^'\"\s]*)", cff, re.M)
    assert match, "CITATION.cff has no version field"
    assert match.group(1) == __version__, (
        f"CITATION.cff says {match.group(1)}, code says {__version__}"
    )

    if os.path.exists(zenodo_path):
        with open(zenodo_path, encoding="utf-8") as fh:
            payload = json.load(fh)  # must parse; Zenodo fails silently otherwise
        if "version" in payload:
            assert payload["version"] == __version__, (
                f".zenodo.json says {payload['version']}, code says {__version__}"
            )
        assert payload["license"] == "gpl-2.0-only"
        assert payload["upload_type"] == "software"
        assert payload["creators"], "a record with no creators lists the GitHub account"

        # The VHP4Safety grant, matching the project's other Zenodo deposits.
        # Note the identifier is NOT the NWA-ORC award number — searching Zenodo
        # for "1292.19.272" returns zero hits, while OpenAIRE indexes this
        # project under NWO grant code 36952. An unresolvable grant id makes
        # archiving fail silently, so the readable number stays in the
        # description and the resolvable code goes in `grants`.
        grant_ids = [g.get("id") for g in payload.get("grants", [])]
        assert "10.13039/501100003246::36952" in grant_ids, (
            "the VHP4Safety NWO grant must be declared; use funder DOI "
            "10.13039/501100003246 (NWO) with code 36952, not the NWA-ORC number"
        )
        assert {c.get("identifier") for c in payload.get("communities", [])} >= {"vhp4safety"}
