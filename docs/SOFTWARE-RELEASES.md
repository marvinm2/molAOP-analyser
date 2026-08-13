# Releasing the software to Zenodo

Runbook for tagging this repository and archiving the source under a citable DOI.

> **The analyser releases after the Builder.** It consumes the Builder's API and its
> `.zenodo.json` declares a `requires` relation to the Builder's software DOI, so that DOI has
> to exist first. The Builder also has a second, separate Zenodo deposit for its curated
> **dataset** — see `docs/RELEASES.md` in that repository. This repository has one deposit
> only: the software.

## How the archive is triggered

Zenodo's GitHub integration is **already enabled** on this repository. There is no script to
run and no token to supply.

- **Publishing a GitHub Release** fires the webhook. Zenodo downloads the source tarball of
  that tag and mints a DOI.
- **Pushing a tag does not.** `git push origin v5.0.0` archives nothing on its own.
- **Pre-releases are archived too**, and mint a real DOI. Never use one as a rehearsal.

## What Zenodo reads

`.zenodo.json` from the **tagged commit** — which is why the release-prep commit lands before
the tag. `CITATION.cff` is present for GitHub's "Cite this repository" button but Zenodo
ignores it whenever `.zenodo.json` exists. Keep their `version` fields in agreement; a test
asserts it.

A minted version DOI cannot be withdrawn, only superseded.

## Checklist

1. Everything that must be true in the archive is on `main`: citation metadata, licence and
   copyright, `data/README.md` provenance, a README a stranger can follow — including the
   system libraries WeasyPrint needs, which `pip install` does not supply.
2. One release-prep commit: `CHANGELOG.md` `[Unreleased]` → the version heading, the version
   constant, and `CITATION.cff` `version` + `date-released`.
3. Merge; wait for CI green **on that exact commit**.
4. Annotated tag, push, confirm `ghcr.io/marvinm2/molaop-analyser:<tag>` exists.
5. Publish a **full** GitHub Release with hand-written notes.
6. Verify the Zenodo record: creators (no `dependabot[bot]`), licence, version, files,
   community, the VHP4Safety grant, and that the `requires` relation points at the Builder's
   real software DOI.
7. Backfill this repository's own DOI into the README badge, `CITATION.cff` **and**
   `.zenodo.json` — the last is re-read from every future tag, so skipping it silently
   reverts the record's cross-links at the next release.

## If the record does not appear

Check Zenodo's **Errors** tab and the webhook delivery log. An invalid `.zenodo.json` fails
silently. Recovery is to delete the GitHub Release, fix, and publish again. Rehearse changes
to the metadata on **sandbox.zenodo.org** first — well-formed JSON is not evidence that Zenodo
accepts the licence identifier, community, or relation types.

## A note on the test suite

The suite runs with outbound network access refused — a `conftest.py` guard raises on any
`AF_INET`/`AF_INET6` connect, so a release cannot go red because the Builder's public API
(100 requests/hour/IP) or the AOP-Wiki SPARQL endpoint was busy. Nothing in the suite needs
the network; the guard made it about 25% faster by removing calls that were being attempted
and silently falling back.

A test whose subject genuinely is the network call opts out with
`@pytest.mark.allow_network`. If you add one, remember it reintroduces exactly the dependency
the guard removes, so the release becomes sensitive to a third party again.

Independently, `pytest.ini` sets a 300-second per-test timeout with `timeout_method = thread`
(issue #104), which dumps every thread's stack. A stalled test therefore fails with a stack
trace naming the culprit rather than hanging.
