"""Tests for scripts/backup_db.py (issue #114).

Each guard in the backup script exists because the obvious version of the check
passes on the failure it is meant to catch. These tests assert the failing
behaviour, not just the happy path:

- PRAGMA integrity_check returns "ok" on a zero-length file, so it cannot tell a
  finished backup from an aborted one — hence the byte floor.
- A structurally valid database with no rows is not a backup of this database —
  hence the content check.
- Connection.backup() retries a locked source forever, so the wait for a
  concurrent writer has to be bounded outside it.
"""

import os
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "backup_db.py"


def make_db(path, experiments=3, batches=2):
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE experiments (id INTEGER PRIMARY KEY)")
    conn.execute("CREATE TABLE batches (id INTEGER PRIMARY KEY)")
    conn.executemany("INSERT INTO experiments (id) VALUES (?)",
                     [(i,) for i in range(experiments)])
    conn.executemany("INSERT INTO batches (id) VALUES (?)",
                     [(i,) for i in range(batches)])
    conn.commit()
    conn.close()
    return path


def run_backup(db_path, backup_dir, **env_overrides):
    env = {**os.environ, "DB_PATH": str(db_path), "BACKUP_DIR": str(backup_dir),
           "MIN_BACKUP_BYTES": "1024"}
    env.update({k: str(v) for k, v in env_overrides.items()})
    return subprocess.run([sys.executable, str(SCRIPT)], env=env,
                          capture_output=True, text=True, timeout=120)


def backups_in(backup_dir):
    return sorted(Path(backup_dir).glob("molAOP_analyser_*.db"))


def test_backup_is_restorable(tmp_path):
    db = make_db(tmp_path / "src.db", experiments=7)
    out = tmp_path / "backups"

    result = run_backup(db, out)

    assert result.returncode == 0, result.stderr
    copies = backups_in(out)
    assert len(copies) == 1
    restored = sqlite3.connect(copies[0])
    assert restored.execute("SELECT COUNT(*) FROM experiments").fetchone()[0] == 7
    restored.close()


def test_integrity_check_alone_would_accept_an_empty_file(tmp_path):
    """The premise of the byte floor, asserted rather than assumed."""
    empty = tmp_path / "zero.db"
    empty.touch()
    conn = sqlite3.connect(empty)
    assert conn.execute("PRAGMA integrity_check;").fetchone()[0] == "ok"
    conn.close()


def test_undersized_backup_is_rejected_and_removed(tmp_path):
    db = make_db(tmp_path / "src.db")
    out = tmp_path / "backups"

    result = run_backup(db, out, MIN_BACKUP_BYTES=50_000_000)

    assert result.returncode == 1
    assert "under the" in result.stderr
    assert backups_in(out) == []


def test_database_without_rows_is_rejected(tmp_path):
    db = make_db(tmp_path / "src.db", experiments=0)
    out = tmp_path / "backups"

    result = run_backup(db, out)

    assert result.returncode == 1
    assert "holds no rows in experiments" in result.stderr
    assert backups_in(out) == []


def test_missing_database_fails_loudly(tmp_path):
    result = run_backup(tmp_path / "absent.db", tmp_path / "backups")

    assert result.returncode == 1
    assert "no database" in result.stderr


@pytest.mark.parametrize("patience_ms,lock_seconds,expect_success", [
    (1000, 6, False),   # writer outlasts our patience — fail, do not hang
    (20000, 2, True),   # writer releases in time — wait it out, then back up
])
def test_bounded_wait_for_a_concurrent_writer(tmp_path, patience_ms,
                                              lock_seconds, expect_success):
    db = make_db(tmp_path / "src.db")
    out = tmp_path / "backups"

    holder = subprocess.Popen(
        [sys.executable, "-c",
         f"import sqlite3,time\n"
         f"c=sqlite3.connect({str(db)!r},isolation_level=None)\n"
         f"c.execute('BEGIN EXCLUSIVE')\n"
         f"time.sleep({lock_seconds})\n"
         f"c.execute('ROLLBACK')"])
    time.sleep(0.5)
    try:
        started = time.monotonic()
        result = run_backup(db, out, BUSY_TIMEOUT_MS=patience_ms)
        elapsed = time.monotonic() - started
    finally:
        holder.wait(timeout=30)

    if expect_success:
        assert result.returncode == 0, result.stderr
        assert len(backups_in(out)) == 1
    else:
        assert result.returncode == 1
        assert "database is locked" in result.stderr
        assert backups_in(out) == []
        # The point of the guard: it returns rather than retrying forever.
        assert elapsed < lock_seconds


def test_retention_prunes_aged_backups_and_stray_sidecars(tmp_path):
    db = make_db(tmp_path / "src.db")
    out = tmp_path / "backups"
    out.mkdir()
    aged = out / "molAOP_analyser_20260101_020000.db"
    sidecar = out / "molAOP_analyser_20260101_020000.db-wal"
    recent = out / "molAOP_analyser_20260803_020000.db"
    for path in (aged, sidecar, recent):
        path.write_bytes(b"x" * 2048)
    old = time.time() - 30 * 86400
    os.utime(aged, (old, old))
    os.utime(sidecar, (old, old))

    result = run_backup(db, out, RETENTION_DAYS=7)

    assert result.returncode == 0, result.stderr
    assert not aged.exists()
    assert not sidecar.exists()
    assert recent.exists()


def test_empty_backups_are_pruned_regardless_of_age(tmp_path):
    db = make_db(tmp_path / "src.db")
    out = tmp_path / "backups"
    out.mkdir()
    # What a pre-fix run left behind: 0 bytes, today's date, sorting newest.
    stale = out / "molAOP_analyser_20260804_010000.db"
    stale.touch()

    result = run_backup(db, out, RETENTION_DAYS=7)

    assert result.returncode == 0, result.stderr
    assert not stale.exists()
