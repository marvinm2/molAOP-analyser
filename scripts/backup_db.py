#!/usr/bin/env python3
"""Back up the analyser's SQLite database using the Online Backup API.

Scheduled daily as a Swarm cron job (`cronjobs_molaop-analyser-backup`), which
runs this script in a one-shot container off this same image with
/mnt/gluster/docker/molaop-analyser/data bind-mounted at /data. It reads the
database file directly on GlusterFS rather than exec'ing into the running
container, so it works on whichever node the scheduler picks and regardless of
where molaop-analyser itself is running.

Run it by hand before a schema migration:
    ssh tgx1 "docker exec \\$(docker ps -qf name=molaop-analyser) \\
        python /app/scripts/backup_db.py"

sqlite3's backup API (https://sqlite.org/backup.html) is safe while the app is
writing; copying the .db file alone is not, because it leaves the -wal and -shm
sidecars behind. Every check below removes the output on failure, so a run that
dies never leaves a file that sorts as the newest backup.
"""
import os
import sqlite3
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

DB_PATH = Path(os.environ.get("DB_PATH", "/data/molAOP_analyser.db"))
BACKUP_DIR = Path(os.environ.get("BACKUP_DIR", "/data/backups"))
RETENTION_DAYS = int(os.environ.get("RETENTION_DAYS", "7"))
# How long to wait for the live writer to release its lock before giving up.
BUSY_TIMEOUT_S = float(os.environ.get("BUSY_TIMEOUT_MS", "60000")) / 1000
# A real backup of this database is >100 MB. Anything under a few pages is a
# failed one, whatever PRAGMA integrity_check says about it — SQLite treats a
# zero-length file as a valid, empty database, so integrity_check returns "ok"
# on precisely the failure it would need to catch. Check the size first.
MIN_BACKUP_BYTES = int(os.environ.get("MIN_BACKUP_BYTES", "65536"))
# Tables that must hold rows for the file to be a backup *of this database*.
# shared_results is deliberately not here: links expire after 30 days, so an
# empty one is a normal state rather than evidence of a truncated copy.
REQUIRED_NONEMPTY = ("experiments", "batches")


def fail(backup_file, message):
    backup_file.unlink(missing_ok=True)
    print(f"[BACKUP ERROR] {message}", file=sys.stderr)
    sys.exit(1)


def main():
    if not DB_PATH.exists():
        print(f"[BACKUP ERROR] no database at {DB_PATH}", file=sys.stderr)
        sys.exit(1)

    BACKUP_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    backup_file = BACKUP_DIR / f"molAOP_analyser_{timestamp}.db"

    source = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True,
                             timeout=BUSY_TIMEOUT_S)
    try:
        # Take the read lock ourselves, under the connection's busy timeout,
        # before handing the connection to backup(). Connection.backup() retries
        # a locked source *forever* — it does not honour the timeout — so
        # calling it first would turn a stuck writer into a job that never
        # returns, and with skip-running=true a job that never returns silently
        # suppresses every run after it. Failing loudly is the safer end state.
        source.execute("BEGIN")
        source.execute("SELECT COUNT(*) FROM sqlite_master")
    except sqlite3.Error as exc:
        source.close()
        fail(backup_file, f"could not read {DB_PATH} within "
                          f"{BUSY_TIMEOUT_S:g}s: {exc}")

    try:
        target = sqlite3.connect(backup_file)
        with target:
            # pages=-1 copies in a single step, and the read transaction above
            # is still open, so a concurrent write can neither block nor
            # restart the copy partway through.
            source.backup(target, pages=-1)
        target.close()
        source.rollback()
    except sqlite3.Error as exc:
        fail(backup_file, f"sqlite backup failed for {DB_PATH}: {exc}")
    finally:
        source.close()

    backup_bytes = backup_file.stat().st_size
    if backup_bytes < MIN_BACKUP_BYTES:
        fail(backup_file, f"{backup_file} is {backup_bytes} bytes, "
                          f"under the {MIN_BACKUP_BYTES}-byte floor")

    check = sqlite3.connect(backup_file)
    try:
        result = check.execute("PRAGMA integrity_check;").fetchone()[0]
        if result != "ok":
            fail(backup_file, f"integrity check failed for {backup_file}: {result}")

        counts = {}
        for table in REQUIRED_NONEMPTY:
            try:
                counts[table] = check.execute(
                    f'SELECT COUNT(*) FROM "{table}"').fetchone()[0]
            except sqlite3.Error as exc:
                fail(backup_file, f"{backup_file} has no usable {table} table: {exc}")
            if counts[table] < 1:
                fail(backup_file, f"{backup_file} holds no rows in {table}")
    finally:
        check.close()

    summary = ", ".join(f"{n} {t}" for t, n in counts.items())
    print(f"[BACKUP OK] {backup_file} — {backup_bytes} bytes, "
          f"integrity ok, {summary}")

    # Prune by age. The glob covers "-wal"/"-shm" sidecars too: a clean close
    # removes them, but one left behind by an interrupted read would otherwise
    # match no prune pattern and stay forever.
    cutoff = datetime.now() - timedelta(days=RETENTION_DAYS)
    for old in BACKUP_DIR.glob("molAOP_analyser_*.db*"):
        if datetime.fromtimestamp(old.stat().st_mtime) < cutoff:
            old.unlink(missing_ok=True)
            print(f"[BACKUP PRUNE] removed {old.name} (older than "
                  f"{RETENTION_DAYS} days)")
    # Prune empty backups regardless of age — one left by a run predating the
    # checks above would otherwise sit there for a full retention period,
    # sorting as the newest.
    for empty in BACKUP_DIR.glob("molAOP_analyser_*.db"):
        if empty.stat().st_size == 0:
            empty.unlink(missing_ok=True)
            print(f"[BACKUP PRUNE] removed empty {empty.name}")


if __name__ == "__main__":
    main()
