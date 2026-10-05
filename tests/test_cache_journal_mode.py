"""The reference cache must not use SQLite WAL mode (#134).

In production CACHE_DIR is on GlusterFS. WAL keeps a shared-memory index in the
-shm file, which SQLite documents as unsafe on network filesystems.
"""
import os

import diskcache

from cache_manager import get_reference_cache

WAL = 2  # SQLite header byte 18: 1 = rollback journal, 2 = WAL


def _journal_bytes(directory):
    return {
        open(os.path.join(directory, f"{i:03d}", "cache.db"), "rb").read(20)[18]
        for i in range(8)
    }


def test_reference_cache_uses_a_rollback_journal(tmp_path):
    cache = get_reference_cache(str(tmp_path))
    cache.set("k", 1)
    cache.close()
    assert WAL not in _journal_bytes(tmp_path)
    assert not any((tmp_path / f"{i:03d}" / "cache.db-shm").exists() for i in range(8))


def test_shards_created_in_wal_mode_are_converted_and_keep_their_data(tmp_path):
    old = diskcache.FanoutCache(directory=str(tmp_path), shards=8)  # diskcache default: WAL
    old.set("k", "kept")
    old.close()
    assert _journal_bytes(tmp_path) == {WAL}

    cache = get_reference_cache(str(tmp_path))
    assert cache.get("k") == "kept"
    cache.close()
    assert WAL not in _journal_bytes(tmp_path)
