"""
Tests for SECRET_KEY resolution in config.

The deployed service supplies the key as a Docker swarm secret via
SECRET_KEY_FILE, so it never appears in the service environment.

These cases call ``config.resolve_secret_key()`` directly. They used to
``importlib.reload(config)`` instead, which left every module that had done
``from config import Config`` holding a stale class — so unrelated fixtures
patching ``config.Config`` silently patched an object nothing read, and the
damage showed up as order-dependent failures elsewhere in the suite (#73).
Nothing here reloads a module any more; keep it that way.
"""

import warnings

import pytest

from config import resolve_secret_key


def _resolve(monkeypatch, *, secret_file=None, secret=None):
    """Resolve the key under a controlled environment; return (key, warnings)."""
    monkeypatch.delenv('SECRET_KEY_FILE', raising=False)
    monkeypatch.delenv('SECRET_KEY', raising=False)
    if secret_file is not None:
        monkeypatch.setenv('SECRET_KEY_FILE', str(secret_file))
    if secret is not None:
        monkeypatch.setenv('SECRET_KEY', secret)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        key = resolve_secret_key()
    return key, [str(w.message) for w in caught]


@pytest.mark.unit
class TestSecretKeyResolution:

    def test_reads_key_from_secret_file(self, monkeypatch, tmp_path):
        f = tmp_path / 'flask_secret'
        f.write_text('s3cret-from-swarm\n')
        key, warns = _resolve(monkeypatch, secret_file=f)
        # Trailing newline stripped — `docker secret create` from a file or a
        # here-string commonly leaves one, and Flask would treat it as part of
        # the key, silently changing it between creation methods.
        assert key == 's3cret-from-swarm'
        assert not any('not set' in w for w in warns)

    def test_secret_file_wins_over_env(self, monkeypatch, tmp_path):
        """The file is the deployed path; env is the dev fallback."""
        f = tmp_path / 'flask_secret'
        f.write_text('from-file')
        key, _ = _resolve(monkeypatch, secret_file=f, secret='from-env')
        assert key == 'from-file'

    def test_falls_back_to_env_var(self, monkeypatch):
        key, warns = _resolve(monkeypatch, secret='from-env')
        assert key == 'from-env'
        assert not any('not set' in w for w in warns)

    def test_missing_file_warns_and_falls_back(self, monkeypatch, tmp_path):
        """An unreadable secret must not take the app down."""
        key, warns = _resolve(
            monkeypatch, secret_file=tmp_path / 'nope', secret='from-env')
        assert key == 'from-env'
        assert any('could not be read' in w for w in warns)

    def test_empty_file_warns_and_falls_back(self, monkeypatch, tmp_path):
        f = tmp_path / 'flask_secret'
        f.write_text('   \n')
        key, warns = _resolve(monkeypatch, secret_file=f, secret='from-env')
        assert key == 'from-env'
        assert any('is empty' in w for w in warns)

    def test_no_source_generates_random_key_and_warns(self, monkeypatch):
        key, warns = _resolve(monkeypatch)
        assert key
        assert len(key) >= 32
        assert any('not set' in w for w in warns)


def test_no_test_module_reloads_config():
    """Guard for #73: module reloads must not come back.

    importlib.reload(config) rebinds config.Config to a new class object while
    every module that did `from config import Config` keeps the old one. The
    result is not a failure here — it is an unrelated fixture elsewhere in the
    suite quietly patching an object nothing reads, surfacing later as a test
    that passes alone and fails in the full suite.

    Parses each test module rather than grepping it, so prose about the trap
    (this docstring included) does not trip the guard.
    """
    import ast
    import pathlib as _pathlib

    def reloads(path):
        for node in ast.walk(ast.parse(path.read_text())):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if isinstance(func, ast.Attribute) and func.attr == 'reload':
                return True
            if isinstance(func, ast.Name) and func.id == 'reload':
                return True
        return False

    tests_dir = _pathlib.Path(__file__).parent
    offenders = [p.name for p in sorted(tests_dir.glob('test_*.py')) if reloads(p)]
    assert offenders == [], (
        f"{offenders} reload a module. Call the resolver directly instead "
        f"(see resolve_secret_key), or the patch-the-wrong-Config trap returns."
    )
