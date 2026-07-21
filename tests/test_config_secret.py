"""
Tests for SECRET_KEY resolution in config.

The deployed service supplies the key as a Docker swarm secret via
SECRET_KEY_FILE, so it never appears in the service environment. Config
resolves it at class-body import time, so each case reloads the module.
"""

import importlib
import warnings

import pytest


def _load_config(monkeypatch, *, secret_file=None, secret=None):
    """Reload config with a controlled environment and return its Config."""
    monkeypatch.delenv('SECRET_KEY_FILE', raising=False)
    monkeypatch.delenv('SECRET_KEY', raising=False)
    if secret_file is not None:
        monkeypatch.setenv('SECRET_KEY_FILE', str(secret_file))
    if secret is not None:
        monkeypatch.setenv('SECRET_KEY', secret)

    import config as config_module
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        importlib.reload(config_module)
    return config_module.Config, [str(w.message) for w in caught]


@pytest.fixture(autouse=True)
def _restore_config():
    """Leave the imported config module in its normal state for other tests."""
    yield
    import config as config_module
    importlib.reload(config_module)


@pytest.mark.unit
class TestSecretKeyResolution:

    def test_reads_key_from_secret_file(self, monkeypatch, tmp_path):
        f = tmp_path / 'flask_secret'
        f.write_text('s3cret-from-swarm\n')
        cfg, warns = _load_config(monkeypatch, secret_file=f)
        # Trailing newline stripped — `docker secret create` from a file or a
        # here-string commonly leaves one, and Flask would treat it as part of
        # the key, silently changing it between creation methods.
        assert cfg.SECRET_KEY == 's3cret-from-swarm'
        assert not any('not set' in w for w in warns)

    def test_secret_file_wins_over_env(self, monkeypatch, tmp_path):
        """The file is the deployed path; env is the dev fallback."""
        f = tmp_path / 'flask_secret'
        f.write_text('from-file')
        cfg, _ = _load_config(monkeypatch, secret_file=f, secret='from-env')
        assert cfg.SECRET_KEY == 'from-file'

    def test_falls_back_to_env_var(self, monkeypatch):
        cfg, warns = _load_config(monkeypatch, secret='from-env')
        assert cfg.SECRET_KEY == 'from-env'
        assert not any('not set' in w for w in warns)

    def test_missing_file_warns_and_falls_back(self, monkeypatch, tmp_path):
        """An unreadable secret must not take the app down."""
        cfg, warns = _load_config(
            monkeypatch, secret_file=tmp_path / 'nope', secret='from-env')
        assert cfg.SECRET_KEY == 'from-env'
        assert any('could not be read' in w for w in warns)

    def test_empty_file_warns_and_falls_back(self, monkeypatch, tmp_path):
        f = tmp_path / 'flask_secret'
        f.write_text('   \n')
        cfg, warns = _load_config(monkeypatch, secret_file=f, secret='from-env')
        assert cfg.SECRET_KEY == 'from-env'
        assert any('is empty' in w for w in warns)

    def test_no_source_generates_random_key_and_warns(self, monkeypatch):
        cfg, warns = _load_config(monkeypatch)
        assert cfg.SECRET_KEY
        assert len(cfg.SECRET_KEY) >= 32
        assert any('not set' in w for w in warns)
