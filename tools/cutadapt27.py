#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Helper to resolve which cutadapt executable to use.
#
# Priority:
#   1) user-specified --cutadapt_path
#   2) $VAULT_CUTADAPT27 environment variable
#   3) "cutadapt" found in current PATH
#
# This module is imported by vault_main.py:
#     from tools import cutadapt27
#     cutadapt_bin = cutadapt27.resolve_cutadapt(getattr(args, "cutadapt_path", None))
#
# The resolved path is then exported to the environment as VAULT_CUTADAPT_BIN
# so that shell scripts (e.g. check_umi.sh) can reuse it.

import os
import shutil
import logging
import subprocess


def _is_executable(path):
    """Return True if path exists and is executable."""
    return bool(path) and os.path.isfile(path) and os.access(path, os.X_OK)


def _is_runnable_cutadapt(path):
    if not _is_executable(path):
        return False
    try:
        result = subprocess.run(
            [path, "--version"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
            timeout=10,
        )
        return result.returncode == 0
    except Exception:
        return False


def resolve_cutadapt(user_path=None, logger=None):
    """
    Decide which cutadapt binary to use.

    Search order:
      1) user_path (CLI option --cutadapt_path)
      2) $VAULT_CUTADAPT27
      3) "cutadapt" in PATH

    Parameters
    ----------
    user_path : str or None
        Explicit path passed from CLI. If not None, it must exist and be executable.
    logger : logging.Logger or None
        Optional logger. If None, use logging.getLogger(__name__).

    Returns
    -------
    str
        Path to the selected cutadapt executable.

    Raises
    ------
    FileNotFoundError
        If no usable cutadapt binary can be found.
    """
    log = logger or logging.getLogger(__name__)

    # 1. User-specified path via CLI
    if user_path:
        if _is_runnable_cutadapt(user_path):
            log.info("[cutadapt] Use user-specified cutadapt: %s", user_path)
            return user_path
        raise FileNotFoundError(
            "[cutadapt] User-specified cutadapt not runnable: %s" % user_path
        )

    # 2. Environment variable VAULT_CUTADAPT27
    env_path = os.environ.get("VAULT_CUTADAPT27")
    if _is_runnable_cutadapt(env_path):
        log.info("[cutadapt] Use cutadapt from $VAULT_CUTADAPT27: %s", env_path)
        return env_path

    # 3. Fallback: cutadapt in current PATH.
    path_in_path = shutil.which("cutadapt")
    if _is_runnable_cutadapt(path_in_path):
        log.info("[cutadapt] Use cutadapt from PATH: %s", path_in_path)
        return path_in_path

    # Nothing worked
    raise FileNotFoundError(
        "[cutadapt] No usable cutadapt executable found.\n"
        "Tried:\n"
        "  1) --cutadapt_path\n"
        "  2) $VAULT_CUTADAPT27\n"
        "  3) cutadapt in PATH\n"
        "Please either:\n"
        "  - install cutadapt in your VAULT conda environment/PATH, or\n"
        "  - provide --cutadapt_path /path/to/cutadapt, or\n"
        "  - set $VAULT_CUTADAPT27 to a valid cutadapt binary."
    )
