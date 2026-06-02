#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Resolve the Sniffles binary used by VAULT."""

import logging
import os
import shutil
import subprocess


def _is_executable(path):
    return bool(path) and os.path.isfile(path) and os.access(path, os.X_OK)


def _run_probe(path, args):
    if not _is_executable(path):
        return None
    try:
        result = subprocess.run(
            [path] + list(args),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=10,
        )
    except Exception as e:
        return None, "probe failed: %s" % e

    output = (result.stdout or "") + "\n" + (result.stderr or "")
    return result.returncode, output


def _probe_output(path, args):
    result = _run_probe(path, args)
    if result is None:
        return None

    returncode, output = result
    if returncode != 0:
        return None
    return output


def _is_runnable(path, args):
    return _probe_output(path, args) is not None


def _candidate_identity(path):
    """Return (identity, detail) for a Sniffles candidate."""
    if not path:
        return None, "not provided"
    if not os.path.exists(path):
        return None, "path does not exist: %s" % path
    if not os.path.isfile(path):
        return None, "not a file: %s" % path
    if not os.access(path, os.X_OK):
        return None, "not executable: %s" % path

    outputs = []
    for args in (["--version"], ["-h"], ["--help"]):
        result = _run_probe(path, args)
        if result is None:
            continue
        _returncode, output = result
        if output:
            outputs.append(output)

    text = "\n".join(outputs)
    lowered = text.lower()

    if "sniffles2" in lowered or "version 2." in lowered:
        return "sniffles2", "appears to be Sniffles2"

    if (
        "version: 1." in lowered
        or "sniffles version 1." in lowered
        or "-m <sorted.bam>" in lowered
    ):
        return "sniffles1", "appears to be Sniffles1"

    if text:
        one_line = " ".join(text.split())[:160]
        return None, "runnable, but version could not be identified: %s" % one_line

    return None, "no probe output from: %s" % path


def _looks_like_requested_caller(path, caller):
    identity, _detail = _candidate_identity(path)
    return identity == caller


def _describe_candidate(label, path, caller, missing_message):
    if not path:
        return "%s: %s" % (label, missing_message)

    identity, detail = _candidate_identity(path)
    if identity == caller:
        return "%s: found %s (%s)" % (label, path, detail)

    if identity in ("sniffles1", "sniffles2"):
        return (
            "%s: found %s, but it is %s; requested --sv_caller %s"
            % (label, path, identity, caller)
        )

    return "%s: %s" % (label, detail)


def _format_failure(caller, details, env_name):
    selected = "--sv_caller %s" % caller
    if caller == "sniffles1":
        selected += " (default unless you changed it)"

    message = [
        "[sniffles] No usable %s binary found." % caller,
        "Currently selected SV caller: %s" % selected,
        "Search order:",
    ]
    message.extend("  - " + detail for detail in details)

    if caller == "sniffles1":
        message.append(
            "Hint: PATH may contain Sniffles2. For Sniffles1, pass "
            "--sniffles_path, or run: export %s=/path/to/sniffles1"
            % env_name
        )
        message.append(
            "If you intentionally want to use Sniffles2, set "
            "--sv_caller sniffles2. Note that Sniffles2 is supported, but "
            "Sniffles1 is usually a better match for VAULT's small "
            "UMI-grouped reads."
        )
    else:
        message.append(
            "Hint: pass --sniffles_path, or run: "
            "export %s=/path/to/sniffles2" % env_name
        )

    return "\n".join(message)


def resolve_sniffles(user_path=None, caller="sniffles1", logger=None):
    """Return the Sniffles executable path for the selected caller style.

    Search order for sniffles1:
      1) --sniffles_path
      2) $VAULT_SNIFFLES1
      3) sniffles in PATH

    Search order for sniffles2:
      1) --sniffles_path
      2) $VAULT_SNIFFLES2
      3) sniffles in PATH
    """
    log = logger or logging.getLogger(__name__)

    if user_path and _looks_like_requested_caller(user_path, caller):
        log.info("[sniffles] Use user-specified binary: %s", user_path)
        return user_path

    env_name = "VAULT_SNIFFLES1" if caller == "sniffles1" else "VAULT_SNIFFLES2"
    env_path = os.environ.get(env_name)
    if env_path and _looks_like_requested_caller(env_path, caller):
        log.info("[sniffles] Use $%s: %s", env_name, env_path)
        return env_path

    path_in_path = shutil.which("sniffles")
    if path_in_path and _looks_like_requested_caller(path_in_path, caller):
        log.info("[sniffles] Use sniffles from PATH: %s", path_in_path)
        return path_in_path

    details = []
    details.append(_describe_candidate(
        "--sniffles_path",
        user_path,
        caller,
        "not provided",
    ))
    details.append(_describe_candidate(
        "$" + env_name,
        env_path,
        caller,
        (
            "not set/exported in this Python process. If it is only a shell "
            "variable, run: export %s=/path/to/sniffles"
        ) % env_name,
    ))
    details.append(_describe_candidate(
        "sniffles in PATH",
        path_in_path,
        caller,
        "not found in PATH",
    ))

    raise FileNotFoundError(_format_failure(caller, details, env_name))
