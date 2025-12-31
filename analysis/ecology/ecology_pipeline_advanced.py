#!/usr/bin/env python3
"""
Wrapper: Advanced ecology pipeline (authoritative).

Delegates to p4_meta_postmerge.py which contains the full reviewer-revision
feature set.
"""
from __future__ import annotations

import subprocess
import sys


def main():
    cmd = [sys.executable or "python3", "p4_meta_postmerge.py"] + sys.argv[1:]
    raise SystemExit(subprocess.call(cmd))


if __name__ == "__main__":
    main()
