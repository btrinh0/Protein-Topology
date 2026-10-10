"""Update the Form A/B preregistration's pinned-file hash table.

Run after the final change to any pinned file. Only the table between the marker
comments is replaced. The command prints the preregistration SHA-256 for the
decision log.
"""
from __future__ import annotations

import re

import form_ab_lib as lib

START, END = "<!-- frozen-hashes:start -->", "<!-- frozen-hashes:end -->"


def table() -> str:
    rows = ["| File | SHA-256 |", "| --- | --- |"]
    rows += [f"| `{name}` | `{digest}` |" for name, digest in lib.frozen_file_hashes().items()]
    return "\n".join(rows)


def main() -> None:
    text = lib.PREREG.read_text(encoding="utf-8")
    if START not in text or END not in text:
        raise SystemExit("Marker comments not found in the preregistration")
    pattern = re.compile(re.escape(START) + r".*?" + re.escape(END), re.S)
    text = pattern.sub(lambda _: f"{START}\n{table()}\n{END}", text)
    lib.PREREG.write_text(text, encoding="utf-8", newline="\n")
    print(f"Preregistration SHA-256: {lib.P.sha256(lib.PREREG)}")
    print(f"Pinned files listed: {len(lib.FROZEN_FILES)}; drift: {lib.frozen_file_drift(text) or 'none'}")


if __name__ == "__main__":
    main()
