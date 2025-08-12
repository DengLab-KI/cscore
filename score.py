"""
Backwards-compatible shim for the previous script interface.

Use the installed CLI instead:
    cscore -i testdata -a comp1.tsv -b comp2.tsv -o out.tsv
"""

from cscore.cli import main

if __name__ == "__main__":
    main()
