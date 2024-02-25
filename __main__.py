# __main__.py

from CompareProteinSeq import *

def main():
    compare_set=("P02189", "P68082", "P02144")

    # new a CompareProteinSeq and compare the protein sequence data from streaming
    c = CompareProteinSeq(compare_set)
    c.main()

if __name__ == "__main__":
    main()