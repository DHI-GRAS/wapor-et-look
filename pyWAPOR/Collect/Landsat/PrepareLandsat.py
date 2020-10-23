import sys
from pathlib import Path
from pyWAPOR.Collect.Landsat.PreprocessLandsat import PreprocessLandsat


def main(data_dir):

    PreprocessLandsat(data_dir)

if __name__ == '__main__':
    main(sys.argv)
