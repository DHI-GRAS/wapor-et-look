import sys
from pyWAPOR.Collect.Landsat.PreprocessLandsat import PreprocessLandsat


def main(data_dir, output_dir):

    PreprocessLandsat(data_dir, output_dir)

if __name__ == '__main__':
    main(sys.argv)
