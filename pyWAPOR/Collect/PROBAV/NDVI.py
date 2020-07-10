import sys
from pyWAPOR.Collect.PROBAV.DataAccess import download_data


def main(download_dir, start_date, end_date, latitude_extent, longitude_extent, username, password):
    """

    """

    download_data(download_dir, start_date, end_date, latitude_extent, longitude_extent, username, password)

if __name__ == '__main__':
    main(sys.argv)