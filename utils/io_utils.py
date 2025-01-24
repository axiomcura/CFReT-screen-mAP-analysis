"""

This module will contain functions pertaining to loading and writing files.
"""
import pathlib

import yaml


def load_config(fpath: str | pathlib.Path ) -> dict:
    """Loads in configuration file if specified path

    Parameters
    ----------
    fpath : str | pathlib.Path
        path pointing to config file

    Returns
    -------
    dict
        contents of the configuration file
    """

    # type checking
    if not isinstance(fpath, (str | pathlib.Path)):
        raise TypeError("'fpath' must be a string or pathlib.Path object")
    if isinstance(fpath, str):
        fpath = pathlib.Path(fpath)

    # check if the path exists and returns the full path
    fpath = fpath.resolve(strict=True)

    # next is to load the yaml file
    with open(fpath) as content:
        return yaml.safe_load(content)
