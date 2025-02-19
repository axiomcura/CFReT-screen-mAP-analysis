"""

This module will contain functions pertaining to loading and writing files.
"""
import json
import pathlib

import yaml


def load_config(fpath: str | pathlib.Path) -> dict:
    """ Loads a configuration file from the specified path.

    Parameters
    ----------
    fpath : str | pathlib.Path
        Path pointing to the configuration file. It can be a string or a pathlib.Path
        object.

    Returns
    -------
    dict
        Contents of the configuration file as a dictionary.

    Raises
    ------
    TypeError
        If 'fpath' is not a string or pathlib.Path object.
    ValueError
        If the file extension is not .json or .yaml.
    FileNotFoundError
        If the file does not exist.
    """

    # Type checking
    if not isinstance(fpath, (str, pathlib.Path)):
        raise TypeError("'fpath' must be a string or pathlib.Path object")
    if isinstance(fpath, str):
        fpath = pathlib.Path(fpath)

    # Resolve the path and check if it exists
    fpath = fpath.resolve(strict=True)

    # Load the configuration file based on its extension
    if fpath.suffix == ".json":
        with open(fpath) as content:
            return json.load(content)
    elif fpath.suffix == ".yaml":
        with open(fpath) as content:
            return yaml.safe_load(content)
    else:
        raise ValueError("Only .json and .yaml files are supported")
