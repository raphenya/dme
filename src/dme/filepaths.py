""" Imports. """
import sys
import os


def determine_path():
    """ Get path."""
    try:
        root = __file__
        if os.path.islink(root):
            root = os.path.realpath(root)
        return os.path.dirname(os.path.abspath(root))
    except ValueError:
        print(f"I'm sorry, but something is wrong. {ValueError}")
        print("There is no __file__ variable. Please contact the author.")
        sys.exit()
