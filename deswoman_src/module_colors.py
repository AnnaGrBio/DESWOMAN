BLACK = '\033[30m'
RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[93m' # orange on some systems
BLUE = '\033[34m'
MAGENTA = '\033[35m'
CYAN = '\033[36m'
LIGHT_GRAY = '\033[37m'
DARK_GRAY = '\033[90m'
BRIGHT_RED = '\033[91m'
BRIGHT_GREEN = '\033[92m'
BRIGHT_YELLOW = '\033[93m'
BRIGHT_BLUE = '\033[94m'
BRIGHT_MAGENTA = '\033[95m'
BRIGHT_CYAN = '\033[96m'
WHITE = '\033[97m'
ORANGE = '\033[33m'

RESET = '\033[0m' # called to return to standard terminal text color


def openFile(NameFile):
    """
    Opens a file, reads its contents, and returns the lines as a list.

    This function opens a file specified by the `NameFile` parameter in read mode, 
    reads all lines from the file, and returns them as a list of strings where each string 
    represents a line from the file.

    Args:
        NameFile (str): The path to the file that should be opened.

    Returns:
        list: A list of strings, each string corresponding to a line in the file.

    Raises:
        FileNotFoundError: If the specified file does not exist or cannot be found.
        IOError: If there are issues with reading the file (e.g., permission denied).
    """
    F = open(NameFile, "r")
    L = F.readlines()
    return L
