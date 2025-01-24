# Global log file
from datetime import datetime


LOG_FILE = f"../logs/fetch_log_{datetime.now().strftime('%Y%m%d%H%M%S')}.log"

def log(message, level="INFO", to_file=True):
    """
    Log a message with a timestamp.

    Parameters:
        message (str): The message to log.
        level (str): Log level (INFO, WARNING, ERROR).
        to_file (bool): Whether to write the log to a file.
    """
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    formatted_message = f"[{timestamp}] [{level}] {message}"
    print(formatted_message)
    if to_file:
        with open(LOG_FILE, "a") as f:
            f.write(formatted_message + "\n")