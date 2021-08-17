class NotDirectory(Exception):
    def __init__(self):
        self.message = "The provided file path is not a directory."

    def __str__(self):
        return self.message