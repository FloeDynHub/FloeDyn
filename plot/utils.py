import time
import os

bcolors = {
    "HEADER" : '\033[95m',
    "OKBLUE" : '\033[94m',
    "OKGREEN" : '\033[92m',
    "WARNING" : '\033[93m',
    "FAIL" : '\033[91m',
    "ENDC" : '\033[0m',
    "BOLD" : '\033[1m',
    "UNDERLINE" : '\033[4m'
}

def colored(s, key):
    return bcolors.get(key) + s + bcolors.get("ENDC")

def timeit(func):
    """decorator to measure function execution time"""
    def timed(*args, **kw):
        ts = time.time()
        result = func(*args, **kw)
        te = time.time()
        print(colored('%r %2.2f s' % (func.__name__, te-ts), "OKGREEN"))
        return result

    return timed


def filename_without_extension(path):
    return os.path.splitext(os.path.basename(path))[0]

def get_unused_path(path):
    base, ext = os.path.splitext(path)
    i = 1
    while os.path.exists(path):
        path = "{}_{}{}".format(base, i, ext)
        i+=1
    return path

def check_path_existence(path):
    return os.path.isdir(path)    

def mkdir_path(path):
    return os.mkdir(path)    