from time import time as _time
from sys import stdout as _stdout
from os.path import split as _split
from pkgutil import extend_path
__path__ = extend_path(__path__, "ome")
del extend_path
import settings


def timing(function):
    def wrapper(*args, **kwargs):
        arg_str = str(args)
        if arg_str[-2] == ",": # trailing comma
            arg_str = arg_str[:-2] + ")"
        print "starting %s%s ...\r" % (function.func_name, arg_str),
        _stdout.flush()
        l = len(function.func_name)
        start = _time()
        res = function(*args, **kwargs)
        print "%s%s complete (%.2f sec)"% (function.func_name, arg_str,
            _time() - start)
        return res
    wrapper.func_doc = function.func_doc
    return wrapper

def notiming(function):
    return function


