import json
import base64
import numpy
import datetime
import os
import cPickle as pickle
import multiprocessing

mutex = multiprocessing.Lock()

"""
Create a cache for a function using

Similar to https://docs.python.org/3/library/functools.html#functools.lru_cache

Use numpy_key where x is a numpy array.

~~Using json instead of pickle appears to be faster, but requires custom encoder and decoder.~~
Pickle with protocol=2 is fast.

michael.james.asher@gmail.com
"""

# Useful key function may want to replace with something more intuitive if your args are short enough for file names
def join_key(*args, **kwargs):
    return str(hash(pickle.dumps({
            "args": args,    
            "kwargs": kwargs   
    }, protocol=pickle.HIGHEST_PROTOCOL)))


def money(func, key_func=join_key, cache_directory=None):
        if cache_directory is None:
            local_cache_directory = os.path.join(os.environ.get('CACHE_DIR', '/tmp'),  func.__module__ + '.' + func.__name__)
        else:
            local_cache_directory = cache_directory

        # @wraps(func) if you want a decorator
        def f_with_cache(*args, **kwargs):
            key = key_func(*args, **kwargs)

            mutex.acquire()
            if not os.path.exists(local_cache_directory):
                os.makedirs(local_cache_directory)
            mutex.release()

            cache_file = os.path.join(local_cache_directory, key)

            if os.path.exists(cache_file + "_input") and os.path.exists(cache_file + "_output"):
                # deals with broken cached files and hash collisions
                try:
                    with open(cache_file + "_input", "rb") as file:
                        y_input = file.read()
                    
                    assert y_input == pickle.dumps({"args": args, "kwargs": kwargs}, protocol=pickle.HIGHEST_PROTOCOL)
                    
                    with open(cache_file + "_output", "rb") as file:
                        return pickle.load(file)
                except:
                    import traceback, sys
                    print "=" * 50
                    traceback.print_exc(file=sys.stdout)
                    print "=" * 50
                    pass
            
            with open(cache_file  + "_input", "wb") as file:
                pickle.dump({"args": args, "kwargs": kwargs}, file, protocol=pickle.HIGHEST_PROTOCOL)
            
            output = func(*args, **kwargs)
            
            with open(cache_file  + "_output", "wb") as file:
                pickle.dump(output, file, protocol=pickle.HIGHEST_PROTOCOL)
            
            return output

        return f_with_cache


class DateNumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.ndarray):
            return {'type': 'numpy.ndarray', 'data_type': str(obj.dtype), 'data_shape': obj.shape, 'data': base64.b64encode(obj.astype(dtype=obj.dtype, order="C_CONTIGUOUS"))}
        elif isinstance(obj, datetime.datetime):
            return {'type': 'datetime.datetime', 'data': obj.strftime("%Y-%m-%dT%H:%M:%S.%fZ")}
        elif isinstance(obj, datetime.date):
            return {'type': 'datetime.date', 'data': obj.strftime("%Y-%m-%d")}
        elif isinstance(obj, datetime.timedelta):
            return (datetime.datetime.min + obj).date().isoformat()
        else:
            return super(DateNumpyEncoder, self).default(obj)


def  DateNumpyDecoder(dct):
    if isinstance(dct, dict) and 'type' in dct and dct['type'] == 'numpy.ndarray':
        a = numpy.frombuffer(base64.b64decode(dct['data']), dtype=dct['data_type']).reshape(dct['data_shape'])
        a.setflags(write=1)
        return a
    elif isinstance(dct, dict) and 'type' in dct and dct['type'] == 'datetime.datetime':
        return datetime.datetime.strptime(dct['data'], "%Y-%m-%dT%H:%M:%S.%fZ")
    elif isinstance(dct, dict) and 'type' in dct and dct['type'] == 'datetime.date':
        return datetime.datetime.strptime(dct['data'], "%Y-%m-%d").date()
    return dct



if __name__ == '__main__':
    import time

    def multiply_today(a, factor):
        time.sleep(1.)
        return factor * a, datetime.date.today()

    a = numpy.random.random((300,700))

    t0 = time.time()
    b = money(multiply_today)(a, 7.)
    print time.time() - t0

    t0 = time.time()
    b2 = money(multiply_today)(a, 7.)
    print time.time() - t0

    assert numpy.all(b[0] == b2[0])
    assert b[1] == b2[1]
