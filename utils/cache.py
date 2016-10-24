import os
import simplejson as json  # import json
# import cPickle as pickle  # import pickle
from functools import wraps
import numpy
import datetime
import base64

"""
Create a cache for a function using

```
@money(key_func, cache_directory)
def f(x):
    return y
```

Similar to https://docs.python.org/3/library/functools.html#functools.lru_cache

Use numpy_key where x is a numpy array.

Using json instead of pickle appears to be faster, but requires custom encoder and decoder.

michael.james.asher@gmail.com
"""

"""
Custom json encoder allows sevearal data types to be saved as json
"""

class DateNumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.ndarray):
            return {'type': 'numpy.ndarray', 'data_type': str(obj.dtype), 'data_shape': obj.shape, 'data': base64.b64encode(obj)}
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
        return numpy.frombuffer(base64.b64decode(dct['data']), dtype=dct['data_type']).reshape(dct['data_shape'])
    elif isinstance(dct, dict) and 'type' in dct and dct['type'] == 'datetime.datetime':
        return datetime.datetime.strptime(dct['data'], "%Y-%m-%dT%H:%M:%S.%fZ")
    elif isinstance(dct, dict) and 'type' in dct and dct['type'] == 'datetime.date':
        return datetime.datetime.strptime(dct['data'], "%Y-%m-%d").date()
    return dct

"""
Some useful key functions
"""

def numpy_key(original_x, *args, **kwargs):
    x = original_x.copy()
    x.flags.writeable = False
    key = str(hash(x.data))
    return key

def join_key(*args, **kwargs):
    return str(hash(''.join([json.dumps(a, cls=DateNumpyEncoder, ignore_nan=True) for a in args + tuple(kwargs.values())])))

def concat_key(*args, **kwargs):
    return '_'.join(str(a) for a in args + tuple(kwargs.values()))


"""
A decorator to cache a some_function(args) by saving each run in <cache_dir>/some_function/<key_func(args)>

"""


def money(key_func=join_key, cache_directory=None):

    def decorator(func):
        if cache_directory is None:
            local_cache_directory = os.path.join(os.environ.get('CACHE_DIR', '/tmp'), 'cache-' + func.__module__ + '.' + func.__name__)
        else:
            local_cache_directory = cache_directory

        @wraps(func)
        def f_with_cache(*args, **kwargs):

            key = key_func(*args, **kwargs)

            if not os.path.exists(local_cache_directory):
                os.makedirs(local_cache_directory)

            cache_file = os.path.join(local_cache_directory, key)

            if os.path.exists(cache_file):
                print "WARNING: using cached " + func.__module__ + '.' + func.__name__
                with open(cache_file, "rb") as file:
                    y = json.load(file, object_hook=DateNumpyDecoder)
                    # y = pickle.load(file)
                return y
            else:
                y = func(*args, **kwargs)
                with open(cache_file, "wb") as file:
                    json.dump(y, file, cls=DateNumpyEncoder, ignore_nan=True)
                    # pickle.dump(y, file, protocol=pickle.HIGHEST_PROTOCOL)
                return y

        return f_with_cache

    return decorator


"""
demo
"""

if __name__ == '__main__':

    import time
    import shutil

    def slow_function(x, a, b):
        time.sleep(0.003)
        return 2*x

    rng = numpy.random.RandomState( 0 )

    xs = rng.normal(loc=0.0, scale=1.0, size=(1000,1))

    t0 = time.time()
    ys = [slow_function(x, 5, 10) for x in xs]
    print time.time()-t0, " sec. for ", len(xs), " samples *without* cache"

    cache_directory = os.path.join(os.environ.get('CACHE_DIR', '/tmp'), 'test_func')
    @money(key_func=numpy_key, cache_directory=cache_directory)
    def cached_slow_function(x, a, b):
        time.sleep(0.003)
        return 2 * x

    t0 = time.time()
    ys_cached = [cached_slow_function(x, 5, 10) for x in xs]
    print time.time()-t0, " sec. for ", len(xs), " samples filling cache"

    assert numpy.allclose(ys, ys_cached)

    t0 = time.time()
    ys = [cached_slow_function(x, 5, 10) for x in xs]
    print time.time()-t0, " sec. for ", len(xs), " samples with cache"

    # clean up demo
    shutil.rmtree(cache_directory)