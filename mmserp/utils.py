
class Suppressor(object):
    """For suppressing prints occuring within functions."""

    def __enter__(self):
        self.stdout = sys.stdout
        sys.stdout = self

    def __exit__(self, type, value, traceback):
        sys.stdout = self.stdout
        if type is not None:
            # Do normal exception handling
            raise

    def write(self, x):
        pass


## use '@timeit' to decorate a function for timing
def timeit(f):
    def timed(*args, **kw):
        ts = time.time()
        for r in range(10): # calls function 100 times
            result = f(*args, **kw)
        te = time.time()
        print('func: %r took: %2.4f sec for %d runs' % (f.__name__, te-ts, 10) )
        return result
    return timed


