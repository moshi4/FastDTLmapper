import datetime as dt
import functools
import time


def print_runtime(func):
    """Print runtime decorator"""

    @functools.wraps(func)
    def wrap(*args, **kargs):
        start_time = time.time()
        ret = func(*args, **kargs)
        end_time = time.time()

        def unixtime_to_datestr(unixtime: float) -> str:
            return dt.datetime.fromtimestamp(unixtime).strftime("%Y/%m/%y %H:%M:%S")

        elapsed_time_s = end_time - start_time
        elapsed_time_h = elapsed_time_s / 3600
        print(
            f"Finished: {unixtime_to_datestr(end_time)} "
            + f"(Elapsed time = {elapsed_time_h:.3f}[h])"
        )
        return ret

    return wrap
