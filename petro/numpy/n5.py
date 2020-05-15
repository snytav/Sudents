import concurrent.futures
import time


start = time.perf_counter()


def do_something(a,b,c,d,e):
    dx = 0.05
    dy = 0.05
    print(f'Sleeping {a} second(s)... b  {b} c{c} d {d} e {e} dx {dx} dy {dy}')
    seconds = a
    time.sleep(seconds)
    return seconds*seconds

with concurrent.futures.ProcessPoolExecutor() as executor:
    secs = [5,4,3,2,1]
    tail = [-1,-2,-3,-4,-5]
    abc  = ['a','b','c','d','e']
    cap  = ['P','Q','X','Y','Z']
    list5  = ['AA','BB','CC','DD','EE']

    results = [executor.submit(do_something,sec,mec,a,c,l) for sec,mec,a,c,l in zip(secs,tail,abc,cap,list5)]
    for f in concurrent.futures.as_completed(results):
        print(f.result())

finish = time.perf_counter()

print(f'Finished in {round(finish-start,2)} second(s)')
