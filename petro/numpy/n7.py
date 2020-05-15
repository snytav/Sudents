import concurrent.futures
import time


start = time.perf_counter()


def do_something(a,b,c,d,e,dx,dy):
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
    dx = 0.01
    dy = 0.02

    results = [executor.submit(do_something,sec,mec,a,c,l,dx,dy) for sec,mec,a,c,l in zip(secs,tail,abc,cap,list5)]
    for f in concurrent.futures.as_completed(results):
        print(f.result())

finish = time.perf_counter()

print(f'Finished in {round(finish-start,2)} second(s)')
