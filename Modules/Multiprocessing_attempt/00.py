import multiprocessing as mp
import time

def f(x):
    time.sleep(1)
    return x * x


if __name__ == '__main__':  
    with mp.Pool(processes=5) as pool:
        results = pool.map(f, range(10))
        print(results)
    print("Check")