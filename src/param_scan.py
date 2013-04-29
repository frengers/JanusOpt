from scoop import futures
import numpy as np


import objective_function as fx

if __name__ == '__main__':
    

    jobs = []
    for tc in np.linspace(80,120,5):
        for tw in np.linspace(10,30,5):
            cmd_list = []
            cmd_list.append('--tauc='+str(tc))
            cmd_list.append('--taucwepp='+str(tw))
            jobs.append(cmd_list)
      
    res = list(futures.map(fx.run, jobs))
    
    print "LENGTH RES = ", len(res)
    # writes output to final_output
    fx.collect_results(res)
    
